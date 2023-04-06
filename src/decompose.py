import cvxopt
import itertools
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.qhull import QhullError
from scipy.sparse import coo_matrix
from quadprog import solve_qp

def get_palette(rgb_data, target_size):
    try:
        rgb_hull = get_simplified_hull(rgb_data, target_size)
        rgb_hull_vertices = rgb_hull.points[rgb_hull.vertices]
        rgb_hull_indices = get_hull_indices(rgb_hull)
    except QhullError as error:
        print("Encountered Qhull error:")
        print(error)
        return np.empty(0, np.float32), np.empty(0, np.uint32)

    print(f"Found palette:\n{rgb_hull_vertices}")
    return np.float32(np.clip(rgb_hull_vertices, 0, 1)), np.uint32(rgb_hull_indices)

def get_weights(rgbxy_data, palette):
    try:
        rgbxy_hull_vertices = rgbxy_data[ConvexHull(rgbxy_data).vertices]
        w_rgbxy = delaunay_coordinates(rgbxy_hull_vertices, rgbxy_data)
        w_rgb = star_coordinates(palette, rgbxy_hull_vertices[:,:3])
    except QhullError as error:
        print("Encountered Qhull error:")
        print(error)
        return np.empty(0, np.float32)

    weights = w_rgbxy.dot(w_rgb)
    return np.float32(weights)

def delaunay_coordinates(vertices, data): # Adapted from Gareth Rees
    # Compute Delaunay tessellation.
    tri = Delaunay(vertices)
    # Find the tetrahedron containing each target (or -1 if not found).
    simplices = tri.find_simplex(data, tol=1e-6)
    assert (simplices != -1).all() # data contains outside vertices.
    # Affine transformation for simplex containing each datum.
    X = tri.transform[simplices, :data.shape[1]]
    # Offset of each datum from the origin of its simplex.
    Y = data - tri.transform[simplices, data.shape[1]]
    # Compute the barycentric coordinates of each datum in its simplex.
    b = np.einsum('...jk,...k->...j', X, Y)
    barycoords = np.c_[b,1-b.sum(axis=1)]
    # Return the weights as a sparse matrix.
    rows = np.repeat(np.arange(len(data)).reshape((-1,1)), len(tri.simplices[0]), 1).ravel()
    cols = tri.simplices[simplices].ravel()
    vals = barycoords.ravel()
    return coo_matrix((vals,(rows,cols)), shape=(len(data),len(vertices))).tocsr()

def star_coordinates(vertices, data):
    ## Find the star vertex
    star = np.argmin(np.linalg.norm(vertices, axis=1))
    ## Make a mesh for the palette
    hull = ConvexHull(vertices)

    # Project outside data points onto the convex hull.
    tri = Delaunay(vertices) # For checking if data point is in hull.
    for i in range(data.shape[0]):
        if tri.find_simplex(data[i]) == -1: # Data point not in hull.
            data[i] = project_to_hull(data[i], hull.equations)

    ## Star tessellate the faces of the convex hull
    simplices = [[star] + list(face) for face in hull.simplices if star not in face]
    barycoords = -1 * np.ones((data.shape[0], len(vertices)))
    ## Barycentric coordinates for the data in each simplex
    for s in simplices:
        s0 = vertices[s[:1]]

        try:
            b = np.linalg.solve((vertices[s[1:]]-s0).T, (data-s0).T).T
        except np.linalg.LinAlgError:
            b = np.linalg.lstsq((vertices[s[1:]]-s0).T, (data-s0).T, rcond=None)[0].T

        b = np.append(1-b.sum(axis=1)[:,None], b, axis=1)
        ## Update barycoords whenever data is inside the current simplex (with threshold).
        mask = (b>=-1e-4).all(axis=1)
        barycoords[mask] = 0.
        barycoords[np.ix_(mask,s)] = b[mask]
    return barycoords

def project_to_hull(z, equations):
    """
    Project `z` to a convex hull defined by the hyperplane equations of its
    facets.

    Based on this Stack Overflow answer: https://stackoverflow.com/a/57631915

    TODO: Consider using cvxopt.solvers.qp instead since we already use cvxopt
    for convex hull simplification.

    Arguments
        z: array, shape (ndim,)
        equations: array shape (nfacets, ndim + 1)

    Returns
        x: array, shape (ndim,)
    """
    G = np.eye(len(z), dtype=float)
    a = np.array(z, dtype=float)
    C = np.array(-equations[:, :-1], dtype=float)
    b = np.array(equations[:, -1], dtype=float)
    x, f, xu, itr, lag, act = solve_qp(G, a, C.T, b, meq=0, factorized=True)
    return x

def get_hull_indices(hull):
    vertex_count = len(hull.vertices)

    indices = []
    for i in range(vertex_count - 1):
        for j in range(i + 1, vertex_count):
            vertex_i = hull.vertices[i]
            vertex_j = hull.vertices[j]

            # We have connections between vertices within the same simplex.
            for simplex in hull.simplices:
                if (vertex_i in simplex) and (vertex_j in simplex):
                    indices.append(i)
                    indices.append(j)
                    break # Avoid counting the same connection multiple times.

    return indices

def get_hull_edges(hull):
    edges = []

    # Edges are all pairs of vertices sharing a face.
    for face in hull.simplices:
        edges.extend(itertools.combinations(face, 2))

    # Adjacent faces produce duplicate edges which should be removed.
    return np.unique(edges, axis=0)

def get_related_faces(hull, edge):
    related_faces = []

    # Find all faces containing at least one of the edges vertices.
    for face in hull.simplices:
        if edge[0] in face or edge[1] in face:
            p0 = hull.points[face[0]]
            p1 = hull.points[face[1]]
            p2 = hull.points[face[2]]
            related_faces.append([p0, p1, p2])

    return related_faces

def get_solver_input(hull, edge):
    """
    Find A, b and c such that minimizing c*x subject to A*x <= b yields a point
    x which will create a positive volume tetrahedron for each face related to
    `edge`. A face is related to an edge if it contains at least one of the
    edge's endpoints. `hull` should be a convex hull in 3D space such that the
    faces are 2D triangles which can be combined with a new point x to form a 3D
    tetrahedron.
    """
    A = []
    b = []
    c = np.zeros(3)

    # Construct the constraints based on the related faces.
    for face in get_related_faces(hull, edge):
        # Compute the face normal vector.
        normal = np.cross(face[1] - face[0], face[2] - face[0])
        normal = normal / np.sqrt(np.dot(normal, normal))

        A.append(normal)
        b.append(np.dot(normal, face[0]))
        c += normal

    A = -np.asfarray(A)
    b = -np.asfarray(b)
    c = np.asfarray(c)
    return A, b, c

def compute_volume(faces, point):
    volume = 0

    # Sum the volume of the tetrahedra made by the faces and the point.
    for face in faces:
        normal = np.cross(face[1] - face[0], face[2] - face[0])
        volume += np.abs(np.dot(normal, point - face[0])) / 6.0

    return volume

def remove_edge(hull):
    min_volume = float("inf")
    min_point = np.empty(3)
    min_edge = [-1, -1]

    for edge in get_hull_edges(hull):
        A, b, c = get_solver_input(hull, edge)
        cvxopt.solvers.options["show_progress"] = False
        cvxopt.solvers.options["glpk"] = dict(msg_lev="GLP_MSG_OFF")
        solution = cvxopt.solvers.lp(cvxopt.matrix(c), cvxopt.matrix(A), cvxopt.matrix(b), solver="glpk")

        if solution["status"] != "optimal":
            # TODO: Investigate in which cases (if any) non-optimal solutions
            # are still useable. Perhaps as a fallback when there are no optimal
            # solutions at all?
            continue

        # Determine how much volume the new point adds.
        point = np.asfarray(solution["x"]).squeeze()
        volume = compute_volume(get_related_faces(hull, edge), point)

        # Keep track of the point adding the least volume.
        if volume < min_volume:
            min_volume = volume
            min_point = point
            min_edge = edge

    if min_volume == float("inf"):
        print("Unable to remove edge!")
        return hull

    # Remove the determined edge by replacing its endpoints with the new point.
    new_hull_vertices = [vertex for vertex in hull.vertices if vertex not in min_edge]
    new_hull_points = hull.points[new_hull_vertices]
    new_hull_points = np.append(new_hull_points, min_point.reshape((1, 3)), axis=0)

    return ConvexHull(new_hull_points)

def calculate_rmse(hull, points, weights):
    """
    Calculate the Root Mean Square Error (RMSE) of every point in `points`'s
    distance to `hull`. A point inside the hull has distance 0.

    The `weights` array specifies the contribution of each point to the MSE part
    of the error. It can be interpreted as specifying the number of occurrences
    of each point (if the points are unique).
    """

    # Stores square distances.
    distances = np.zeros(len(points))

    # Stores the sum of used weights (not all weights are used).
    weight_sum = 0

    # Delaunay tessellation for checking which points are in hull.
    tri = Delaunay(hull.points[hull.vertices])
    simplices = tri.find_simplex(points)

    for i in range(len(points)):
        # Skip points inside the hull.
        if simplices[i] != -1:
            continue

        projection = project_to_hull(points[i], hull.equations)
        distances[i] = np.dot(points[i] - projection, points[i] - projection)
        distances[i] *= weights[i]
        weight_sum += weights[i]

    mse = np.sum(distances) / weight_sum
    return np.sqrt(mse)

def get_simplified_hull(points, target_size):
    if target_size == 0:
        return get_automatically_simplified_hull(points)
    else:
        return get_manually_simplified_hull(points, target_size)

def get_manually_simplified_hull(points, target_size):
    # A 3D convex hull has at least 4 vertices
    if target_size < 4:
        target_size = 4

    hull = ConvexHull(points)
    previous_size = len(hull.vertices)

    if len(hull.vertices) <= target_size:
        return hull

    while True:
        hull = remove_edge(hull)

        if len(hull.vertices) <= target_size or len(hull.vertices) == previous_size:
            return hull

        previous_size = len(hull.vertices)

def get_automatically_simplified_hull(points, error_threshold=5.0/255.0):
    hull = ConvexHull(points)
    previous_hull = hull
    previous_size = len(hull.vertices)

    unique_points, counts = np.unique(points, axis=0, return_counts=True)

    while True:
        hull = remove_edge(hull)

        if len(hull.vertices) <= 10:
            rmse = calculate_rmse(hull, unique_points, counts)
            if rmse > error_threshold:
                return previous_hull

        if len(hull.vertices) == 4 or len(hull.vertices) == previous_size:
            return hull

        previous_hull = hull
        previous_size = len(hull.vertices)
