import cvxopt
import itertools
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.spatial import ConvexHull, Delaunay, KDTree
from scipy.spatial.qhull import QhullError
from scipy.sparse import coo_matrix

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

    Arguments
        z: array, shape (ndim,)
        equations: array shape (nfacets, ndim + 1)

    Returns
        x: array, shape (ndim,)
    """
    G = np.eye(len(z), dtype=float)
    a = np.array(-z, dtype=float)
    C = np.array(equations[:, :-1], dtype=float)
    b = np.array(-equations[:, -1], dtype=float)

    solution = cvxopt.solvers.qp(cvxopt.matrix(G), cvxopt.matrix(a), cvxopt.matrix(C), cvxopt.matrix(b))
    return np.asfarray(solution["x"]).squeeze()

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
    for index, face in enumerate(hull.simplices):
        if edge[0] in face or edge[1] in face:
            related_faces.append(index)

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
        normal = hull.equations[face][0:3] # Must point outward from hull.
        p0 = hull.points[hull.simplices[face][0]] # Get any vertex of the face.

        A.append(normal)
        b.append(np.dot(normal, p0))
        c += normal

    A = -np.asfarray(A)
    b = -np.asfarray(b)
    c = np.asfarray(c)
    return A, b, c

def compute_volume(hull, faces, point):
    volume = 0

    # Sum the volume of the tetrahedra made by the faces and the point.
    for face in faces:
        normal = hull.equations[face][0:3] # Must point outward from hull.
        p0 = hull.points[hull.simplices[face][0]] # Get any vertex of the face.
        volume += np.abs(np.dot(normal, point - p0)) / 6.0

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
            continue

        # Determine how much volume the new point adds.
        point = np.asfarray(solution["x"]).squeeze()
        volume = compute_volume(hull, get_related_faces(hull, edge), point)

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

    # The RMSE is zero if every point is inside the hull.
    if (weight_sum == 0):
        return 0

    mse = np.sum(distances) / weight_sum
    return np.sqrt(mse)

def project_hull_to_rgb_cube(hull):
    hull.points[hull.vertices] = np.clip(hull.points[hull.vertices], 0, 1)

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
    project_hull_to_rgb_cube(hull)
    previous_size = len(hull.vertices)

    if len(hull.vertices) <= target_size:
        return hull

    while True:
        hull = remove_edge(hull)
        project_hull_to_rgb_cube(hull)

        if len(hull.vertices) <= target_size or len(hull.vertices) == previous_size:
            return hull

        previous_size = len(hull.vertices)

def get_automatically_simplified_hull(points, error_threshold=5.0/255.0):
    hull = ConvexHull(points)
    project_hull_to_rgb_cube(hull)
    previous_hull = hull
    previous_size = len(hull.vertices)

    unique_points, counts = np.unique(points, axis=0, return_counts=True)

    while True:
        hull = remove_edge(hull)
        project_hull_to_rgb_cube(hull)

        if len(hull.vertices) <= 10:
            rmse = calculate_rmse(hull, unique_points, counts)
            if rmse > error_threshold:
                return previous_hull

        if len(hull.vertices) == 4 or len(hull.vertices) == previous_size:
            return hull

        previous_hull = hull
        previous_size = len(hull.vertices)

def compute_neighbor_centers(vertices, points, tree):
    """
    Finds the nearest neighbors of vertices in points and computes the center of
    each of those groupds of neighbors. The tree parameter should be a SciPy
    KDTree constructed using KDTree(points).
    """

    # TODO: In the paper they disallow duplicate neighbors.
    # TODO: Automatically determine a reasonable number of neighbors.
    _, neighbors = tree.query(vertices, 3)
    return np.mean(points[neighbors], axis=-2)

def reconstruction_loss(hull, points, weights):
    # Stores Euclidean distances.
    distances = np.zeros(len(points))

    # Delaunay tessellation for checking which points are in hull.
    tri = Delaunay(hull.points[hull.vertices])
    simplices = tri.find_simplex(points)

    for i in range(len(points)):
        # Skip points inside the hull.
        if simplices[i] != -1:
            continue

        projection = project_to_hull(points[i], hull.equations)
        distances[i] = np.sqrt(np.dot(points[i] - projection, points[i] - projection))
        distances[i] *= weights[i]

    # TODO: Note that we normalize over ALL distances, even those that are 0
    # from points inside the hull. This is different from compute_rmse() but it
    # does match the paper. Do we want to change this?
    return np.sum(distances) / np.sum(weights)

def representative_loss(hull, points, tree):
    vertices = hull.points[hull.vertices]
    centers = compute_neighbor_centers(vertices, points, tree)
    distances = np.linalg.norm(vertices - centers, axis=1)
    return np.sum(distances) / len(vertices)

def energy_function(hull, points, tree, factor):
    unique_points, counts = np.unique(points, axis=0, return_counts=True)
    rec = reconstruction_loss(hull, unique_points, counts)
    rep = representative_loss(hull, points, tree)
    return factor * rec + rep

def optimize_vertex_positions(vertices, points):
    """
    Optimizes a palette's vertex positions following the method described in
    doi.org/10.1111/cgf.13812.
    """
    unique_points, counts = np.unique(points, axis=0, return_counts=True)
    tree = KDTree(points) # Structure for quickly finding nearest neighbors.

    # We only do one cycle.
    for vertex in range(len(vertices)):
        # Locate the current vertex and the center of its nearest neighbors.
        position = vertices[vertex].copy()
        center = compute_neighbor_centers(vertices, points, tree)[vertex]

        # Don't optimize if vertex is perceptually identical to neighbor center.
        distance = np.linalg.norm(position - center)
        if (distance < 1.0 / 255.0):
            continue

        # The energy function to optimize. Specifies a trade-off between
        # reconstruction loss and representation loss.
        def energy_function(t):
            new_position = (1 - t) * position + t * center
            vertices[vertex] = new_position
            hull = ConvexHull(vertices)

            rec = reconstruction_loss(hull, unique_points, counts)
            rep = representative_loss(hull, points, tree)
            return 5 * rec + rep

        # Start with a coarse, global search.
        t_count = 10 # "k" in the paper.
        ts = np.linspace(-0.5, 1.0, t_count)
        energies = [energy_function(t) for t in ts]
        min_index  = np.argmin(energies)

        # Follow up with a fine, local search.
        bounds = (ts[max(min_index - 1, 0)], ts[min(min_index + 1, t_count - 1)])
        result = minimize_scalar(energy_function, bounds=bounds, method="bounded")

        # Update the palette vertex position.
        final_t = result["x"]
        final_position = (1 - final_t) * position + final_t * center
        vertices[vertex] = final_position

    return vertices
