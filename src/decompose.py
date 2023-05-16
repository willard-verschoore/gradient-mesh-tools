import cvxopt
import itertools
import numpy as np
from enum import Enum
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

class WeightType(Enum):
    RGBXY = 0
    MVC = 1

def get_weights(rgbxy_data, palette, weight_type=WeightType.RGBXY):
    if weight_type == WeightType.RGBXY.value:
        weights = get_rgbxy_weights(rgbxy_data, palette)
    elif weight_type == WeightType.MVC.value:
        weights = get_mvc_weights(rgbxy_data[:, :3], palette)
    else:
        print(f"Unknown weight type: {weight_type}")
        weights = np.empty(0, np.float32)

    return weights

def get_rgbxy_weights(rgbxy_data, palette):
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

def get_mvc_weights(rgb_data, palette):
    """
    Determines the Mean Value Coordinates (MVCs) with respect to the palette for
    each input color. The palette colors are interpreted as the vertices of a
    convex hull. This yields a triangular mesh, which allows us to use the
    techique described in doi.org/10.1145/1073204.1073229.
    """
    hull = ConvexHull(palette)
    vertices = hull.points[hull.vertices]
    simplices = consistent_winding_order_simplices(hull)
    weights = []

    # Project outside data points onto the convex hull.
    rgb_data = project_outside_points_to_hull(rgb_data, hull)

    for color in rgb_data:
        weight_vector = np.zeros(len(vertices))

        ds = np.linalg.norm(vertices - color, axis=1) # Distances to vertices.

        # If a color is very close to a vertex we set that vertex's weight to 1.
        min_index = np.argmin(ds)
        if ds[min_index] < 1e-6:
            weight_vector[min_index] = 1.0
            weights.append(weight_vector)
            continue

        us = (vertices - color) / ds[:, np.newaxis] # Unit directions to vertices.

        # Loop over every triangle of the hull.
        for simplex in simplices:
            d_1 = ds[simplex[0]]
            d_2 = ds[simplex[1]]
            d_3 = ds[simplex[2]]

            u_1 = us[simplex[0]]
            u_2 = us[simplex[1]]
            u_3 = us[simplex[2]]

            l_1 = np.linalg.norm(u_2 - u_3)
            l_2 = np.linalg.norm(u_3 - u_1)
            l_3 = np.linalg.norm(u_1 - u_2)

            theta_1 = 2 * np.arcsin(l_1 / 2)
            theta_2 = 2 * np.arcsin(l_2 / 2)
            theta_3 = 2 * np.arcsin(l_3 / 2)

            # Use barycentric coordinates if the color lies on the hull triangle.
            h = (theta_1 + theta_2 + theta_3) / 2
            if np.pi - h < 1e-6:
                coordinates = barycentric_coordinates(hull.points[simplex], color)
                weight_vector = np.zeros(len(vertices))
                weight_vector[simplex] = coordinates
                break

            c_1 = (2 * np.sin(h) * np.sin(h - theta_1)) / (np.sin(theta_2) * np.sin(theta_3)) - 1
            c_2 = (2 * np.sin(h) * np.sin(h - theta_2)) / (np.sin(theta_3) * np.sin(theta_1)) - 1
            c_3 = (2 * np.sin(h) * np.sin(h - theta_3)) / (np.sin(theta_1) * np.sin(theta_2)) - 1

            sign = np.sign(np.linalg.det(np.vstack((u_1, u_2, u_3))))
            s_1 = sign * np.sqrt(1 - c_1 * c_1)
            s_2 = sign * np.sqrt(1 - c_2 * c_2)
            s_3 = sign * np.sqrt(1 - c_3 * c_3)

            # Skip if the color lies outside but on the same plane as the hull triangle.
            if np.abs(s_1) < 1e-6 or np.abs(s_2) < 1e-6 or np.abs(s_3) < 1e-6:
                continue

            w_1 = (theta_1 - c_2 * theta_3 - c_3 * theta_2) / (d_1 * np.sin(theta_2) * s_3)
            w_2 = (theta_2 - c_3 * theta_1 - c_1 * theta_3) / (d_2 * np.sin(theta_3) * s_1)
            w_3 = (theta_3 - c_1 * theta_2 - c_2 * theta_1) / (d_3 * np.sin(theta_1) * s_2)

            weight_vector[simplex[0]] += w_1
            weight_vector[simplex[1]] += w_2
            weight_vector[simplex[2]] += w_3

        weight_vector /= np.sum(weight_vector)
        weights.append(weight_vector)

    return np.float32(weights)

def barycentric_coordinates(triangle, point):
    edge_1 = triangle[2] - triangle[1]
    edge_2 = triangle[0] - triangle[1]
    edge_3 = point - triangle[1]
    w_0 = np.linalg.norm(np.cross(edge_1, edge_3)) / np.linalg.norm(np.cross(edge_1, edge_2))

    edge_1 = triangle[0] - triangle[2]
    edge_2 = triangle[1] - triangle[2]
    edge_3 = point - triangle[2]
    w_1 = np.linalg.norm(np.cross(edge_1, edge_3)) / np.linalg.norm(np.cross(edge_1, edge_2))

    w_2 = 1.0 - w_0 - w_1
    return np.array([w_0, w_1, w_2])

def consistent_winding_order_simplices(hull):
    simplices = np.copy(hull.simplices)

    # Find two edges for each triangle
    edges_1 = hull.points[hull.simplices[:, 1]] - hull.points[hull.simplices[:, 0]]
    edges_2 = hull.points[hull.simplices[:, 2]] - hull.points[hull.simplices[:, 0]]

    # Infer the triangle normals from the cross product of two edges.
    inferred_normals = np.cross(edges_1, edges_2)
    inferred_normals /= np.linalg.norm(inferred_normals, axis=1)[:, np.newaxis]

    # Compare the inferred normals with the true, outward pointing ones.
    true_normals = hull.equations[:, 0:3] # Outward pointing normals from Qhull.
    differences = true_normals - inferred_normals
    magnitudes = np.einsum("...ij,...ij->...i", differences, differences) # Squared magnitude.

    # If the normals point in opposite directions the squared difference is 4 (which is > 1).
    inconsistencies = np.argwhere(magnitudes > 1).flatten()

    # Swap two indices to reverse the winding order for inconsistent triangles.
    for index in inconsistencies:
        simplices[index, [0, 1]] = simplices[index, [1, 0]]

    return simplices


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
    data = project_outside_points_to_hull(data, hull)

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

def project_points_to_triangles(points, triangles):
    """
    Projects an array of points to each triangle in an array of triangles.
    Copied from: https://stackoverflow.com/a/32529589.
    """
    with np.errstate(all='ignore'):
        # Unpack triangle points.
        p0,p1,p2 = np.asarray(triangles).swapaxes(0,1)

        # Calculate triangle edges.
        e0 = p1-p0
        e1 = p2-p0
        a = np.einsum('...i,...i', e0, e0)
        b = np.einsum('...i,...i', e0, e1)
        c = np.einsum('...i,...i', e1, e1)

        # Calculate determinant and denominator.
        det = a*c - b*b
        inv_det = 1. / det
        denom = a-2*b+c

        # Project to the edges.
        p  = p0-points[:,np.newaxis]
        d = np.einsum('...i,...i', e0, p)
        e = np.einsum('...i,...i', e1, p)
        u = b*e - c*d
        v = b*d - a*e

        # Calculate numerators.
        bd = b+d
        ce = c+e
        numer0 = (ce - bd) / denom
        numer1 = (c+e-b-d) / denom
        da = -d/a
        ec = -e/c

        # Vectorize test conditions.
        m0 = u + v < det
        m1 = u < 0
        m2 = v < 0
        m3 = d < 0
        m4 = (a+d > b+e)
        m5 = ce > bd

        t0 =  m0 &  m1 &  m2 &  m3
        t1 =  m0 &  m1 &  m2 & ~m3
        t2 =  m0 &  m1 & ~m2
        t3 =  m0 & ~m1 &  m2
        t4 =  m0 & ~m1 & ~m2
        t5 = ~m0 &  m1 &  m5
        t6 = ~m0 &  m1 & ~m5
        t7 = ~m0 &  m2 &  m4
        t8 = ~m0 &  m2 & ~m4
        t9 = ~m0 & ~m1 & ~m2

        u = np.where(t0, np.clip(da, 0, 1), u)
        v = np.where(t0, 0, v)
        u = np.where(t1, 0, u)
        v = np.where(t1, 0, v)
        u = np.where(t2, 0, u)
        v = np.where(t2, np.clip(ec, 0, 1), v)
        u = np.where(t3, np.clip(da, 0, 1), u)
        v = np.where(t3, 0, v)
        u *= np.where(t4, inv_det, 1)
        v *= np.where(t4, inv_det, 1)
        u = np.where(t5, np.clip(numer0, 0, 1), u)
        v = np.where(t5, 1 - u, v)
        u = np.where(t6, 0, u)
        v = np.where(t6, 1, v)
        u = np.where(t7, np.clip(numer1, 0, 1), u)
        v = np.where(t7, 1-u, v)
        u = np.where(t8, 1, u)
        v = np.where(t8, 0, v)
        u = np.where(t9, np.clip(numer1, 0, 1), u)
        v = np.where(t9, 1-u, v)

        # Return closest points.
        return (p0.T +  u[:, np.newaxis] * e0.T + v[:, np.newaxis] * e1.T).swapaxes(2,1)

def project_to_hull(points, hull):
    """
    Projects an array of points onto a convex hull. Returns the projections and
    their squared Euclidean distances to the original point.
    """
    # Project the points onto each triangular facet of the hull.
    triangles = hull.points[hull.simplices]
    projections = project_points_to_triangles(points, triangles)

    # Determine the squared Euclidean distance of each point to its projections.
    differences = projections - points[:, np.newaxis, :]
    distances = np.einsum("...ij,...ij->...i", differences, differences)

    # Determine for each point which triangle's projection is closest.
    min_indices = np.argmin(distances, axis=-1)
    min_indices = (np.arange(len(points)), min_indices)

    return projections[min_indices], distances[min_indices]

def project_outside_points_to_hull(points, hull):
    """
    Applies project_to_hull() only to the points outside the hull. Points already
    inside the hull are unaltered.
    """
    projections = np.copy(points)

    # Use Delaunay tessellation to check which points are in the hull.
    tri = Delaunay(hull.points[hull.vertices])
    simplices = tri.find_simplex(points)
    outside_indices = np.argwhere(simplices == -1).flatten()

    # Project outside data points onto the convex hull.
    projections[outside_indices], _ = project_to_hull(projections[outside_indices], hull)

    return projections

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

    # Use Delaunay tessellation to check which points are in the hull.
    tri = Delaunay(hull.points[hull.vertices])
    simplices = tri.find_simplex(points)
    outside_indices = np.argwhere(simplices == -1).flatten()

    # There is no error if all points are inside the hull.
    if (outside_indices.size == 0):
        return 0

    outside_points = points[outside_indices]
    outside_weights = weights[outside_indices]

    _, distances = project_to_hull(outside_points, hull) # Squared distances.
    distances *= outside_weights # Weight the contribution of each point.

    mse = np.sum(distances) / np.sum(outside_weights)
    return np.sqrt(mse)

def project_hull_to_rgb_cube(hull):
    hull.points[hull.vertices] = np.clip(hull.points[hull.vertices], 0, 1)

def simplify_hull(vertices):
    hull = ConvexHull(vertices)
    hull = remove_edge(hull)

    hull_vertices = np.clip(hull.points[hull.vertices], 0, 1)
    hull_indices = get_hull_indices(hull)

    print(f"Found palette:\n{hull_vertices}")
    return np.float32(hull_vertices), np.uint32(hull_indices)

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

    # A 3D convex hull cannot have fewer than 4 vertices
    if len(hull.vertices) <= 4:
        return hull

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
    neighbor_count = len(points) // 1000 # Use the paper recommended 0.1% of points.
    neighbor_count = max(neighbor_count, 3) # Always use at least 3 points.
    _, neighbors = tree.query(vertices, neighbor_count)
    return np.mean(points[neighbors], axis=-2)

def reconstruction_loss(hull, points, weights):
    # Use Delaunay tessellation to check which points are in the hull.
    tri = Delaunay(hull.points[hull.vertices])
    simplices = tri.find_simplex(points)
    outside_indices = np.argwhere(simplices == -1).flatten()

    # There is no error if all points are inside the hull.
    if (outside_indices.size == 0):
        return 0

    outside_points = points[outside_indices]
    outside_weights = weights[outside_indices]

    _, distances = project_to_hull(outside_points, hull) # Squared distances.
    distances = np.sqrt(distances) # Euclidean distances.
    distances *= outside_weights # Weight the contribution of each point.

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

def optimize_vertex_positions(vertices, points, lam=10.0):
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
            return lam * rec + rep

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
