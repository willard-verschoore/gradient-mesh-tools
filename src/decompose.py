import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.qhull import QhullError
from scipy.sparse import coo_matrix

def get_palette(rgb_data):
    try:
        rgb_hull = ConvexHull(rgb_data)
        rgb_hull_vertices = rgb_hull.points[rgb_hull.vertices]
    except QhullError as error:
        print("Encountered Qhull error:")
        print(error)
        return np.empty(0, np.float32)

    print(f"Found palette:\n{rgb_hull_vertices}")
    return np.float32(rgb_hull_vertices)

def get_palette_indices(rgb_data):
    try:
        rgb_hull = ConvexHull(rgb_data)
    except QhullError as error:
        print("Encountered Qhull error:")
        print(error)
        return np.empty(0, np.uint32)

    vertex_count = len(rgb_hull.vertices)

    rgb_hull_indices = []
    for i in range(vertex_count - 1):
        for j in range(i + 1, vertex_count):
            vertex_i = rgb_hull.vertices[i]
            vertex_j = rgb_hull.vertices[j]

            # We have connections between vertices within the same simplex.
            for simplex in rgb_hull.simplices:
                if (vertex_i in simplex) and (vertex_j in simplex):
                    rgb_hull_indices.append(i)
                    rgb_hull_indices.append(j)
                    break # Avoid counting the same connection multiple times.

    return np.uint32(rgb_hull_indices)

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
        mask = (b>=-1e-8).all(axis=1)
        barycoords[mask] = 0.
        barycoords[np.ix_(mask,s)] = b[mask]
    return barycoords
