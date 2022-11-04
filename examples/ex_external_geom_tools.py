"""Script that shows how to use two external packages, `shapely` and `scipy.spatial`
to deal with geometric shapes and related computations (lengths, areas, volumes).

Recall: `dualGPy` can deal with standard shapes such as parallelograms and triangles
in 2D and tetra, orthogonal hexahedra and wedges in 3D. On the contrary, those two
library, whenever they can be used, they do not have such constraint.

`shapely` only deals with 2D shapes. Even though it admits points with 3 coordinates,
it discards the z-component, thus it cannot be used in 3D.

We take advantage of the function to compute convex hulls of points in
`scipy.spatial`. It works with the hypothesis that the faces are their own convex
hull. It cannot deal with lines, so no segment support in 2D. When using 3D points,
it fails if the hull is actually 2D. This means that it cannot be used for faces.

That being said, even with these two packages, we lack support for computing face
surfaces in 3D, that is not that bad since `dualGPy` does a honest job with this kind
of computation.
"""

from dualGPy.Mesh import Mesh2D, Mesh3D
import numpy as np

def max_error(a,b): return np.max(np.abs(a - b))

n = 6
anisotropic = False

# SET REFERENCE
m = Mesh2D(n, anisotropic)
m.get_boundary_faces()
m.ComputeGeometry()


# USING SHAPELY FOR 2D SHAPES
from shapely.geometry import Polygon, LineString
sh_areas   = np.array([ Polygon(m.mesh.points[c,:]).area for c in m.cells ],
                   dtype = np.float64)
sh_lengths = np.array([ LineString(m.mesh.points[s,:]).length for f in m.faces.values() for s in f ],
                   dtype = np.float64)

print("[2D - Shapely] Max error on surfaces:", max_error(sh_areas, m.volume))
print("[2D - Shapely] Max error on segments:", max_error(sh_lengths, m.area))
print("")

# USING (PARTIALLY) SCIPY.SPATIAL FOR 2D SHAPES
# We take advantage of ConvexHull. This works because we suppose that the shapes are convex
# Since it cannot deal with length, the edges cannot be dealt with
from scipy.spatial import ConvexHull
ch_areas = np.array([ ConvexHull(m.mesh.points[c,:]).area for c in m.cells ],
                    dtype = np.float64)
print("[2D - SciPy] Max error on surfaces:", max_error(sh_areas, m.volume))
print("")

# USING (PARTIALLY) SCIPY.SPATIAL FOR 3D SHAPES
del(m)
m = Mesh3D(n, anisotropic)
m.get_boundary_faces()
m.ComputeGeometry()

ch_volumes = np.array([ ConvexHull(m.mesh.points[c,:]).volume for c in m.cells ],
                    dtype = np.float64)
# The following fails...
try:
    ch_areas   = np.array([ ConvexHull(m.mesh.points[s,:]).length for f in m.faces.values() for s in f ],
                       dtype = np.float64)
except:
    print('{ WARNING: As expected, the computation of a (D-1)-dimensional hull fails }')
print("[3D - Cartesian - SciPy] Max error on volumes:",  max_error(ch_volumes, m.volume))
