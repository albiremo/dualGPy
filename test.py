import meshio
import dualGPy as dm
# For plotting both the mesh and dual mesh
import matplotlib.pyplot as plt
import numpy as np
# two triangles and one quad
points = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [2.0, 1.0],
    [3.0, 0.0],
    [3.0, 1.0],
    [3.0, 2.0],
    [2.0, 2.0],
    [1.0, 2.0],
    [0.0, 2.0],
    [0.0, -1.0],
    [1.0, -1.0],
    [2.0, -1.0],
    [3.0, -1.0],
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3],[4,6,7,5],[5,7,8,9],[3,5,9,10],[3,10,11,2],[14,15,6,4],[13,14,4,1],[12,13,1,0]]),
]

mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    # Each item in cell data must match the cells array
    cell_sets = {
    "Boundary": {"line": [[0, 2],[2,11],[11,10],[10,9],[9,8],[8,7],[7,6],[15,6],[14,15],[13,14],[12,13],[12,0]]},
}
)
#print(dir(mesh.cells))
#print(mesh.cell_sets["Boundary"])
for key, value in mesh.cell_sets.items():
 for key1,value1 in mesh.cell_sets[key].items():
   print(mesh.cell_sets[key][key1])

#for i in range(len(mesh.cell_sets)):
#	print(mesh.cell_sets[i])

mesh.write(
    "foo.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)
dm.get_dual(mesh, order=False)

