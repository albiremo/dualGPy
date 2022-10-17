"""This script takes a pseudo-3D mesh, meaning a 2D mesh extruded in one direction with just one layer of cells, and transforms it in a truly 2D mesh"""
from dualGPy.Utils import dualize_mesh
import meshio as mio

mshname = 'raebis'
in_ext = ".vtu"
# WARNING: some extensions, like vtu, do not fully support 2D, meaning that
# they will still use 3D arrays and, when the mesh is read with meshio, it will
# think it's dealing with 3D and not 2D.
# * Formats that support true-2D: dat (Tecplot), med, msh (Ansys, gmsh as well,
#   although more complex call if multi-type element mesh)
# * Formats than do NOT support true-2D: vtu/vtk
out_ext, out_format = (".dat", None) # (".med", None) # (".msh", "ansys") # (".vtu", None) #
# Extrusion dimension, the one to "delete"
del_dir = 1

msh_3D = mio.read(mshname + in_ext)

# dualize_mesh need a list of list of cells as input meshio.cells is a list of
# cellblocks. A cellblock has a cell_type (e.g. 'triangle') and data (a numpy
# matrix)
cells_3D = [c for cellblock in msh_3D.cells for c in cellblock.data.tolist()]

pts_2D, cells_2D = dualize_mesh(msh_3D.points.tolist(), cells_3D, del_dir)

# Although the following should be accepted, a warning is issued by numpy which
# is used by meshio under the hood, hence we change the call and use a
# dictionary msh_2D = mio.Mesh(points = pts_2D, cells = cells_2D)
#
cells_2D_d = { cell_type:nodes for cell_type,nodes in cells_2D }
msh_2D = mio.Mesh(points = pts_2D, cells = cells_2D_d)

out_name = '2D_' + mshname + out_ext
if out_format is not None:
    msh_2D.write(out_name, file_format = out_format)
else:
    msh_2D.write(out_name)

# Testing if output is really 2D
msh_test = mio.read(out_name)
print("Is the mesh \"truly\" 2D in output format?",
      msh_test.points.shape[-1] == 2)
