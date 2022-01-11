import meshio
from dualGPy.CGNS_Adapter import CGNS_Adapter
from dualGPy.Graph import Graph2D
from dualGPy.Mesh import Mesh2D
import pandas as pd

name = "rae.cgns"
adapter = CGNS_Adapter(name)
points,cells =  adapter.prepare_cgns_for_meshio(0)


mesh = meshio.Mesh(
    points,
    cells,
)

mesh.write("input.vtk")

