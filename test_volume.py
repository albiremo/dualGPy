from dualGPy.Utils import *
from dualGPy.Mesh import Mesh2D, Mesh3D
import numpy as np

Mesh =  Mesh3D(3,False)
Mesh.get_boundary_faces()
Mesh.ComputeGeometry()

