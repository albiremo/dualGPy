# content of test_class.py
import sys
import os
sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *
from dualGPy.Mesh import Mesh2D, Mesh3D

Mesh =  Mesh3D(3,False)
Mesh.get_boundary_faces()
Mesh.ComputeGeometry()
