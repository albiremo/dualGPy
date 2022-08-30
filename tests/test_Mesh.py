import sys
import os
sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *
from dualGPy.Mesh import Mesh2D, Mesh3D
import numpy as np

class Test2D:
    def test_one(self):
        Mesh = Mesh2D(2,False)
        assert(Mesh.mesh.points[0][0]==0)
class Test3D:
    def test_Volume_Hexa(self):
        Mesh =  Mesh3D(3,False)
        Mesh.get_boundary_faces()
        Mesh.ComputeGeometry()
        assert(np.all(Mesh.area))
        assert(np.all(Mesh.volume))

