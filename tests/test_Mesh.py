# content of test_class.py
import sys
import os
sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *
from dualGPy.Mesh import Mesh2D

class Test2D:
    def test_one(self):
        Mesh = Mesh2D(2,False)
        assert(Mesh.mesh.points[0][0]==0)

