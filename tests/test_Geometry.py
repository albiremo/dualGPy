import sys
import os
sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *
from dualGPy.Geometry import Face3D
import numpy as np
class Test3D:
    def test_square(self):
        points =  np.array([0,0,0,1,0,0,1,1,0,0,1,0])
        Face = Face3D(points)
        Face.ComputeArea()
        assert(Face.area==1)
    def test_triangle(self):
        points =  np.array([0,0,0,1,0,0,0.5,1,0])
        Face = Face3D(points)
        Face.ComputeArea()
        assert(Face.area==0.5)

