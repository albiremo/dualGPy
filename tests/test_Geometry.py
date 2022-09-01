import sys
import os
sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *
from dualGPy.Geometry import Face3D,Hexa,Tetra,Wedge
import numpy as np
class Test3D:
    def test_square(self):
        points =  np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]])
        Face = Face3D(points)
        Face.ComputeArea()
        assert(Face.area==1)
    def test_triangle(self):
        points =  np.array([[0,0,0],[1,0,0],[0.5,1,0]])
        Face = Face3D(points)
        Face.ComputeArea()
        assert(Face.area==0.5)
    def test_cube(self):
        g_points =  np.array([[0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0]]
        )
        c_points = g_points
        faces = [[1,5,6,2],[4,5,6,7],[3,2,6,7]]
        cella = Hexa(c_points,g_points,faces)
        cella.ComputeVolume()
        assert(cella.volume==1)
    def test_tetra(self):
        c_points = np.array([[0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.5, 0.5, 0.5]])
        faces = [[0,1,2],[1,2,3],[0,2,3]]
        cella = Tetra(c_points,c_points,faces)
        cella.ComputeVolume()
        assert(cella.volume==0.08333333333333333) #1/3*A0*h
    def test_wedge(self):
        c_points = np.array([[0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0]])
        faces = [[0,1,2],[1,2,4,5],[0,2,3,4],[1,0,3,5],[3,4,5]]
        cella = Wedge(c_points,c_points,faces)
        cella.ComputeVolume()
        assert(cella.volume==0.5) #1/3*A0*h
        
