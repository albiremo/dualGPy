import sys
import os
sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *
from dualGPy.Mesh import Mesh2D, Mesh3D
import numpy as np
from meshio import Mesh

class Test2D:
    def test_one(self):
        Mesh = Mesh2D(2,False)
        assert(Mesh.mesh.points[0][0]==0)

    def test_multi_zone(self):
        """Testing on this mesh composed of 2 rectangles and two triangles (| and _ = 0.5*L)
    6
   /|\
  / | \
  | | |
5/__|__\3
 |  |4 |
 |  |  |
 |  |  |
0|__|__|2
    1
"""
        L = 1
        points = L * np.array([
            [0, 0], # idx = 0
            [1, 0], # idx = 1
            [2, 0], # idx = 2
            [2, 2], # idx = 3
            [1, 2], # idx = 4
            [0, 2], # idx = 5
            [1, 4]  # idx = 6
            ])
        cells = {
            'quad':     [[0,1,4,5], [1,2,3,4]],
            'triangle': [[4,5,6],   [3,4,6]]
            }
        base_mesh = Mesh(points = points, cells = cells)
        mesh = Mesh2D(base_mesh)
        mesh.get_boundary_faces()
        mesh.ComputeGeometry()
        ref_vol = L*L * np.array([2,2,1,1], dtype = float)
        assert(np.all(mesh.volume == ref_vol))

class Test3D:
    def test_Volume_Hexa(self):
        Mesh =  Mesh3D(3,False)
        Mesh.get_boundary_faces()
        Mesh.ComputeGeometry()
        assert(np.all(Mesh.area))
        assert(np.all(Mesh.volume))

    def test_multi_zone(self):
        """A cube of 2*L edge divided into 2 on dimension z, mounted by a pyramid dived into 2 tetrahedra"""
        L = 1
        points = L * np.array([
            [0, 0, 0], # idx = 0
            [2, 0, 0], # idx = 1
            [2, 2, 0], # idx = 2
            [0, 2, 0], # idx = 3
            [0, 0, 1], # idx = 4
            [2, 0, 1], # idx = 5
            [2, 2, 1], # idx = 6
            [0, 2, 1], # idx = 7
            [0, 0, 2], # idx = 8
            [2, 0, 2], # idx = 9
            [2, 2, 2], # idx = 10
            [0, 2, 2], # idx = 11
            [1, 1, 5]  # idx = 12
        ])
        cells = {
        'hexahedron': [[0,1,2,3,4,5,6,7], [4,5,6,7,8,9,10,11]],
        'tetra'     : [[8,9,11,12], [9,10,11,12]]
        }
        base_mesh = Mesh(points = points, cells = cells)
        mesh = Mesh3D(base_mesh)
        mesh.get_boundary_faces()
        # Since the algorithm store only internal faces, we should manually add some faces for the parallelograms
        mesh.faces[0].extend([[0,1,5,4],[1,2,6,5]])
        mesh.faces[1].extend([[8,9,5,4],[4,7,11,8]])
        mesh.ComputeGeometry()
        ref_vol = L*L*L * np.array([4,4,2,2], dtype = float)
        assert(np.max(np.abs(mesh.volume - ref_vol)) < 1e-14)
