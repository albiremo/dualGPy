import sys
import os
import numpy as np
import meshio

sys.path.append(os.path.abspath('../'))
from dualGPy.Utils import *  # noqa: E402
from dualGPy.Mesh import Mesh2D, Mesh3D  # noqa: E402

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
        base_mesh = meshio.Mesh(points = points, cells = cells)
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

    def test_Volume_non_rectangle_Hexa(self):
        Lx, Ly, Lz = 2, 2, 2
        dx, dy = .5*Lx, .25*Ly
        points = np.array([
            [0     , 0     , 0] , # idx = 0
            [Lx    , dy    , 0] , # idx = 1
            [Lx    , Ly+dy , 0] , # idx = 2
            [0     , Ly    , 0] , # idx = 3
            [dx    , 0     , Lz], # idx = 4
            [Lx+dx , dy    , Lz], # idx = 5
            [Lx+dx , Ly+dy , Lz], # idx = 6
            [dx    , Ly    , Lz], # idx = 7
            ], dtype = float)
        cells = [('hexahedron', [[0,1,2,3,4,5,6,7]])]
        base_mesh = meshio.Mesh(points = points, cells = cells)
        mesh = Mesh3D(base_mesh)
        # Since we are working with only one cell, we have to give it the faces ourselves
        mesh.faces[0] = [ [0,1,2,3], [0,1,5,4], [0,3,4,7] ]
        mesh.ComputeGeometry()
        ref_vol = Lx*Ly*Lz
        assert(abs(mesh.volume[0] - ref_vol) < 1e-14)

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
        base_mesh = meshio.Mesh(points = points, cells = cells)
        mesh = Mesh3D(base_mesh)
        mesh.get_boundary_faces()
        # Since the algorithm store only internal faces, we should manually add some faces for the parallelograms
        mesh.faces[0].extend([[0,1,5,4],[1,2,6,5]])
        mesh.faces[1].extend([[8,9,5,4],[4,7,11,8]])
        mesh.ComputeGeometry()
        ref_vol = L*L*L * np.array([4,4,2,2], dtype = float)
        assert(np.max(np.abs(mesh.volume - ref_vol)) < 1e-14)
