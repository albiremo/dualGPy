import numpy as np
from numba import njit, prange



def get_dual_points(compliant_cells : int, index : int) -> int:
      """ Function that returns the points of the dual mesh nearest to the point in the mesh given by the index.
          Parameters:
          mesh:       meshio.Mesh object
              Input mesh.
          index:      int
              Index of the point in the input mesh for which to calculate the nearest points of the dual mesh.
          Return:
              compliant 
        # Find the cells where the given index appears, REMEMBER the where statement gives you immediately back the indexof the compliant cell
        # building the compliant cells list
      """
      compliant = [e for i,e in enumerate(compliant_cells) if e is index]
      return compliant


@njit(parallel=True)
def parallel_nonzero_count(arr):
    flattened = arr.ravel()
    sum_ = 0
    for i in prange(flattened.size):
        sum_ += flattened[i] != 0
    return sum_
