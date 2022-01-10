import numpy as np
from numba import njit, prange

# http://www.katrin-affolter.ch/Python/performance_of_all_and_any
@njit
def any_non_zero_values(elements,index):

    for el in elements:

        if el == index:

            return True

    return False




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
      compliant = [i for i,e in enumerate(compliant_cells) if any_non_zero_values(e,index)]
 
      return compliant

