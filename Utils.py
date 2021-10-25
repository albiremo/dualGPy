def get_dual_points(compliant_cells : int, index : int) -> int:
      """ Static method that returns the points of the dual mesh nearest to the point in the mesh given by the index.
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
      compliant=[]
      for i in range(len(compliant_cells)):
      # compress with the use of any to determine the compliant cells
         if any(compliant_cells[i]==index):
           compliant.append(i) 
      # Find the centers of all the cells
      return compliant
 