import meshio
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


# WIP
def get_area(points):
    """Returns the area of a polygon given the vertices.
    Parameters:
        points:     numpy.ndarray
            Vertices of the polygon
    Warning: the points need to be ordered clockwise or anticlockwise."""

    # shift all the points by one
    shifted = np.roll(points, 1, axis=0)

    # Use the shoelace formula
    area = 0.5 * np.sum((shifted[:, 0] + points[:, 0])*(shifted[:, 1] - points[:, 1]))

    return np.abs(area)


def get_dual_points(compliant_cells, index):
    """Returns the points of the dual mesh nearest to the point in the mesh given by the index.
    Parameters:
        mesh:       meshio.Mesh object
            Input mesh.
        index:      int
            Index of the point in the input mesh for which to calculate the nearest points of the dual mesh.
    """
    # Find the cells where the given index appears, REMEMBER the where statement gives you immediately back the indexof the compliant cell
    # building the compliant cells list
    _compliant = []
    for i in range(len(compliant_cells)):
    # compress with the use of any to determine the compliant cells
      if any(compliant_cells[i]==index):
        _compliant.append(i) 
    # Find the centers of all the cells
    return  _compliant

def compute_cell_centers(mesh,cells):
    """Returns the center points of the mesh elements, given the list of the different elements that are
   given as ordered vectors of points.
    Parameters:
        cells:       list of arrays
            list of arrays of point representing the cells. The cell Id is the position in the list.
    """
    x_centerpoint=0
    y_centerpoint=0
    #initialization of the vector of all the centerpoints of the cells
    centro=[]
    for elemento in cells:
        # cycle on the point of the elements
        for i in range(len(elemento)):
        # compute the centerpoint (accumulating)
         x_centerpoint += mesh.points[elemento[i],0]
         y_centerpoint += mesh.points[elemento[i],1]
         # define the center point
        x_centerpoint/=len(elemento)
        y_centerpoint/=len(elemento)
      # append to the centerpoints vector the element computed
        centro.append([x_centerpoint,y_centerpoint])
      # re-initialize the accumulation vectors
        x_centerpoint=0
        y_centerpoint=0
    return(centro) 

def get_dual(mesh,  order=False):
    """Returns the dual mesh held in a dictionary with dual["points"] giving the coordinates and
    dual["cells"] giving the indicies of all the cells of the dual mesh.
    Parameters:
        mesh:       meshio.Mesh object
            Input mesh.
        order:      boolean
            Whether to reorder the indices of each cell, such that they are in anticlockwise order.
    """
    
    assert isinstance(mesh, meshio.Mesh)
    # parce cells and boundaries
    #the alpha list is the list of cells (joined all the types, a todo should be to consider cells of different 
    #types
    alpha = []
    beta = []
    for i in range(len(mesh.cells)):
       for j in range(len(mesh.cells[i][1])):
        	alpha.append(mesh.cells[i][1][j])
    # to parse the boundaries we make use of the dictionary iterators
    for key, value in mesh.cell_sets.items():
        for key1,value1 in mesh.cell_sets[key].items():
             for i in range(len(value1)):
               beta.append(value1[i])

    print(alpha)
    print(beta)
    # create a dictionary for the graph
    graph ={}
    # Get the first set of points of the dual mesh
    compliant_cells = get_dual_points(alpha, 0)    
   # define a new key of the dictionary for each cell of the mesh
    graph={}
    for i in range(len(alpha)):
     graph.update({i :[]})
    # cycle on the points
    for idx in range(1, len(mesh.points)):
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        compliant_cells = get_dual_points(alpha, idx)
        print(compliant_cells)
        # in this part we build the graph: for each point of the mesh we have the compliant cells
        # and we cycle over the compliant cells (two nested loop, with an if that avoids to inspect the same cell)
        # me create the list inter that check the common point between two vectors (that can have also different 
        # dimension, considering that they can represent cells of completely different shape. 
        # checked that we have more than two vertex in common (WE ARE IN 2D HERE), and that the node is not already
        # connected with the analysed cell, we add it to the respective dictionary key.
        for i in compliant_cells:
          for j in compliant_cells:
             if i!=j:
               inter = list(set(alpha[i]).intersection(alpha[j]))
               if ((len(inter)>=2) and (j not in graph[i])):
                 graph[i].append(j)     
    print(graph)
    # Define the boundary cells in a similar way we did for the graph
    bnd = []       
    for i in range(len(alpha)):
    # We cycle on the boundary cells marked
        for j in beta:
           inter = list(set(alpha[i]).intersection(j))
           if ((len(inter)>=2) and (i not in bnd)):
             bnd.append(i) 
    print(bnd)        
    # definition of the graph from the dictionary as sugegsted by networkx
    g = nx.Graph(graph)
    # initialization of the figure
    plt.figure()
    # cycle on the elements
    for elemento in alpha:
      # cycle on the point of the elements
      for i in range(len(elemento)):
       # plot the grid
       x_value = [mesh.points[elemento[i-1],0],mesh.points[elemento[i],0]]
       y_value = [mesh.points[elemento[i-1],1],mesh.points[elemento[i],1]]
       plt.plot(x_value,y_value)  
    # compute the centerpoints
    cp = compute_cell_centers(mesh,alpha)
    # draw the adjacency graph
    nx.draw(g,pos=cp, with_labels=True)
    plt.savefig('foo.png')
