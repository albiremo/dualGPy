import meshio
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def array_intersection(a, b):
    """Returns a boolean array of where b array's elements appear in the a array"""
    # Source: https://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
    if not isinstance(a, np.ndarray):
        a = np.array(a)
    if not isinstance(b, np.ndarray):
        b = np.array(b)
    tmp = np.prod(np.swapaxes(a[:, :, None], 1, 2) == b, axis=2)
    return np.sum(np.cumsum(tmp, axis=0) * tmp == 1, axis=1).astype(bool)



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


def get_dual_points(mesh, compliant_cells, index):
    """Returns the points of the dual mesh nearest to the point in the mesh given by the index.
    Parameters:
        mesh:       meshio.Mesh object
            Input mesh.
        index:      int
            Index of the point in the input mesh for which to calculate the nearest points of the dual mesh.
    """
    ## For each type of cell do the following
    # Find the cells where the given index appears, REMEMBER the where statement gives you immediately back the indexof the compliant cell
    _idxs = [np.where(x[1] == index)[0] for x in mesh.cells]
# building the compliant cells list
    _compliant = []
    for i in range(len(compliant_cells)):
# compress with the use of any to determine the compliant cells
      if any(compliant_cells[i]==index):
        _compliant.append(i) 
#    	_compliant. = [np.where(x == index)[0] for x,i in enumerate(compliant_cells)]
    # Find the centers of all the cells
    _vs = [mesh.points[x[1][_idxs[i]]].mean(axis=1) for i,x in enumerate(mesh.cells)]
    return np.concatenate(_vs, axis=0), _compliant

#def get_neigh(compliant_cells, analyded_node, id2D):
#    """Returns the points of the dual mesh nearest to the point in the mesh given by the index.
#    Parameters:
#        mesh:       meshio.Mesh object
#            Input mesh.
#        index:      int
#            Index of the point in the input mesh for which to calculate the nearest points of the dual mesh.
#    """
#    assert isinstance(mesh, meshio.Mesh)
#    ## For each type of cell do the following
#    # Find the cells where the given index appears
#    # Find the centers of all the cells    
#    for x in compliant_cells:
#	for i in range(len(x)):
#		for j in range(len(analysed_node))
#    			if (len(np.where(mesh.cells[1][i] ==  analysed_node[1][j])[0])>=2):
#                           neigh = np.concatenate(x, axis=0)
#    return neigh

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
    new_points, compliant_cells = get_dual_points(mesh, alpha, 0)    
    # Create the containers for the points and the polygons of the dual mesh
    dual_points = new_points
    # define a new key of the dictionary for each cell of the mesh
    graph={}
    for i in range(len(alpha)):
     graph.update({i :[]})
    # cycle on the points
    for idx in range(1, len(mesh.points)):
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        new_points, compliant_cells = get_dual_points(mesh, alpha, idx)
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
        for j in beta:
           inter = list(set(alpha[i]).intersection(j))
           if ((len(inter)>=2) and (i not in bnd)):
             bnd.append(i) 
    print(bnd)        
    g = nx.Graph(graph)
    plt.figure()
    for elemento in alpha:
      for i in range(len(elemento)):
       x_value = [mesh.points[elemento[i-1],0],mesh.points[elemento[i],0]]
       y_value = [mesh.points[elemento[i-1],1],mesh.points[elemento[i],1]]
       plt.plot(x_value,y_value)  
    nx.draw(g, with_labels=True)
    plt.savefig('foo.png')
