from CoMMA import *
import meshio
from dualGPy.CGNS_Adapter import CGNS_Adapter
from dualGPy.Graph import Graph2D
from dualGPy.Mesh import Mesh2D
import pandas as pd

name = "rae.cgns"
adapter = CGNS_Adapter(name)
points,cells =  adapter.prepare_cgns_for_meshio(0)


mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    # Each item in cell data must match the cells array
)

mesh.write("input.vtk")
Mesh1 = Mesh2D(mesh)
print("setup done")
Mesh1.get_boundary_faces()
Dfaces = pd.DataFrame.from_dict(Mesh1.Dfaces, orient="index")
Dfaces.to_csv("Dfaces.csv")
faces = pd.DataFrame.from_dict(Mesh1.faces, orient="index")
faces.tp_csv("faces.csv")
print("area computed")
Mesh1.ComputeVolume()
print("volume comuted")
Mesh1.ComputeArea()
rint("bf computed")
Mesh1.boundary_detection()
print("boundary detection performed")
# graph
Graph1 = Graph2D(Mesh1)
Graph1.get_adj_matrix()
print("adjacency matrix written")
adjacency = np.asarrray(Graph1.adj)
np.savetxt('adj.csv', adjacency, delimiter=',')
Graph1.adj_to_csr()
edges = np.asarrray(Graph1.edges)
vertex = np.asarrray(Graph1.vertex)
np.savetxt('edges.csv', edges, delimiter=',')
np.savetxt('vertex.csv', vertex, delimiter=',')
print("shift to CSR done")

nb_fc = len(Graph1.vertex)-1
adjMatrix_row_ptr= np.array(Graph1.vertex , dtype='long')
adjMatrix_col_ind= np.array(Graph1.edges ,dtype='long')
adjMatrix_areaValues=np.array(Mesh1.area,dtype='double')
volumes = np.array(Mesh1.volume,dtype='double')
isOnBnd = np.array(Mesh1.boundary_cells,dtype='long')
array_isOnRidge=np.array(Mesh1.onRidge,dtype='long')
array_isOnValley=np.array(Mesh1.onValley, dtype='long')
array_isOnCorner=np.array(Mesh1.onCorner, dtype='long')
fc_to_cc = np.full(nb_fc, -1,dtype='long')
indCoarseCell = 0
minCard = -1
goalCard = -1
maxCard = -1
verbose = 0
arrayOfFineAnisotropicCompliantCells = np.arange(nb_fc,dtype='long')
agglomerationLines_Idx = np.zeros(nb_fc,dtype='long')
agglomerationLines = np.zeros(nb_fc,dtype='long')
isFirstAgglomeration = 1
isAnisotropic = 1
dimension = 2
is_basic_or_triconnected = 1

print("CoMMA call")
fc_to_cc_res,agglomerationLines_Idx_res,agglomerationLines_res=agglomerate_one_level(adjMatrix_row_ptr, adjMatrix_col_ind, adjMatrix_areaValues, volumes,arrayOfFineAnisotropicCompliantCells,isOnBnd,array_isOnValley,array_isOnRidge,array_isOnCorner,isFirstAgglomeration,isAnisotropic,fc_to_cc,agglomerationLines_Idx,agglomerationLines,is_basic_or_triconnected,dimension,goalCard,minCard,maxCard,verbose)

print("finalizing")
fine_cells = []

for j in range(len(fc_to_cc_res)):
        fine_cells.append(np.float64(fc_to_cc_res[j]))
meshOUT = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    # Each item in cell data must match the cells array
    cell_data = {"agglomerate":[fine_cells] },
)

print("writing")
meshOUT.write(".vtk")

