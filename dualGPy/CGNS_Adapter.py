import CGNS.MAP as CGM
import CGNS.PAT.cgnsclass as CGC
import CGNS.PAT.cgnskeywords as CGK
import CGNS.PAT.cgnsutils as CGU
import CGNS.PAT.cgnslib as CGL
import numpy as np

class CGNS_Adapter:
 """ Class that adapts the CGNS file to 
 be processed with DualGPy
 """ 
 def __init__(self, name):
  self.name = name
 
 @staticmethod
 def get_NGon_Node(n_zone):
  """ Method that retrives the NGon method """
  # Search the element types in the childs of the zone
  # in order to retrive the node relative to the ngon elements
  for n_elts in n_zone.nextChild(sidstype=CGK.Elements_ts):
      # NGON NODE
   if n_elts.data[0] == CGK.NGON_n:
    return n_elts


 def get_ngon_a_from_zone_node(self,n_zone):
  """ Retrive the element range and the connectivity ARRAY of the current
        zone. remember that it is composed as follows:
        elementConnectivityArray = [num_el_face0, coord[0].... coord[num_el], num_el_face1, coord[0].....coord[num_el]]"""
  nGonNode = self.get_NGon_Node(n_zone).node
  if not nGonNode:
        raise ValueError("No NGon Node in the tree")
  elementConnectivityArray = CGU.hasChildName(nGonNode, CGK.ElementConnectivity_s)
  elemRange = CGU.hasChildName(nGonNode, CGK.ElementRange_s)
  return elemRange[1], elementConnectivityArray[1]


 @staticmethod
 def initialize_coordinates(n_zone):
  """ Retrives the coordinates of the points composing the mesh """
  tpath_coord = [CGK.Zone_ts, CGK.GridCoordinates_ts]
  zpl_coord = CGU.getAllNodesByTypeList(n_zone, tpath_coord)
  # For tree t_initial, we get the nodes of types tpath
  for zp_coord in zpl_coord:
            parent_node = CGU.getNodeByPath(n_zone, zp_coord)
            x = CGU.hasChildName(parent_node, "CoordinateX")[1]
            y = CGU.hasChildName(parent_node, "CoordinateY")[1]
            z = CGU.hasChildName(parent_node, "CoordinateZ")[1]
  return(x,y,z) 

 def prepare_cgns_for_meshio(self,trid):
  """ extraxt the faces and the points for the meshio elaboration. The flag 3D defines
         if we treat with a 3D or 2D mesh""" 
  (t_initial, t_ini_lk, t_ini_path) = CGM.load(self.name)
  n_tree = CGC.CGNSPython(t_initial)
  for n_base in n_tree.nextChild(sidstype=CGK.CGNSBase_ts):
        for n_zone in n_base.nextChild(sidstype=CGK.Zone_ts):
               x,y,z = self.initialize_coordinates(n_zone.node)
               ngon_node = self.get_NGon_Node(n_zone).node
               elem_range, element_connectivity_array = self.get_ngon_a_from_zone_node(n_zone)
               number_of_faces=[]
               number_of_faces.append(elem_range[1] - elem_range[0] + 1)  # Why a +1 is necessary?
               reshaped = np.reshape(element_connectivity_array,(number_of_faces[0],5))
               modified = []
               for actual_l in reshaped:
                deleted = np.delete(actual_l,0)
                minus_one = [number - 1 for number in deleted]
                modified.append(minus_one)
               print("ended")
  temp = []
  coord = []
  for num in range(len(x)):
        temp.append(x[num])
       # temp.append(y[num])
        temp.append(z[num])
        coord.append(temp)
        temp=[]
  points = coord
  cells =[("quad",modified)] 
  return(points,cells) 

