import CGNS.MAP as CGM
import CGNS.PAT.cgnsclass as CGC
import CGNS.PAT.cgnskeywords as CGK
import CGNS.PAT.cgnsutils as CGU
import CGNS.PAT.cgnslib as CGL

class CGNS_Adapter:
  """ Class that adapts the CGNS file to be processed
      with DualGPy"""
    def get_NGon_Node(self,n_zone):
    """ Method that retrives the NGon method """
       for n_elts in n_zone.nextChild(sidstype=CGK.Elements_ts):
        # NGON NODE
        if n_elts.data[0] == CGK.NGON_n:
            return n_elts
    def get_NFace_ElementRange_Value_From_Zone_Node(n_zone):
     nFaceNode = None
     for n_elts in n_zone.nextChild(sidstype=CGK.Elements_ts):
        # NGON NODE
        if n_elts.data[0] == CGK.NFACE_n:
            nFaceNode = n_elts.node
     if nFaceNode:
        elemRange = CGU.hasChildName(nFaceNode, CGK.ElementRange_s)
        return elemRange[1]

    def get_ngon_a_from_zone_node(n_zone):
     nGonNode = get_NGon_Node(n_zone).node
     if not nGonNode:
        raise ValueError("No NGon Node in the tree")
     elementConnectivityArray = CGU.hasChildName(nGonNode, CGK.ElementConnectivity_s)
     elemRange = CGU.hasChildName(nGonNode, CGK.ElementRange_s)
     return elemRange[1], elementConnectivityArray[1]

     def initialize_coordinates(n_zone):
        tpath_coord = [CGK.Zone_ts, CGK.GridCoordinates_ts]
        zpl_coord = CGU.getAllNodesByTypeList(n_zone, tpath_coord)
        # For tree t_initial, we get the nodes of types tpath
        for zp_coord in zpl_coord:
            parent_node = CGU.getNodeByPath(n_zone, zp_coord)
            x = CGU.hasChildName(parent_node, "CoordinateX")[1]
            y = CGU.hasChildName(parent_node, "CoordinateY")[1]
            z = CGU.hasChildName(parent_node, "CoordinateZ")[1]
        return(x,y,z) 
     def prepare_cgns_for_meshio(self,3D):
# We load the full CGNS tree
      (t_initial, t_ini_lk, t_ini_path) = CGM.load("ls59.cgns")
      n_tree = CGC.CGNSPython(t_initial)
      for n_base in n_tree.nextChild(sidstype=CGK.CGNSBase_ts):
        for n_zone in n_base.nextChild(sidstype=CGK.Zone_ts):
               i+=1
               print(i)
               x,y,z = initialize_coordinates(n_zone.node)
               print(len(x),len(y),len(z)) 
               ngon_node = get_NGon_Node(n_zone).node
               elem_range, element_connectivity_array = get_ngon_a_from_zone_node(n_zone)
               number_of_faces=[]
               number_of_faces.append(elem_range[1] - elem_range[0] + 1)  # Why a +1 is necessary?
               reshaped = np.reshape(element_connectivity_array,(number_of_faces[0],5))
               modified = []
               for actual_l in reshaped:
                modified.append(np.delete(actual_l,0))
               n_face_element_range = get_NFace_ElementRange_Value_From_Zone_Node(n_zone)
      temp = []
      coord = []
      for num in range(len(x)):
        temp.append(x[num])
        temp.append(y[num])
        if 3D:
         temp.append(z[num])
        coord.append(temp)
        temp=[]
      points = coord
      cells =[("quad",modified)] 
      return(points,cells) 

