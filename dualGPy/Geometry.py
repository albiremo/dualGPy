import abc 
import numpy as np
from numba import njit, prange

class Face(abc.ABC):
    def __init__(self,points):
        self.points = points
    @property
    def area(self):
        return self._area

    @area.setter
    def area(self, new_area):
        if (new_area>0):
            self._area = new_area
        else:
            print("Error: the area of this solid is 0")
    @area.deleter
    def area(self):
        del self._area
    @property
    def n_vertices(self):
        return self._n_vertices
    @n_vertices.setter
    def n_vertices(self, new_n_vertices):
        if (new_n_vertices > 0):
           self._n_vertices = new_n_vertices
        else:
           print("Error: the n_segments of a face is always > 0")
    @n_vertices.deleter
    def n_vertices(self):
        del self._n_vertices
    @abc.abstractmethod
    def ComputeArea(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError

class Face2D(Face) :
    def __init__(self, segments,*args, **kwargs):
        self.segments = segments
        self.leng_segments = []
        super().__init__(*args, **kwargs)
        self.n_vertices = int(len(self.points)/2)


    def ComputeArea(self):
        cell_points_reshaped=np.reshape(self.points,(self.n_vertices,2))
        shifted = np.roll(cell_points_reshaped, 1, axis=0)
        area_0 = 0.5 * np.sum((shifted[:, 0] + cell_points_reshaped[:, 0])*(shifted[:, 1] - cell_points_reshaped[:, 1]))
        self.area = abs(area_0)
                  
    def ComputeLength(self,global_points):
        for segment in self.segments:
             leng = np.sqrt((global_points[segment[1]][1]-global_points[segment[0]][1])**2+(global_points[segment[1]][0]-global_points[segment[0]][0])**2)
             self.leng_segments.append(leng)




#class Solid(abc.ABC):
#    """ Interface class to compute the different characteristics of an element of
#        a 3D mesh """
#    def __init__(self,Faces):
#        self.Faces = Faces
#        self.volume = 0
#
#    @property
#    def volume(self):
#        return self._volume
#
#    @volume.setter
#    def volume(self, new_volume):
#	if new_volume > 0:
#	    self._volume = new_volume
#	else:
#	    print("Error: the volume of this solid is 0")
#    @volume.deleter
#    def volume(self):
#	del self._volume
#
