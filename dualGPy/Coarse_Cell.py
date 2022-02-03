import numpy as np
class Coarse_cell:
    """ Interface class to compute all the geometrical characteristics of the mesh """
    def __init__(self, fc_to_cc):
        self.fc_to_cc = fc_to_cc

