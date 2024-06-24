import itertools
import math
import numpy as np
import torch
from pyar.afir.restraints import get_covalent_radius

hartree2kcalmol = 627.509
bohr2angstroms = 0.529177
hartree2kjmol = 2625.5


class AFIRPotential:
    def __init__(self, **kwarg):
        self.config = kwarg
        self.hartree2kcalmol = hartree2kcalmol
        self.bohr2angstroms = bohr2angstroms
        self.hartree2kjmol = hartree2kjmol
        self.element_list = self.config["element_list"]
        self.num_atoms = len(self.element_list)
        self.num_fragm_1 = len(self.config["AFIR_Fragm_1"])
        self.num_fragm_2 = len(self.config["AFIR_Fragm_2"])
        self.num_fragm = self.num_fragm_1 + self.num_fragm_2
        self.num_fragm_1_list = self.config["AFIR_Fragm_1"]
        self.num_fragm_2_list = self.config["AFIR_Fragm_2"]
        self.gamma = self.config["AFIR_gamma"]
        self.alpha = 0.0
        self.A = 0.0
        self.B = 0.0
        self.p = 6.0
        self.R_0 = 3.8164/self.bohr2angstroms
        self.EPSIRON = 1.0061/self.hartree2kjmol
        self.omega = 0.0
        self.energy = 0.0
        return
    def calc_energy(self, geom_num_list):
        """
        # required variables: self.config["AFIR_gamma"], 
                             self.config["AFIR_Fragm_1"], 
                             self.config["AFIR_Fragm_2"],
                             self.config["element_list"]
        """
        """
        ###  Reference  ###
            Chem. Rec., 2016, 16, 2232
            J. Comput. Chem., 2018, 39, 233
            WIREs Comput. Mol. Sci., 2021, 11, e1538
        """
        R_0 = 3.8164/self.bohr2angstroms #ang.→bohr
        EPSIRON = 1.0061/self.hartree2kjmol #kj/mol→hartree
        if self.config["AFIR_gamma"] > 0.0 or self.config["AFIR_gamma"] < 0.0:
            alpha = (self.config["AFIR_gamma"]/self.hartree2kjmol) / ((2 ** (-1/6) - (1 + math.sqrt(1 + (abs(self.config["AFIR_gamma"]/self.hartree2kjmol) / EPSIRON))) ** (-1/6))*R_0) #hartree/Bohr
        else:
            alpha = 0.0
        A = 0.0
        B = 0.0
        
        p = 6.0

        for i, j in itertools.product(self.config["AFIR_Fragm_1"], self.config["AFIR_Fragm_2"]):
            R_i = get_covalent_radius(self.config["element_list"][i-1])
            R_j = get_covalent_radius(self.config["element_list"][j-1])
            vector = torch.linalg.norm(geom_num_list[i-1] - geom_num_list[j-1], ord=2) #bohr
            omega = ((R_i + R_j) / vector) ** p #no unit
            A += omega * vector
            B += omega
        
        energy = alpha*(A/B)#A/B:Bohr
        return energy #hartree


def maain():
    pass

if __name__ == "__main__":
    pass