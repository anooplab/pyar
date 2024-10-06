

import logging
import os
import subprocess as subp

import numpy as np
import re
import glob

from pyar.interface import SF, write_xyz, which

orca_logger = logging.getLogger('pyar.orca')
Bohr2Angstrom =  0.52917721092

class Orca(SF):
    def __init__(self, molecule, qc_params):

        super(Orca, self).__init__(molecule)

        self.start_coords = molecule.coordinates
        self.inp_file = 'trial_' + self.job_name + '.inp'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.optimized_coordinates = []
        self.energy = 0.0
        # print(custom_keyword)
        keyword = f"! {qc_params['method']} {qc_params['basis']}"

        if any(x >= 21 for x in molecule.atomic_number):
            keyword += ' def2-ECP'
        keyword += '  def2/J'
        if self.scftype == 'uks':
            keyword += ' UKS'
        nprocs = qc_params['nprocs']
        # if custom_keyword is not None:
        #     keyword += custom_keyword
        keyword += f"\n%pal nprocs {nprocs} end\n"
        keyword += f"%scf maxiter {qc_params['scf_cycles']} end\n"
        self.keyword = keyword

    def prepare_input(self):
        keyword = self.keyword
        coords = self.start_coords
        with open(self.inp_file, "w") as f1:
            f1.write(keyword + "\n")
            f1.write("*xyz {0} {1}\n".format(str(self.charge), str(self.multiplicity)))
            for i in range(self.number_of_atoms):
                f1.write(
                    " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
            f1.write("*")

    def optimize(self):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        # TODO: Add a return 'CycleExceeded'

        # max_cycles = options['opt_cycles']  # noqa: F841
        # gamma = options['gamma']  # noqa: F841
        # convergence = options['opt_threshold']  # noqa: F841

        self.keyword = self.keyword + '!engrad'
        self.prepare_input()

        with open(self.out_file, 'w') as fopt:
            out = subp.Popen([which("orca"), self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            f = open(self.out_file, "r")
            line = f.readlines()
            if "****ORCA TERMINATED NORMALLY****" in line[-2]:
                self.energy = self.get_energy()
                self.optimized_coordinates = np.loadtxt(self.inp_file[:-4] + ".xyz", dtype=float, skiprows=2,
                                                        usecols=(1, 2, 3))
                write_xyz(self.atoms_list,
                          self.optimized_coordinates,
                          self.result_xyz_file, energy=self.energy)
                f.close()
                return True
            else:
                print("Error: OPTIMIZATION PROBABLY FAILED. "
                      "CHECK THE .out FILE FOR PARTIAL OPTIMIZTION ")
                print("Check for partial optimization.")
                print("Location: {}".format(os.getcwd()))
                return False

    def get_energy(self):
        """
        :return:This object will return energy from an orca calculation. It will return Hartree units.
        """
        try:
            with open(self.out_file, "r") as out:
                line = out.readlines()
                en_steps = [item for item in line if
                            "FINAL SINGLE POINT ENERGY" in item]
                if en_steps:
                    energy_in_hartrees = float((en_steps[-1].strip().split())[-1])
                else:
                    energy_in_hartrees = 0.0
            return energy_in_hartrees
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")
    


    def get_gradient(self, path):
        results = {}
        engrad_fn = glob.glob(os.path.join(path, "*.engrad"))
        if not engrad_fn:
            raise Exception("ORCA calculation failed.")

        assert len(engrad_fn) == 1
        engrad_fn = engrad_fn[0]
        with open(engrad_fn) as handle:
            engrad = handle.read()
        engrad = re.findall(r"([\d\-\.]+)", engrad)
        atoms = int(engrad.pop(0))
        energy = float(engrad.pop(0))
        force = -np.array(engrad[: 3 * atoms], dtype=float)
        results["energy"] = energy
        results["forces"] = force

        return results
    
    def parse_orca_output(self, molecule):
        natom = len(molecule.atoms)
        
        if self.calculate_energy:
            if not self.calculate_energy_gradients:
                with open(f'{self.inpfile}_property.txt', 'r') as orcaout:
                    orcaout_lines = orcaout.readlines()
                    for ii in range(len(orcaout_lines)):
                        if 'Total Energy' in orcaout_lines[ii]: # ? check SCf energy or total enenrgy
                            molecule.energy = float(orcaout_lines[ii].split()[-1]) 
            else:
                with open(f'{self.inpfile}.engrad', 'r') as orcaout:
                    orcaout_lines = orcaout.readlines()
                    for ii in range(len(orcaout_lines)):
                        if 'total energy' in orcaout_lines[ii]:
                            molecule.energy = float(orcaout_lines[ii+2])

        if self.calculate_energy_gradients:
            with open(f'{self.inpfile}.engrad', 'r') as orcaout:
                orcaout_lines = orcaout.readlines()
                for ii in range(len(orcaout_lines)):
                    if 'gradient' in orcaout_lines[ii]:
                        grad = orcaout_lines[ii+2: ii+2+natom*3]
                        grad = [float(g.split()[0])/Bohr2Angstrom for g in grad]
            for ii in range(natom):
                a = molecule.atoms[ii]
                a.energy_gradients = grad[3*ii: 3*ii+3] 
        
        if self.calculate_hessian:

            with open(f'{self.inpfile}.hess', 'r') as orcaout:
                orcaout_lines = orcaout.readlines()
                hessian_index = orcaout_lines.index('$hessian\n') # start from $hessian
                ncoordinate = int(orcaout_lines[hessian_index+1]) # second line is #coordinates
                ncircle = int((ncoordinate-0.5) / 5)+1 # blocks

                hessian_matrix = np.zeros((ncoordinate, ncoordinate))

                for ii in range(ncircle):
                    start_index = hessian_index+2+(ncoordinate+1)*ii
                    cols_index_list = orcaout_lines[start_index].split()
                    for jj in range(ncoordinate):
                        hessian_line = orcaout_lines[start_index+1+jj].split()
                        for kk in range(len(cols_index_list)):
                            col_index = cols_index_list[kk]
                            hessian_matrix[int(hessian_line[0])][int(col_index)] = float(hessian_line[kk+1])/Bohr2Angstrom**2
            molecule.hessian = hessian_matrix

        if not self.calculate_energy and self.output_keywords:
            for keyword in self.output_keywords:
                keyword_name = self.input_file+'_'+keyword.lower().replace(' ', '_')
                with open(f'{self.inpfile}_property.txt', 'r') as orcaout:
                    orcaout_lines = orcaout.readlines()
                    for ii in range(len(orcaout_lines)):
                        if keyword in orcaout_lines[ii]:
                            molecule.__dict__[keyword_name] = float(orcaout_lines[ii].split()[-1]) 


def main():
    pass


if __name__ == "__main__":
    main()
