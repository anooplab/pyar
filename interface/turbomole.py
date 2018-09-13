"""
turbomole.py - interface to turbomole program

Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""
import os
import shutil
import subprocess as subp
import sys
import glob
import datetime
import numpy as np
import socket
import time
import re
import interface
import interface.babel
from afir import restraints
from interface import SF
from interface import which
from units import angstrom2bohr, bohr2angstrom

import logging

turbomole_logger = logging.getLogger('pyar.turbomole')


class Turbomole(SF):
    def __init__(self, molecule, method):
        if which('define') is None:
            print('set Turbomole path')
            turbomole_logger.error('set Turbomole path')
            sys.exit()

        super(Turbomole, self).__init__(molecule)

        self.charge = method['charge']
        self.scf_type = method['scftype']
        self.multiplicity = method['multiplicity']
        self.start_coords = angstrom2bohr(molecule.coordinates)
        self.atoms_in_fragments = molecule.fragments
        self.energy = None
        self.optimized_coordinates = None

    def optimize(self, max_cycles=350, gamma=0.0):
        """This is the python implementation  of jobex of turbomole"""

        make_coord(self.atoms_list, self.start_coords)
        prepare_control()

        if gamma == 0.0:
            return self.run_turbomole_jobex()

        # Test for Turbomole input
        if not os.path.isfile('control'):
            turbomole_logger.critical("Turbomole control file should exist")
            return False

        remove('converged')
        with open('not.converged', 'w') as fp:
            fp.write('$convcrit\n')
            fp.write('   energy              6\n')
            fp.write('   cartesian gradient  3\n')
            fp.write('   basis set gradient  3\n')
            fp.write('   cycles              1\n')
            fp.write('$convergence not reached\n')
            fp.write('$end\n')

        remove('abend.*')
        remove('statistics')
        remove('dscf_problem')
        remove('job.[0-9]*')
        remove('job.last')
        remove('GEO_OPT_CONVERGED')
        remove('GEO_OPT_FAILED')
        remove('GEO_OPT_RUNNING')
        with open('GEO_OPT_RUNNING', 'w') as fs:
            fs.write("Job started at  %s\n" % datetime.datetime.now())
            fs.write("  on machine %s\n"% socket.gethostname())
            fs.write("  by user %s\n" % os.getlogin())
            fs.write("  in directory %s\n" % os.getcwd())
            fs.write(" CAUTION: this file will not be removed if user or the system decides \n")
            fs.write("          to kill the job !!!\n")
            fs.write("          To let jobex stop smoothly, simple create\n")
            fs.write("          an empty file called \'stop\' in this directory\n")
            fs.write("          jobex will stop at the next possible step\n")

        if os.path.isfile('nextstep'):
            with open('nextstep') as ff:
                first_step = ff.read().strip()
        else:
            first_step = 'dscf'
        turbomole_logger.debug('the first step is %s' % first_step)

        initial_status = False
        initial_message = ""
        initial_energy = None
        if first_step == 'statpt':
            initial_status, initial_message, initial_energy = self.calc_energy
        elif first_step == 'relax':
            update_coord()
            initial_status, initial_message, initial_energy = self.calc_energy
        elif first_step == 'grad':
            initial_status, initial_message, initial_energy = get_energy()
        elif first_step == 'dscf':
            initial_status, initial_message, initial_energy = self.calc_energy
        else:
            turbomole_logger.critical('first step is not defined')
        turbomole_logger.debug('First step: %s, %s, %f' %(initial_status, initial_message, initial_energy))

        for cycle in range(max_cycles):
            # Calculate Gradients
            status, message, _ = self.calc_gradients
            if status is False:
                turbomole_logger.critical('Gradient evaluation failed')
                turbomole_logger.critical('Check the job in %s' % os.getcwd())
                turbomole_logger.error(message)
                return False

            grad_line = ""
            for line in open('gradient').readlines():
                if 'cycle' in line:
                    grad_line = line.strip()
            turbomole_logger.debug(grad_line)

            # Calculate afir gradient if gamma is greater than zero
            if gamma > 0.0:
                restraint_energy, trg, rg = restraints.isotropic(self.atoms_in_fragments, self.atoms_list, get_coords(), gamma)
                rewrite_turbomole_energy_and_gradient_files(self.number_of_atoms, rg, restraint_energy, trg)
                turbomole_logger.debug('restraint energy = %f\ntotal restraint gradients %f' % (restraint_energy, trg))

            # Update coordinates and check convergence.
            status = update_coord()
            if status is True:
                turbomole_logger.info('converged at %d' % cycle)
                self.energy = get_energy()
                self.optimized_coordinates = bohr2angstrom(get_coords())
                interface.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                                    self.job_name,
                                    energy=self.energy)
                return True

            # Calculate energy
            status, message, energy = self.calc_energy
            if status is False:
                turbomole_logger.critical('Energy evaluation failed')
                turbomole_logger.critical(message)
                turbomole_logger.error('check %s' %os.getcwd())
                return False

            with open('energy.dat','a') as fe:
                if gamma > 0.0:
                    fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy+restraint_energy))
                else:
                    fe.writelines("{:3d} {:15.8f}\n".format(cycle, energy))
        else:
            turbomole_logger.error("OPTIMIZATION DID NOT CONVERGE WITHIN $cycles CYCLES\n")
            turbomole_logger.error("  Restarting it after checking the gradient "
                                   "norms might be a good idea... "
                                   "(grep cycle gradient)")
            remove('GEO_OPT_RUNNING')

            with open('GEO_OPT_FAILED') as fe:
                fe.write("OPTIMIZATION DID NOT CONVERGE WITHIN $cycles CYCLES\n"
                         "   Restarting it after checking the gradient norms\n "
                         "   might be a good idea... (grep cycle gradient)")

            status = 'cycle_exceeded'
            self.energy = get_energy()
            self.optimized_coordinates = bohr2angstrom(get_coords())
            return status

    def run_turbomole_jobex(self):
        """
        return one of the following:
            True: converged
            SCFFailed
            GradFailed
            UpdateFailed
            CycleExceeded
            False: unknown error or jebex excution error
        :rtype: string or boolean
        """
        with open('jobex.out', 'w') as fj:
            try:
                subp.check_call(['jobex', '-ri', '-c', '350'], stdout=fj,
                            stderr=fj)
            except subp.CalledProcessError as e:
                turbomole_logger.error('jobex failed, check %s/jobex.out'
                                       % os.getcwd())
                return False

        if os.path.isfile('GEO_OPT_FAILED'):
            message = open('GEO_OPT_FAILED').read()
            if 'ERROR: Module' in message:
                turbomole_logger.error('Error in module!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return False
            elif 'ERROR in statpt step,' in message:
                turbomole_logger.error('Statpt failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'UpdateFailed'
            elif 'ERROR in relax step,' in message:
                turbomole_logger.error('Relax failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'UpdateFailed'
            elif 'ERROR: your energy calculation did not converge !!,' in message:
                turbomole_logger.error('SCF failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'SCFFailed'
            elif 'ERROR in dscf step' in message:
                turbomole_logger.error('SCF failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'SCFFailed'
            elif 'ERROR in gradient step' in message:
                turbomole_logger.error('Gradient failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'GradientFailed'
            elif 'OPTIMIZATION DID NOT CONVERGE' in message:
                turbomole_logger.error('Geometry did not converge!\n check %s' % os.getcwd())
                return 'CycleExceeded'
            else:
                turbomole_logger.error('Unknown Error!\n chcek the files in %s' % os.getcwd())
                return False

        elif os.path.isfile('GEO_OPT_CONVERGED'):

            turbomole_logger.info('jobex done')
            self.energy = get_energy()
            self.optimized_coordinates = bohr2angstrom(get_coords())
            interface.write_xyz(self.atoms_list, self.optimized_coordinates,
                                self.result_xyz_file,
                                self.job_name,
                                energy=self.energy)
            return True

    @property
    def calc_energy(self):

        status = run_turbomole_module('ridft')
        msg = [line for line in open('ridft.out').readlines() if 'ended' in line]
        if status is False or 'abnormally' in msg:
            return False, msg, None
        elif os.path.isfile('dscf_problem'):
            msg = "SCF Failure. Check files in"+os.getcwd()
            remove('GEO_OPT_RUNNING')
            return False, msg, None
        else:
            with open('job.last', 'a') as fp:
                fp.write(open('ridft.out').read())
            return True, msg, get_energy()

    @property
    def calc_gradients(self):

        run_turbomole_module('rdgrad')
        msg = [line for line in open('rdgrad.out').readlines() if 'ended' in line]
        if 'abnormally' in msg:
            return False, msg, None
        else:
            gradients = get_gradients(self.number_of_atoms)
            with open('job.last','a') as fp:
                fp.write(open('rdgrad.out').read())
            return True, msg, gradients


def abend(energy_program='ridft', gradient_program='rdgrad', relax_program='statpt'):
    msg = subp.check_output(['actual', '-e', energy_program, "-g", gradient_program, "-c", relax_program])
    turbomole_logger.info(msg)
    for failed in glob.glob('abend.*'):
        failed_program = failed.split('.')[1]
        turbomole_logger.critical('Error in '+failed_program)
        remove(failed)
        remove('GEO_OPT_RUNNING')
        err_msg = "ERROR: Module " + failed_program + " failed to run " \
                  "properly - please check output files for the reason"
        turbomole_logger.critical(err_msg)
        with open('GEO_OPT_FAILED', 'w') as fp:
            fp.write(err_msg)
        return False


def remove(file_name):
    files = glob.glob(file_name)
    for file in files:
        if os.path.exists(file_name):
            os.remove(file_name)
    return


def statpt_cart():
    statpt_status = run_turbomole_module('statpt')
    soft_abend('statpt')
    return statpt_status


def make_coord(atoms_list, coordinates, outfile='coord'):
    """
    """

    with open(outfile, 'w') as fcoord:
        fcoord.write('$coord\n')
        for s, c in zip(atoms_list, coordinates):
            fcoord.write("{:20.13} {:20.13} {:20.13} {:3s}\n".format(c[0], c[1], c[2], s))
        fcoord.write('$user-defined bonds\n')
        fcoord.write('$end')
    return


def prepare_control(basis="def2-SVP", func="b-p", ri="on",
                    memory=6000, grid="m4", scf_conv=7, scf_maxiter=500,
                    dftd3="yes", charge=0, multiplicity=1,
                    read_define_input="no", coordinates='internal'):

    """This function will prepare the control file. Most of the arguments have
       default values(all have in its current form). The read_define_input var
       if "yes", then define_input_file should provide. The define will run
       from this input file."""

    # Turbomole 'coord' file should exist
    if not os.path.isfile('coord'):
        turbomole_logger.critical("Turbomole coordinate file, 'coord', should exist")
        return False

    define_input_file = "define.inp"
    if read_define_input == "no":
        with open(define_input_file, "w") as fdefine:
            fdefine.write("\n\n")
            fdefine.write("a coord\n")
            if coordinates is 'internal':
                fdefine.write("ired\n*\n")
            elif coordinates is 'cartesian':
                fdefine.write("*\nno\n")
            fdefine.write("bb all %s\n*\n" % basis)
            fdefine.write("eht\n")
            fdefine.write("y\n")
            fdefine.write("%d\n" % charge)
            if multiplicity > 2:
                fdefine.write("n\n")
                fdefine.write("u %d\n\n" % multiplicity)
            else:
                fdefine.write("y\n\n")
            fdefine.write("\n\n\n")
            fdefine.write("dft\n")
            fdefine.write("on\n")
            fdefine.write("func %s\n" % func)
            fdefine.write("grid %s\n" % grid)
            fdefine.write("\n")
            if ri == "on":
                fdefine.write("ri\n on\n m\n")
                fdefine.write("%d\n" % memory)
                fdefine.write("\n")
            if dftd3=="yes":
                fdefine.write("dsp\n on\n")
                fdefine.write("\n")
            fdefine.write("scf\nconv\n")
            fdefine.write("%d\n" % scf_conv)
            fdefine.write("iter\n")
            fdefine.write("%d\n" % scf_maxiter)
            fdefine.write("\n*")

        define_log_file = "define.log"

        with open(define_log_file, 'w') as fout, open('define.inp') as fin:
            try:
                subp.check_call(['define'], stdin=fin, stdout=fout,
                                stderr=fout, universal_newlines=False)
            except subp.CalledProcessError as e:
                turbomole_logger.critical('Error in define')
                turbomole_logger.critical(e.output, e.returncode)
                return False
    return True


def soft_abend(program):
    with open('job.last', 'a') as fout:
        subp.check_call(['actual', '-e', 'ridft', '-g', 'rdgrad', '-c', program],
                        stdout=fout, stderr=fout)


def add_tmole_key():
    control = open('control').read()

    if '$tmole' not in control:
        print('no tmole')
        new_control = re.sub('end', 'tmole\n$end', control)

        with open('control', 'w') as fc:
            fc.write(new_control)


def get_metric_value():
    control = open('control').read()
    for line in control.split('\n'):
        if 'metric' in line:
            return int(line.split()[1])
    else:
        return 3


def change_or_add_metfic_in_redund_inp(metric):
    control = open('control').read()

    for i in control.split('$'):
        if 'redund_inp' in i:
            new_control = re.sub('metric '+str(metric), 'metric '+str(metric-1), control)
            break
    else:
        new_control = re.sub('end', 'redund_inp\n    metric '+str(metric-1)+'\n$end', control)

    with open('control', 'w') as fc:
        fc.write(new_control)


def generate_internal_coordinates():

    if not os.path.isfile('coord'):
        turbomole_logger.error("Turbomole coordinate file, 'coord', should exist")
        return False

    define_input_file = "def_inp"
    define_log_file = "define.log"

    add_tmole_key()
    current_metric = get_metric_value()

    for metric in range(current_metric, -3, -1):
        if metric == 0:
            continue
        change_or_add_metfic_in_redund_inp(metric)

        with open(define_log_file, 'w') as fout, open(define_input_file) as fin:
            try:
                subp.check_call(['define'], stdin=fin, stdout=fout,
                                stderr=fout, universal_newlines=False)
            except subp.CalledProcessError as e:
                pass

        if 'ATOMIC ATTRIBUTE DATA' in open(define_log_file).read():
            if 'intdef' in open('control').read():
                subp.check_call(['kdg', 'cartesianstep'])
                subp.check_call(['kdg', 'redund_inp'])
                remove('hessapprox')
                break
            else:
                subp.check_call(['kdg', 'redund_inp'])
                break
    else:
        subp.check_call(['kdg', 'redund_inp'])
        turbomole_logger.error('error in define step while generating new internal coordiantes')

    remove('tmp.input')
    return True

# TODO: insert redund_inp and metric

def get_energy():
    with open('energy') as fp:
        return float(fp.readlines()[-2].split()[1])


def get_gradients(natoms):
    start =natoms+1
    grad = [line.replace('D', 'E').split() for line in open('gradient').readlines()[-start:-1]]
    return grad


def get_coords():
    with open('coord') as fp:
        f = fp.readlines()
        if f[0].strip() != '$coord':
            turbomole_logger.critical("'$coord' not found. Check the file.")
            sys.exit()
        geometry_section = []
        for line in f[1:]:
            if line.strip() == '$user-defined bonds':
                break
            geometry_section.append(line.split())
        coordinates = []
        for i, c in enumerate(geometry_section):
            try:
                x_coord = float(c[0])
                y_coord = float(c[1])
                z_coord = float(c[2])
                symbol = c[3].capitalize()
            except:
                turbomole_logger.error("Something wrong in line: ", i + 1)
                turbomole_logger.error(c)
                sys.exit()
            coordinates.append([x_coord, y_coord, z_coord])

        return np.array(coordinates)


def update_coord():
    """ This function is kept, in case, we want to implement
    relax function and choose which one to use"""
    return statpt_step()


def statpt_step():
    if os.path.exists('def_inp'):
        os.remove('def_inp')

    statpt_status = run_turbomole_module('statpt')
    if statpt_status is False:
        turbomole_logger.critical('Statpt  failed')
        turbomole_logger.error('check %s' % os.getcwd())
        return False


    if os.path.exists('def_inp'):
        define_status = generate_internal_coordinates()
        if define_status is False:
            turbomole_logger.critical('Statpt  failed')
            turbomole_logger.error('check %s' % os.getcwd())
        else:
            turbomole_logger.debug('generated internal coordinates')
        os.remove('def_inp')

    soft_abend('statpt')

    if os.path.isfile('abend.statpt'):
        remove('abend.statpt')
        statpt_status = run_turbomole_module('statpt')
        if statpt_status is False:
            turbomole_logger.critical('Statpt  failed')
            turbomole_logger.error('check %s' % os.getcwd())
            return False
        soft_abend('statpt')
        if os.path.isfile('abend.statpt'):
            remove('abend.statpt')
            statpt_cart()
            if os.path.exists('def_inp'):
                generate_internal_coordinates()
                turbomole_logger.debug('generated internal coordinates')
                os.remove('def_inp')
        remove('relax_problem')

    return check_geometry_convergence()


def run_turbomole_module(module):
    output_file = module+'.out'
    with open(output_file, 'w') as fc:
        try:
            subp.check_call([module], stdout=fc, stderr=fc)
        except subp.CalledProcessError as e:
            time.sleep(1)
            return check_output_for_fail(module)
    return True


def check_output_for_fail(module):
    key_phrase = module+' : all done'
    log_file = module+'.out'
    if key_phrase not in open(log_file).read():
        turbomole_logger.critical('Error in %s step' % module)
        with open('GEO_OPT_FAILED', 'w') as ferr:
            ferr.write("ERROR in %s step, check the output" % module)
            remove('GEO_OPT_RUNNING')
        return False
    return True


def check_geometry_convergence():
    convergence_status = False

    run_turbomole_module('convgrep')

    if os.path.isfile('converged'):
        remove('nextstep')
        remove('optinfo')
        turbomole_logger.info("THE OPTIMIZATION IS CONVERGED.\n")
        convergence_status = True
        remove("not.converged")
        remove('GEO_OPT_RUNNING')
    return convergence_status


def rewrite_turbomole_gradient(number_of_atoms, total_restraint_energy, total_restraint_gradients,
                               restraint_gradient):
    """Rewrite new gradient file in turbomole format"""
    import re
    with open('gradient', 'r') as gradient_file:
        gradient_contents = gradient_file.read()
        gradients = re.split("cycle", gradient_contents)
        prelog = "cycle".join(gradients[:-1])
        last_gradient = gradients[-1]
        list_from_last_gradient = last_gradient.split('\n')

        cycle_number = int(list_from_last_gradient[0].split()[1])
        scf_energy = float(list_from_last_gradient[0].split()[5])
        total_gradients = float(list_from_last_gradient[0].split(' =')[-1])
        new_energy = scf_energy + total_restraint_energy
        new_total_gradients = total_gradients + total_restraint_gradients
        new_gradients = "cycle = %6d    SCF energy =%20.10f   |dE/dxyz| =%10.6f \n" \
                        % (cycle_number, new_energy, new_total_gradients)

        for line in list_from_last_gradient[1:number_of_atoms + 1]:
            new_gradients = new_gradients + line + "\n"

        i = 0
        for line in list_from_last_gradient[number_of_atoms + 1:number_of_atoms * 2 + 1]:
            dx, dy, dz = line.split()[0], line.split()[1], line.split()[2]
            dx = float(re.sub('D', 'E', dx)) - restraint_gradient[i][0]
            dy = float(re.sub('D', 'E', dy)) - restraint_gradient[i][1]
            dz = float(re.sub('D', 'E', dz)) - restraint_gradient[i][2]
            temp_line = "{:22.13f} {:22.13f} {:22.13f}".format(dx, dy, dz)
            temp_line = re.sub('E', 'D', temp_line)
            new_gradients = new_gradients + str(temp_line) + "\n"
            i += 1
    with open('gradient', 'w') as g:
        g.write(prelog + new_gradients + '$end')


def rewrite_turbomole_energy(total_restraint_energy):
    import re
    with open('energy', 'r') as energy_file:
        energy = energy_file.readlines()
        till_last_energy = energy[:-2]
        old_energy = energy[-2].split()[1]
        new_energy = '{:15.10f}'.format(float(old_energy) + total_restraint_energy)
    with open('energy', 'w') as new_energy_file:
        new_energy_file.write(''.join(till_last_energy) + re.sub(old_energy, str(new_energy), energy[-2]) + '$end\n')


def rewrite_turbomole_energy_and_gradient_files(number_of_atoms,
                                                restraint_gradients,
                                                total_restraint_energy,
                                                total_restraint_gradients):

    rewrite_turbomole_gradient(number_of_atoms, total_restraint_energy,
                               total_restraint_gradients,
                               restraint_gradients)
    rewrite_turbomole_energy(total_restraint_energy)



def main():

    turbomole_logger = logging.getLogger('pyar.turbomole')
    handler = logging.FileHandler('turbomole.log', 'w')
    turbomole_logger.addHandler(handler)
    turbomole_logger.setLevel(logging.DEBUG)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--charge', type=int, default=0,
                        help='charge')
    parser.add_argument('-m', '--multiplicity', type=int, default=1,
                        help='multiplicity')
    parser.add_argument('--scftype', type=str, default='rhf',
                        choices=['rhf', 'uhf'],
                        help='SCF type (rhf/uhf)')
    parser.add_argument("input_file", metavar='file',
                        type=str,
                        help='input coordinate file')
    parser.add_argument('-o', '--opt', action='store_true',
                        help='optimize')
    parser.add_argument('-g', '--gamma', type=float, default=0.0,
                        help='optimize')
    parser.add_argument('--cycle', type=int, default=350,
                        help='maximum number of optimization cycles')
    args = parser.parse_args()
    from Molecule import Molecule
    mol = Molecule.from_xyz(args.input_file)
    mol.fragments = [[0], [3]]

    method_args = {
        'charge': args.charge,
        'multiplicity': args.multiplicity,
        'scftype': args.scftype,
        'software': 'turbomole'
    }
    geometry = Turbomole(mol, method_args)
    if args.opt:
        import optimiser

        turbomole_logger.info('optimising')
        optimiser.optimise(mol, method_args, args.gamma)
    else:
        make_coord(mol.atoms_list, mol.coordinates)
        prepare_control()
        turbomole_logger.info('created input file')


if __name__ == "__main__":
    main()
