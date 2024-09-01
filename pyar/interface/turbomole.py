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
import datetime
import glob
import logging
import os
import getpass
import re
import shutil
import socket
import subprocess as subp
import sys
import tempfile
import time

import numpy as np

from pyar import interface
from pyar.afir import restraints
from pyar.data.units import angstrom2bohr, bohr2angstrom
from pyar.interface import SF
from pyar.interface import which

turbomole_logger = logging.getLogger('pyar.turbomole')


def remove(file_name):
    files = glob.glob(file_name)
    for file in files:
        if os.path.exists(file_name):
            os.remove(file_name)
    return


def safe_rewrite_file(modified_data_groups, file_name):
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        tmp_file.write(modified_data_groups)
    shutil.copystat(file_name, tmp_file.name)
    shutil.move(tmp_file.name, file_name)


def sed_inplace(pattern, replace_to, filename='control'):
    """
    Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    `sed -i -e 's/'${pattern}'/'${replace_to}' "${filename}"`.
    https://stackoverflow.com/questions/4427542/how-to-do-sed-like-text-replace-with-python#4427835

    """
    pattern_compiled = re.compile(pattern)

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(replace_to, line))

    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)


class Turbomole(SF):

    def __init__(self, molecule, qc_params):
        if which('define') is None:
            turbomole_logger.error('set Turbomole path')
            sys.exit('Set turbomole path')

        super(Turbomole, self).__init__(molecule)

        self.start_coords = angstrom2bohr(molecule.coordinates)
        self.atoms_in_fragments = molecule.fragments
        self.energy = None
        self.optimized_coordinates = None
        self.basis = qc_params['basis']
        self.method = qc_params['method']

    def optimize(self):
        """This is the python implementation  of jobex of turbomole

        :returns: Union(True,
                  'SCFFailed',
                  'GradFailed',
                  'UpdateFailed',
                  'CycleExceeded',
                  False)

        """
        # max_cycles = options['opt_cycles']
        # gamma = options['gamma']
        # convergence = options['opt_threshold']

        # if convergence == 'loose':
        #     scf_conv = 6
        # elif convergence == 'tight':
        #     scf_conv = 8
        # else:
        #     scf_conv = 7
        scf_conv = 6
        gamma = 100.0
        max_cycles = 200
        basis_set = self.basis
        functional = 'bp'
        if self.method.lower() == 'bp86':
            functional = 'bp'
        if self.method.lower() == 'b3lyp':
            functional = 'b3-lyp'

        make_coord(self.atoms_list, self.start_coords)
        define_status = prepare_control(scf_conv=scf_conv,
                                        charge=self.charge,
                                        multiplicity=self.multiplicity,
                                        basis=basis_set,
                                        func=functional)
        if define_status is False:
            turbomole_logger.debug('Initial Define failed, converting to cartesian coordinate')
            remove('control')
            define_status = prepare_control(scf_conv=scf_conv, coordinates='cartesian', charge=self.charge,
                                            multiplicity=self.multiplicity)
        if define_status is False:
            turbomole_logger.debug('Initial Define failed again. Quit')
            return 'UpdateFailed'

        if gamma is None:
            return self.run_turbomole_jobex(max_cycles=max_cycles)

        # Test for Turbomole input
        if not os.path.isfile('control'):
            turbomole_logger.critical("Turbomole control file should exist")
            return False

        with open('not.converged', 'w') as fp:
            fp.write('$convcrit\n')
            fp.write('   energy              6\n')
            fp.write('   cartesian gradient  3\n')
            fp.write('   basis set gradient  3\n')
            fp.write('   cycles              1\n')
            fp.write('$convergence not reached\n')
            fp.write('$end\n')

        remove('statistics')
        remove('dscf_problem')

        turbomole_logger.info("Turbomole Optimization started\n")
        turbomole_logger.info("  at  %s\n" % datetime.datetime.now())
        turbomole_logger.info("  on machine %s\n" % socket.gethostname())
        turbomole_logger.info("  by user %s\n" % getpass.getuser())
        turbomole_logger.info("  in directory %s\n" % os.getcwd())

        initial_status, initial_energy = calc_energy()
        if initial_status is False:
            turbomole_logger.warning('Initial energy evaluation failed.')
            return 'SCFFailed'

        turbomole_logger.debug('First step: %s, %f' % (initial_status, initial_energy))

        for cycle in range(max_cycles):
            # Calculate Gradients
            gradient_status = calc_gradients()
            if gradient_status is False:
                turbomole_logger.error('Gradient evaluation failed in cycle %d' % cycle)
                return 'GradFailed'

            for line in open('gradient').readlines():
                if 'cycle' in line:
                    turbomole_logger.debug(line.strip())

            # Calculate afir gradient if gamma is greater than zero
            # if gamma > 0.0:
            afir_energy, afir_gradients = restraints.isotropic(self.atoms_in_fragments, self.atoms_list, get_coords(),
                                                               gamma)
            rewrite_turbomole_energy_and_gradient_files(self.number_of_atoms, afir_energy, afir_gradients)
            turbomole_logger.debug(f'restraint energy = {afir_energy:f}')

            # Update coordinates and check convergence.
            update_status = update_coord()
            if update_status is False:
                turbomole_logger.critical('Coordinate update failed in cycle %d' % cycle)
                turbomole_logger.critical('Check the job in %s' % os.getcwd())
                return 'UpdateFailed'

            convergence_status = check_geometry_convergence()
            if convergence_status is True:
                turbomole_logger.info('converged at %d' % cycle)
                self.energy = get_energy()
                self.optimized_coordinates = bohr2angstrom(get_coords())
                interface.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                                    self.job_name,
                                    energy=self.energy)
                return True

            # Calculate energy
            scf_status, energy = calc_energy()
            if scf_status is False:
                turbomole_logger.critical('Energy evaluation failed in cycle %d' % cycle)
                return 'SCFFailed'

            with open('energy.dat', 'a') as fe:
                # if gamma > 0.0:
                fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy + afir_energy))
                # else:
                #     fe.writelines("{:3d} {:15.8f}\n".format(cycle, energy))
        else:
            turbomole_logger.info("OPTIMIZATION DID NOT CONVERGE WITHIN "
                                  "%d CYCLES\n Restarting it after checking "
                                  "the gradient norms\n might be a good idea..."
                                  " (grep cycle gradient)" % max_cycles)

            self.energy = get_energy()
            self.optimized_coordinates = bohr2angstrom(get_coords())
            return 'CycleExceeded'

    def run_turbomole_jobex(self, max_cycles=350):
        """
        Run a turbomole optimisation.

        return one of the following:
        True: converged
        SCFFailed
        GradFailed
        UpdateFailed
        CycleExceeded
        False: unknown error or jebex execution error

        :rtype: string or boolean
        """
        with open('jobex.out', 'w') as fj:
            try:
                subp.check_call(['jobex', '-ri', '-c', str(max_cycles)], stdout=fj,
                                stderr=fj)
            except subp.CalledProcessError as e:
                turbomole_logger.debug('jobex failed, check %s/jobex.out'
                                       % os.getcwd())
                return False

        if os.path.isfile('GEO_OPT_FAILED'):
            message = open('GEO_OPT_FAILED').read()
            if 'ERROR: Module' in message:
                turbomole_logger.debug('Error in module!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return False
            elif 'ERROR in statpt step,' in message:
                turbomole_logger.debug('Statpt failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'UpdateFailed'
            elif 'ERROR in relax step,' in message:
                turbomole_logger.debug('Relax failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'UpdateFailed'
            elif 'ERROR: your energy calculation did not converge !!,' in message:
                turbomole_logger.debug('SCF failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'SCFFailed'
            elif 'ERROR in dscf step' in message:
                turbomole_logger.debug('SCF failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'SCFFailed'
            elif 'ERROR in gradient step' in message:
                turbomole_logger.debug('Gradient failed!\n chcek %s/GEO_OPT_FAILED' % os.getcwd())
                return 'GradientFailed'
            elif 'OPTIMIZATION DID NOT CONVERGE' in message:
                turbomole_logger.debug('Geometry did not converge in %d cycles.\ncheck %s' % (max_cycles, os.getcwd()))
                self.energy = get_energy()
                self.optimized_coordinates = bohr2angstrom(get_coords())
                return 'CycleExceeded'
            else:
                turbomole_logger.debug('Unknown Error!\n chcek the files in %s' % os.getcwd())
                return False

        elif os.path.isfile('GEO_OPT_CONVERGED'):

            turbomole_logger.debug('jobex done')
            self.energy = get_energy()
            self.optimized_coordinates = bohr2angstrom(get_coords())
            interface.write_xyz(self.atoms_list, self.optimized_coordinates,
                                self.result_xyz_file,
                                self.job_name,
                                energy=self.energy)
            return True


def run_single_point():
    """
    Run a single point calculation

    return one of the following:
    True: SCF Converged
    False: error

    :rtype: string or boolean
    """

    with open('ridft.out', 'w') as fj:
        try:
            subp.check_call(['ridft'], stdout=fj, stderr=fj)
        except subp.CalledProcessError as e:
            turbomole_logger.debug('rdift failed, check %s/rdift.out'
                                   % os.getcwd())
            return False

    return True


def calc_energy():
    run_status = run_turbomole_module('ridft')
    msg = {
        line.strip()
        for line in open('ridft.out').readlines()
        if 'ended' in line
    }

    if run_status is False or 'abnormally' in msg:
        turbomole_logger.debug(msg)
        turbomole_logger.debug('Check the the files in %s' % os.getcwd())
        return False, None
    elif os.path.isfile('dscf_problem'):
        turbomole_logger.debug("SCF Failure. Check files in" + os.getcwd())
        return False, None
    else:
        return True, get_energy()


def calc_gradients():
    run_status = run_turbomole_module('rdgrad')
    msg = [line for line in open('rdgrad.out').readlines() if 'ended' in line]
    if run_status is False or 'abnormally' in msg:
        turbomole_logger.debug(msg)
        turbomole_logger.debug('Gradient calculation failed!')
        turbomole_logger.debug('Chcek files in %s' % os.getcwd())
        return False
    else:
        return True


def choose_coordinate_system(choice):
    if choice == 'int':
        sed_inplace('internal   off', 'internal   on')
        sed_inplace('redundant  off', 'redundant  on')
        sed_inplace('cartesian  on', 'cartesian  off')
        turbomole_logger.debug('changed to internal')
    elif choice == 'crt':
        sed_inplace('internal   on', 'internal   off')
        sed_inplace('redundant  on', 'redundant  off')
        sed_inplace('cartesian  off', 'cartesian  on')
        turbomole_logger.debug('changed to cartesian')
    else:
        turbomole_logger.error('bad choice for coordinate system %s' % choice)


def statpt_cart():
    turbomole_logger.debug('entering cartesian step')
    choose_coordinate_system('crt')
    statpt_status = run_turbomole_module('statpt')
    choose_coordinate_system('int')
    return statpt_status


def make_coord(atoms_list, coordinates, output_file='coord'):
    """
    Create a turbomole 'coord' file from the provided list of atoms and
    cartesian coordinates.

    :type atoms_list: a list of atoms in the coordinates of the molecule.
                      e.g. [H, C, N, O]
    :type coordinates: an array of cartesian coordinates of the molecule,
                       in the same order of atoms.
                       e.g.:
                       [[0.00000 0.00000 0.000000]
                       [1.00000 0.00000 0.000000]
                       [0.00000 1.00000 0.000000]]
    :type output_file: Filename of the output file, usually 'coord'
    """

    with open(output_file, 'w') as fcoord:
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
    """
    This function will prepare the control file. Most of the arguments have
    default values(all have in its current form). The read_define_input var
    if "yes", then define_input_file should provide. The define will run
    from this input file.
    """

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
            fdefine.write("eht\n\n")
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
            if dftd3 == "yes":
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
                turbomole_logger.debug('Error in define')
                turbomole_logger.debug("%s: %s\n" % (e.output, e.returncode))
                return False
        with open(define_log_file) as fout:
            if "define : all done" not in fout.read() or 'abnormally' in fout.read():
                turbomole_logger.critical('Error in define')
                return False
    return True


def actual(*args):
    control = open('control').read()
    if 'reset' in args:
        modified_control = re.sub('actual step', 'last step', control)
        turbomole_logger.error('actual step found in control file, indicating error termination')
        safe_rewrite_file(modified_control, 'control')
    else:
        match = re.search('actual step', control)
        if match:
            turbomole_logger.error('Error in %s' % match.group().split()[-1])
            return False
    return True


def add_tmole_key():
    control = open('control').read()

    if '$tmole' not in control:
        turbomole_logger.debug("Adding '$tmole' to control file")
        new_control = re.sub('end', 'tmole\n$end', control)

        safe_rewrite_file(new_control, 'control')


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
            new_control = re.sub('metric ' + str(metric), 'metric ' + str(metric - 1), control)
            break
    else:
        new_control = re.sub('end', 'redund_inp\n    metric ' + str(metric - 1) + '\n$end', control)

    with open('control', 'w') as fc:
        fc.write(new_control)


def generate_internal_coordinates():
    if not os.path.isfile('coord'):
        turbomole_logger.error("Turbomole file, 'coord', should exist")
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
                if e.output:
                    fout.write(e.output)
                else:
                    turbomole_logger.error('Unknown error in define')

        if 'ATOMIC ATTRIBUTE DATA' in open(define_log_file).read():
            if 'intdef' in open('control').read():
                remove_data_group('cartesian_step')
                remove('hessapprox')

            remove_data_group('redund_inp')
            break
    else:
        remove_data_group('redund_inp')
        turbomole_logger.error('error in define step while generating '
                               'new internal coordinates')
        return False

    remove('tmp.input')
    return True


def remove_data_group(group):
    """similar to turbomole kdg script"""
    data_groups = open('control', 'r').read().split('$')

    for each_group in data_groups[1:]:
        if each_group.split()[0] == group:
            data_groups.remove(each_group)

    modified_data_groups = '$'.join(data_groups)

    safe_rewrite_file(modified_data_groups, 'control')


def get_energy():
    with open('energy') as fp:
        return float(fp.readlines()[-2].split()[1])


def get_gradients():
    return [
        line.replace('D', 'E').split()
        for line in open('gradient').read().split('cycle')[-1].split('\n')
        if len(line.split()) == 3
    ]


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
    if os.path.exists('def_inp'):
        os.remove('def_inp')

    statpt_status = run_turbomole_module('statpt')
    if statpt_status is False:
        turbomole_logger.error('Statpt  failed once')

    if os.path.exists('def_inp'):
        turbomole_logger.info('found def_inp')
        gen_int_status = generate_internal_coordinates()
        if gen_int_status is False:
            turbomole_logger.error('Generation of internal coordinates failed')
        else:
            turbomole_logger.debug('Generated internal coordinates')
        os.remove('def_inp')

    if actual() is False:
        turbomole_logger.debug('trying again')
        statpt_status = run_turbomole_module('statpt')
        if statpt_status is False:
            turbomole_logger.error('Statpt  failed again')
        if actual() is False:
            if statpt_cart() is False:
                turbomole_logger.error('Statpt with cartesian step also failed! Quit')
                turbomole_logger.error('check %s' % os.getcwd())
                return False
            if os.path.exists('def_inp'):
                gen_int_status = generate_internal_coordinates()
                if gen_int_status is False:
                    turbomole_logger.error('Generation of internal coordinates failed')
                else:
                    turbomole_logger.debug('generated internal coordinates')
                os.remove('def_inp')
        remove('relax_problem')

    return True


def run_turbomole_module(module):
    output_file = module + '.out'
    with open(output_file, 'w') as fc:
        try:
            subp.check_call([module], stdout=fc, stderr=fc)
        except subp.CalledProcessError as e:
            time.sleep(1)
            turbomole_logger.info('Error in %s' % module)
            if e.output:
                fc.write(e.output)

    key_phrase = module + ' : all done'
    if key_phrase not in open(output_file).read():
        turbomole_logger.critical('Error in %s' % module)
        return False

    return True


def check_geometry_convergence():
    convergence_status = False

    with open('convgrep.out', 'w') as fc:
        try:
            subp.check_call(['convgrep'], stdout=fc, stderr=fc)
        except subp.CalledProcessError as e:
            fc.write(e.output)
            turbomole_logger.error('Convgrep failed')
            return False
    if os.path.isfile('converged'):
        remove('optinfo')
        turbomole_logger.info("THE OPTIMIZATION IS CONVERGED.\n")
        convergence_status = True
        remove("not.converged")
        remove('GEO_OPT_RUNNING')
    return convergence_status


def rewrite_turbomole_gradient(number_of_atoms, restraint_energy,
                               restraint_gradient):
    """Rewrite new gradient file in turbomole format"""

    dft_grad = np.array(get_gradients(), dtype='float')
    afir_grad = restraint_gradient
    combined_gradients = dft_grad + afir_grad
    combined_total_gradients = np.sqrt(np.sum(combined_gradients ** 2))

    with open('gradient', 'r') as gradient_file:
        gradient_file_contents = gradient_file.read()
    gradients = re.split("cycle", gradient_file_contents)
    contents_upto_current_cycle = "cycle".join(gradients[:-1])
    current_gradient = gradients[-1]
    lines_in_last_gradient = current_gradient.split('\n')

    first_line = lines_in_last_gradient[0].split()
    cycle_number = int(first_line[1])
    scf_energy = float(first_line[5])

    combined_energy = scf_energy + restraint_energy

    new_gradients = "cycle = %6d    SCF energy =%20.10f   |dE/dxyz| =%10.6f \n" \
                    % (cycle_number, combined_energy, combined_total_gradients)

    number_of_atoms = len(dft_grad)
    for line in lines_in_last_gradient[1:number_of_atoms + 1]:
        new_gradients = new_gradients + line + "\n"

    i = 0
    for line in lines_in_last_gradient[number_of_atoms + 1:number_of_atoms * 2 + 1]:
        dx, dy, dz = line.split()[0], line.split()[1], line.split()[2]
        dx = float(re.sub('D', 'E', dx)) - restraint_gradient[i][0]
        dy = float(re.sub('D', 'E', dy)) - restraint_gradient[i][1]
        dz = float(re.sub('D', 'E', dz)) - restraint_gradient[i][2]
        temp_line = "{:22.13f} {:22.13f} {:22.13f}".format(dx, dy, dz)
        temp_line = re.sub('E', 'D', temp_line)
        new_gradients = new_gradients + str(temp_line) + "\n"
        i += 1
    new_string_for_gradients = contents_upto_current_cycle + new_gradients + '$end'
    safe_rewrite_file(new_string_for_gradients, 'gradient')


def rewrite_turbomole_energy(total_restraint_energy):
    import re
    with open('energy', 'r') as energy_file:
        energy = energy_file.readlines()
        till_last_energy = energy[:-2]
        old_energy = energy[-2].split()[1]
        new_energy = '{:15.10f}'.format(float(old_energy) +
                                        total_restraint_energy)

    str_of_new_energy_file = ''.join(till_last_energy) + re.sub(old_energy, str(new_energy),
                                                                energy[-2]) + '$end\n'
    safe_rewrite_file(str_of_new_energy_file, 'energy')


def rewrite_turbomole_energy_and_gradient_files(number_of_atoms,
                                                restraint_energy,
                                                restraint_gradients):
    rewrite_turbomole_gradient(number_of_atoms, restraint_energy,
                               restraint_gradients)
    rewrite_turbomole_energy(restraint_energy)


def movie_maker():
    grad = bohr2angstrom(np.array(get_gradients(), dtype=float))
    coordinates = bohr2angstrom(get_coords())
    energy = get_energy()
    atoms_list = [c.split()[-1] for c in open('coord').read().split('$')[1].split('\n')[1:-1]]

    with open('movie.xyz', 'w') as fp:
        fp.write("%3d\n" % len(coordinates))
        fp.write('original:' + str(energy) + '\n')
        for a, c in zip(atoms_list, coordinates):
            fp.write("{:<2}{:12.5f}{:12.5f}{:12.5f}\n".format(a.upper(),
                                                              c[0], c[1], c[2]))

    steps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
             -1 - 2, -3 - 4 - 5, -6, -7, -8. - 9, -10, -9, -8, -7, -6, -5, -4,
             -3, -2, -1, 0]

    for i in steps:
        with open('movie.xyz', 'a') as fp:
            fp.write("%3d\n" % len(coordinates))
            fp.write(str(i) + ':' + str(energy) + '\n')
            for a, c, g in zip(atoms_list, coordinates, grad):
                fp.write("{:<2}{:12.5f}{:12.5f}{:12.5f}\n".format(a.upper(),
                                                                  c[0] + g[0] * i,
                                                                  c[1] + g[1] * i,
                                                                  c[2] + g[2] * i))


def plot_energy(params):
    import matplotlib.pyplot as plt
    for i in params:
        energies = np.loadtxt(i + '/energy', comments='$', usecols=(0, 1))
        plt.plot(energies[:, 0], energies[:, 1], label=i)
    plt.xlabel('cycles')
    plt.ylabel('energy')
    plt.title('SCF Convergence')
    plt.legend(bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure)
    plt.show()


def main():
    pass


if __name__ == "__main__":
    main()
