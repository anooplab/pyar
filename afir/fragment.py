"""
fragment_make.py - to make 'fragment' file from 2 xyz files in PyAR
"""
import sys
from functools import reduce

'''
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
'''


def get_no_of_atoms(files):
    """This function will test the xyz file format.
       if it finds it alright, then it will return
       the number of atoms.
    """
    with open(files) as xyzfile:
        filecont = xyzfile.readlines()
    no_of_atoms = len(filecont) - 2
    if no_of_atoms != int(filecont[0]):
        raise IncorrectXyzFileError
    return no_of_atoms


def get_no_of_atoms_list(*files):
    """This function is for making a list of all the atom numbers
       in all the xyz files (fragments) provided
    """
    all_file_atoms = []
    for ifile in files:
        no_of_atom = get_no_of_atoms(ifile)
        all_file_atoms.append(no_of_atom)
    return all_file_atoms


def get_fragment_string(all_xyz_file_atom_number):
    """This function will make the lines specified by the dftd3 code
       (modified by A. Anoop & M. Waller) for fragment file.
    """
    init = 1
    final = 0
    line = []
    for i in range(len(all_xyz_file_atom_number)):
        final += int(all_xyz_file_atom_number[i])
        string = str(init) + "-" + str(final)
        init = final + 1
        line.append(string)
    return line


def frag_lines_mod(lines_list):
    """This function will check if for a line the number in both side of
       the 'hiphen' is same or not. If they same, then it will replace the
       line with one number.
    """
    return_line = []
    for lines in lines_list:
        start, end = map(int, lines.split('-'))
        if start == end:
            lines = start
            return_line.append(lines)
        else:
            return_line.append(lines)
    return return_line


def make_fragment_file(number_of_atoms_in_fragement_1, number_of_atoms_in_fragment_2):
    """This is the functional which by help of the above functions
       create the fragment file
    """
    no_of_atoms_list = [number_of_atoms_in_fragement_1, number_of_atoms_in_fragment_2]
    lines_raw = get_fragment_string(no_of_atoms_list)
    all_lines_for_fragment = frag_lines_mod(lines_raw)
    with open('fragment', 'w') as f:
        for i in all_lines_for_fragment:
            f.write("%s\n" % i)
    return


def main():
    pass


if __name__ == "__main__()":
    main()


def read_fragments(n):
    import operator
    try:
        fp = open('./fragment').readlines()
    except IOError:
        print("./fragment file should exist!")
        sys.exit()
    lines = [lines for lines in fp if lines.strip()]
    atoms_in_fragments = []

    for line in lines:
        for delimitter in ',;':
            line = line.replace(delimitter, ' ')
        a = []
        for things in line.split():
            if "-" in things:
                start, end = map(int, things.split('-'))
                a += [i for i in range(start, end + 1)]
            else:
                a += [int(things)]
        atoms_in_fragments.append(a)
    collected_atoms = set([item for sublist in atoms_in_fragments for item in sublist])
    if len(collected_atoms) < n:
        all_atoms_in_molecules = set([i for i in range(1, n + 1)])
        remaining = list(all_atoms_in_molecules - collected_atoms)
        atoms_in_fragments.append(remaining)
    common = reduce(operator.iand, map(set, atoms_in_fragments))
    if common:
        print("Fragments contain common atoms")
        sys.exit()
    n_fragments = len(atoms_in_fragments)
    return n_fragments, atoms_in_fragments