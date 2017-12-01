"""
files_dirs.py - handling files and directories in PyAR
"""
'''
Copyright (C) 2016 by AnoopLab at Indian Institute of Technology Kharagpur, India
Authro: Surajit Nandi

file_manager.py is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

import glob
import os
import shutil
import sys


def important_file_check(*arguments):
    """This function will check if the files given as argument exists
       or not
    """
    for i in arguments:
        try:
            open(i).readlines()
        except:
            raise IOError("file <%s> does not exist" % i)
    return


def file_exists_check_with_return_status(file_path):
    """This function will check if files in a location exists or not
       If exists, it will return True else, it will return False.
    """
    if os.path.exists(file_path):
        status = True
    else:
        status = False
    return status


def make_directories(*arguments):
    """This function will generate the directories given as arguments.
       It will check if the directory exists. If not, it will 
       make the directory. If the directory exists, it will first delete
       it and make new directory.
    """
    for dirs in arguments:
        if not os.path.isdir(dirs):
            os.makedirs(dirs)
        else:
            shutil.rmtree(dirs)
            os.makedirs(dirs)
    return


def delete_files(*arguments):
    """This function will delete files which ar supplied as
       arguments.
    """
    for files in arguments:
        if os.path.isfile(files):
            os.remove(files)
        else:
            print("file: ", files, "does not exists.")
    return


def delete_directories(*arguments):
    """This function will delete directories which are provided with
       arguments if they exists
    """
    for directory in arguments:
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        else:
            print("directory: ", directory, "does not exists.")
    return


def get_files(extension, destination="./"):
    """This function will return all file names as a list from a directory
    'destination'. By default it is current directory. Extension should be
    simply "xyz", "out" , "log" etc.
    """
    file_list = []
    files_needed = destination + "*." + extension
    print(files_needed)
    files = glob.glob(files_needed)
    for ifile in files:
        file_list.append(os.path.basename(ifile))
    return file_list


def get_dirs_files(destination="./", wildcard="*"):
    """This function will return all the subdirectory of the parent directory
       destination. Note, it will use glob. So, will not test for files"""
    dirs = glob.glob(destination + "/" + wildcard)
    subdirs = []
    for i in dirs:
        last_name = os.path.basename(i)
        subdirs.append(last_name)
    subdirs.sort()
    return subdirs


def scopy(present, destination):
    """This function will copy a file to a location. destination_name should not
       contain the pathname starting with slash '/'. For details, see os.path.join
       documentation in www.python.org website.
    """
    shutil.copy(present, destination)
    return


def smove(current_file, moved_file, path_to_move):
    """This function will move a file current_file to a destination directory
       path_to_move with a name moved_file
    """
    final_file = os.path.join(path_to_move, moved_file)
    shutil.move(current_file, final_file)
    return


def mmove(*args):
    """
    :param args: This is the argument lists. It will take tha last name as the directory name.
    If any of the files match with the strings it provided, it will move this file to the
    directory
    :return: It will return nothing.
    """
    dirname = args[-1]
    print("dirname:", dirname)
    for items in args[:-1]:
        for name in glob.glob('./*' + str(items) + '*'):
            smove(name, name, dirname)
    return


def mcopy(*args):
    """
    :param args: This are the argument lists. It will take the last name as directory name.
     if any of the files match with the strings it provided, it will copy this files to the
     directory. Similar to mmove function
    :return: Nothing
    """
    dirname = args[-1]
    print("dirname:", dirname)
    for items in args[:-1]:
        for name in glob.glob("./*" + str(items) + "*"):
            scopy(name, dirname)
    return


def write_result(filename, **kwargs):
    """This function will write the results as final files.
    """
    mode = kwargs.pop("mode", 'w')
    result_file = open(filename, mode)
    for (name, value) in kwargs.items():
        result_file.write("%s = %s   " % (name, value))
    result_file.write("\n")
    result_file.close()


def check_stop(filename):
    """ This function is for checking by presence of a file. If this file
        present, The program will stop after completion of a cycle. This
        function is like the 'STOP' file concept in turbomole.
    """
    if os.path.isfile(filename):
        print("file", filename, "is present.")
        print("The program will stop.")
        sys.exit(0)


def main():
    pass


if __name__ == "__main__":
    main()
