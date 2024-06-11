"""
Handle files and directories in pyar

Copyright (C) 2016 by AnoopLab at Indian Institute of Technology Kharagpur, India
Author: Surajit Nandi

file_manager.py is part of the pyar project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

"""

import glob
import logging
import os
import shutil
import sys

file_manager_logger = logging.getLogger('pyar.file_manager')


def important_file_check(*arguments):
    """Checks if the files given as argument exists

    """

    for i in arguments:
        try:
            open(i).readlines()
        except IOError:
            print(f"file <{i}> does not exist")
    return


def file_exists_check_with_return_status(file_path):
    """This function will check if files in a location exists or not
       If exists, it will return True else, it will return False.

    """

    return bool(os.path.exists(file_path))


def make_directories(*arguments):
    """This function will generate the directories given as arguments.
       It will check if the directory exists. If not, it will 
       make the directory. If the directory exists, it will first delete
       it and make new directory.

    """
    for dirs in arguments:
        if not os.path.isdir(dirs):
            file_manager_logger.debug(f'making directory: {dirs}')
        else:
            file_manager_logger.debug(f'{dirs} exists! remvoing and recreating the directory')

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
    files_needed = f"{destination}*.{extension}"
    print(files_needed)
    files = glob.glob(files_needed)
    return [os.path.basename(each_file) for each_file in files]


def get_dirs_files(destination="./", wildcard="*"):
    """This function will return all the subdirectory of the parent directory
       destination. Note, it will use glob. So, will not test for files

    """
    dirs = glob.glob(f"{destination}/{wildcard}")
    sub_directories = []
    for i in dirs:
        last_name = os.path.basename(i)
        sub_directories.append(last_name)
    sub_directories.sort()
    return sub_directories


def move_file(current_file, moved_file, path_to_move):
    """This function will move a file current_file to a destination directory
       path_to_move with a name moved_file

    """

    final_file = os.path.join(path_to_move, moved_file)
    shutil.move(current_file, final_file)
    return


def bulk_move(*args):
    """

    :param args: This is the argument lists. It will take tha last name as the directory name.
        If any of the files match with the strings it provided, it will move this file to the
        directory
    :return: It will return nothing.

    """
    directory_name = args[-1]
    print("directory_name:", directory_name)
    for items in args[:-1]:
        for name in glob.glob(f'./*{str(items)}*'):
            move_file(name, name, directory_name)
    return


def bulk_copy(*args):
    """
    :param args: This are the argument lists. It will take the last name as directory name.
     if any of the files match with the strings it provided, it will copy this files to the
     directory. Similar to bulk_move function
    :return: Nothing

    """
    directory_name = args[-1]
    print("directory_name:", directory_name)
    for items in args[:-1]:
        for name in glob.glob(f"./*{str(items)}*"):
            shutil.copy(name, directory_name)
    return


def write_result(filename, **kwargs):
    """This function will write the results as final files.

    """
    mode = kwargs.pop("mode", 'w')
    with open(filename, mode) as result_file:
        for name, value in kwargs.items():
            result_file.write(f"{name} = {value}   ")
        result_file.write("\n")


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
