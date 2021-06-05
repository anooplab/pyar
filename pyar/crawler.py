def make_formula(at_ls):
    freq = {items: at_ls.count(items) for items in at_ls}
    return ''.join(f"{key}{value}" for key, value in freq.items())


def collect_files(start, exclude_pattern, pattern):
    import os
    paths = []
    for root, directory, files in os.walk(start):
        for file in files:
            if pattern in file and \
                    'selected' in root and \
                    os.path.splitext(file)[-1] == '.xyz' and \
                    exclude_pattern not in root:
                paths.append(os.path.join(root, file))
    return paths


def collect_data(xyz_files):
    import pandas as pd
    from pyar import Molecule
    from pyar.interface import babel
    ne = []
    for xyz_file in xyz_files:
        inchi_string = babel.make_inchi_string_from_xyz(xyz_file)
        smile_string = babel.make_smile_string_from_xyz(xyz_file)
        atoms_list, mol_coordinates, name, title, energy = Molecule.read_xyz(xyz_file)
        nat = len(atoms_list)
        formula = make_formula(atoms_list)
        ne.append([nat, formula, xyz_file, atoms_list, mol_coordinates, energy, smile_string, inchi_string])
    df = pd.DataFrame(ne, columns=['n_atoms', 'formula', 'Name', 'Atoms', 'coordinates', 'Energy', 'SMILE', 'InChi'])
    df.to_csv('data.csv')
    return df


def find_best_geometry():
    try:
        df = pd.read_csv('data.csv')
    except IOError:
        print(f'data.csv is not found')
        return
    for i in df.loc[df.groupby('formula').Energy.agg('idxmin')].Name:
        print(i)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Program to search for \
             relevant files from a pyar run, extract relevant data, pick \
             global minima, or analyse important properties.')
    parser.add_argument('--starting-directory', default='./aggregates',
                        help='I will search for .xyz files recursively in the inner directories starting from the provided path. default="./aggregates"')
    parser.add_argument('--exclude', metavar='EXCLUDE-PATTERN', default='tmp',
                        help='exclude files that contain <EXCLUE-PATTERN>. defualt = "tmp"')
    parser.add_argument('--pattern', default='result', help='pattern of file names. default = "result_"')
    parser.add_argument('--mine', action='store_true', help='Explore the data')
    parser.add_argument('--find-gm', action='store_true', help='Find global \
                         minimum geometry from the files.')
    args = parser.parse_args()

    if args.mine:
        xyz_files = collect_files(args.starting_directory, args.exclude, args.pattern)
        collect_data(xyz_files)

    if args.find_best_geometries:
        find_best_geometry()


if __name__ == "__main__":
    main()
