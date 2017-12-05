import csv


def read_data(filename):
    pass

def main():
    inputs = sys.argv[1:]

    with open('energy.csv', 'a+') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Name", "Energy"])
        for each_file in inputs:
            name, energy = read_data(each_file)
            if name and energy:
                writer.writerow([name, energy])
            else:
                writer.writerow([name, None])


if __name__ == '__main__':
    main()