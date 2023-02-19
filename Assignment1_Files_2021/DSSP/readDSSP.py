#!/usr/bin/python3

# Author: Nick Bounatsos
# Student Number: 2768686

from sys import argv
import os


def read_AccUnfold():
    """ Read file with maximal possible surface accessibility
    per residue. Returns the values in a dictionary """
    # open file
    f = open('Data/AccUnfold.data', 'r')
    
    # create dictionary
    unfolded_acc = {}
    for line in f:
        splitline = line.strip().split()
        # unfolded_acc[amino_acid] = area
        unfolded_acc[splitline[0]] = float(splitline[2])
    f.close()
    return unfolded_acc

def decide_if_buried(aa, acc, unfolded_acc):
    """ Calculate the fraction of buried surface area. If this
    fraction is less than seven percent, the amino acid is
    considered buried """

    """
    The fraction of buried surface in a protein can be calculated using the following formula:
    Fraction of buried surface = (Total solvent accessible surface area - Total exposed surface area) / Total solvent accessible surface area
    Here, "solvent accessible surface area" (SASA) refers to the area of the protein surface that is accessible to solvent molecules, 
    and "exposed surface area" (ESA) refers to the area of the protein surface that is exposed to solvent molecules. 
    The difference between the total solvent accessible surface area and the total exposed surface area gives the total buried surface area.
    """

    buried = False
    ### START CODING HERE
    # normal amino acids
    # print(aa)
    # print(acc)
    # print(unfolded_acc)
    # print(unfolded_acc[aa])
    fraction = acc/unfolded_acc[aa]
    # print(fraction)
    if fraction < 0.07: buried = True
    # print(buried)
    ### END CODING HERE
    return buried

def read_dir(d, unfolded_acc):
    """ For each DSSP file in the input directory, extract the
    total number of each amino acid type and the number of buried
    amino acids per amino acid type
    """
    # initialize dictionaries
    # count all amino acids
    #   all_aa_count[aa] = count
    all_aa_count = {}
    # count buried amino acids
    #   buried_aa_count = count
    buried_aa_count = {}

    # parse all DSSP files in the directory
    for filename in [fn for fn in os.listdir(d) if fn.endswith('.dssp')]:
        f = open(os.path.join(d, filename), 'r')
        
        # start reading at the first line starting with '  # '
        start_reading = False
        for line in f:
            if line.startswith('  # '):
                start_reading = True
            elif start_reading:
                line = line.rstrip()

                # amino acid type
                aa_type = line[13]

                # skip amino acids marked '!'
                if aa_type == '!':
                    continue

                if aa_type.islower():
                    print(aa_type)


                ### START CODING HERE
                # write conditional statements indicating what
                # needs to happen with 'special' residues

                if aa_type.islower():
                    aa_type = aa_type.upper()

                # skip unknown amino acids
                if aa_type == 'X' or aa_type == 'B' or aa_type == 'J' or aa_type == 'O' or aa_type == 'U' or aa_type == 'Z':
                    continue

                ### END CODING HERE

                # residue number
                res_num = int(line[5:10])
                # chain ID
                chain = line[11]
                # accessible surface area
                acc = float(line[34:38])

                # if the amino acid type is in the dictionaries,
                # add 1 to the total count, else create a new
                # entry in both dictionaries
                if aa_type in all_aa_count:
                    all_aa_count[aa_type] += 1
                else:
                    all_aa_count[aa_type] = 1
                    buried_aa_count[aa_type] = 0

                # check if the amino acid is buried
                buried = decide_if_buried(aa_type, acc, unfolded_acc)
                if buried:
                    buried_aa_count[aa_type] += 1

        # close file
        f.close()
    # return dictionaries
    return all_aa_count, buried_aa_count
# print(read_dir('/home/banoffee/Documents/sb-project/Assignment1_Files_2021/DSSP/Data/DSSP_files_small_lib', read_AccUnfold()))

def print_propensities(all_aa_count, buried_aa_count, outfile):
    """ For each amino acid, calculate the propensity to be
    buried and write it into an output file """
    f = open(outfile, 'w')
    list_aa = sorted(all_aa_count.keys())

    ### START CODING HERE
    # you should calculate the propemsity for each amino acid type to be buried
    # you can use all_aa_count and buried_aa_count[aa]
    # you can use the following loop structure over all amino acids:

    # print(all_aa_count)
    # print(buried_aa_count)

    # n_aa_s = []
    # n_aa = []
    n_total_s = sum(buried_aa_count.values())
    n_total = sum(all_aa_count.values())

    y = n_total_s / n_total

    # print(list_aa)

    for aa in list_aa:
        x = buried_aa_count[aa] / all_aa_count[aa]

        propensity_buried = x/y

        # print(propensity_buried)

        # pass
    
    # print(n_aa, n_aa_s)

    # to print to the output file you can use:
        print(aa, propensity_buried, file=f)

    ### END CODING HERE

def main():
    # check input directory
    d_in = argv[1]
    if not os.path.isdir(d_in):
        print(d_in, 'is not a directory')
        exit(1)
    outfile = 'Output/propensity_buried.txt'
    unfolded_acc = read_AccUnfold()
    all_aa_count, buried_aa_count = read_dir(d_in, unfolded_acc)
    print_propensities(all_aa_count, buried_aa_count, outfile)
    # print('Output to file ', outfile)

if __name__ == '__main__':
    main()