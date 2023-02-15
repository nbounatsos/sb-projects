#!/usr/bin/python3

# Author: Nick Bounatsos
# Student Number: 2768686

from sys import argv
import os
from math import sqrt, atan2, degrees

# Vector functions that we need to calculate the angles
def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1]*v2[2] - v1[2]*v2[1]
    j = v1[2]*v2[0] - v1[0]*v2[2]
    k = v1[0]*v2[1] - v1[1]*v2[0]
    return [i,j,k]

def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

# PDB file parser
def readPDB(PDB_file):
    """ Reads a PDB file and stores the atom
    coordinates and amino acid types of the protein """
    # open the file
    f = open(PDB_file, 'r')

    # dictionaries to store the output
    # pdb atom coordinates:
    #     pdbcoord[chain][residue_number][atom_type] = coordinates
    pdbcoord = {}
    # residue type per chain and residue number (i.e. store the sequence)
    #     pdbseq[chain][resnum] = restype
    pdbseq = {}

    # parse each line in the file
    for line in f:
        # remove whitespace at the end of the line
        line = line.strip()
        # only parse the lines containing atom coordinates
        if line[:4] == 'ATOM':
            # ATOM type (e.g. C-alpha)
            atom_type = line[12:16].strip()
            # AMINO ACID type (e.g. alanine)
            aa_type = line[17:20].strip()
            # residue number
            res_num = int(line[22:26])
            # Protein chain
            chain = line[21]
            # coordinates
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])

            # if chain does not exists create new entry
            if not chain in pdbcoord:
                pdbcoord[chain] = {}
                pdbseq[chain] = {}
            # if resnum does not exists create new entry
            if not res_num in pdbcoord[chain]:
                pdbcoord[chain][res_num] = {}

            # store coordinates as a vector
            pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]
            # store sequence
            pdbseq[chain][res_num] = aa_type

    # close file
    f.close()

    # return dictionaries
    return pdbcoord, pdbseq

def calculateDihedral(a1, a2, a3, a4):

    ### START CODING HERE

    con_vector1 = [a2[i] - a1[i] for i in range(len(a1))]
    con_vector2 = [a3[i] - a2[i] for i in range(len(a1))]
    con_vector3 = [a4[i] - a3[i] for i in range(len(a1))]

    sin_angle = magnitude(con_vector2) * dot_product(con_vector1, cross_product(con_vector2, con_vector3))

    cos_angle = dot_product(cross_product(con_vector1, con_vector2), cross_product(con_vector2, con_vector3))

    angle = atan2(sin_angle, cos_angle)

    dihedral = degrees(angle)

    ### END CODING HERE

    return dihedral

print(calculateDihedral([1, 9, 2], [3, 2, 1], [2, 4, 7], [8, 2, 5]))

def assign_ss(phi, psi):
    """ Assign a secondary structure type based on the phi
    and psi angles of a residue """
    ### START CODING HERE
    if phi <= 0:
        if psi >= 0:
            secondary_structure = 'beta'
        elif psi <=0:
            secondary_structure = "alpha"
    else:
        secondary_structure = "loop"

    ### END CODING HERE
    return secondary_structure
# print(assign_ss(55, 25))

def print_phi_psi(pdbcoord, pdbseq, outfile):
    """ given the PDB coordinates, calculate the dihedral
    angles of all the residues, assign secondary structure
    types and write them into an output file """
    f = open(outfile, 'w')

    # get the chains from the PDB file
    list_chains = sorted(pdbcoord.keys())
    for chain in list_chains:
        # get the sorted residue numbers from the pdbcoord dictionary
        list_residue_numbers = sorted(pdbcoord[chain].keys())
        for res_num in list_residue_numbers:
            try:
                ### START CODING HERE

                n_spot = pdbcoord[chain][res_num]['N']
                pn_spot = pdbcoord[chain][res_num + 1]['N']
                c_spot = pdbcoord[chain][res_num]['C']
                nc_spot = pdbcoord[chain][res_num - 1]['C']
                ca_spot = pdbcoord[chain][res_num]['CA']

                phi = calculateDihedral(nc_spot, n_spot, ca_spot, c_spot)
                psi = calculateDihedral(n_spot, ca_spot, c_spot, pn_spot)
                ss = assign_ss(phi,psi)
                ### END CODING HERE

            except KeyError:
                print('WARNING: KeyError:', KeyError, 'in residue', chain, res_num)
                phi = 0 
                psi = 0
                ss = 'NaN'

            # get amino acid
            aa_type = pdbseq[chain][res_num]
            # write into output file
            print(chain, res_num, aa_type, phi, psi, ss, file=f)
    f.close()
    print('written:', outfile)

def main():
    # input PDB file
    f_in = argv[1]
    f_out = 'Output/phi_psi.txt'
    # read PDB file
    pdbcoord, pdbseq = readPDB(f_in)
    print_phi_psi(pdbcoord, pdbseq, f_out)
    # for testing
    # for i in ['1TIM', '3PG8']:
    #     f_in = '{}.pdb'.format(i)
    #     print(f_in)
    #     f_out = 'student/output/phi_psi_{}.txt'.format(i)
    
    #     # read PDB file
    #     pdbcoord, pdbseq = readPDB(f_in)
    #     print_phi_psi(pdbcoord, pdbseq, f_out)

if __name__ == '__main__':
    main()
