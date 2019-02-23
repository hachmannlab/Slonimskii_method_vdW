'''
# python
# 
# Mohammad Atif Faiz Afzal
# April 10, 2018
# 
# A python routine to calculate the van der Waals volume of molecules using Slonimskii's method.
# Note: This funtion only works for molecules containing C, N, O, F, Cl, Br, I, P, S, Si, Se, and Te atoms. 

'''

import sys
#sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel
#import subprocess
#import time
import argparse
from mpi4py import MPI
import glob, os

def print_l(sentence):
    if rank == 0:
        print sentence
        logfile.write(str(sentence)+"\n")

def print_le(sentence, msg="Aborting the run"):
    if rank == 0:
        print sentence
        #logfile.write(sentence+"\n")
        #error_file.write(sentence+"\n")
        sys.exit(msg)
    else:
        sys.exit()
def print_e(sentence):
    if rank == 0:
        print sentence
        #error_file.write(sentence+"\n")


def check_if_mol(smiles, line=None, file_name=None):

    '''
    This function returns True if input string is either a valid SMILES or InChI string. Input is SMILES/InChI string.
    '''

    if check_if_inchi(smiles) == True:            
        this_mol = pybel.readstring("inchi",smiles)
        smiles = str(this_mol)
        smiles = smiles.strip()
    elif check_if_smiles(smiles) == False:    
        ## check if smiles
        if line == None:
            tmp_str = 'Error: The SMILES/InChI string(\'{}\') provided is not valid. Please provide correct SMILES/InChI.'.format(smiles)
        else:
            tmp_str = 'Error: The SMILES/InChI string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct SMILES/InChI.'.format(smiles, line, file_name)
        print_le(tmp_str, "Aborting due to wrong molecule description.")
    return smiles

## This function is to check if the provided SMILES string are valid    
def check_if_smiles(smiles):

    '''
    This function returns True if input string is a valid SMILES string. Input is SMILES string.
    '''

    data = True
    if rank == 0:
        try:
            mol = pybel.readstring("smi", smiles)
        except:
            data = False
    data = comm.bcast(data, root=0)
    
    return data

## This function is to check if the provided InChI string are valid    
def check_if_inchi(inchi):

    '''
    This function returns True if input string is a valid InChI string. Input is InChI string.
    '''

    try:
        mol = pybel.readstring("inchi", inchi)
    except:
        return False
    return True

def get_vdw(smiles):

    '''
    This function returns vdW of a molecule. Input is SMILES string
    '''

    ## reading SMILES

    mol =  pybel.readstring("smi", smiles)
    mol.OBMol.AddHydrogens()

    no_atoms = len(mol.atoms)

    bonds = mol.OBMol.NumBonds()

    no_ar = 0
    no_non_ar = 0

    ## calculate the no of aromatic and non-aromatic rings
    for r in mol.OBMol.GetSSSR():
        if r.IsAromatic():
            no_ar = no_ar+1
        else:
            no_non_ar = no_non_ar+1

    def get_num_struc(smarts):
        smarts = pybel.Smarts(smarts)        
        num_unique_matches = len(smarts.findall(mol))
        return num_unique_matches

    ## Calculating no.of fusions

    no_f_ring_AlAr = get_num_struc('[A]~@[*](~@[a])~@[*](~@[a])~@[A]')
    no_f_ring_AlAl = get_num_struc('[A]~@[*](~@[A])~@[*](~@[A])~@[A]')
    no_f_ring_ArAr = get_num_struc('[a]~@[*](~@[a])~@[*](~@[a])~@[a]')
    no_f_ring_S = get_num_struc('[s]1~@[c](~@[c])~@[c](~@[c])~@[c](~@[c])~@[c]1(~@[c])')


    ## Calculating no.of atoms

    no_of_C = smiles.count('c') + smiles.count('C')
    no_of_N = smiles.count('n') + smiles.count('N')
    no_of_O = smiles.count('o') + smiles.count('O')
    no_of_F = smiles.count('f') + smiles.count('F')
    no_of_Cl = smiles.count('Cl')
    no_of_Br = smiles.count('Br')
    no_of_I = smiles.count('i') + smiles.count('I')
    no_of_P = smiles.count('p') + smiles.count('P')
    no_of_S = smiles.count('s') + smiles.count('S')
    no_of_Si = smiles.count('si') + smiles.count('Si')
    no_of_Se = smiles.count('se') + smiles.count('Se')
    no_of_Te = smiles.count('te') + smiles.count('Te')

    no_of_H = no_atoms - (no_of_C + no_of_N + no_of_O + no_of_F + no_of_Cl + no_of_Br + no_of_I + no_of_P + no_of_S + no_of_Si + no_of_Se + no_of_Te)

    no_of_C = no_of_C - no_of_Cl
    no_of_S = no_of_S - no_of_Si - no_of_Se

    V_vdw = (no_of_H)*7.24 + (no_of_C)*20.58 + (no_of_N)*15.6 + (no_of_O)*14.71 + (no_of_F)*13.31 + (no_of_Cl)*22.45 + (no_of_Br)*26.52 + (no_of_I)*32.52 + (no_of_P)*24.43 + (no_of_S)*24.43 + (no_of_Si)*38.79 + (no_of_Se)*28.73 + (no_of_Te)*36.62

    V_vdw = V_vdw - 5.92*(bonds) - 14.7*(no_ar) - 3.8*(no_non_ar) + 5*(no_f_ring_ArAr) + 3*(no_f_ring_AlAr) + 1*(no_f_ring_AlAl) - 5*(no_f_ring_S) 
    #print 'V_vdw',V_vdw
    return V_vdw 






if __name__ == "__main__":

    ## initializing MPI to time, to check the MPI efficiency
    wt1 = MPI.Wtime()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpisize = comm.Get_size()

    ## Argument parser desription
    parser = argparse.ArgumentParser(description='This code calculates van der Waals volume of molecules using Slonimskii\'s method.\
# Note: This funtion only works for molecules containing C, N, O, F, Cl, Br, I, P, S, Si, Se, and Te atoms.\nInput can be provided as SMILES or any other format supported by OpenBabel. \
Requirements include pybel and mpi4py.\n\n\
Install pybel using the below command or follow instructions on https://pypi.python.org/pypi/openbabel\n\n\
pip install openbabel\n\n\
Install mpi4py using the below command\n\n\
pip install mpi4py', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', "--input", action='store', dest='file_name', default='.dat', help="\
Provide SMILES/InChI string of the molecule. If input is a file, then provide the file name. Default is None. The file can contain multiple SMILES/InChI strings. Do not provide this option if the calculations are to be performed on mutiple files in a batch form. In that case, only provide the molecule_type option. If the input file is a smiles file, provide the file extension as '.smi'.")

    parser.add_argument('-t', "--molecule_type", action='store', dest='mol_type', default='smi', help="Mention the molecule type in this option. Default is smi.")

    parser.add_argument('-o', "--output_file_name", action='store', dest='output_file', default='vdw.out', help="Mention the output file name. Default is vdw.out")

    parser.add_argument('-p', "--print_screen", action='store', dest='print_screen', default='True', help="If print output to the command line.")

    parser.add_argument('-b', "--batch_calc", action='store', dest='batch_cal', default='False', help="Mention if calculations are batch type. True only if the batch calculations are for mutiple files. Please provide molecule type option for batch calculations. Do not provide -i option when performing batch calculations.")

    #parser.add_argument('-p', "--print_screen", action='store', dest='print_screen', default='True', help="If the outout to display on the command line.")

    args = parser.parse_args()
    file_name = args.file_name
    mol_type = args.mol_type.lower()
    logfile = open(args.output_file,'w')
    print_screen = args.print_screen.lower()
    batch_cal = args.batch_cal.lower()

    print_l('Molecule\tvdW_volume')

    if file_name[-4:] == '.smi' or  file_name[-7:] == '.inchi':
        for i,line in enumerate(open(file_name,'r')):
            smiles = line.strip()
            if smiles.isspace() or len(smiles) == 0 or smiles[0] == '#':
                continue
                
            smiles = check_if_mol(smiles,i+1,args.file_name)
            vdw = get_vdw(smiles)
            print_l(smiles+'\t'+str(vdw))
        sys.exit()

    if batch_cal == 'false':
        
        if mol_type == 'smi' or mol_type == 'inchi':
            smiles = check_if_mol(file_name)
            vdw = get_vdw(smiles)
            
            print_l(smiles+'\t'+str(vdw))
            sys.exit()

        mol = pybel.readfile(mol_type, file_name).next()
        vdw = get_vdw(mol.write("smi"))
        print_l(file_name+'\t'+str(vdw))
        sys.exit()
        
    print file_name

    
    if batch_cal == 't' or batch_cal == 'true':
        for infile in glob.glob("*."+mol_type):

            mol = pybel.readfile(mol_type, infile).next()
            
            vdw = get_vdw(mol.write("smi"))
            print_l(infile+'\t'+str(vdw))
        sys.exit()
    
