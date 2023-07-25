#!/usr/bin/env python

'''
This module is for writing info into a particular format. 
'''

import numpy as np

class Writers(object):
    '''
    '''

    def write_xyz(coords, atoms, filename, open_type, i=False):
        xyz_file = open(filename, open_type)
        for molecule in range(len(coords)):
            count = molecule
            if i:
                count = i
            xyz_file.write('{}\n{}\n'.format(len(atoms), count+1))
            for atom in range(len(atoms)):
                x = coords[molecule][atom][0]
                y = coords[molecule][atom][1]
                z = coords[molecule][atom][2]
                xyz_file.write('{:4} {:11.6f} {:11.6f} {:11.6f}\n'.format(
                        atoms[atom], x, y, z))
        xyz_file.close()

    def write_pdb(a, b, molecule, coords, filename, open_type, frame):
        pdb_file = open(filename, open_type)
        #pdb_file.write('MODEL {}\n'.format(frame)) #chimera
        for i, s in zip(range(len(coords)), molecule.solute_atoms):
            pdb_file.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s}' \
                    ' {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}' \
                    '{:6.3f}{:6.3f}          {:>2s}{:2s}'.format(
                    'ATOM', int(str(s)[0:5]), 
                    molecule.atom_names[s], ' ', molecule.resnames[s][0:3], 
                    'A', int(str(molecule.resids[s])[0:4]), ' ', 
                    float(coords[i][0]), 
                    float(coords[i][1]), 
                    float(coords[i][2]), 
                    a[i], b[i], ' ', ' ')]) + '\n')
        pdb_file.write('END\n') #vmd
        #pdb_file.write('ENDMDL\n') #chimera
        #pdb_file.write('TER\n') #amber
        pdb_file.close()


    def write_csv(list_, filename, header):
        np.savetxt('{}.csv'.format(filename), np.column_stack(list_), 
                delimiter=',', header=header, fmt='%.3f')

