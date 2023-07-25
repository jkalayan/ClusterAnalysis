#!/usr/bin/env python

'''
This module is used to store information from MD simulations run with AMBER. 
'''

import math
import numpy as np
import sys

class Molecule(object):
    '''
    Base class for atom info in f frames for N atoms.
    '''
    def __str__(self):
        return ('\nInput files: %s, \nN atoms: %s, ' \
                '\nN solute atoms: %s, \nN frames: %s\n' % 
                (self.filenames, len(self.atom_names), 
                len(self.solute_atoms), len(self.solute_coords)))

    def populate_molecule(self, other):
        self.filenames = other.filenames
        self.atom_names = other.atom_names
        self.masses = other.masses
        self.resids = other.resids
        self.resnames = other.resnames
        self.charges = self.get_2D_array(other.charges).reshape(
                -1,len(self.atom_names))
        self.solute_atoms = other.solute_atoms
        self.solute_coords = self.get_3D_array(other.solute_coords)
        #self.dimensions = np.array(other.dimensions).reshape(-1,6)
        self.atoms_neighbours = other.atoms_neighbours
        del other #remove other to clear memory

    def get_3D_array(self, np_list):
        np_list = np.reshape(np.vstack(np_list), 
                (-1,len(self.solute_atoms),3))
        return np_list

    def get_2D_array(self, np_list):
        return np.vstack(np_list)
