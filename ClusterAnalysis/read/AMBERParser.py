#!/usr/bin/env python

"""
This module is for reading in forces, coordinates and energies from 
AMBER topolgy and trajectory files. 
"""

import MDAnalysis       
from MDAnalysis import *
import numpy as np
import sys
from ..calculate.Cluster import Cluster


class AMBERParser(object):
    """
    A class used to parse Amber topology and trajectory files.
    """

    def __init__(self, amb_top, amb_coords, molecule, 
            start, end, step, excluded, hydrogen_bonding, max_cutoff):
        """
        Parameters
        ----------
        amb_top : str
            the path and name of the Amber topology file
        amb_coords : str
            the path and name of the Amber coordinate file
        molecule : class
            a class for saving coord, neighbour, atom and residue info
        start : int
            starting frame number from which analysis starts (default is 0)
        end : int
            starting frame number from which analysis starts (default is 1)
        step : int
            step between frames for analysis (default is 1)
        excluded : list of strings
            list of solute molecules to exclude from clustering (solvent 
            included by default)
        hydrogen_bonding : boolean
            to include analysis of water molecules near solute atoms (default 
            is False)
        maxcut_off : float
            maximum distance for inclusion of atoms in neighbour analysis 
            (default is 5)
        """
        self.amb_top = amb_top
        self.amb_coords = amb_coords
        self.filenames = [amb_top, amb_coords]
        self.solute_coords = []
        #self.dimensions = []
        self.atoms_neighbours = []
        self.read_files(start, end, step, excluded, hydrogen_bonding, 
                max_cutoff)
        molecule.populate_molecule(self) #populate molecule class

    def __str__(self):
        return ('\nMD files: %s, \natoms: %s, N atoms: %s, ' \
                '\nN frames: %s, \ncoords: %s' % 
                (self.filenames,
                ' '.join(map(str, self.atoms)), 
                len(self.atoms), len(self.coords), self.coords.shape))

    def read_files(self, start, end, step, excluded, hydrogen_bonding, 
            max_cutoff):
        """Reads in Amber topology and trajectory files and saves nearest 
        neighbour list.

        Parameters
        ----------
        start : int
            starting frame number from which analysis starts (default is 0)
        end : int
            starting frame number from which analysis starts (default is 1)
        step : int
            step between frames for analysis (default is 1)
        excluded : list of strings
            list of solute molecules to exclude from clustering (solvent 
            included by default)
        hydrogen_bonding : boolean
            to include analysis of water molecules near solute atoms (default 
            is False)
        maxcut_off : float
            maximum distance for inclusion of atoms in neighbour analysis 
            (default is 5)
        """

        #read_coords = Universe(self.amb_top, self.amb_coords, format='TRJ')
        read_coords = Universe(self.amb_top, self.amb_coords, format='NC')
        self.atom_names = np.array([tp.name for tp in read_coords.atoms])
        self.masses = np.array([tp.mass for tp in read_coords.atoms])
        self.resids = np.array([tp.resid for tp in read_coords.atoms])
        self.resnames = np.array([tp.resname for tp in read_coords.atoms])
        self.charges = np.array([tp.charge for tp in read_coords.atoms])                     

        self.solute_atoms = [i for i, s in enumerate(self.resnames) 
                if s not in excluded]

        print('n_frames: {}'.format(len(read_coords.trajectory)))

        for f in range(start, end, step):
            print(read_coords.trajectory[f])
            coords = np.array(read_coords.trajectory[f].positions)
            solute_coords = coords[self.solute_atoms]
            self.solute_coords.append(solute_coords)
            dimensions = np.array(read_coords.trajectory[f].dimensions)
            ## find closest neighbours within cutoff, solute or solvent
            neighbour_coords = solute_coords
            if hydrogen_bonding:
                neighbour_coords = coords
            neighbours_list = []
            for a in self.solute_atoms:
                atom_coords = coords[a]
                neighbours = Cluster.get_neighbours(atom_coords, 
                        neighbour_coords, dimensions, max_cutoff)
                neighbours = list(neighbours)
                # if solvent is not included, then find correct atom indices
                # for solvent atoms
                if hydrogen_bonding == False:
                    new_neighbours = []
                    for n in neighbours:
                        new_neighbours.append(self.solute_atoms.index(n))
                    neighbours = new_neighbours
                neighbours_list.append(neighbours)
            self.atoms_neighbours.append(neighbours_list)
            sys.stdout.flush()
 
    def get_3D_array(self, np_list, n_structures):
        """
        Method used to reshape a numpy array to 3D
        """
        return np.reshape(np.vstack(np_list), (n_structures,-1,3))

    def get_2D_array(self, np_list):
        """
        Method used to reshape a numpy array to 2D
        """
        return np.vstack(np_list)

       
