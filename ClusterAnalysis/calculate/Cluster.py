#!/usr/bin/env python

'''
This module is used to store information from MD simulations run with AMBER. 
'''

import math
import numpy as np
from numpy import linalg as LA
import sys
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)
        #create nested dict in one go
# from sklearn.cluster import AgglomerativeClustering
from ..write.Writers import Writers

import MDAnalysis       
from MDAnalysis import *

class Cluster(object):
    '''
    Base class for cluster analyses.
    '''

    def __init__(self):
        self.res_name_id_dict = nested_dict() #solute molecule count
        self.cluster_dict = nested_dict()
        self.per_frame_cluster_sizes = []
        self.per_frame_cluster_max_size = []
        self.per_frame_cluster_counts = []
        self.solute_water_hb_dict = nested_dict()
        self.solute_solute_dict = nested_dict()

    def get_clusters(self, molecule, solvent, excluded, hydrogen_bonding):
        '''From neighbours within a cutoff distance, find all atoms that
        are in the same cluster. For each atom, get the neighbour list, then 
        for each atom in the neighbour list, append the neighbours neighbours 
        to the same cluster, skipping repeats.
        '''

        resid_list = []
        for resname, resid in zip(molecule.resnames[molecule.solute_atoms], 
                molecule.resids[molecule.solute_atoms]):
            if resname not in self.res_name_id_dict:
                self.res_name_id_dict[resname] = 0
            if resid not in resid_list:
                self.res_name_id_dict[resname] += 1
                resid_list.append(resid)

        n_atoms = len(molecule.solute_atoms)
        n_frames = len(molecule.atoms_neighbours)
        filename = '{}.pdb'.format('clusters')
        open(filename, 'w')
        for f in range(n_frames):

            '''
            ## get neighbour list for each atom within the max_cutoff r
            dimensions = molecule.dimensions[f]
            neighbour_coords = molecule.coords[f]
            neighbours_list = []
            for a in range(n_atoms):
                atom_coords = molecule.coords[f][a]
                neighbours = Cluster.get_neighbours(atom_coords, 
                        neighbour_coords, dimensions, max_cutoff)
                neighbours_list.append(neighbours)
            '''
            neighbours_list = molecule.atoms_neighbours[f]

            if hydrogen_bonding:
                self.solute_water_hb_dict = Cluster.get_hb_counts(
                        self, self.solute_water_hb_dict, f, n_frames, 
                        neighbours_list, solvent, molecule)

            self.solute_solute_dict = Cluster.get_sol_sol_counts(self, 
                    self.solute_solute_dict, f, n_frames, neighbours_list,
                    solvent, excluded, molecule)

            ## now find clusters, skipping those atoms that are already 
            ## assigned to a cluster. Need to used indices rather than
            ## atom nums so that we can go back into neighbour_list
            clusters_list = []
            index_list = []
            for atom_i in molecule.solute_atoms:
                index_i = molecule.solute_atoms.index(atom_i)
                if index_i not in index_list:
                    i_neighbours = []
                    for index_k in neighbours_list[index_i]:
                        if index_k in molecule.solute_atoms:
                            new_index_k = molecule.solute_atoms.index(
                                    index_k)
                            i_neighbours.append(new_index_k)
                        else:
                            continue

                    for index_m in i_neighbours:
                        if index_m not in index_list:
                            index_list.append(index_m)
                        for index_j in neighbours_list[index_m]:
                            if index_j in molecule.solute_atoms:
                                new_index_j = \
                                        molecule.solute_atoms.index(
                                        index_j)
                                if new_index_j not in i_neighbours:
                                    i_neighbours.append(new_index_j)
                                else:
                                    continue
                            else:
                                continue


                    #create a new list of separate clusters only
                    if len(i_neighbours) > 0:
                        clusters_list.append(sorted(i_neighbours))

            ## label each cluster with a number between 0 and 1
            n_clusters = len(clusters_list)
            ave_cluster_size = 0
            max_cluster_size = 0
            count = 0
            cluster_indices = np.zeros((n_atoms))
            for i in clusters_list:
                ave_cluster_size += len(i)
                if len(i) > max_cluster_size:
                    max_cluster_size = len(i)
                self.cluster_dict = Cluster.get_cluster_populations(
                        self.cluster_dict, molecule, i)
                cluster_size = len(i)
                count += 1
                for j in i:
                    cluster_indices[j] = \
                            count / n_clusters
            ave_cluster_size = ave_cluster_size / n_clusters

            #print(f, n_clusters, ave_cluster_size, max_cluster_size)
            self.per_frame_cluster_sizes.append(ave_cluster_size)
            self.per_frame_cluster_max_size.append(max_cluster_size)
            self.per_frame_cluster_counts.append(n_clusters)

            ## write to pdb file
            Writers.write_pdb(cluster_indices, cluster_indices, 
                    molecule, molecule.solute_coords[f], filename, 'a', f)

        Writers.write_csv([range(n_frames), self.per_frame_cluster_sizes,
                self.per_frame_cluster_counts, 
                self.per_frame_cluster_max_size], 
                'per-frame-clusters', 
                header='frame,ave_cluster_size,n_clusters,max_cluster_size')

        sizes = []
        counts = []
        for info, atoms_key in self.cluster_dict.items():
            if info == 'size':
                for size, count in sorted(list(atoms_key.items())):
                    #print(size, count)
                    sizes.append(size)
                    counts.append(count)
            '''
            else:
                print()
                for size2, atoms2_key in sorted(list(atoms_key.items())):
                    print(size2, self.cluster_dict['size'][size2])
                    for atoms, count2 in sorted(list(atoms2_key.items())):
                        print('\t', atoms, count2)
            '''

        Writers.write_csv([sizes, counts], 
                'cluster-size-distribution', 
                header='cluster_size,frequency')

        if len(self.solute_water_hb_dict) > 0:
            Writers.write_csv([range(n_frames)] + 
                        list(self.solute_water_hb_dict.values()), 
                    'solute-water-contact-counts', 
                    header='%s' % (','.join(map(str, 
                    ['frame'] + list(self.solute_water_hb_dict.keys())
                    )))
                    )

        if len(self.solute_solute_dict) > 0:
            Writers.write_csv([range(n_frames)] + 
                        list(self.solute_solute_dict.values()), 
                    'solute-solute-contact-counts', 
                    header='%s' % (','.join(map(str, 
                    ['frame'] + list(self.solute_solute_dict.keys())
                    )))
                    )

    def get_hb_counts(self, hb_dict, f, n_frames, neighbours_list,
            solvent, molecule):
        '''If a solute atom neighbours a water molecule, count this on
        a per-frame basis.'''
        for atom_i in molecule.solute_atoms:
            n_mols = self.res_name_id_dict[molecule.resnames[atom_i]]
            res_atom = '{}_{}_{}'.format(molecule.resnames[atom_i], 
                    molecule.atom_names[atom_i], n_mols)
            index_i = molecule.solute_atoms.index(atom_i)
            for atom_j in neighbours_list[index_i]:
                if molecule.resnames[atom_j] in solvent:
                    if res_atom not in hb_dict:
                        hb_dict[res_atom] = np.zeros((n_frames))
                    hb_dict[res_atom][f] += 1 / n_mols
                else:
                    continue
        return hb_dict

    def get_sol_sol_counts(self, sol_sol_dict, f, n_frames, neighbours_list,
            solvent, excluded, molecule):
        '''If a solute atom neighbours any other solute atom not with the
        same resid, count this on a per-frame basis.'''
        for atom_i in molecule.solute_atoms:
            n_mols = self.res_name_id_dict[molecule.resnames[atom_i]]
            res_atom = '{}_{}_{}'.format(molecule.resnames[atom_i], 
                    molecule.atom_names[atom_i], n_mols)
            index_i = molecule.solute_atoms.index(atom_i)
            for atom_j in neighbours_list[index_i]:
                if molecule.resids[atom_i] != molecule.resids[atom_j]:
                    if molecule.resnames[atom_j] not in solvent+excluded:
                        if res_atom not in sol_sol_dict:
                            sol_sol_dict[res_atom] = np.zeros((n_frames))
                        sol_sol_dict[res_atom][f] += 1 / n_mols
                    else:
                        continue
                else:
                    continue
        return sol_sol_dict

    def get_cluster_populations(cluster_dict, molecule, cluster):
        cluster_size = len(cluster)
        if cluster_size not in cluster_dict['size']:
            cluster_dict['size'][cluster_size] = 0
        cluster_dict['size'][cluster_size] += 1
        for index_i in cluster:
            atom_i = molecule.solute_atoms[index_i]
            atom = molecule.atom_names[atom_i]
            if atom not in cluster_dict['atoms'][cluster_size]:
                cluster_dict['atoms'][cluster_size][atom] = 0
            cluster_dict['atoms'][cluster_size][atom] += 1
        return cluster_dict


    def get_neighbours(atom_coords, neighbour_coords, dimensions, 
            max_cutoff):
        '''Get atom indices for neighbours within a distance cutoff'''
        pairs, dists = MDAnalysis.lib.distances.capped_distance(atom_coords, 
                neighbour_coords, max_cutoff=max_cutoff, min_cutoff=None, 
                box=dimensions, method=None, return_distances=True)
        nearest = np.argsort(np.array(dists))
        ## sort neighbours from closest to furthest
        neighbours = np.array([pair[1] for pair in pairs])[nearest]
        #neighbours = sorted([pair[1] for pair in pairs])
        return neighbours

    def find_bonded_atoms(atoms, coords, max_cutoff):
        '''Get adjacency matrix'''
        n_atoms = len(atoms)
        A = np.zeros((n_atoms, n_atoms))
        for i in range(n_atoms):
            for j in range(i):
                r = Cluster.get_r(coords[i], coords[j])
                #print(i+1,j+1,r)
                if r <= max_cutoff:
                    A[i,j] = 1
                    A[j,i] = 1
                else:
                    continue
        return A

    def get_r(coordsA, coordsB):
        return np.linalg.norm(coordsA-coordsB)
