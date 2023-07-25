#!/usr/bin/env python
"""
A script for running solute cluster analysis. 
"""

from datetime import datetime
import argparse
import sys
import numpy as np
from collections import Counter
from ClusterAnalysis import Molecule, AMBERParser, Cluster, Writers

def run_cluster_analysis(top_file='top_file', traj_file='traj_file', 
        start='start', end='end', step='step', max_cutoff='max_cutoff', 
        solvent='solvent', excluded='excluded', 
        hydrogen_bonding='hydrogen_bonding'):

    startTime = datetime.now()
    print(startTime)
    # create list of residues excluded from solute cluster analysis
    excluded = excluded + solvent 
    # initiate molecule class
    molecule = Molecule()
    # parse Amber topology and trajectory files
    AMBERParser(top_file, traj_file, molecule, start, end, step, excluded,
            hydrogen_bonding, max_cutoff)
    # print out information for user about solute and solvent in the system
    print('HB: {}'.format(hydrogen_bonding))
    print(molecule)
    resnames_count = Counter(molecule.resnames)
    print('residues in system:', resnames_count.keys())
    print('solvent:', solvent)
    print('excluded from cluster analysis:', excluded)
    # find solute clusters and optional analysis of waters around solute atoms
    cluster = Cluster()
    print('running cluster and HB analysis')
    Cluster.get_clusters(cluster, molecule, solvent, excluded, 
            hydrogen_bonding)


    print(datetime.now() - startTime)


def main():

    try:
        usage = 'runClusterAnalysis.py [-h]'
        parser = argparse.ArgumentParser(description='Program for '\
                'analysing solute clusters from MD simulations.', usage=usage, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('Options')
        group = parser.add_argument('-top', '--top_file', metavar='file', 
                help='path to and name of topology file.')
        group = parser.add_argument('-traj', '--traj_file', metavar='file',
                help='path to and name of coordinate trajectory file.')
        group = parser.add_argument('-s', '--start', action='store', type=int,
                default=0, help='start frame number for analysis.')
        group = parser.add_argument('-e', '--end', default=1, action='store',
                type=int, help='end frame number for analysis.')
        group = parser.add_argument('-df', '--step', default=1, type=int,
                action='store', help='step between frames for analysis.')
        group = parser.add_argument('-r', '--max_cutoff', action='store', 
                default=5, type=float, 
                help='cutoff value for neighbour distance.')
        group = parser.add_argument('-solv', '--solvent', nargs='+', 
                default=['WAT'], help='name of water molecules in topology.')
        group = parser.add_argument('-ex', '--excluded', nargs='+', 
                default=[], help='names of solute molecules to exclude '\
                'from cluster analysis.')
        group = parser.add_argument('-hb', '--hydrogen_bonding', 
                action='store_true', 
                help='include hydrogen-bond analysis (includes solvent). '\
                'Use with caution as this is slower to run')

        op = parser.parse_args()
    except argparse.ArgumentError:
        logging.error('Command line arguments are ill-defined, '
        'please check the arguments.')
        raise
        sys.exit(1)

    run_cluster_analysis(top_file=op.top_file, traj_file=op.traj_file, 
            start=op.start, end=op.end, step=op.step, 
            max_cutoff=op.max_cutoff, solvent=op.solvent, 
            excluded=op.excluded, hydrogen_bonding=op.hydrogen_bonding)

if __name__ == '__main__':
    main()


