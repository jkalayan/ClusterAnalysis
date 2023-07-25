
# Description
ClusterAnalysis is a tool for finding neighbouring molecules within a specified radial distance. Molecules are defined to be in the same cluster if any atoms in a given molecule are within a user specified cutoff radius from atoms in another molecule. Additionally, users can also use ClusterAnalysis to find which solute atoms are within a radial distance from atoms in water molecules by including the `hydrogen_bonding` flag. Please note however, the inclusion of water molecules in distance calculations will be slower to run. 

AMBER topology and trajectory files are read in using MDAnalysis and information pertaining to Cartesian coordinates, atom types, atom indices, residue types and residue indices are used to analyse molecules with ClusterAnalysis. 

For each simulation frame, nearest atomic neighbour distances within a given cutoff radius are determined for each atom by applying brute force distance calculations to find relative distances via MDAnalysis, specifically using the method `MDAnalysis.lib.distances.capped_distance`. For atoms at periodic boundaries, the minimum image convention is applied when computing distances (checks to see if the distance between two atoms is greater than half the length of the box).

For all neighbouring atoms within the cutoff radius of an initial atom, if a neighbouring atom belongs to another molecule in the system (it will have a different residue index), then that other molecule is considered to be in the same cluster as the molecule of the initial atom and the molecules are assigned to the same cluster. This process is repeated for solute molecules that haven't yet been assigned a cluster. 

Once all simulation frames are processed and each molecule is assigned a cluster in each frame, calculations such as average cluster size, number of clusters and maximum cluster size per frame are outputted to the per-frame-clusters.csv file.

If water is included in neighbour analysis, then the average number of surrounding water molecules within the cutoff radius in each frame for each solute atom type are counted. If there are multiple solute molecules of the same type, then the number of water molecules counted around an atom type are averaged over the total number of that type of solute molecule in the system. For example, the average number of water molecules around the hydroxyl oxygen atom over all ethanol solute molecules in a system is outputted. 

# Usage
An example for running ClusterAnalysis with bash scripting:

```
EXE="$HOME/bin/ClusterAnalysis/runClusterAnalysis.py" 
    # path to folder containing executable script

DIR="/path/to/files" # path to where topology and trajectory files are

TOPOLOGY="molecule.top" # topology file name
TRAJECTORY="molecule.mdcrd" # trajectory file name

START=0 # starting frame for analysis
END=100 # end frame for analysis
STEP=1 # step in-between frames to analyse

CUTOFF=3 # cutoff distance in Angstrom
SOLVENT="WAT" # name of solvent molecules as in topology file
EXCLUDE="Cl- Na+ DMS" # names of molecules to exclude from cluster analysis

## setup instructions for running code and execute, inclusion of -hb 
## means hydrogen bond analysis is included.
EXE2="${EXE} -top ${DIR}/${TOPOLOGY} 
        -traj ${DIR}/${TRAJECTORY} 
        -s ${START} -e ${END} -df ${STEP} 
        -r ${CUTOFF} 
        -solv ${SOLVENT} 
        -ex ${EXCLUDE} 
        -hb > file.log"
$(eval ${EXE2} )
```

# Output files
- `file.log`: Default log file to check progress of analysis that prints out if hydrogen bond analysis is used, number of frames in the trajectory analysis and an output of each frame that is processed for neighbour calculations. Also included is information about the number of atoms in the system, number of solute atoms, solute molecules, solvent molecules, molecules excluded from cluster analysis and the time taken to perform analysis. An example:

```
2021-09-02 10:39:22.800962
HB: True
n_frames: 23428
< Timestep 0 with unit cell dimensions [130.606   130.606   130.606   109.47122 109.47122 109.47122] >
< Timestep 10 with unit cell dimensions [118.895226 118.895226 118.895226 109.47122  109.47122  109.47122 ] >
< Timestep 20 with unit cell dimensions [118.647484 118.647484 118.647484 109.47122  109.47122  109.47122 ] >
.
.
.
< Timestep 22990 with unit cell dimensions [118.6785  118.6785  118.6785  109.47122 109.47122 109.47122] >

Input files: ['test_directory/noWATFnoION.Mol.top', 'test_directory/final_strip.nc'], 
N atoms: 396, 
N solute atoms: 396, 
N frames: 2300

residues in system: dict_keys(['*0'])
solvent: ['WAT']
excluded from cluster analysis: ['Cl-', 'Na+', 'DMS', 'WAT']
running cluster and HB analysis
0:10:23.500780
```

- `cluster-size-distribution.csv`: Output of the distribution of cluster sizes across all frames analysed, the first column is the number of atoms in a cluster and the second column is how often a cluster of that size was found in the frames analysed. An example:
```
# cluster_size,frequency
36.000,5115.000
72.000,980.000
108.000,607.000
144.000,1409.000
180.000,574.000
216.000,408.000
324.000,10.000
360.000,536.000
```

- `per-frame-clusters.csv`: Output of the number of clusters found per-frame analysed. The first column is the frame number analysed, the second column is the average number of atoms per cluster in a frame and the third column is the number of clusters in a given frame.
```
# frame,ave_cluster_size,n_clusters
0.000,39.600,10.000
1.000,39.600,10.000
2.000,39.600,10.000
.
.
.
2297.000,198.000,2.000
2298.000,198.000,2.000
2299.000,198.000,2.000
```

- `solute-water-contact-counts.csv`: If `hydrogen_bonding=True`, output of the average number of hydrogen-bonds between a solute atom (that is not in the excluded molecules list) and a water atom. A hydrogen-bond is defined if an atom on water is `max_cutoff` distance or less from a solute atom. The first column is the frame number analysed and the rest of the columns contain the average number of hydrogen-bonds with water per solute atom. Each header for the solute atoms is formatted as molecule name_atom name_number of molecules in the simulation.
```
frame,*0_O_11,*0_O1_11,*0_H2_11 ...
0,1.091	0.909,0.455
1,0.909	1.091,1.091
2,1.273	0.545,0.364
3,1.091	1.636,0.273
.
.
.
```

- `solute-solute-contact-counts.csv`: Output of the average number of contacts between a solute atom (that is not in the excluded molecules list) and any solute atom in a different molecule. A contact is defined if a solute atom in another molecule is `max_cutoff` distance or less from a solute atom. The first column is the frame number analysed and the rest of the columns contain the average number of solute-solute contacts per solute atom. Each header for the solute atoms is formatted as molecule name_atom name_number of molecules in the simulation.
```
frame,*0_O1_11,*0_H10_11,*0_H9_11 ...
0,0.091,0.091,0
1,0,0,0.273
2,0.091,0.091,0.091
3,0.091,0,0
.
.
.
```

- `clusters.pdb`: this pdb file contains the solute atoms only for the frames analysed, if opened in VMD, the `pdbbfactor.tcl` script can be used to colour molecules by what cluster they are in. Note the cluster number is given in the beta-factor column and cluster numbers have values between 0-1. An example:
```
ATOM      0  C   *0  A   1      66.863  66.669  63.916 0.100 0.100              
ATOM      1  C1  *0  A   1      66.545  68.089  63.994 0.100 0.100              
ATOM      2  C2  *0  A   1      66.121  68.724  65.137 0.100 0.100              
ATOM      3  O   *0  A   1      65.940  70.025  65.100 0.100 0.100                          
.
.
.
ATOM    392  H7  *0  A  11      72.799  58.695  60.141 0.500 0.500              
ATOM    393  H8  *0  A  11      70.436  59.109  60.358 0.500 0.500              
ATOM    394  H9  *0  A  11      69.996  55.787  63.110 0.500 0.500              
ATOM    395 H10  *0  A  11      72.499  55.421  62.847 0.500 0.500              
END  
```


