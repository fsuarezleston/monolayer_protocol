#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Fabián Suárez Lestón
"""

import os
import re
import argparse
import MDAnalysis
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(description =
    '''PullAndPush, v1.0''')
parser.add_argument( "-f", "--folder", type=str,
    help= ''' Folder with the system.''' )
parser.add_argument( "-a", "--apl_range", type = float, nargs='+', default = [ 45, 105, 5 ],
    help= '''Initial, final and increment of area per lipid in the isotherm.
    Default: %(default)s''')
args = parser.parse_args()


if __name__ == "__main__":
    # Change working directory
    os.chdir( args.folder )
    
    # The target values of APL
    target_apl = np.arange(*args.apl_range)
    
    # Files with the values of the APL along the expansion and compression processes
    data_c = pd.read_csv( "area_com.dat", header=None, skiprows=1, sep="\s+" )
    data_e = pd.read_csv( "area_exp.dat", header=None, skiprows=1, sep="\s+" )
    
    # Read the topol.top file
    with open( "topol.top", "r+" ) as file:
        top = file.readlines()
    
    # Iterate over the desired APLs
    for apl in target_apl:
        # Difference between scaned and target APL in the compression and expansion processes
        dif_c = abs( data_c[1] - apl )
        dif_e = abs( data_e[1] - apl )
        
        # Find the time at which the APL was closest to the target and the
        # process by which this APL was achieved
        if min( dif_c ) < min( dif_e ):
            time = data_c[0].iloc[ dif_c.argmin() ]
            proc = "com"
        else:
            time = data_e[0].iloc[ dif_e.argmin() ]
            proc = "exp"
            
        # Create a folder for the point in the isotherm
        try:
            os.mkdir( f"APL_{apl}" )
        except:
            pass
        
        # Extract the state at the target APL
        os.system( f"echo 0 | gmx trjconv -s NPT_{proc}.tpr -f NPT_{proc}.trr -vel -o APL_{apl}/INITIAL.trr -b {time} -dump {time}" )
        os.system( f"echo 0 | gmx trjconv -s NPT_{proc}.tpr -f APL_{apl}/INITIAL.trr -vel -o APL_{apl}/INITIAL.pdb" )
        
        
        # Load the system
        u = MDAnalysis.Universe( f"APL_{apl}/INITIAL.pdb" )
        
        # The non-water resiudes
        system = u.select_atoms("not resname SOL TIP3 OPC W")
        
        # Increase the Z dimension
        Zsys = 1.25 * ( max( system.positions[:,2] ) - min( system.positions[:,2] ) )
        
        # Select all the atoms whithin the limits of the non-water resiudes
        # This is done to remove water molecules which may have been adsorbed into the vacuum region
        atoms = u.select_atoms( f"same resid as ( prop z <= {max(system.positions[:,2])} and prop z >= {min(system.positions[:,2])} )" )
        
        # Center in the box
        atoms.positions -= np.array( [0, 0, atoms.center_of_geometry()[2] - Zsys/2] )
        
        # Redefine the dimensions of the system
        dims = np.copy(u.dimensions)
        dims[2] = Zsys
        u.dimensions = dims
        
        # Save the changes
        atoms.write( f"APL_{apl}/INITIAL.pdb" )
        
        # The system composition
        comp = { res: len( atoms.select_atoms( f"resname {res}" ).residues )for res in np.unique( atoms.resnames ) }
        
        # Generate a new topology file with the correct composition
        with open( f"APL_{apl}/topol.top", "w+" ) as newtop:
            for line in top:
                try:
                    if line.split()[0] in comp.keys():
                        # Write the residue and its abundance
                        newtop.write( f"{line.split()[0]:<7s}    {comp[line.split()[0]]}\n" )
                    
                    elif "toppar" in line:
                        # Modify relative path to toppar folder
                        newtop.write( re.sub( "toppar", "../toppar", line ) ) 
                        
                    else:
                        newtop.write( line )
                except:
                    newtop.write( line )
        
        # Remove checkpoints
        try:
            os.system( f"rm APL_{apl}/*#" )
        except:
            pass
    
        # Generate a new index file
        u = MDAnalysis.Universe( f"APL_{apl}/INITIAL.pdb" )
        with MDAnalysis.selections.gromacs.SelectionWriter( f'APL_{apl}/index.ndx', mode = 'w' ) as ndx:
            ndx.write( u.atoms,
                      name = 'SYSTEM' )
            ndx.write( u.select_atoms("resname SOL TIP3 OPC W ION CLA CL Cl- NA SOD POT Na+"),
                      name = 'SOLV' )
            ndx.write( u.select_atoms("not resname SOL TIP3 OPC W ION CLA CL Cl- NA SOD POT Na+"),
                      name = 'MEMB' )
        
        # Generate a TPR file for the NVT simulation
        os.system( f"gmx grompp -f NVT.mdp -o APL_{apl}/NVT.tpr -c APL_{apl}/INITIAL.pdb -p APL_{apl}/topol.top -n APL_{apl}/index.ndx " )
