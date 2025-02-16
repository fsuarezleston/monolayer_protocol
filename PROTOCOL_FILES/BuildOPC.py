#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Fabián Suárez Lestón
"""

import argparse
import MDAnalysis
import numpy as np


parser = argparse.ArgumentParser(description =
    '''Water to OPC, v1.0''')
parser.add_argument( "-f", "--file", type= str, default = 'system.gro',
    help= ''' File with the system.\n
    Default: %(default)s ''' )
parser.add_argument( "-wr", "--water_resname", type = str, default = "SOL",
    help = ''' Name of the current water model residue.\n
    Default: %(default)s ''' )
parser.add_argument( "-wn", "--water_names", type = str, nargs='+', default = [ "OW", "HW1", "HW2" ],
    help = ''' Name of the atoms in the current water model.\n
    Default: %(default)s ''' )
parser.add_argument( "-or", "--opc_resname", type = str, default = "SOL",
    help = ''' Name of the OPC model residue.\n
    Default: %(default)s ''' )
parser.add_argument( "-on", "--opc_names", type = str, nargs='+', default = [ "OH2", "H1", "H2", "MW" ],
    help = ''' Name of the atoms in the OPC model.\n
    Default: %(default)s ''' )
parser.add_argument( "-ir", "--ions_resname", nargs='+', default = ["SOD", "CLA", "POT", "CAL"],
    help = ''' Name of the ions residues.\n
    Default: %(default)s ''' ) 
 
parser.add_argument( "-o", "--output", type = str, default="system_out.gro",
    help= ''' New system with OPC water model.\n
    Default: %(default)s ''' )
args = parser.parse_args()


def DefVS( water ) -> np.ndarray:
    '''
    Define the position of the virtual site in a 4-point model

    Parameters
    ----------
    water : MDAnalysis.core.AtomGroup
        A water molecule.

    Returns
    -------
    numpy.ndarray
        Coordinates of the virtual site.

    '''
    
    OW, H1, H2 = water.atoms
    
    vector = (H1.position - OW.position) + (H2.position - OW.position)
    vector /= np.linalg.norm( vector )
    
    return vector * 0.1594 + OW.position
    

if __name__ == "__main__":
    u = MDAnalysis.Universe( args.file )
    
    # Non-water residues
    NON_WATER = u.select_atoms( f"not resname {args.water_resname}" )
    
    # Water residues
    WATER = u.select_atoms( f"resname {args.water_resname}" ).residues
    
    OPC_pos = []
    for mol in WATER:
        # Position of the 4 points in the OPC water model
        OPC_pos.append( np.vstack( [ mol.atoms.positions, DefVS( mol ) ] ) )
    OPC_pos = np.vstack( OPC_pos )
    
    
    # Parameters of the new system
    NAtoms = len(NON_WATER.atoms) + len(OPC_pos)
    NResidues = len(u.residues)
    ResIndices = np.hstack( [ NON_WATER.resindices, [ i for i in np.unique( WATER.atoms.resindices ) for _ in args.opc_names ] ] )
    SegIndices = [0] * NResidues
    
    # Generate a new system
    U = MDAnalysis.Universe.empty( NAtoms,
                                   n_residues = NResidues,
                                   atom_resindex = ResIndices,
                                   residue_segindex = SegIndices,
                                   trajectory = True )
        
    # Add atributes
    U.add_TopologyAttr( 'name', np.hstack( [ NON_WATER.names, args.opc_names*len(WATER) ] ) )
    U.add_TopologyAttr( 'type', np.hstack( [ NON_WATER.types, ["O","H","H",""]*len(WATER) ] ) )
    U.add_TopologyAttr( 'resname', np.hstack( [ NON_WATER.residues.resnames, [ args.opc_resname ]*len(WATER) ]) )
    U.add_TopologyAttr( 'resids', list(range(1,NResidues+1) ) )
    U.atoms.positions = np.vstack( [ NON_WATER.positions, OPC_pos ] )
    U.dimensions = u.dimensions

    # Save the system
    U.atoms.write( args.output )
    

    with MDAnalysis.selections.gromacs.SelectionWriter( 'index.ndx', mode = 'w' ) as ndx:
        ndx.write( U.atoms,
                  name = 'SYSTEM' )
        ndx.write( U.select_atoms( f"resname {args.opc_resname} {' '.join(args.ions_resname)}" ),
                  name = 'SOLV' )
        ndx.write( U.select_atoms( f"not resname {args.opc_resname} {' '.join(args.ions_resname)}" ),
                  name = 'MEMB' )
