#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Fabián Suárez Lestón
"""

import argparse
import MDAnalysis

parser = argparse.ArgumentParser(description =
    '''Prepare PLUMED, v1.0''')
parser.add_argument( "-f", "--file", type= str, default = 'step7_1.gro',
    help= ''' File with the system.\n
    Default: %(default)s ''' )
parser.add_argument( "-lr", "--lipid_resname", nargs='+',
    help = ''' Name of the lipids residues.\n
    Default: %(default)s ''' )
parser.add_argument( "-c", "--comp_max", type=float, default = 40,
    help= ''' Maximum compression APL (in Å^2).\n
    Default: %(default)s ''' )
parser.add_argument( "-e", "--exp_max", type=float, default = 140,
    help= ''' Maximum expansion APL (in Å^2).\n
    Default: %(default)s ''' )

args = parser.parse_args()

# Load the equilibrated monolayer system
u = MDAnalysis.Universe( args.file, args.file )

# The number of lipids per leaflet
NLeaflet =  int( len( u.select_atoms( f'resname {" ".join( args.lipid_resname )}' ).residues ) / 2 )

# The APL in the final state
apl = u.dimensions[0] * u.dimensions[1] / NLeaflet

# Generate a plumed_XXX.dat file for the compression and the expansion
for proc, end in [ ("com", args.comp_max), ("exp", args.exp_max) ]:
    with open( f"plumed_{proc}.dat", "w+" ) as FILE:
        FILE.write(  "cell: CELL\n" )
        FILE.write(  "x2: COMBINE ARG=cell.ax,cell.ay,cell.az POWERS=2,2,2 PERIODIC=NO\n" )
        FILE.write(  "y2: COMBINE ARG=cell.ax,cell.ay,cell.az POWERS=2,2,2 PERIODIC=NO\n" )
        FILE.write( f"area: MATHEVAL ARG=x2,y2 FUNC=(x^0.5*y^0.5)/{NLeaflet:d}*100 PERIODIC=NO\n" )
        FILE.write(  "MOVINGRESTRAINT ...\n" )
        FILE.write(  " ARG=area\n" )
        FILE.write( f" STEP0=0        AT0={apl:.2f} KAPPA0=100\n" )
        FILE.write( f" STEP1=1000000  AT1={end}\n" )
        FILE.write( f" STEP2=1100000  AT2={end}\n" )
        FILE.write(  "...\n" )
        FILE.write( f"PRINT ARG=area STRIDE=3000 FILE=area_{proc}.dat\n" )
