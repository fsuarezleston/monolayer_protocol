#!/bin/bash

# Copy files in new folder
echo "Copying files"
cp PROTOCOL_FILES/CHARMMGUI_README ${1}/README
cp PROTOCOL_FILES/NVT.mdp ${1}
cp PROTOCOL_FILES/toppar/* ${1}/toppar/
cd ${1}/

# Add correct ITPs
echo "Modifying topology"
sed -i"" 's/\#include "toppar\/forcefield.itp"/\#include "toppar\/forcefield.itp"\n\#include "toppar\/opc_types.itp"/' topol.top
sed -i"" 's/\#include "toppar\/TIP3.itp"/\#include "toppar\/opc4.itp"/' topol.top
sed -i"" 's/TIP3/SOL/' topol.top

# Change water model to OPC
echo "Changing water model to OPC"
mv step5_input.gro ORIGINALstep5_input.gro
python ../PROTOCOL_FILES/BuildOPC.py -f ORIGINALstep5_input.gro -wr TIP3 -wn OH2 H1 H2 -o step5_input.gro
cd &>/dev/null

echo "Done!"
