# Monolayer Simulation Protocol
A protocol for the generation of π-A isotherms of lipid monolayers in GROMACS

It is designed for the use of monolayer systems generated with _CHARMM-GUI_, 
specifically in CHARMM36. If other method for the assembly of the monolayers 
or other forcefield is used modifications are required.

## Steps

1) Generate the monolayer system using the [Monolayer Builder tool of CHARMM-GUI](https://www.charmm-gui.org/?doc=input/membrane.monolayer).
   A reasonable number of water molecules per lipid must be provided (see comment [^1]).

2) Copy the folder with the GROMACS inputs into the monolayer protocol
   folder `./monolayer_protocol/PROTOCOL_FILES/MONOLAYER_FOLDER`

3) Execute the `ToOPC.sh` script to substitute the TIP3 water model by
   the 4-point OPC model:
   
   ```
   $ bash PROTOCOL_FILES/ToOPC.sh MONOLAYER_FOLDER
   ```
   
   Or, if you have your own method to convert 3-point into 4-point water
   molecules, you can use it instead.

5) Follow CHARMM-GUI's minimization and equilibration protocol as usual.
   The CHARMM-GUI's README file has been modified in the previous step
   to execute only one production simulation. You can execute the whole
   protocol as
   
   ```
   $ cd MONOLAYER_FOLDER
   $ csh README
   ```
   
   making sure that the commands are correct for the execution in 
   your machine (see comment [^2]).

7) Prepare the `plumed_XXX.dat` files for the expansion (exp) and compression 
   (com) processess.

   ```
   $ python ../PROTOCOL_FILES/PreparePLUMED.py -f step7_1.gro -lr LIPIDS -c MIN_APL -e MAX_APL
   ```

   `LIPIDS` refers to the list of lipid resiudes separated by spaces, e.g. `POPC POPE DPPC`; while
   `MIN_APL` and `MAX_APL` are the target APL values of the compression and
   expansion processes, respectively. Please, note that the expansion and
   compression times are preset, and match the time defined in `NPT.mdp`. If
   you want to change the expansion or compression velocities, you must modify
   accodordingly these files.

8) Execute the expansion and compression simulations, starting from the 
   final stage of the previous equilibration protocol.

   ```
   $ gmx grompp -f ../PROTOCOL_FILES/NPT.mdp -o NPT_exp.tpr -c step7_1.gro -t step7_1.cpt -p topol.top -n index.ndx
   $ gmx mdrun -v -deffnm NPT_exp -plumed plumed_exp.dat
   $ gmx grompp -f ../PROTOCOL_FILES/NPT.mdp -o NPT_com.tpr -c step7_1.gro -t step7_1.cpt -p topol.top -n index.ndx
   $ gmx mdrun -v -deffnm NPT_com -plumed plumed_com.dat
   ```

   For moderate size systems (~100 lipids and ~80 water molecules per lipid) these processes, if
   executed for the default time, are not expensive, and can be executed in an HPC within an hour.

9) Execute the `PullAndPush.py` script to generate the initial structures of
   the points in the isotherm.

   ```
   $ python ../PROTOCOL_FILES/PullAndPush.py -f . -a INITIAL_APL FINAL_APL INCREMENT_APL
   ```

   This will create a series of `APL_XXX`, where `XXX` is the value of the ApL of the system ranging
   from `INITIAL_APL` to `FINAL_APL` (not included), in `INCREMENT_APL` increments. Each of them contains
   a `NVT.tpr` which corresponds to a NVT simulation at its value of ApL. The simulations must be carried out
   at a HPC facility.

   Feel free to regulate the simulation by properly modifying `NVT.mdp` provided by us prior
   to the execution of this step.


## Useful information

**Example of CHARMM-GUI's README modification**
```
#!/bin/csh
#SBATCH --time=0-02:30:00  
#SBATCH --nodes=1  
#SBATCH --tasks-per-node=2 --cpus-per-task=4  
#SBATCH --gpus-per-node=1

# ...
# ...

srun --ntasks=1 gmx_mpi grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr \  
       -c ${init}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx -maxwarn 0  
srun gmx_mpi mdrun -v -deffnm ${mini_prefix}

# ...
# ...
```

## Comments

[^1]:  Number of water molecules per lipid should be provided with the following consideration: During the expansion your water layer will be contracted in the same extent as the area is extended. Let's assume you want to scan from 50 Å2 to 120 Å2, and with CHARMM-GUI, you create a system with ApL = 60 Å2. Then, as at the maximum expansion is twice the initial area, the water layer will be twice as thin, so if you think a water layer 4 nm is enough, then you must generate an initial 8 nm layer.

[^2]: You can create a slurm-script out of README csh script, but then you need to add `SBATCH` definitions, load the environment, and then run all `gmx` using `srun`. Please, see the example.
