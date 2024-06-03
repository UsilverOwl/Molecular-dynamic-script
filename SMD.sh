#!/bin/bash
#PDB FILE, will be replace in the code later by $PDB
FILE="DR4t" #<<<<<<<<<<<<<<<<<<<<<<<<< PUT THE PDB NAME HERE (without the extension)

#---------  SIMU SETUP  -----------
BOXSIZE=1.5 #cubic simulation boxsize
BOXTYPE=cubic #Box type
NT=60 #Number of core.
WATER=tip3p #Water type
NUMBEROFREPLICAS=3 #Number of replica
FF=charmm36 #Force field
#Simulation time in nanosec. Will be converted in fs and modified in the mdp file.
minim=minim
nvt=nvt
npt=npt
md=md

#---------  HPC SETUP  -----------
MPI="mpirun -np 16" #If you have to submit jobs with MPI softwares like "mpirun -np 10". Add the command here
GMX=gmx #GMX command (can be "$GMX_mpi" sometimes. Just change it here
#THOSE COMMANDS
GPU0="-gpu_id 0,1 -nt 60  -ntmpi 2   -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 1 "
GPU1="-gpu_id 1 -nt  20  -ntmpi  -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 1"
export GMX_GPU_PME_PP_COMMS=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_GPU_DD_COMMS=1
MDRUN_CPU="$GMX mdrun -nt $NT"
MDRUN_GPU="$GMX mdrun $GPU0"


# Create replica directories
#for i in {1..3}
#do
#  replica_dir="replica$i"
#  mkdir -p $replica_dir/graph $replica_dir/gro $replica_dir/results
#done

# Run GROMACS for each replica
#for i in {1..3}
#do
#  replica_dir="replica$i"
#  $gmx_cmd -v -deffnm $replica_dir -nt 4 &> $replica_dir/gmx.log &



#for ((i=0; i<$NUMBEROFREPLICAS; i++))
#	do
#	echo ">>>>> replica_"$((i+1))
#	mkdir "replica_"$((i+1))
#	cd "replica_"$((i+1))
#	cp -R ../mdp .
#	cp ../$PDB"_solv_ions.gro" .
#	cp ../topol.top .
#	cp ../*.itp .
  #  cp ../index.ndx . 2> /dev/null


# Create replica directories
#for i in {1..3}
#do
#  replica_dir="replica$i"
#  mkdir -p $replica_dir/graph $replica_dir/gro $replica_dir/results
#done

# Create results subdirectories
#for replica_dir in replica{1..3}
#do
#  results_dir="$replica_dir/results"
#  mkdir -p $results_dir/mini $results_dir/nvt $results_dir/npt $results_dir/prod
#done


#set "PDB" name (all the simulation filenames are based on this variable).
	PDB=$FILE
	########################
	##   TOPOLOGIE Creation
	#######################
	echo "6 6" | $GMX pdb2gmx -f $PDB".pdb" -o $PDB"_processed.gro" -water $WATER  -ignh -ff $FF -ter
	########################
	##   Box Creation
	#######################
	#Default editconf, changer if you want...
	$GMX editconf -f  $PDB"_processed.gro" -o $PDB"_newbox.gro" -d $BOXSIZE -bt $BOXTYPE


######################
###SOLVATATION
#####################
$GMX solvate -cp $PDB"_newbox.gro" -cs tip3p.gro -o $PDB"_solv.gro" -p topol.top

cp topol.top topol2.top

#######################
## ADDING IONS
#######################
$GMX grompp -f mdp/ions.mdp -c $PDB"_solv.gro" -p topol.top -o ions.tpr --maxwarn 1

echo "7" | $GMX genion -s ions.tpr -o $PDB"_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral

# Run GROMACS for each replica
	for ((i=0; i<$NUMBEROFREPLICAS; i++)); do
		echo ">>>>> replica_"$((i+1))
		mkdir -p "replica_"$((i+1))/dmso/LF
		mkdir -p "replica_"$((i+1))/gro
		cp $PDB"_solv.gro" -t "replica_"$((i+1))/dmso
		cp dmso.pdb -t "replica_"$((i+1))/dmso
		cp dmso.itp -t "replica_"$((i+1))/dmso
		cp -r charmm36.ff/ -t "replica_"$((i+1))/dmso
				cp -r charmm36ff/ -t "replica_"$((i+1))
		cp topol2.top -t "replica_"$((i+1))/dmso
		cp -r mdp2 -t "replica_"$((i+1))/dmso

		cp -r charmm36.ff -t "replica_"$((i+1))
		cd "replica_"$((i+1))
		export GMXLIB=$GMXLIB:~MD/"replica_"$((i+1))/charmm36.ff
		cp -R ../mdp .
		cp ../$PDB"_solv_ions.gro" .
		cp ../topol.top .
		cp ../*.itp .
		cp ../index.ndx . 2> /dev/null
		mkdir -p results/mini
		mkdir graph
		mv DR4t_solv_ions.gro -t gro/


#######################
	## MINIMIZATION
	#######################
		$GMX grompp -f mdp/$minim.mdp -c gro/$PDB"_solv_ions.gro" -p topol.top -o em.tpr $INDEX
		$GMX mdrun -v -deffnm em  -c  -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2

####Cleaning


	mv em* -t results/mini/
	mv mdout.mdp -t results/mini/
	mv results/mini/*.gro -t gro/

	#potential energy graph
		echo "11 0" | $GMX energy -f results/mini/em.edr -o mini_"$PDB"_pot.xvg

		mv mini_"$PDB"_pot.xvg -t graph/
##################
###EQUILIBRATION
####################

###TEMPERATURE###
		$GMX grompp -f mdp/nvt.mdp -c gro/em.gro -r gro/em.gro -p topol.top -o nvt.tpr
		$GMX mdrun -v -deffnm nvt -c -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2

	#cleaning


		mkdir -p results/nvt
		mv nvt* results/nvt/ 2> /dev/null
		mv mdout.mdp results/nvt/
		mv results/nvt/*.gro -t gro/
		#temperature_graph
		echo "16 0" | $GMX energy -f results/nvt/nvt.edr -o temperature.xvg


###PRESSURE###
		$GMX grompp -f mdp/npt.mdp -c gro/em.gro -r gro/em.gro -t results/nvt/nvt.cpt -p topol.top -o npt.tpr

			$GMX mdrun -v -deffnm npt  -c -gpu_id 0 #-pin on -pinstride 1 -nt $NT

		#cleaning
		mkdir -p results/npt
		mv npt* results/npt/ 2> /dev/null
		mv mdout.mdp results/npt/
		mv results/npt/*.gro -t gro/

#Pression and density graph

		echo "17 0" | $GMX energy -f results/npt/npt.edr -o pressure.xvg
		echo "22 0" | $GMX energy -f results/npt/npt.edr -o density.xvg

mv *.xvg -t graph/
#######################
	## Production
#######################


		$GMX grompp -f mdp/md.mdp -c gro/npt.gro -t results/npt/npt.cpt -p topol.top -o "md_"$PDB"_prod.tpr"  $INDEX
		$MDRUN_GPU -deffnm "md_"$PDB"_prod"  -v -c

		mv *xvg -t graph/


#################################################################################################################################################################################################
		echo "1 1" | $GMX rms -s "results/prod/md_"$PDB"_prod.tpr" -f "results/prod/md_"$PDB"_prod.trr" -o graph/rmsd.xvg -tu ns
		echo "1 0" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_DR4t_prod.xtc" -o "md_"$PDB"_LF1.pdb" -pbc mol  -center 1 -dump -1
		echo "1 0" |$GMX trjconv -pbc mol -f md_DR4t_prod.xtc -o nopbc.xtc -center 1 -s md_DR4t_prod.tpr

		cd results/
		mkdir prod
		cd prod/
		echo "1 0" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_prod.xtc" -o "md_"$PDB"_clean_temp.xtc" -pbc mol -center 1

		echo "Other System" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_prod.trr" -o "md_"$PDB"_clean_temp.xtc" -pbc nojump -ur compact -center

		echo "Other System" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_full.xtc" -fit rot+trans

		echo "Other non-Water" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_nowat.xtc" -fit rot+trans

		echo "Other non-Water" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_clean_nowat.pdb" -pbc nojump -ur compact -center -b 0 -e 0

 	#extract last frame in PDB only for the protein.
		echo "1 1" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_"$PDB"_clean_temp.xtc" -o "md_"$PDB"_LF2.pdb" -pbc mol  -center 1 -dump -1

		rm "md_"$PDB"_clean_temp.pdb"

		echo "non-Water" | $GMX convert-tpr -s "md_"$PDB"_prod.tpr" -o tpr_nowat.tpr

	# Create a smooth trajectory
		echo "1" | $GMX filter -s tpr_nowat.tpr -f "md_"$PDB"_clean_nowat.xtc" -ol "md_"$PDB"_clean_nowat_filtered.xtc" -all -fit


cd ..
cd ..
#################################################################################################################################################################################################
##DMSO
#################################################################################################################################################################################################
cp -r charmm36ff/ -t "replica_"$((i+1))/dmso/
cp -r mdp/ -t "replica_"$((i+1))/dmso/
cp -r dmsip.itp/ -t "replica_"$((i+1))/dmso/
cp -r dmso.prm/ -t "replica_"$((i+1))/dmso/

cd "replica_"$((i+1))/dmso
export GMXLIB=$GMXLIB:~/MD/replica_1/dmso/charmm36.ff

export GMXLIB=$GMXLIB:~/MD/replica_1/dmso/dmso/charmm36.ff
	# Set the topology file and the DMSO PDB file



# Calculate the number of DMSO molecules to add
n_water=$(grep -oP '(?<=SOL\s).*' "$topol_file" | awk '{sum+=$1} END {print sum}')
ndmso=$(printf "%d" $(echo "scale=0; $n_water * 0.006780473358264" | bc))


# Add the DMSO molecules using GROMACS
		echo "7" | $GMX insert-molecules -f DR4t_solv.gro -ci dmso.pdb -o dmso_5%.gro -nmol $ndmso -replace -try 100
# Add the DMSO molecules using GROMACS
		echo "7" | $GMX insert-molecules -f md_DR4t_LF1.pdb -ci dmso.pdb -o dmso_5%_LF1.gro -nmol $ndmso -replace -try 100



# Define the line to add
LINE="#include \"dmso.itp\""
LINE="#include \"dmso.prm\""


		mv *LF* -t dmso/LF
	##################################################################################################

	echo "2 3" | $GMX pdb2gmx -f dmso_5%.gro -o dmso"_processed.gro" -water $WATER -ignh -ff $FF -ter



# Open the file in read-only mode
filename="topol.top"
tempfile=$(mktemp)

# Read the file line by line
while IFS= read -r line; do
  # Check if the line matches the target line
  if [ "$line" = "LINE=\"#include \\\"dmso.itp\\\"\"" ]; then
    # Add the new lines after the target line
    echo "$line" >> $tempfile
    echo "; additional params for the molecule" >> $tempfile
    echo "#include \"dmso.prm\"" >> $tempfile
  else
    # Just copy the original line
    echo "$line" >> $tempfile
  fi
done < "$filename"

# Replace the original file with the modified one
mv $tempfile $filename
#########################################################################################
# Read the topol.top file line by line
	while IFS= read -r line
	do
	# Check if the line contains "#include "topol_RNA.itp""
		if [[ $line == *"#include "topol_RNA.itp""* ]]
		then
		# If it does, replace it with "#include "topol_dmso.itp""
		line="#include \"topol_dmso.itp\""
		fi
	# Print the modified line
	echo "$line"

	# Check if the line contains "#include "topol_RNA.itp""
		if [[ $line == *"#include \charmm36.ff\ions.itp"* ]]
		then
		# If it does, replace it with "#include "topol_dmso.itp""
		line="#include dmso.itp"
		fi
	# Print the modified line
	echo "$line"


	# Check if the line contains "#include "topol_RNA.itp""
		if [[ $line == *"RNA                 1"* ]]
		then
		# If it does, replace it with "#include "topol_dmso.itp""
		line="DMSO              $ndmso"
		fi
	# Print the modified line
	echo "$line"

	done < "topol.top" > "modified_topol.top"

		mkdir -p results/mini
		mkdir graph


	#######################
	## MINIMIZATION
	#######################
		$GMX grompp -f mdp2/$minim.mdp -c dmso"_processed.gro" -p topol.top -o em.tpr $INDEX
		$GMX mdrun -v -deffnm em  -c  -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2


####Cleaning


	mv em* -t results/mini/
	mv mdout.mdp -t results/mini/
	mv results/mini/*.gro -t gro/

	#potential energy graph
		echo "11 0" | $GMX energy -f results/mini/em.edr -o mini_"$PDB"_pot.xvg


#######################
## ADDING IONS
#######################
		$GMX grompp -f mdp2/ions.mdp -c dmso_%5.gro -p topol.top -o ions.tpr --maxwarn 1

		echo "7" | $GMX genion -s ions.tpr -o $PDB"_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral


#######################
	## MINIMIZATION
	#######################
		$GMX grompp -f mdp2/$minim.mdp -c gro/$PDB"_solv_ions.gro" -p topol.top -o em_dmso.tpr $INDEX
		$GMX mdrun -v -deffnm em_dmso -gpu_id 0,1 -c -nt $NT -pin on -pinstride 1

####Cleaning


	mv em* -t results/mini/
	mv mdout.mdp -t results/mini/
	mv results/mini/*.gro -t gro/

	#potential energy graph
		echo "11 0" | $GMX energy -f results/mini/em.edr -o mini_"$PDB"_pot.xvg

##################
###EQUILIBRATION
####################

###TEMPERATURE###
		$GMX grompp -f mdp2/nvt.mdp -c em_dmso.gro -r em_dmso.gro -p topol.top -o nvt_dmso.tpr
		$GMX mdrun -v -deffnm nvt_dmso -c -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2

	#cleaning


		mkdir -p results/nvt
		mv nvt* results/nvt/ 2> /dev/null
		mv mdout.mdp results/nvt/
		mv results/nvt/*.gro -t gro/
		#temperature_graph
		echo "16 0" | $GMX energy -f results/nvt/nvt.edr -o temperature.xvg


###PRESSURE###
		$GMX grompp -f mdp2/npt.mdp -c nvt_dmso.gro -r nvt_dmso.gro -t nvt_dmso.cpt -p topol.top -o npt_dmso.tpr

			$GMX mdrun -v -deffnm npt  -c -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2

		#cleaning
		mkdir -p results/npt
		mv npt* results/npt/ 2> /dev/null
		mv mdout.mdp results/npt/
		mv results/npt/*.gro -t gro/

#Pression and density graph

		echo "17 0" | $GMX energy -f results/npt/npt.edr -o pressure.xvg
		echo "22 0" | $GMX energy -f results/npt/npt.edr -o density.xvg

mv *.xvg -t graph/
#######################
	## Production
#######################


		$GMX grompp -f mdp2/md.mdp -c npt_dmso.gro -t npt_dmso.cpt -p topol.top -o "md_dmso_prod.tpr"  $INDEX
		$MDRUN_GPU -deffnm "md_dmso_prod"  -v -c

		mv *xvg -t graph/


#################################################################################################################################################################################################
		echo "1 1" | $GMX rms -s "results/prod/md_"$PDB"_prod.tpr" -f "results/prod/md_"$PDB"_prod.trr" -o graph/rmsd.xvg -tu ns
		echo "1 0" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_DR4t_prod.xtc" -o "md_"$PDB"_LF1.pdb" -pbc mol  -center 1 -dump -1
		echo "1 0" |$GMX trjconv -pbc mol -f md_DR4t_prod.xtc -o nopbc.xtc -center 1 -s md_0_1.tpr

		cd results/prod

		echo "1 0" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_prod.xtc" -o "md_dmso_clean_temp.xtc" -pbc mol -center 1

		echo "Other System" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_prod.trr" -o "md_dmso_clean_temp.xtc" -pbc nojump -ur compact -center

		echo "Other System" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_clean_full.xtc" -fit rot+trans

		echo "Other non-Water" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_clean_nowat.xtc" -fit rot+trans

		echo "Other non-Water" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_clean_nowat.pdb" -pbc nojump -ur compact -center -b 0 -e 0

 	#extract last frame in PDB only for the protein.
		echo "1 1" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_LF2.pdb" -pbc mol  -center 1 -dump -1

		rm "md_dmso_clean_temp.pdb"

		echo "non-Water" | $GMX convert-tpr -s "md_dmso_prod.tpr" -o tpr_nowat.tpr

	# Create a smooth trajectory
		echo "1" | $GMX filter -s tpr_nowat.tpr -f "md_dmso_clean_nowat.xtc" -ol "md_dmso_clean_nowat_filtered.xtc" -all -fit


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

scp -r mdp/ -t LF/
cp topol.top -t LF/

cp -r charmm36.ff/ -t LF/
cp -r dmsip.itp/ -t dmso/LF
cp -r dmso.prm/ -t dmso/LF
cd LF/
#################################################################################################################################################################################################







	echo "2 3" | $GMX pdb2gmx -f dmso_5%_LF1.gro -o dmso"_processed.gro" -water $WATER -ignh -ff $FF -ter



# Open the file in read-only mode
filename="topol.top"
tempfile=$(mktemp)

# Read the file line by line
while IFS= read -r line; do
  # Check if the line matches the target line
  if [ "$line" = "LINE=\"#include \\\"dmso.itp\\\"\"" ]; then
    # Add the new lines after the target line
    echo "$line" >> $tempfile
    echo "; additional params for the molecule" >> $tempfile
    echo "#include \"dmso.prm\"" >> $tempfile
  else
    # Just copy the original line
    echo "$line" >> $tempfile
  fi
done < "$filename"

# Replace the original file with the modified one
mv $tempfile $filename
#########################################################################################
# Read the topol.top file line by line
	while IFS= read -r line
	do
	# Check if the line contains "#include "topol_RNA.itp""
		if [[ $line == *"#include "topol_RNA.itp""* ]]
		then
		# If it does, replace it with "#include "topol_dmso.itp""
		line="#include \"topol_dmso.itp\""
		fi
	# Print the modified line
	echo "$line"

	# Check if the line contains "#include "topol_RNA.itp""
		if [[ $line == *"#include \charmm36.ff\ions.itp"* ]]
		then
		# If it does, replace it with "#include "topol_dmso.itp""
		line="#include \charmm36.ff\dmso.itp"
		fi
	# Print the modified line
	echo "$line"


	# Check if the line contains "#include "topol_RNA.itp""
		if [[ $line == *"RNA                 1"* ]]
		then
		# If it does, replace it with "#include "topol_dmso.itp""
		line="DMSO              $ndmso"
		fi
	# Print the modified line
	echo "$line"

	done < "topol.top" > "modified_topol.top"

		mkdir -p results/mini
		mkdir graph


	#######################
	## MINIMIZATION
	#######################
		$GMX grompp -f mdp2/$minim.mdp -c dmso"_processed.gro" -p topol.top -o em.tpr $INDEX
		$GMX mdrun -v -deffnm em  -c  -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2


####Cleaning


	mv em* -t results/mini/
	mv mdout.mdp -t results/mini/
	mv results/mini/*.gro -t gro/

	#potential energy graph
		echo "11 0" | $GMX energy -f results/mini/em.edr -o mini_"$PDB"_pot.xvg


#######################
## ADDING IONS
#######################
		$GMX grompp -f mdp2/ions.mdp -c dmso_%5.gro -p topol.top -o ions.tpr --maxwarn 1

		echo "7" | $GMX genion -s ions.tpr -o $PDB"_solv_ions.gro" -p topol.top -pname NA -nname CL -neutral


#######################
	## MINIMIZATION
	#######################
		$GMX grompp -f mdp2/$minim.mdp -c gro/$PDB"_solv_ions.gro" -p topol.top -o em_dmso.tpr $INDEX
		$GMX mdrun -v -deffnm em_dmso -gpu_id 0,1 -c -nt $NT -pin on -pinstride 1

####Cleaning


	mv em* -t results/mini/
	mv mdout.mdp -t results/mini/
	mv results/mini/*.gro -t gro/

	#potential energy graph
		echo "11 0" | $GMX energy -f results/mini/em.edr -o mini_"$PDB"_pot.xvg

##################
###EQUILIBRATION
####################

###TEMPERATURE###
		$GMX grompp -f mdp2/nvt.mdp -c em_dmso.gro -r em_dmso.gro -p topol.top -o nvt_dmso.tpr
		$GMX mdrun -v -deffnm nvt_dmso -c -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2

	#cleaning


		mkdir -p results/nvt
		mv nvt* results/nvt/ 2> /dev/null
		mv mdout.mdp results/nvt/
		mv results/nvt/*.gro -t gro/
		#temperature_graph
		echo "16 0" | $GMX energy -f results/nvt/nvt.edr -o temperature.xvg


###PRESSURE###
		$GMX grompp -f mdp2/npt.mdp -c nvt_dmso.gro -r nvt_dmso.gro -t nvt_dmso.cpt -p topol.top -o npt_dmso.tpr

			$GMX mdrun -v -deffnm npt  -c -gpu_id 0 -pin on -pinstride 1 -nt 60   -ntmpi 2

		#cleaning
		mkdir -p results/npt
		mv npt* results/npt/ 2> /dev/null
		mv mdout.mdp results/npt/
		mv results/npt/*.gro -t gro/

#Pression and density graph

		echo "17 0" | $GMX energy -f results/npt/npt.edr -o pressure.xvg
		echo "22 0" | $GMX energy -f results/npt/npt.edr -o density.xvg

mv *.xvg -t graph/
#######################
	## Production
#######################


		$GMX grompp -f mdp2/md.mdp -c npt_dmso.gro -t npt_dmso.cpt -p topol.top -o "md_dmso_prod.tpr"  $INDEX
		$MDRUN_GPU -deffnm "md_dmso_prod"  -v -c

		mv *xvg -t graph/


#################################################################################################################################################################################################
		echo "1 1" | $GMX rms -s "results/prod/md_"$PDB"_prod.tpr" -f "results/prod/md_"$PDB"_prod.trr" -o graph/rmsd.xvg -tu ns
		echo "1 0" | $GMX trjconv -s "md_"$PDB"_prod.tpr" -f "md_DR4t_prod.xtc" -o "md_"$PDB"_LF1.pdb" -pbc mol  -center 1 -dump -1
		echo "1 0" |$GMX trjconv -pbc mol -f md_DR4t_prod.xtc -o nopbc.xtc -center 1 -s md_0_1.tpr

		cd results/prod

		echo "1 0" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_prod.xtc" -o "md_dmso_clean_temp.xtc" -pbc mol -center 1

		echo "Other System" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_prod.trr" -o "md_dmso_clean_temp.xtc" -pbc nojump -ur compact -center

		echo "Other System" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_clean_full.xtc" -fit rot+trans

		echo "Other non-Water" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_clean_nowat.xtc" -fit rot+trans

		echo "Other non-Water" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_clean_nowat.pdb" -pbc nojump -ur compact -center -b 0 -e 0

 	#extract last frame in PDB only for the protein.
		echo "1 1" | $GMX trjconv -s "md_dmso_prod.tpr" -f "md_dmso_clean_temp.xtc" -o "md_dmso_LF2.pdb" -pbc mol  -center 1 -dump -1

		rm "md_dmso_clean_temp.pdb"

		echo "non-Water" | $GMX convert-tpr -s "md_dmso_prod.tpr" -o tpr_nowat.tpr

	# Create a smooth trajectory
		echo "1" | $GMX filter -s tpr_nowat.tpr -f "md_dmso_clean_nowat.xtc" -ol "md_dmso_clean_nowat_filtered.xtc" -all -fit



cd ..

		done
	done
