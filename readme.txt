 ____       _          _____ 
/ ___|  ___| | ___  __|_   _|  
\___ \ / _ \ |/ _ \/ __|| |   
 ___) |  __/ |  __/ (__ | |   
|____/ \___|_|\___|\___||_|

Created by RidgeLab Group, BYU Bioinformatics

Readme file for understanding input for all phases of SelecT pipeline
Primary instructions are created for use on a supercomputer. Modifications for individual setup my be necessary.
This readme is for using .jar executables only. Src code is for reference only.

________________________________________________________________________________________
PHASE 1

Requred Arguments:
[1]	Data Directory		Directory should contain all phased VCF or HAP/LEGEND file required for selection analysis
				Files names must contain proper flags and file extensions
				Includes optional (but highly recommended) Ancestral data embedded in VCF files or as separate LEGEND/EMF file

[2]	Map Directory		Directory that contains all required genetic map files for SelecT analysis
				File names must contain proper chromosome flags

[3]	Start Chromosome	Must be a number between 1-22; sex chromosomes not yet supported

[4]	End Chromosome		Must be a number between 1-22 and greater or equal to Start Chromosome; sex chromosomes not yet supported

[5]	Target Population	Population identifier for experimental population 	
				TST can be used if no standard indentifier exists

[6]	Cross Population	Population identifier for cross population
				TST can be used if no standard indentifier exists
				Cross Population cannont be the same as Target Population

Optional Arguments:
--out_pop	Outgroup Population	Population identifier for outgroup population
					TST can be used if no standard indentifier exists
					Outgroup Population cannont be the same as Target Population

--working_dir	Working Directory	Defines the directory where SelecT will create a new working directory
					Default is current directory

--win_size	Window Size		For changing SelecT analysis window size (in megabases)
					Default is 0.5Mb

Examples:
java -Xmx[MB]m —jar EnviSetup.jar [1] [2] [3] [4] [5] [6] --working_dir=path/to/directory
java -Xmx3000m -jar EnviSetup.jar clean_data/data clean_data/map 21 21 CEU YRI --working_dir=clean_data



________________________________________________________________________________________
PHASE 2

Possible BASH script for supercomputer running SLURM

for i in `seq 0 3`;
do
	sbatch selection.slurm 21 $i
done



________________________________________________________________________________________
PHASE 3

Required Arguments:
[1]	Working Directory	SelecT working directory created in Phase 1
				Working Directory name can be changed but subdirectory names must be unchanged

[2]	Simulation Directory	Directory where simulations can be found
				Must contain simulation file neutral_simulation.tsv and selection_simulation.tsv
				These can be found here: https://github.com/hsmith36/SelecT/tree/master/clean_data/sim

[3]	Chromosome		Chromosome number where window can be found

[4]	Window Number		Window index number as defined by SelecT evironment setup
				See SelecT_workspace/envi_files/all_wins for window ranges

Optional Arguments:
-inon		Non-absolute iHS	Runs iHS score probabilities where large negative scores ONLY are associated with selection (replicate CMS_local)
					Defualt is large positive AND negative iHS scores equate to greater selection (Voight)

--prior_prob	Prior Probability	Set Prior Probability to a custom value between 0.0 and 1.0
					Defaults to (1 / actual number of variants within window)

-pp		Prior Probability Flag	Sets Prior Probability to 1/10,000 (replicate CMS_local)

--daf_cutoff	DAF Cutoff		Defines the derived dllele frequency cutoff for compose score calculation
					Defaults to a DAF value of 0.2
					Special case: DAF value of 0.0 indicates incomplete MoP score calculation, PoP is unchanged

Exmaples:
java -Xmx[MB]m -jar StatsCalc.jar [1] [1] [3] [4]
java -jar StatsCalc.jar clean_data/SelecT_workspace clean_data/sim 21 2



________________________________________________________________________________________
PHASE 4

Required Arguments:
[1] 	Working Directory	SelecT working directory created in Phase 1
				Working Directory name can be changed but subdirectory names must be unchanged

[2]	Chromosome		Chromosome number where window can be found

Optional Arguments:
-co		Combine Only	Only runs first half of analysis where windows are combined into one file

--combine_fltr	Combine Filter	Uses a specific filter for printing specific combination of stats in output file
				Can only be used in conjunction with the -co flag
				i=iHS, x=XPEHH, h=iHH, dd=dDAF, d=DAF, f=Fst, up=unstd_PoP, um=unstd_MoP, p=PoP, m=MoP
				Each tag should be separated by a colon (i.e. i:x:h:dd:d:f:up:um:p:m)

--p_value	p_Value		Sets the p-value cutoff for significance check on composite scores
				Defaults to 0.01

-wc		Write Combine	Similar to -co, but also runs significance filtering

-rn		Normalization	Runs normalization step across the entire dataset/chromosome
				Normalizes by standardization (mean 0; standard deviation 1)

-ui		Use Incomplete	Use incomplete data when analyzing MoP scores

-im		Ignore MoP	Ignores all MoP scores and finds significance based upon PoP only

-ip		Ignore PoP 	Ignores all PoP scores and finds significance based upon MoP only
				If both -im and -ip flags are present significance is found by looking at either PoP or MoP
				If neither -im or -ip flag is present significance is found by looking at both PoP and MoP	








