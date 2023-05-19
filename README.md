# SEM_popsim

## Population simulator with phenotypic recursive effects

### Considerations

The code is free to use and modification, parameters loaded by default are used by the authors.

In case of use for research, please cite us.

*THERE IS NO WARRANTY THAT THE PROGRAM IS ERROR FREE. COMMENTS ARE WELCOME.*

The code simulate the following causal network, path coefficients are easy to change, for different causal networks major modifications are required, including tables and dataframes.

<img src='images/network.jpg' width='200'>

### Use of BLUPF90
During simulations RENUMF90 and BLUPF90 executables are required (Animal breeding and genetic group, University of Georgia).   
Can be downloaded from http://nce.ads.uga.edu/html/projects/programs/.  
The files names mut be:
- `renumf90.exe` and `blupf90.exe` if you run in windows
- `renumf90` and `blupf90` if you run in linux
  
Parameter files are used in simulation:  
- `RENUMF90_for_diagonal.par` for diagonal covariances  
- `RENUMF90_for_nodiagonal.par` for no diagonal covariances
  
Code available:  
- `SEM-Sim.single - Windows.R` for simulation a single replicate with a known simulation seed using windows store values for single individuals

## User customization - Single replicate - SEM-Sim.single - Windows.R

**If you are running on Linux system make the following changes**
In section `RUN RENUMF90 AND BLUPF90` in   
*lines 518 - 519; lines 714 - 715; lines 802 - 803 and lines 1006 - 1007*:    
- Replace `system(command = "renumf90.exe", input = MTM.par)` by `system(command = "./renumf90", input = MTM.par)` 
- Replace `system(command = "blupf90.exe", input = MTM.par)` by `system(command = "echo renf90.par | ./blupf90")`

OUTPUT of execution:
- **REGISTROS_base.txt** - Generations without selection, by default 6 (base + 5 random mating generations)
- **REGISTROS_MTM_L** and **REGISTROS_SEM_L**: animals simulated in long format, each phenotipic record in a single row
      - Columns: ID; SIRE (as PAD); DAM (as MAD); GENERATION (as GEN); SEX (as SEX); HERD (as ROD); TRAIT (as Trait); TRUE BREEDING VALUE (as BV); Y; EBV. 
- **REGISTROS_MTM_W** and **REGISTROS_SEM_W**: animals simulated in wide format, each animal in a single row
      - Columns: ID; SIRE (as PAD); DAM (as MAD); GENERATION (as GEN); SEX (as SEX); HERD (as ROD); five columns for TRUE BREEDING VALUE (as BV_#); five columns for Y (as Y_#); five columns for EBV (as EBV_#). 

**PARAMETERS TO CHECK BEFORE RUNNING**
L50: Define working directory, BLUPF90 executables and files must be in this folder.

L53: Define storage directory. Output files will be saved in this folder.

L56: Parameter file for RENUMF90 used for MTM

L60: Parameter file for RENUMF90 for equivalent MTM for SEM estimations

*Parameter files could be the same if the variance parameters are the same between models.*  
*If you want to use estimated parameters these parameters could be different (by error of estimation).*

L67: Simulation seed.

L72 to L80: Population structure. $NVAC$ * $PREPOSH$ and $NTOR$ * $PREPOSM$ must be integer.

L86: S2U: Additive covariance structure for simulation in matrix format.

L94: S2E: Residual covariance structure for simmulation in matrix format

L102: LAMBDA: Causal structure for simulation in matrix format. For non reciprocal relations the matrix resultant is lower diagonal. Read as row is caused by column.

L128: LAMBDA.e: Causal structure for estimation in matrix format. Used for backsolving EBVs, if the process is executed with all known parameters $LAMBDA=LAMBDA.e$.

*After run BLOCK I the matrix for the equivalent model are printed in the console, if you are running the estimations using known parameters these matrix must be used for BLUPF90. Check in RENUMF90 parameter files.*

L547; L584; L838; L876: in `slice_max(EBV#` state the EBV used for selection in format "`EBV#`" (i.e. EBV1, EBV2, ....., EBV5)

## User customization - Multiple replicate - SEM-Sim.multi - Windows.R

**If you are running on Linux system make the following changes**
In section `RUN RENUMF90 AND BLUPF90` in   
*lines 533 - 534; lines 729 - 730; lines 824 - 825 and lines 1028 - 1029*:    
- Replace `system(command = "renumf90.exe", input = MTM.par)` by `system(command = "./renumf90", input = MTM.par)` 
- Replace `system(command = "blupf90.exe", input = MTM.par)` by `system(command = "echo renf90.par | ./blupf90")`

OUTPUT of execution:
- **REGISTROS_base.txt** - Generations without selection, by default 6 (base + 5 random mating generations)
- **REGISTROS_MTM_L** and **REGISTROS_SEM_L**: animals simulated in long format, each phenotipic record in a single row for last simulation
      - Columns: ID; SIRE (as PAD); DAM (as MAD); GENERATION (as GEN); SEX (as SEX); HERD (as ROD); TRAIT (as Trait); TRUE BREEDING VALUE (as BV); Y; EBV. 
- **REGISTROS_MTM_W** and **REGISTROS_SEM_W**: animals simulated in wide format, each animal in a single row for last simulation
      - Columns: ID; SIRE (as PAD); DAM (as MAD); GENERATION (as GEN); SEX (as SEX); HERD (as ROD); five columns for TRUE BREEDING VALUE (as BV_#); five columns for Y (as Y_#); five columns for EBV (as EBV_#).
- **RTA_SEL**: means for each trait by generation for all replicates
      - Columns: METHOD (as MET); SEED; TRAIT UNDER SELECTION (as TR.SEL); GENERATION (as GEN); RESPONSE TRAIT (as TR.RTA); Y; TRUE BREEDING VALUE (as TBV); EBV. 

**PARAMETERS TO CHECK BEFORE RUNNING**
L50: Define working directory, BLUPF90 executables and files must be in this folder.

L53: Define storage directory. Output files will be saved in this folder.

L56: Parameter file for RENUMF90 used for MTM

L60: Parameter file for RENUMF90 for equivalent MTM for SEM estimations

*Parameter files could be the same if the variance parameters are the same between models.*  
*If you want to use estimated parameters these parameters could be different (by error of estimation).*

L67: Simulation seed.

L72 to L80: Population structure. $NVAC$ * $PREPOSH$ and $NTOR$ * $PREPOSM$ must be integer.

L86: S2U: Additive covariance structure for simulation in matrix format.

L94: S2E: Residual covariance structure for simmulation in matrix format

L102: LAMBDA: Causal structure for simulation in matrix format. For non reciprocal relations the matrix resultant is lower diagonal. Read as row is caused by column.

L128: LAMBDA.e: Causal structure for estimation in matrix format. Used for backsolving EBVs, if the process is executed with all known parameters $LAMBDA=LAMBDA.e$.

*After run BLOCK I the matrix for the equivalent model are printed in the console, if you are running the estimations using known parameters these matrix must be used for BLUPF90. Check in RENUMF90 parameter files.*

L162: Number of replicates

L562; L599; L860; L898: in `slice_max(EBV#` state the EBV used for selection in format "`EBV#`" (i.e. EBV1, EBV2, ....., EBV5)

l772; L1084: the number in **3rd column** must be the number of trait under selection
