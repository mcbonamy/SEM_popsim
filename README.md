# SEM_popsim

## Population simulator with phenotypic recursive effects

### Considerations

The code is free to use and modification.

In case of use for research, please cite us.

*THERE IS NO WARRANTY THAT THE PROGRAM IS ERROR FREE. COMMENTS ARE WELCOME.*

The code simulate the following causal network, path coefficients are easy to change, for different causal networks major modifications are required, including tables and dataframes.

<img src='images/network.jpg' width='200'>

### Use of BLUPF90
During simulations RENUMF90 and BLUPF90 executables are required (Animal breeding and genetic group, University of Georgia).   
Can be downloaded from http://nce.ads.uga.edu/html/projects/programs/
  
Parameter files are used in simulation:  
- `RENUMF90_for_diagonal.par` for diagonal covariances  
- `RENUMF90_for_nodiagonal.par` for no diagonal covariances

## [User customization]

L49: Define working directory, BLUPF90 executables and files must be in this folder.

L52: Define storage directory. Output files will be saved in this folder.

L55: Parameter file for RENUMF90 used for MTM

L59: Parameter file for RENUMF90 for equivalent MTM for SEM estimations

*Parameter files could be the same if the variance parameters are the same between models. If you want to use estimated parameters these parameters could be different (by error of estimation).*

L66: Simulation seed.

L71 to L79: Population structure. $NVAC*PREPOSH$ and $NTOR*PREPOSM$ must be integer.

L85: S2U: Additive covariance structure for simulation in matrix format.

L93: S2E: Residual covariance structure for simmulation in matrix format

L101: LAMBDA: Causal structure for simulation in matrix format. For non reciprocal relations the matrix resultant is lower diagonal. Read as row is caused by column.

L127: LAMBDA.e: Causal structure for estimation in matrix format. Used for backsolving EBVs, if the process is executed with all known parameters $LAMBDA=LAMBDA.e$.

*After run BLOCK I the matrix for the equivalent model are printed in the console, if you are running the estimations using known parameters these matrix must be used for BLUPF90. Check in RENUMF90 parameter files.*

L546; L583; L837; L875: in `slice_max(EBV#` state the EBV used for selection in format "`EBV#`" (i.e. EBV1, EBV2, ....., EBV5)
