#----------------------------------------------------#
#   Population simulator for causal related traits   #
#   Last UPDATE: May 2023                            #
#   Bonamy Martin - La Plata University              #
#----------------------------------------------------#

# CONSIDERATIONS
# ----------------------------------------------------------------------------------
#   The code is free to use and modification.
#   In case of use for research, please cite us.  
#   
# 
#   The code simulate the following causal network.  
#   
#       Y1 --> Y2 --> Y4 <-- Y3
#                     Y4 --> Y5
#
#   For different causal networks, modifications of Lambda matrix and recursions on
#   phenotypic values are required.  
#
#   THERE IS NO WARRANTY THAT THE PROGRAM IS ERROR FREE. COMMENTS ARE WELCOME.
#   
#   For user documentation please go to readme file
#
# ------------------------------------------------------------------------------------

#   LIBRARIES USED
# ---------------------------

    library(MASS)
    library(Matrix)
    library(dplyr)
    library(pedigree)
    library(AGHmatrix)
    library(tidyr)
    library(cowplot)

# ----------------------------------------------------------------------------------------
# BLOCK I
# Definitions
# ----------------------------------------------------------------------------------------

# Files and storage  
# ------------------------- 

  {     # Run all definitions.

#     WORKING DIRECTORY - RENUMF90 and BLUPF90 and parameter file must be in this folder
      setwd("~/Folder/Subfolder")
      
#     STORAGE DIRECTORY - where the outfiles will be saved
      store <- "/Folder/Subfolder/Something/"

#     Parameter file for RENUMF90 for Multiple trait model estimations
      MTM.par <- "Parameter_file_MTM.txt"
      
#     Parameter file for RENUMF90 for Structural equation model estimations
#     for BLUPF90 documentation go to http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all2.pdf
      SEM.par <- "Parameter_file_SEM.txt"

# Simulations Parameters
# -------------------------------   
      
#     Random number generator
#     -----------------------
      semilla <- 34145     # Seed number
      set.seed(semilla)

#     Population strucutre
#     --------------------
      NROD = 4      # Number of herds
      NVAC = 250    # Females by herd [NVAC/NTOR must be an integer]
      NTOR = 50     # Total of males 
      NTRAIT = 5    # Number of traits
      PREPOSH = 0.3 # Proportion for FEMALE replace by generation [NVAC*PREPOSH must be an integer]
      PREPOSM = 0.5 # Proportion for MALE replace by generation [NTOR*PREPOSM must be an integer]
      NGEN.R  = 5   # Number of initial random generations
      NGEN.MTM = 5  # Number of generations under MTM
      NGEN.SEM = 5  # Number of generations under SEM

#     Covariance structure
#     ---------------------------------

#     Genetic co-variance matrix 
      S2U <- matrix(c(0.5,  0,    0,    0,    0,
                      0,    0.5,  0,    0,    0,
                      0,    0,    0.5,  0,    0,
                      0,    0,    0,    0.5,  0,
                      0,    0,    0,    0,    0.5), 
                    NTRAIT, NTRAIT , byrow = T)

#     Residual variance            
      S2E <- matrix(c(1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0,
                      0, 0, 1, 0, 0,
                      0, 0, 0, 1, 0,
                      0, 0, 0, 0, 1), 
                    NTRAIT, NTRAIT, byrow = T)
  
#     Causal structure - read as row is caused by column    
      LAMBDA = matrix(c(0,0,0,0,0,           
                        0.4,0,0,0,0,         
                        0,0,0,0,0,           
                        0,-0.3,0.4,0,0,      
                        0,0,0,0.3,0),        
                      NTRAIT, NTRAIT, byrow = T)

#     Aditive variance 
#     -----------------
      I5 <- diag(1, NTRAIT)  
      IL <- I5 - LAMBDA
      ILi <- solve(IL)
      G0s <- ILi %*% S2U %*% t(ILi)
      print("G0s")
      print(G0s)
      print("This matrix must be loaded in parameter file for RENUMF90 if you work with known parameters")

#     Residual variance under equivalen model
#     ---------------------------------------
      R0s <- ILi %*% S2E %*% t(ILi)
      print("R0s")
      print(R0s)
      print("This matrix must be loaded in parameter file for RENUMF90 if you work with known parameters")

#     Estimated causal coefficients - used by backsolving EBV for SEM
#     In known parameters scenario this matrix is equal to LAMBDA
      LAMBDA.e = matrix(c(0,0,0,0,0,           
                          0.4,0,0,0,0,         
                          0,0,0,0,0,           
                          0,-0.3,0.4,0,0,      
                          0,0,0,0.3,0),        
                        NTRAIT, NTRAIT, byrow = T)
      
      IL.e <- I5 - LAMBDA.e
      ILi.e <- solve(IL.e)      
      
# ----------------------------------------------------------------------------------------

# OBJECTS STORED
# --------------

#   Actual geration, actualized in each iteration.  

    #   VACA [COW] - Store Herd#, order#, trait, TrueBV, Phenotype, ID
    VACA <- array(data = NA, dim = c(NROD, NVAC, NTRAIT, 3))
      
    #   TORO [BULL] - Store order#, trait, TrueBV, Phenotype, ID
    TORO <- array(data = NA, dim = c(NTOR, NTRAIT, 3))
      
    #   CRIA [CALF] - Store Herd#, order#, trait, TrueBV, Phenotype, ID, sex, sire, dam 
    CRIA <- array(data = NA, dim = c(NROD, NVAC, NTRAIT, 6))  

  }

# ----------------------------------------------------------------------------------------
# BLOCK II 
# Random mating generations
# ----------------------------------------------------------------------------------------
  
# --------------------  
# BASE GENERATION (G1)
# --------------------
  
  {     # Create and store G1

#   --------------      
#   CREATE ANIMALS
#   --------------

    #    BULLS
    for (tor in 1:NTOR) {
      
      BV <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2U)    # True Breeding Values Sampling 
      
      for (tr in 1:NTRAIT) {
        #Store True Breeding Value
        
        TORO[tor, tr,1] = BV[tr]
      }
      
      res <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2E)  # Residual sampling
      
      # Phenotype - For differents causal networks reassign the recursive effects
      TORO[tor,1,2] = TORO[tor,1,1] + res[1]
      TORO[tor,2,2] = TORO[tor,2,1] + res[2] + 
        TORO[tor,1,2]*LAMBDA[2,1]                     # 1 --> 2
      TORO[tor,3,2] = TORO[tor,3,1] + res[3]
      TORO[tor,4,2] = TORO[tor,4,1] + res[4] + 
        TORO[tor,2,2]*LAMBDA[4,2] +                   # 2 --> 4 
        TORO[tor,3,2]*LAMBDA[4,3]                     # 3 --> 4
      TORO[tor,5,2] = TORO[tor,5,1] + res[5] + 
        TORO[tor,4,2]*LAMBDA[5,4]                     # 4 --> 5 
    
      # ID  
      TORO[tor,,3] = tor
    }
    
    #    COWS
    for (rod in 1:NROD) {
      for (vac in 1:NVAC) {
        
        BV <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2U)             # True Breeding Values Sampling
        
        for (tr in 1:NTRAIT) {
          #True Breeding Value
          
          VACA[rod,vac,tr,1] = BV[tr]
          
        }
        
        res <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2E)            # Residual sampling
        
        # Phenotype - For differents causal networks redefine the recursive effects
        VACA[rod,vac,1,2] = VACA[rod,vac,1,1] + res[1]
        VACA[rod,vac,2,2] = VACA[rod,vac,2,1] + res[2] + 
          VACA[rod,vac,1,2]*LAMBDA[2,1]                   # 1 --> 2
        VACA[rod,vac,3,2] = VACA[rod,vac,3,1] + res[3]
        VACA[rod,vac,4,2] = VACA[rod,vac,4,1] + res[4] + 
          VACA[rod,vac,2,2]*LAMBDA[4,2] +                 # 2 --> 4
          VACA[rod,vac,3,2]*LAMBDA[4,3]                   # 3 --> 4
        VACA[rod,vac,5,2] = VACA[rod,vac,5,1] + res[5] + 
          VACA[rod,vac,4,2]*LAMBDA[5,4]                   # 4 --> 5
        
        # ID
        VACA[rod,vac,,3] = (NTOR + (rod-1)*NVAC) + vac
      }
    }
    
#   ---------------------------
#   STORE ALL DATA IN DATAFRAME
#   ---------------------------
    
#   LAYOUT: ID, SIRE, DAM, GENERATION, SEX, HERDE, BV1 .... BVn, PHENO1 ..... PHENOn

#   Add MALES
    REGISTROS <- data.frame(cbind(TORO[,1,3],   # ID
                                  0,            # SIRE
                                  0,            # DAM
                                  1,            # GENERATION
                                  1,            # SEX
                                  0,            # HERD
                                  TORO[,1,1],   # BV1
                                  TORO[,2,1],   # BV2 
                                  TORO[,3,1],   # BV3 
                                  TORO[,4,1],   # BV4 
                                  TORO[,5,1],   # BV5
                                  TORO[,1,2],   # PHENOTYPE1
                                  TORO[,2,2],   # PHENOTYPE2 
                                  TORO[,3,2],   # PHENOTYPE3 
                                  TORO[,4,2],   # PHENOTYPE4 
                                  TORO[,5,2])   # PHENOTYPE5
                            )

#   Add FEMALES
    REGISTROS <- rbind(REGISTROS,
                       data.frame(cbind(as.vector(t(VACA[,,1,3])),   # ID
                                        0,                           # SIRE
                                        0,                           # DAM
                                        1,                           # GENERATION
                                        2,                           # SEX
                                        rep(1:NROD, each = NVAC),    # HERD
                                        as.vector(t(VACA[,,1,1])),   # BV1
                                        as.vector(t(VACA[,,2,1])),   # BV2 
                                        as.vector(t(VACA[,,3,1])),   # BV3 
                                        as.vector(t(VACA[,,4,1])),   # BV4 
                                        as.vector(t(VACA[,,5,1])),   # BV5
                                        as.vector(t(VACA[,,1,2])),   # PHENOTYPE1
                                        as.vector(t(VACA[,,2,2])),   # PHENOTYPE2 
                                        as.vector(t(VACA[,,3,2])),   # PHENOTYPE3 
                                        as.vector(t(VACA[,,4,2])),   # PHENOTYPE4 
                                        as.vector(t(VACA[,,5,2]))    # PHENOTYPE5
                                        )    
                                  )
                        )    
    
    } # End base generation
    
# ----------------------------------------------------------------------------------------  
    
# ---------------------------------  
# INITIAL RANDOM MATING GENERATIONS 
# ---------------------------------      
  
  { 

#   NEWBORN ANIMALS
#   ---------------    
    
    for (GEN in 1:NGEN.R) {
      

#     Permutation vector to randomize mating    
#     Herd are in rows, females in columns    

      perm = matrix(NA, nrow = NROD, ncol=NVAC)
      
      for (i in 1:nrow(perm)) {
        perm[i,] = rep(sample(seq(1:NTOR)), NVAC/NTOR)
        }
      
#     COMPUTE INBREEDING
#     ------------------  
      
      Fi <- calcInbreeding(REGISTROS[,c(1:3)])    # use later for residual medelian sampling
      
      
#     GENERATE ANIMALS
#     ----------------
      for (rod in 1:NROD) {                           
        for (vac in 1:NVAC) {
          
          # Mendelian residual
          ResMen <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), 
                            Sigma = 0.5*(1-Fi[VACA[rod,vac,1,3]]/2 + Fi[TORO[perm[rod,vac],1,3]]/2)*S2U) 
          
          for (tr in 1:NTRAIT) {
            # True breeding value
            
            CRIA[rod,vac,tr,1] <- 0.5*(VACA[rod,vac,tr,1] + TORO[perm[rod,vac],tr,1]) +    # Sire + Dam
              ResMen[tr]                                                                   # residal
          } 
          
          res <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2E)
          
          # PHENOTYPE
          CRIA[rod,vac,1,2] = CRIA[rod,vac,1,1] + res[1]
          CRIA[rod,vac,2,2] = CRIA[rod,vac,2,1] + res[2] + 
            CRIA[rod,vac,1,2]*LAMBDA[2,1]
          CRIA[rod,vac,3,2] = CRIA[rod,vac,3,1] + res[3]
          CRIA[rod,vac,4,2] = CRIA[rod,vac,4,1] + res[4] + 
            CRIA[rod,vac,2,2]*LAMBDA[4,2] + 
            CRIA[rod,vac,3,2]*LAMBDA[4,3]
          CRIA[rod,vac,5,2] = CRIA[rod,vac,5,1] + res[5] + 
            CRIA[rod,vac,4,2]*LAMBDA[5,4]
          
          
          # ID
          CRIA[rod,vac,,3] = (NTOR + (NVAC * NROD)) + (NVAC*NROD* (GEN - 1)) + (((rod - 1)*NVAC) + vac)
          
          # SEX
          CRIA[rod,vac,,4] = rbinom(1,1,0.5)  # Sampled as 0 and 1
          
          # SIRE
          CRIA[rod,vac,,5] = TORO[perm[rod,vac],1,3]
          
          # DAM
          CRIA[rod,vac,,6] = VACA[rod,vac,1,3]
          
        }
        
      }
      
#       ADD NEWBORN ANIMALS
#       -------------------
              REGISTROS <- rbind(REGISTROS,
                           data.frame(cbind(as.vector(t(CRIA[,,1,3])),   # ID
                                            as.vector(t(CRIA[,,1,5])),   # SIRE
                                            as.vector(t(CRIA[,,1,6])),   # DAM
                                            GEN[]+1,                     # GENERATION
                                            as.vector(t(ifelse(CRIA[,,1,4] == 1, 2, 1))),   # SEX (Change 0 and 1, by 1 and 2)
                                            rep(1:NROD, each = NVAC),    # HERD
                                            as.vector(t(CRIA[,,1,1])),   # BV1
                                            as.vector(t(CRIA[,,2,1])),   # BV2 
                                            as.vector(t(CRIA[,,3,1])),   # BV3 
                                            as.vector(t(CRIA[,,4,1])),   # BV4 
                                            as.vector(t(CRIA[,,5,1])),   # BV5
                                            as.vector(t(CRIA[,,1,2])),   # PHENOTYPE1
                                            as.vector(t(CRIA[,,2,2])),   # PHENOTYPE2 
                                            as.vector(t(CRIA[,,3,2])),   # PHENOTYPE3 
                                            as.vector(t(CRIA[,,4,2])),   # PHENOTYPE4 
                                            as.vector(t(CRIA[,,5,2]))    # PHENOTYPE5
                                            )
                                      )
                           )
      
#     MALES FOR CULLING [RANDOM]
#     --------------------------
      M.DESC <- sample(seq(1:NTOR),size = NTOR*PREPOSM, replace = FALSE)
      
#     REPOSITION MALES [RANDOM]
#     -------------------------
      M.REPOS <- REGISTROS %>% filter(X4 == GEN+1 & X5 == 1) %>% slice_sample(n = length(M.DESC))

#     REPLACE MALES
#     -------------
      for (i in 1:length(M.DESC)) {
        # Replace ID
        TORO[M.DESC[i],,3] <- M.REPOS[i,1]
        
        # Replace BV and phenotype for trait 1
        TORO[M.DESC[i],1,1] <- M.REPOS[i,7]
        TORO[M.DESC[i],1,2] <- M.REPOS[i,12]
        
        # Replace BV and phenotype for trait 2
        TORO[M.DESC[i],2,1] <- M.REPOS[i,8]
        TORO[M.DESC[i],2,2] <- M.REPOS[i,13]
        
        # Replace BV and phenotype for trait 3
        TORO[M.DESC[i],3,1] <- M.REPOS[i,9]
        TORO[M.DESC[i],3,2] <- M.REPOS[i,14]
        
        # Replace BV and phenotype for trait 4
        TORO[M.DESC[i],4,1] <- M.REPOS[i,10]
        TORO[M.DESC[i],4,2] <- M.REPOS[i,15]
        
        # Replace BV and phenotype for trait 5
        TORO[M.DESC[i],5,1] <- M.REPOS[i,11]
        TORO[M.DESC[i],5,2] <- M.REPOS[i,16]
        }
      
#     FEMALES FOR CULLING
#     -------------------
      
#     Same positions for females to cull in all herds
      H.DESC <- sample(seq(1:NVAC),size = NVAC*PREPOSH,replace = FALSE) 
      
#     FEMALES FOR REPOSITION
#     ----------------------
      
      H.REPOS <- REGISTROS %>% filter(X4 == GEN+1 & X5 == 2) %>% group_by(X6) %>% slice_sample(n = length(H.DESC)) %>% ungroup()
      H.REPOS <- as.matrix(H.REPOS)
      
#     REPLACE FEMALES
#     ---------------
      
      for (rod in 1:NROD) {
        for (i in 1:length(H.DESC)) {
          # Reemplazar ID
          VACA[rod,H.DESC[i],,3] <- H.REPOS[(rod-1)*length(H.DESC)+i,1]
          
          # Replace BV and phenotype for trait 1
          VACA[rod,H.DESC[i],1,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,7]
          VACA[rod,H.DESC[i],1,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,12]
          
          # Replace BV and phenotype for trait 2
          VACA[rod,H.DESC[i],2,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,8]
          VACA[rod,H.DESC[i],2,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,13]
          
          # Replace BV and phenotype for trait 3
          VACA[rod,H.DESC[i],3,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,9]
          VACA[rod,H.DESC[i],3,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,14]
          
          # Replace BV and phenotype for trait 4
          VACA[rod,H.DESC[i],4,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,10]
          VACA[rod,H.DESC[i],4,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,15]
          
          # Replace BV and phenotype for trait 5
          VACA[rod,H.DESC[i],5,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,11]
          VACA[rod,H.DESC[i],5,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,16]
          }   
        }   
      } 

#   STORE DATA IN PHISICAL MEMORY
#   -----------------------------
    write.table(REGISTROS, file = paste0(store,"REGISTROS_Base.txt"), 
                sep = " ", row.names = F, col.names = T)
  }

# ----------------------------------------------------------------------------------------      
# END OF INITIAL RANDOm MATING
# ----------------------------------------------------------------------------------------

###########################################################################################
###########################################################################################

#   START OF SELECTION BASED ON EBVs

###########################################################################################
###########################################################################################

  {   # RUN ALL
    
    start_time_tot <- Sys.time()
    
# ------------------------------------------------------------------------------
# 
# BLOCK III  
# SELECTION BASED ON MTM ESTIMATED EBVs
# 
# ------------------------------------------------------------------------------
    
  {
#   UPDATE SEED
#   -----------
    set.seed(semilla)
    
#   INPUT BASE POPULATION
#   ---------------------
    REGISTROS <- read.csv(paste0(store, "REGISTROS_Base.txt"), 
                          sep = " ", header = T)    

#   MAKE A COPY FROM PREVIOUS ARRAYS [base population]
#   -------------------------------------------------- 
    REGISTROS.MTM <- REGISTROS
    TORO.MTM <- TORO
    VACA.MTM <- VACA
    CRIA.MTM <- CRIA    
    
#   GENERATE NEWBORN ANIMALS
#   ------------------------    
    
    {
      start_time <- Sys.time()        
      for (GEN in (NGEN.R+1):((NGEN.R+NGEN.MTM)+1)) {   
        
        
#     WRITE DATA AND PED FILE FOR BLUPF90
#     -----------------------------------
      write.table(REGISTROS.MTM[,c(1:3)], file = "PED.txt", sep = " ", row.names = F, col.names = F) # Pedigree
      write.table(REGISTROS.MTM[,c(1,4:16)], file = "DATOS.txt", sep = " ", row.names = F, col.names = F) # Data
        
#     RUN RENUMF90 AND BLUPF90        
#     ------------------------
        
#     WINDOWS      
      system(command = "renumf90.exe", input = MTM.par)
      system(command = "blupf90.exe", input = "renf90.par")
      
#     IMPORT SOLUTIONS FROM BLUPF90
#     -----------------------------      
        
      sol <- read.table("solutions", header = F, sep = "", skip = 1, dec = ".")
      colnames(sol) <- c("Trait","Efecto","Nivel","Sol")
        
#     TABLE FORMAT FOR ANIMALS SOLUTIONS
#     ----------------------------------
      EBV = sol %>% filter(Efecto == 4) %>% pivot_wider(names_from = Trait, values_from = Sol)

      renadd <- read.table("renadd04.ped", colClasses = c("integer", rep("NULL", 8), "integer"))
        
      EBV <- merge(renadd, EBV, by.x = "V1", by.y = "Nivel")
      colnames(EBV) <- c("r.ID", "ID", "Eff", "EBV1", "EBV2", "EBV3", "EBV4", "EBV5")
      # V10 is original ID, V1 ID given by RENUMF90
        
#     CREATE A DF FOR SELECTION
#     --------------------------------
      MTM.SEL <- merge(REGISTROS.MTM, EBV[,c(2,4:8)], by.x = "X1", by.y = "ID")

#     MALES CULLING [Random]
#     ----------------------
      M.DESC <- sample(seq(1:NTOR),size = NTOR*PREPOSM,replace = FALSE)
      
#     SELECT MALES FOR REPOSITION
#     ---------------------------        
      M.REPOS <- MTM.SEL %>% filter(X4 == GEN & X5 == 1) %>% slice_max(EBV1, n = length(M.DESC))
        
#     REPLACE MALES
#     -------------
      
      for (i in 1:length(M.DESC)) {
          # Replace ID
          TORO.MTM[M.DESC[i],,3] <- M.REPOS[i,1]
          
          # Replace BV and phenotype for trait 1
          TORO.MTM[M.DESC[i],1,1] <- M.REPOS[i,7]
          TORO.MTM[M.DESC[i],1,2] <- M.REPOS[i,12]
          
          # Replace BV and phenotype for trait 2
          TORO.MTM[M.DESC[i],2,1] <- M.REPOS[i,8]
          TORO.MTM[M.DESC[i],2,2] <- M.REPOS[i,13]
          
          # Replace BV and phenotype for trait 3
          TORO.MTM[M.DESC[i],3,1] <- M.REPOS[i,9]
          TORO.MTM[M.DESC[i],3,2] <- M.REPOS[i,14]
          
          # Replace BV and phenotype for trait 4
          TORO.MTM[M.DESC[i],4,1] <- M.REPOS[i,10]
          TORO.MTM[M.DESC[i],4,2] <- M.REPOS[i,15]
          
          # Replace BV and phenotype for trait 5
          TORO.MTM[M.DESC[i],5,1] <- M.REPOS[i,11]
          TORO.MTM[M.DESC[i],5,2] <- M.REPOS[i,16]
          }   # Cierra el for de reemplazo de machos
        
        
#     FEMALE CULLING [Random]
#     -----------------------
      H.DESC <- sample(seq(1:NVAC),size = NVAC*PREPOSH,replace = FALSE) 
        
#     SELECT FEMALES FOR REPOSITION
#     -----------------------------      
      H.REPOS <- MTM.SEL %>% filter(X4 == GEN & X5 == 2) %>% group_by(X6) %>% slice_max(EBV1, n = length(H.DESC)) %>% ungroup()
        
      H.REPOS <- as.matrix(H.REPOS)

#     REPLACE FEMALES
#     ---------------
      for (rod in 1:NROD) {
        for (i in 1:length(H.DESC)) {
            # Replace ID
            VACA.MTM[rod,H.DESC[i],,3] <- H.REPOS[(rod-1)*length(H.DESC)+i,1]
            
            # Replace BV and phenotype for trait 1
            VACA.MTM[rod,H.DESC[i],1,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,7]
            VACA.MTM[rod,H.DESC[i],1,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,12]
            
            # Replace BV and phenotype for trait 2
            VACA.MTM[rod,H.DESC[i],2,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,8]
            VACA.MTM[rod,H.DESC[i],2,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,13]
            
            # Replace BV and phenotype for trait 3
            VACA.MTM[rod,H.DESC[i],3,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,9]
            VACA.MTM[rod,H.DESC[i],3,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,14]
            
            # Replace BV and phenotype for trait 4
            VACA.MTM[rod,H.DESC[i],4,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,10]
            VACA.MTM[rod,H.DESC[i],4,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,15]
            
            # Replace BV and phenotype for trait 5
            VACA.MTM[rod,H.DESC[i],5,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,11]
            VACA.MTM[rod,H.DESC[i],5,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,16]
            }   # Close females within herds
        }     # Close females replace 
        
        
#     GENERATE NEWBORN ANIMALS - MTM SELECTED PARENTS
#     -----------------------------------------------

      #     Permutation vector to randomize mating    
      #     Herd are in rows, females in columns    
      
      perm = matrix(NA, nrow = NROD, ncol=NVAC)
      
      for (i in 1:nrow(perm)) {
        perm[i,] = rep(sample(seq(1:NTOR)), NVAC/NTOR)
        } 
        
#     GET INBREEDING
#     ---------------
      Fi <- calcInbreeding(REGISTROS.MTM[,c(1:3)])

#     PRODUCE NEWBORN ANIMALS
#     -----------------------    
      for (rod in 1:NROD) {
        for (vac in 1:NVAC) {
          
          ResMen <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = 0.5*(1-Fi[VACA.MTM[rod,vac,1,3]]/2 + 
                                                                            Fi[TORO.MTM[perm[rod,vac],1,3]]/2)*S2U) # Mendelian sampling residual
          
          for (tr in 1:NTRAIT) {
            # True breeding value
            CRIA.MTM[rod,vac,tr,1] <- 0.5*(VACA.MTM[rod,vac,tr,1] + TORO.MTM[perm[rod,vac],tr,1]) +    # Parents effect
              ResMen[tr]                                                                               # Residual effect
          } #  Close traits within animals within herd
          
          res <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2E)  # Sampling error
          
          # PHENOTYPES
          CRIA.MTM[rod,vac,1,2] = CRIA.MTM[rod,vac,1,1] + res[1]       # Trait 1
          CRIA.MTM[rod,vac,2,2] = CRIA.MTM[rod,vac,2,1] + res[2] +     # Trait 2
            CRIA.MTM[rod,vac,1,2]*LAMBDA[2,1]                          # 1 --> 2
          CRIA.MTM[rod,vac,3,2] = CRIA.MTM[rod,vac,3,1] + res[3]       # Trait 3
          CRIA.MTM[rod,vac,4,2] = CRIA.MTM[rod,vac,4,1] + res[4] +     # Trait 4
            CRIA.MTM[rod,vac,2,2]*LAMBDA[4,2] +                        # 2 --> 4
            CRIA.MTM[rod,vac,3,2]*LAMBDA[4,3]                          # 3 --> 4
          CRIA.MTM[rod,vac,5,2] = CRIA.MTM[rod,vac,5,1] + res[5] +     # Trait 5
            CRIA.MTM[rod,vac,4,2]*LAMBDA[5,4]                          # 4 --> 5
          
          
          # ID
          CRIA.MTM[rod,vac,,3] = (NTOR + (NVAC * NROD)) + (NVAC*NROD* (GEN-1)) + (((rod - 1)*NVAC) + vac)
          
          # SEX
          CRIA.MTM[rod,vac,,4] = rbinom(1,1,0.5)  # Sampled as 0 and 1
          
          # SIRE
          CRIA.MTM[rod,vac,,5] = TORO.MTM[perm[rod,vac],1,3]
          
          # DAM
          CRIA.MTM[rod,vac,,6] = VACA.MTM[rod,vac,1,3]
          
        }  # Close cows within herds
      }   # close generations
      
#     ADD NEWBORN ANIMALS TO DATABASE
#     -------------------------------
      REGISTROS.MTM <- rbind(REGISTROS.MTM,
                               data.frame(cbind(as.vector(t(CRIA.MTM[,,1,3])),   # ID
                                                as.vector(t(CRIA.MTM[,,1,5])),   # SIRE
                                                as.vector(t(CRIA.MTM[,,1,6])),   # DAM
                                                GEN[]+1,                         # GENERATION
                                                as.vector(t(ifelse(CRIA.MTM[,,1,4] == 1, 2, 1))),   # SEX
                                                rep(1:NROD, each = NVAC),        # HERD
                                                as.vector(t(CRIA.MTM[,,1,1])),   # BV1
                                                as.vector(t(CRIA.MTM[,,2,1])),   # BV2 
                                                as.vector(t(CRIA.MTM[,,3,1])),   # BV3 
                                                as.vector(t(CRIA.MTM[,,4,1])),   # BV4 
                                                as.vector(t(CRIA.MTM[,,5,1])),   # BV5
                                                as.vector(t(CRIA.MTM[,,1,2])),   # PHENO1
                                                as.vector(t(CRIA.MTM[,,2,2])),   # PHENO2 
                                                as.vector(t(CRIA.MTM[,,3,2])),   # PHENO3 
                                                as.vector(t(CRIA.MTM[,,4,2])),   # PHENO4 
                                                as.vector(t(CRIA.MTM[,,5,2]))    # PHENO5
                                                )
                                          )
                             )
        
      } # CLOSE MTM LOOP
      
#     CREATE TABLES FOR EXPORT
#     ------------------------

#     WRITE DATAFILE AND PED FOR BLUPF90 (The last generations do not have EBVs)
#     ----------------------------------------------------------------------------
      write.table(REGISTROS.MTM[,c(1:3)], file = "PED.txt", sep = " ", row.names = F, col.names = F) # Pedigree
      write.table(REGISTROS.MTM[,c(1,4:16)], file = "DATOS.txt", sep = " ", row.names = F, col.names = F) # Data
      
#     CALL BLUPF90        
#     ------------
      
#     WINDOWS      
      system(command = "renumf90.exe", input = MTM.par)
      system(command = "blupf90.exe", input = "renf90.par")
      
#     Import EBVs
      sol <- read.table("solutions", header = F, sep = "", skip = 1, dec = ".")
      colnames(sol) <- c("Trait","Efecto","Nivel","Sol")
      
      EBV = sol %>% filter(Efecto == 4) %>% pivot_wider(names_from = Trait, values_from = Sol)
      
      renadd <- read.table("renadd04.ped", colClasses = c("integer", rep("NULL", 8), "integer"))
      
      EBV <- merge(renadd, EBV, by.x = "V1", by.y = "Nivel")
      colnames(EBV) <- c("r.ID", "ID", "Eff", "EBV1", "EBV2", "EBV3", "EBV4", "EBV5")

#     Create an unique DF
      MTM.END <- merge(REGISTROS.MTM, EBV[,c(2,4:8)], by.x = "X1", by.y = "ID")
      colnames(MTM.END) <- c("ID","PAD","MAD","GEN","SEX","ROD",
                             "BV_1","BV_2","BV_3","BV_4","BV_5",
                             "Y_1","Y_2","Y_3","Y_4","Y_5",
                             "EBV_1","EBV_2","EBV_3","EBV_4","EBV_5")
      
      MTM.END2 = MTM.END %>% pivot_longer(-c(ID, PAD, MAD, GEN, SEX, ROD)) %>% separate(name, "_", into = c("Trait","Valor")) %>% pivot_wider(names_from = Trait, values_from = value)
      colnames(MTM.END2) <- c("ID","PAD","MAD","GEN","SEX","ROD","Trait","BV","Y","EBV")
      
      MTM.END2 <- as.data.frame(MTM.END2)
      
#     WRITE FILES FOR MTM SELECTION PROCESS
#     -------------------------------------
      write.table(MTM.END, file = paste0(store,"REGISTROS_MTM_W.txt"), 
                  sep = " ", row.names = F, col.names = T)
      write.table(MTM.END2, file = paste0(store,"REGISTROS_MTM_L.txt"), 
                  sep = " ", row.names = F, col.names = T)
      
      print("FIN PROCESO MTM")
      
      end_time <- Sys.time()
      
      end_time - start_time
      
      } # CLOSE TIME
    }   # CLOSE MTM
    
# ----------------------------------------------------------------------------------------      
# END FOR MTM SELECTION PROCESS
# ----------------------------------------------------------------------------------------
    
# ------------------------------------------------------------------------------
# 
# SECTION IV  
# SELECTION BASED ON SEM ESTIMATED EBVs
# 
# ------------------------------------------------------------------------------  
  
    {
      
#   UPDATE SEED
#   ------------------
    set.seed(semilla)
    
#   INPUT BASE PUPOLATION
#   ---------------------
    REGISTROS <- read.csv(paste0(store, "REGISTROS_Base.txt"), 
                          sep = " ", header = T)     
    
    #   MAKE A COPY FROM PREVIOUS ARRAYS [base population]
#   ------------------------------------------------------
    REGISTROS.SEM <- REGISTROS
    TORO.SEM <- TORO
    VACA.SEM <- VACA
    CRIA.SEM <- CRIA    
    
    #   GENERATE NEWBORN ANIMALS
#   ----------------------------    
    
    {
      start_time <- Sys.time()        
      
      for (GEN in (NGEN.R+1):((NGEN.R+NGEN.SEM)+1)) {

#     WRITE DATA AND PED FILE FOR BLUPF90
#     -----------------------------------
      write.table(REGISTROS.SEM[,c(1:3)], file = "PED.txt", sep = " ", row.names = F, col.names = F) # Pedigree
      write.table(REGISTROS.SEM[,c(1,4:16)], file = "DATOS.txt", sep = " ", row.names = F, col.names = F) # Datos
        
#     CALL BLUPF90        
#     ------------
      
#     WINDOWS      
      system(command = "renumf90.exe", input = SEM.par)
      system(command = "blupf90.exe", input = "renf90.par")
      

#     IMPORT SOLUTIONS FOR BLUPF90
#     -----------------------------  
      sol <- read.table("solutions", header = F, sep = "", skip = 1, dec = ".")
      colnames(sol) <- c("Trait","Efecto","Nivel","Sol")
        
#     TABLE FORMAT FOR ANIMALS SOLUTIONS
#     ----------------------------------
      EBV = sol %>% filter(Efecto == 4) %>% pivot_wider(names_from = Trait, values_from = Sol)
      
      renadd <- read.table("renadd04.ped", colClasses = c("integer", rep("NULL", 8), "integer"))
      
      EBV <- merge(renadd, EBV, by.x = "V1", by.y = "Nivel")
      colnames(EBV) <- c("r.ID", "ID", "Eff", "EBV1", "EBV2", "EBV3", "EBV4", "EBV5")
      #  V10 is original ID, V1 ID given by RENUMF90
        
#     GET SEM EBVs BY BACKSOLVING
#     ---------------------------
      for (i in 1:nrow(EBV)) {
          EBV[i,c(4:8)] <- t(IL.e%*%t(as.vector(EBV[i,c(4:8)])))
        }

#     CREATE A DF FOR SELECTION
#     -------------------------
      SEM.SEL <- merge(REGISTROS.SEM, EBV[,c(2,4:8)], by.x = "X1", by.y = "ID")
        
        
#     MALES CULLING
#     -------------
      M.DESC <- sample(seq(1:NTOR),size = NTOR*PREPOSM,replace = FALSE)

#     SELECT MALES TO REPLACE
#     -----------------------
      M.REPOS <- SEM.SEL %>% filter(X4 == GEN & X5 == 1) %>% slice_max(EBV1, n = length(M.DESC))
        
        
#     REPLACE MALES
#     -------------
      for (i in 1:length(M.DESC)) {
          # Replace ID
          TORO.SEM[M.DESC[i],,3] <- M.REPOS[i,1]
          
          # Replace BV and phenotype for trait 1
          TORO.SEM[M.DESC[i],1,1] <- M.REPOS[i,7]
          TORO.SEM[M.DESC[i],1,2] <- M.REPOS[i,12]
          
          # Replace BV and phenotype for trait 2
          TORO.SEM[M.DESC[i],2,1] <- M.REPOS[i,8]
          TORO.SEM[M.DESC[i],2,2] <- M.REPOS[i,13]
          
          # Replace BV and phenotype for trait 3
          TORO.SEM[M.DESC[i],3,1] <- M.REPOS[i,9]
          TORO.SEM[M.DESC[i],3,2] <- M.REPOS[i,14]
          
          # Replace BV and phenotype for trait 4
          TORO.SEM[M.DESC[i],4,1] <- M.REPOS[i,10]
          TORO.SEM[M.DESC[i],4,2] <- M.REPOS[i,15]
          
          # Replace BV and phenotype for trait 5
          TORO.SEM[M.DESC[i],5,1] <- M.REPOS[i,11]
          TORO.SEM[M.DESC[i],5,2] <- M.REPOS[i,16]
          
        }   # Close male reposition
        
        
#     FEMALE CULLING
#     --------------
      H.DESC <- sample(seq(1:NVAC),size = NVAC*PREPOSH,replace = FALSE) 
      
#     SELECT FEMALES FOR REPOSITION
#     -----------------------------       
      H.REPOS <- SEM.SEL %>% filter(X4 == GEN & X5 == 2) %>% group_by(X6) %>% slice_max(EBV1, n = length(H.DESC)) %>% ungroup()
        
      H.REPOS <- as.matrix(H.REPOS)
        
        
#     REPLACE FEMALES
#     ---------------
      for (rod in 1:NROD) {
          for (i in 1:length(H.DESC)) {
            # Replace ID
            VACA.SEM[rod,H.DESC[i],,3] <- H.REPOS[(rod-1)*length(H.DESC)+i,1]
            
            # Replace BV and phenotype for trait 1
            VACA.SEM[rod,H.DESC[i],1,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,7]
            VACA.SEM[rod,H.DESC[i],1,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,12]
            
            # Replace BV and phenotype for trait 2
            VACA.SEM[rod,H.DESC[i],2,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,8]
            VACA.SEM[rod,H.DESC[i],2,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,13]
            
            # Replace BV and phenotype for trait 3
            VACA.SEM[rod,H.DESC[i],3,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,9]
            VACA.SEM[rod,H.DESC[i],3,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,14]
            
            # Replace BV and phenotype for trait 4
            VACA.SEM[rod,H.DESC[i],4,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,10]
            VACA.SEM[rod,H.DESC[i],4,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,15]
            
            # Replace BV and phenotype for trait 5
            VACA.SEM[rod,H.DESC[i],5,1] <- H.REPOS[(rod-1)*length(H.DESC)+i,11]
            VACA.SEM[rod,H.DESC[i],5,2] <- H.REPOS[(rod-1)*length(H.DESC)+i,16]
          } # Close females within herd
        }   # Close female reposition      
        
        
#     GENERATE NEWBORN ANIMALS - SEM SELECTED PARENTS
#     -----------------------------------------------
        
      #     Permutation vector to randomize mating    
      #     Herd are in rows, females in columns 
              
      perm = matrix(NA, nrow = NROD, ncol=NVAC)
        
      for (i in 1:nrow(perm)) {
          perm[i,] = rep(sample(seq(1:NTOR)), NVAC/NTOR)
        } 
        
#     GET INBREEDING
#     --------------
      Fi <- calcInbreeding(REGISTROS.SEM[,c(1:3)])    

#     PRODUCE NEWBORN ANIMALS
#     -----------------------
      for (rod in 1:NROD) {
        for (vac in 1:NVAC) {
          
          ResMen <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = 0.5*(1-Fi[VACA.SEM[rod,vac,1,3]]/2 + 
                                                                            Fi[TORO.SEM[perm[rod,vac],1,3]]/2)*S2U) # Medelian sampling residual
          
          for (tr in 1:NTRAIT) {
            # True breeding value
            CRIA.SEM[rod,vac,tr,1] <- 0.5*(VACA.SEM[rod,vac,tr,1] + TORO.SEM[perm[rod,vac],tr,1]) +    # Parents
              ResMen[tr]                                                                               # Residual
          } # Close traits within animals within herd
          
          res <- mvrnorm(n = 1, mu = matrix(c(rep(0,5))), Sigma = S2E)
          
          # PHENOTYPES
          CRIA.SEM[rod,vac,1,2] = CRIA.SEM[rod,vac,1,1] + res[1]
          CRIA.SEM[rod,vac,2,2] = CRIA.SEM[rod,vac,2,1] + res[2] + 
            CRIA.SEM[rod,vac,1,2]*LAMBDA[2,1]
          CRIA.SEM[rod,vac,3,2] = CRIA.SEM[rod,vac,3,1] + res[3]
          CRIA.SEM[rod,vac,4,2] = CRIA.SEM[rod,vac,4,1] + res[4] + 
            CRIA.SEM[rod,vac,2,2]*LAMBDA[4,2] + 
            CRIA.SEM[rod,vac,3,2]*LAMBDA[4,3]
          CRIA.SEM[rod,vac,5,2] = CRIA.SEM[rod,vac,5,1] + res[5] + 
            CRIA.SEM[rod,vac,4,2]*LAMBDA[5,4]
          
          
          # ID
          CRIA.SEM[rod,vac,,3] = (NTOR + (NVAC * NROD)) + (NVAC*NROD* (GEN-1)) + (((rod - 1)*NVAC) + vac)
          
          # SEX
          CRIA.SEM[rod,vac,,4] = rbinom(1,1,0.5)  # Aca muestrea como 0 y 1
          
          # SIRE
          CRIA.SEM[rod,vac,,5] = TORO.SEM[perm[rod,vac],1,3]
          
          # DAM
          CRIA.SEM[rod,vac,,6] = VACA.SEM[rod,vac,1,3]
        } # Close cows within herds
      }   # Close generations
      
#     ADD NEWBORN ANIMALS TO DATABASE
#     -------------------------------
      REGISTROS.SEM <- rbind(REGISTROS.SEM,
                               data.frame(cbind(as.vector(t(CRIA.SEM[,,1,3])),   # ID
                                                as.vector(t(CRIA.SEM[,,1,5])),   # SIRE
                                                as.vector(t(CRIA.SEM[,,1,6])),   # DAM
                                                GEN[]+1,                         # GENERATION
                                                as.vector(t(ifelse(CRIA.SEM[,,1,4] == 1, 2, 1))),   # SEX
                                                rep(1:NROD, each = NVAC),        # HERD
                                                as.vector(t(CRIA.SEM[,,1,1])),   # BV1
                                                as.vector(t(CRIA.SEM[,,2,1])),   # BV2 
                                                as.vector(t(CRIA.SEM[,,3,1])),   # BV3 
                                                as.vector(t(CRIA.SEM[,,4,1])),   # BV4 
                                                as.vector(t(CRIA.SEM[,,5,1])),   # BV5
                                                as.vector(t(CRIA.SEM[,,1,2])),   # PHENO1
                                                as.vector(t(CRIA.SEM[,,2,2])),   # PHENO2 
                                                as.vector(t(CRIA.SEM[,,3,2])),   # PHENO3 
                                                as.vector(t(CRIA.SEM[,,4,2])),   # PHENO4 
                                                as.vector(t(CRIA.SEM[,,5,2]))    # PHENO5
                                                )
                                          )
                             )
        
      } # CLOSE SEM LOOP

#     CREATE TABLES FOR EXPORT
#     ------------------------
      
#     WRITE DATAFILE AND PED FOR BLUPF90 (The last generations do not have EBVs)
#     --------------------------------------------------------------------------
      write.table(REGISTROS.SEM[,c(1:3)], file = "PED.txt", sep = " ", row.names = F, col.names = F) # Pedigree
      write.table(REGISTROS.SEM[,c(1,4:16)], file = "DATOS.txt", sep = " ", row.names = F, col.names = F) # Data
      
#     CALL BLUPF90        
#     ------------
      
#     WINDOWS      
      system(command = "renumf90.exe", input = SEM.par)
      system(command = "blupf90.exe", input = "renf90.par")

#     IMPORT EBVs
      sol <- read.table("solutions", header = F, sep = "", skip = 1, dec = ".")
      colnames(sol) <- c("Trait","Efecto","Nivel","Sol")
      
      EBV = sol %>% filter(Efecto == 4) %>% pivot_wider(names_from = Trait, values_from = Sol)
      
      renadd <- read.table("renadd04.ped", colClasses = c("integer", rep("NULL", 8), "integer"))
      
      EBV <- merge(renadd, EBV, by.x = "V1", by.y = "Nivel")
      colnames(EBV) <- c("r.ID", "ID", "Eff", "EBV1", "EBV2", "EBV3", "EBV4", "EBV5")
      
#     GET SEM EBV BY BACKSOLVING
#     --------------------------
      for (i in 1:nrow(EBV)) {
        EBV[i,c(4:8)] <- t(IL.e%*%t(as.vector(EBV[i,c(4:8)])))
      }
      
#     CREATE AN UNIQUE DF
#     -------------------
      SEM.END <- merge(REGISTROS.SEM, EBV[,c(2,4:8)], by.x = "X1", by.y = "ID")
      colnames(SEM.END) <- c("ID","PAD","MAD","GEN","SEX","ROD",
                             "BV_1","BV_2","BV_3","BV_4","BV_5",
                             "Y_1","Y_2","Y_3","Y_4","Y_5",
                             "EBV_1","EBV_2","EBV_3","EBV_4","EBV_5")
      
      SEM.END2 = SEM.END %>% pivot_longer(-c(ID, PAD, MAD, GEN, SEX, ROD)) %>% separate(name, "_", into = c("Trait","Valor")) %>% pivot_wider(names_from = Trait, values_from = value)
      colnames(SEM.END2) <- c("ID","PAD","MAD","GEN","SEX","ROD","Trait","BV","Y","EBV")
      
      SEM.END2 <- as.data.frame(SEM.END2)
      
#     WRITE FILES FOR SEM PROCESS
#     ---------------------------
      write.table(SEM.END, file = paste0(store,"REGISTROS_SEM_W.txt"), 
                  sep = " ", row.names = F, col.names = T)
      write.table(SEM.END2, file = paste0(store,"REGISTROS_SEM_L.txt"), 
                  sep = " ", row.names = F, col.names = T)      

      print("FIN PROCESO SEM")
      
      end_time <- Sys.time()
      
      end_time - start_time
      } # Cierra el loop del tiempo
    } # Cierra el loop de SEM
    
# ----------------------------------------------------------------------------------------      
# END FOR SEM SELECTION PROCESS
# ----------------------------------------------------------------------------------------
  

    end_time_tot <- Sys.time()
    
    print("FIN PROCESO COMPLETO")
    end_time_tot - start_time_tot
    
  } # CLOSE COMPLETE PROCESS [MTM + SEM]    
    
      
