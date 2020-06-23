
# Simulation model — cultural compression and inflating heritability. 
# for Uchiyama, Spicer & Muthukrishna, "Cultural evolution of genetic heritability"


rm(list = ls())
# =*=*=*=*=*=*=*=*=*=*=*
library(tidyverse)
library(ggplot2)
library(moments)
source("/Users/ryutaro/Desktop/Muthukrishna_Lab/sociality-IQ/simulation run scripts/heritability_sim_functions_02.R")
# ^ set directory containing heritability_sim_functions_xx.R

loops <- 20
repeat_runs <- 1 # number of runs for each generative model configuration (of the coefficient matrix). Simulation is run this number of times for each of the 3 'cult.mode' or tightness levels. So coefficient matrix gets reconfigured every 3*repeat_runs loops    
N <- 500 #number of agents
### cult.mode <- "sel_95tight"  # { "sel_95tight", "sel_95mid", "sel_95loose" } 
cult.steps <- 4
G_vars <- 4
E_vars <- 2
C_vars <- 4
filestring_header <- "0618-test" # set header of file name

interaction.ratio <- 1 # out of possible interactions what proportion should be non-zero 
# interaction.strength <- 1 # baseline is 1: expected sum of all of var x's interactions with other vars is equal to expected size of its main effect. 0.5 halves interaction magnitudes.  
interaction.strength <- G_vars+E_vars+C_vars
ranking_length <- 50 
phenotype.name <- c("smart")
init_alpha_beta <- 12 # initial value for alpha and beta parameters of symmetrical beta distribution


results.variables <- c("step","cult.mode", "interact.ratio", "interact.strength", "init.alphabeta", 
                       "gene.effect","eco.effect","cult.effect",
                       "EandC.effect","GandEandC.effect","interact.effect","full.effect","h2", 
                       "P.IQmean", "P.IQvar","P.IQskewness","P.IQkurtosis","N",
                       "vars_G","vars_E","vars_C")

#WWWWWWWWWWWWWWW OUTER-LOOP WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

for (count in 1:loops){
  print(paste0("loop: ",count, " of ",loops))

  ### step <- 1 # This line is not neccesary. Just use when running the script without the looping through the main loop 
  results <- data.frame(matrix(NA,1:cult.steps,length(results.variables))) 
  names(results) <- results.variables
  GEC_vars <- G_vars + E_vars + C_vars
  ind_G <- 1:G_vars
  ind_E <- (G_vars+1):(G_vars+E_vars)
  ind_C <- (G_vars+E_vars+1):GEC_vars
  ind_GC <- c(1:G_vars, (G_vars+E_vars+1):GEC_vars)
  ind_GE <- 1:(G_vars+E_vars)
  ind_EC <- (G_vars+1):GEC_vars

  # configuring the coefficient matrix, but only once every (repeat_runs * 3) loops
  if (count %% (repeat_runs*3) ==1 ){ #### keeping generative model coefficients constant across cult_mode settings
    print("coefficients updated")
    init_alpha_beta_var <- 1/(4 *(init_alpha_beta*2 +1)) # simplified variance calculation of beta distribution for alpha=beta
    init_alpha_beta_sd <- init_alpha_beta_var^0.5
    maineffects <- rnorm(GEC_vars) /init_alpha_beta_sd 
    interactions <- get_interactions(GEC_vars, init_alpha_beta_sd, interaction.ratio, interaction.strength)
    interactions.double <- (interactions + t(interactions))
    }

  # setting cult.mode, and computing beta distributions
  if (count %% (repeat_runs*3) == 1 ){cult.mode <- "sel_95tight"}
  if (count %% (repeat_runs*3) == (1+repeat_runs) ){cult.mode <- "sel_95mid"}
  if (count %% (repeat_runs*3) == (1+repeat_runs*2) | count %% (repeat_runs*3) == 0 ){cult.mode <- "sel_95loose"}
  cult.alpha <- get_cult.alpha(cult.mode, maineffects, interactions.double, ind_C, C_vars, cult.steps, init_alpha_beta)
  cult.beta <- get_cult.beta(cult.mode, maineffects, interactions.double, ind_C, C_vars, cult.steps, init_alpha_beta)
  
  # prepare tables for logging data
  modelparam_list_full <- data.frame(matrix(NA,cult.steps,GEC_vars+1)) 
  modelparam_list_full[,1] <- 1:cult.steps
  modelparam_list_GEC <- data.frame(matrix(NA,cult.steps,GEC_vars+1)) 
  modelparam_list_GEC[,1] <- 1:cult.steps 
  rankings <- data.frame(matrix(NA, ranking_length *2 *cult.steps, G_vars+E_vars+C_vars+3)) 
  names(rankings)[1:2] <- c("cult.steps","rank")
  rankings$cult.steps <- rep(1:cult.steps, each=ranking_length*2)
  rankings$rank <- rep(c(1:ranking_length,(N-ranking_length+1):N), cult.steps)
  cultparams <- data.frame(matrix(NA,C_vars*4,cult.steps+2))
  cultparams[,1] <- rep(1:C_vars,each=4)
  cultparams[,2] <- c("alpha", "beta","cult.mean","cult.var")
  names(cultparams) <- c("variable","param",as.character(1:cult.steps))
  filestring <- paste0(filestring_header, '-', as.character(count+1000), "--G-", G_vars, "-E-", 
                       E_vars, "-C-", C_vars, "-ratio-", interaction.ratio*10, "-str-", 
                       interaction.strength,"-",cult.mode) 
  
  
  #@@@@@@ INNER LOOP START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  for (step in 1:cult.steps){
   
    pop <- get_pop(N,G_vars,E_vars,C_vars,init_alpha_beta,step,phenotype.name) # generate population scores
    names(rankings)[3:dim(rankings)[2]] <- names(pop)
    
    # each gene has N=(length(eco)+length(cult)) interactions;
    # each eco factor has N=length(cult) interactions (eco X gene is redundant with above); 10 eco, 100 interactions
    interactions.wg.rowsums <- get_interactions.wg.rowsums(N,GEC_vars,interactions)
    maineffects.wg <- sweep(pop[,1:GEC_vars], MARGIN=2, maineffects, '*')
    phenotype <-  rowSums(maineffects.wg) + rowSums(interactions.wg.rowsums) # phenotype scores
    
    if (step==1){  # setting the scale of the phenotype during the first iteration, and continuing to use it to show Flynn effect 
      Pscale.sd <- sd(phenotype)
      Pscale.mean <- mean(phenotype/Pscale.sd *15) # taking mean only after rescaling variability
    }
    phenotype.IQformat <- (phenotype/Pscale.sd *15) - Pscale.mean + 100
    pop[,GEC_vars+1] <- phenotype.IQformat # the phenotype in IQ format, with mean 100 and SD 15
    pop.st <- get_pop.st(pop) # standardized population scores, for input to regression models
    
    #====== GENERATE REGRESSION MODELS =====================
    # generate strings to enter into lm()
    string_G <- character()
    string_E <- character()
    string_C <- character()
    for (i in 1:G_vars){ string_G <- str_c(string_G, paste0("gene_", as.character(i),"+ ")) }
    for (i in 1:E_vars){ string_E <- str_c(string_E, paste0("eco_", as.character(i),"+ ")) }
    for (i in 1:C_vars){ string_C <- str_c(string_C, paste0("cult_", as.character(i),"+ ")) }
    string_G <- substr(string_G, 1, nchar(string_G)-2) # shave off the last "+ " in string
    string_E <- substr(string_E, 1, nchar(string_E)-2) 
    string_C <- substr(string_C, 1, nchar(string_C)-2) 
    
    string_I <- character() #string of all interactions
    for (i in 1:G_vars){ 
      for (j in 1:E_vars){string_I <- str_c(string_I, paste0("gene_",as.character(i),"*eco_",as.character(j),"+ "))}
      for (k in 1:C_vars){string_I <- str_c(string_I, paste0("gene_",as.character(i),"*cult_",as.character(k),"+ "))}}
    for (j in 1:E_vars){
      for (k in 1:C_vars){string_I <- str_c(string_I, paste0("eco_",as.character(j),"*cult_",as.character(k),"+ "))}}
    string_I <- substr(string_I, 1, nchar(string_I)-2) # shave off the last "+ " in string
    
    # generate full string for evaluating model
    modelstring_G <- paste0("smart ~ ", string_G) 
    modelstring_E <- paste0("smart ~ ", string_E)
    modelstring_C <- paste0("smart ~ ", string_C)
    modelstring_EC <- paste0("smart ~ ", string_E,"+",string_C)
    modelstring_GEC <- paste0("smart ~ ", string_G,"+",string_E,"+",string_C)
    modelstring_I <- paste0("smart ~ ", string_I)
    modelstring_full <- paste0("smart ~ ", string_G,"+",string_E,"+",string_C,"+",string_I)
      
    lm_G <- lm( eval(parse(text=modelstring_G)), data=pop.st)
    lm_E <- lm( eval(parse(text=modelstring_E)), data=pop.st)
    lm_C <- lm( eval(parse(text=modelstring_C)), data=pop.st)
    lm_EC <- lm( eval(parse(text=modelstring_EC)), data=pop.st)
    lm_GEC <- lm( eval(parse(text=modelstring_GEC)), data=pop.st)
    lm_I <- lm( eval(parse(text=modelstring_I)), data=pop.st)
    lm_full <- lm( eval(parse(text=modelstring_full)), data=pop.st)
    
    model_G <- summary(lm_G) 
    model_E <- summary(lm_E)
    model_C <- summary(lm_C)
    model_EC <- summary(lm_EC)
    model_GEC <- summary(lm_GEC)
    model_I <- summary(lm_I)
    model_full <- summary(lm_full)
    
    effect_G <- round(model_G$r.squared, digits=3) # extract R^2
    effect_E <- round(model_E$r.squared, digits=3)
    effect_C <- round(model_C$r.squared, digits=3)
    effect_EC <- round(model_EC$r.squared, digits=3)
    effect_GEC <- round(model_GEC$r.squared, digits=3)
    effect_I <- round(model_I$r.squared, digits=3)
    effect_full <- round(model_full$r.squared, digits=3)
    h2 <- round(effect_G/(effect_G + effect_E + effect_C), digits=3)
    cultparams <- get_cultparams(cult.alpha, cult.beta, step)
    
    results[step,] <- c(step, cult.mode, interaction.ratio, interaction.strength, init_alpha_beta,   
                     effect_G, effect_E, effect_C, effect_EC, effect_GEC,
                     effect_I, effect_full,  h2, round(mean(phenotype.IQformat),1),
                     round(var(phenotype.IQformat),1), round(skewness(phenotype.IQformat),3), 
                      round(kurtosis(phenotype.IQformat),3), N, G_vars, E_vars, C_vars)
    
    sortedagents <- cbind(sort(phenotype.IQformat, decreasing=T, index.return=T)$x[c(1:ranking_length,(N-ranking_length+1):N)], 
                    sort(phenotype.IQformat, decreasing=T, index.return=T)$ix[c(1:50,(N-49):N)])
    rankings[which(rankings$cult.steps==step),(3:dim(rankings)[2])] <- round(pop[sortedagents[,2],],3)
  
    # paste model parameters into cumulative matrices
    modelparam_list_full[step, 2:(GEC_vars+1)] <- t(round(model_full$coef,3)[2:(GEC_vars+1)])
    modelparam_list_GEC[step, 2:(GEC_vars+1)] <- t(round(model_GEC$coef,3)[2:(GEC_vars+1)])
  
  } # INNER LOOP END
  
  
  #@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # prepare both ground-truth effects and model estimates for output 
  maineffects.out <- data.frame(t(round(data.frame(maineffects),3))) 
  colnames(maineffects.out) <- names(pop)[1:(dim(pop)[2]-1)]
  alleffects.beta <- round(data.frame(interactions),3) 
  diag(alleffects.beta) <- round(maineffects,3)
  colnames(alleffects.beta) <- names(pop)[1:(dim(pop)[2]-1)]
  alleffects.beta$names <- names(pop)[1:(dim(pop)[2]-1)]
  alleffects.beta$vars_influence <- round(maineffects + colSums(interactions.double), 3)

  # summaries of change in inferred model parameters
  modelparam_list_GEC[cult.steps+1,] <- round(modelparam_list_GEC[cult.steps,]-modelparam_list_GEC[1,], 3)
  modelparam_list_GEC[cult.steps+2,] <- sign(modelparam_list_GEC[cult.steps,]*modelparam_list_GEC[1,])
  modelparam_list_GEC[cult.steps+1,1] <- 0 # change in coefficient from first to last timestep
  modelparam_list_GEC[cult.steps+2,1] <- 0 # change in sign
  
  modelparam_list_GEC
  
  colnames(modelparam_list_full) <- c("cult.step", names(pop)[1:(dim(pop)[2]-1)])
  colnames(modelparam_list_GEC) <- c("cult.step", names(pop)[1:(dim(pop)[2]-1)])

  time <- as.POSIXlt(Sys.time()) #timestamp
  timeB <- paste0("-",time$hour,"h", time$min,"m")
  
  # output files
  write.csv(results, paste0(filestring, timeB, "-summary.csv"))
  write.csv(rankings, paste0(filestring, timeB,"-rankings.csv"))
  write.csv(maineffects.out, paste0(filestring, timeB, "-groundtruth-maineffects.csv"))##
  write.csv(alleffects.beta, paste0(filestring, timeB, "-trueeffects.csv"))
  write.csv(modelparam_list_full, paste0(filestring, timeB, "-maineffects-fullmodel.csv"))
  write.csv(modelparam_list_GEC, paste0(filestring, timeB, "-maineffects-GECmodel.csv"))
  write.csv(cultparams, paste0(filestring, timeB, "-cultparams.csv"))


  } # END OUTER LOOP





