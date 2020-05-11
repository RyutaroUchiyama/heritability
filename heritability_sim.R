


rm(list = ls())
# =*=*=*=*=*=*=*=*=*=*=*
library(tidyverse)
library(ggplot2)
library(moments)


N <- 5000 #number of agents
cult.mode <- "sel_95tight"  # { "sel_95tight", "sel_95mid", "sel_95loose"  }
cult.steps <- 6
G_vars <- 4
E_vars <- 2
C_vars <- 4
P_vars <- 1
interaction.ratio <- 0.2 # out of possible interactions what proportion should be non-zero 
interaction.strength <- 4 # baseline is 1: expected sum of all of var x's interactions with other vars is equal to expected size of its main effect. 0.5 halves interaction magnitudes.  
ranking_length <- 50  
pheno.names <- c("smart")
init_alpha_beta <- 12 # initial value for alpha and beta parameters of beta distribution
filestring <- "0506-matchingacross-tightloose"

results.variables <- c("step","cult.mode", "interact.ratio", "interact.strength", "init.alphabeta", 
                       "gene.effect","eco.effect","cult.effect",
                       "EandC.effect","GandEandC.effect","interact.effect","full.effect","h2", 
                       "P.IQmean", "P.IQvar","P.IQskewness","P.IQkurtosis","N",
                       "vars_G","vars_E","vars_C")

results <- data.frame(matrix(NA,1:cult.steps,length(results.variables))) 
names(results) <- results.variables

step <- 1 # This is not neccesary. Just for running the script without the main loop 
GEC_vars <- G_vars + E_vars + C_vars
ind_G <- 1:G_vars
ind_E <- (G_vars+1):(G_vars+E_vars)
ind_C <- (G_vars+E_vars+1):GEC_vars
ind_GC <- c(1:G_vars, (G_vars+E_vars+1):GEC_vars)
ind_GE <- 1:(G_vars+E_vars)
ind_EC <- (G_vars+1):GEC_vars

# setting up parameters of the generative model
init_alpha_beta_var <- 1/(4 *(init_alpha_beta*2 +1)) # simplified variance calculation of beta distribution for alpha=beta
init_alpha_beta_sd <- init_alpha_beta_var^0.5

maineffects <- rnorm(GEC_vars) /init_alpha_beta_sd 
interactions <- matrix(0, GEC_vars, GEC_vars)
entries <- upper.tri(interactions)
interactions[entries] <- rnorm(sum(entries)) /init_alpha_beta_sd
interaction.removal <- c(rep(0, round( (1-interaction.ratio)*sum(entries) )),
                         rep(1, round( interaction.ratio*sum(entries) )) )
if (length(interaction.removal) < sum(entries)){ # just in case the round() above distorts the number of entries. Only assuming discrepancy of one item
  interaction.removal <- c(interaction.removal, sample(0:1,1,replace=T,prob=c((1-interaction.ratio),interaction.ratio)))}
if (length(interaction.removal) > sum(entries)){ 
  interaction.removal <- interaction.removal[1:(length(interaction.removal)-1)] }
interactions[entries] <-  interactions[entries] *sample(interaction.removal)
interactions <- interactions /((GEC_vars-1) *interaction.ratio) * interaction.strength # first controlling for number of interactions culled by interaction.ratio, then multiplying by interaction.strength 
interactions.double <- (interactions + t(interactions))

sel_95tight <- c(800,42); sel_95mid <- c(80,4.2); sel_95loose <- c(8,0.42);
sel <- eval(parse(text=cult.mode))
cult.sign <- sign( maineffects[ind_C] + colSums(interactions.double[,ind_C])) # which direction do the cult distribuitions need to move
cult.alpha <- matrix(NA, cult.steps, C_vars)
cult.beta <- matrix(NA, cult.steps, C_vars)
for (i in 1:C_vars){
  if (cult.sign[i]== 1){
    cult.alpha[,i] <-  exp(seq(log(init_alpha_beta), log(sel[1]),length.out=cult.steps))
    cult.beta[,i] <-   exp(seq(log(init_alpha_beta), log(sel[2]),length.out=cult.steps))}
  if (cult.sign[i]== -1){
    cult.alpha[,i] <-  exp(seq(log(init_alpha_beta), log(sel[2]),length.out=cult.steps))
    cult.beta[,i] <-   exp(seq(log(init_alpha_beta), log(sel[1]),length.out=cult.steps))}
}
# when alpha=beta=12= exp(2.4849), variance is at 0.01 (initial distribution)
# when alpha=beta=1250 = exp(7.1305) , variance is at 0.0001 (final dist in stabilizing selection)
# when alpha=beta= 2 = exp(0.6931) , variance is at 0.05 (final dist in entropic "selection")
# when alpha=98=exp(4.5849) & beta=1, variance is at 0.0001 (final dists in directional selection)
# when alpha=beta=0.75 , variance is at 0.1 (final dist in divergent selection)

# prepare tables for logging data
modelparam_list_full <- data.frame(matrix(NA,cult.steps,GEC_vars+1)) 
modelparam_list_full[,1] <- 1:cult.steps
modelparam_list_GEC <- data.frame(matrix(NA,cult.steps,GEC_vars+1)) 
modelparam_list_GEC[,1] <- 1:cult.steps 
modelparam_matrix <- matrix(NA , GEC_vars*cult.steps, GEC_vars+1)
modelparam_matrix[,1] <- rep(1:cult.steps, each=GEC_vars)
rankings <- data.frame(matrix(NA, ranking_length *2 *cult.steps, G_vars+E_vars+C_vars+P_vars+2)) 
names(rankings)[1:2] <- c("cult.steps","rank")
rankings$cult.steps <- rep(1:cult.steps, each=ranking_length*2)
rankings$rank <- rep(c(1:ranking_length,(N-ranking_length+1):N), cult.steps)
cultparams <- data.frame(matrix(NA,C_vars*4,cult.steps+2))
cultparams[,1] <- rep(1:C_vars,each=4)
cultparams[,2] <- c("alpha", "beta","cult.mean","cult.var")
names(cultparams) <- c("variable","param",as.character(1:cult.steps))



#@@@@@@ MAIN LOOP START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
for (step in 1:cult.steps){
 
  # initial scores for each genetic, ecological, and cultural factor
  venus_G <- matrix(0, N, G_vars)
  venus_G <- apply(venus_G, 2, function(x) x + rbeta(N, init_alpha_beta, init_alpha_beta)) # eco distribuition is same as initial dist for culture
  venus_E <- matrix(0, N, E_vars)
  venus_E <- apply(venus_E, 2, function(x) x + rbeta(N, init_alpha_beta, init_alpha_beta)) # eco distribuition is same as initial dist for culture
  venus_C <- matrix(0, N, C_vars)
  for (i in 1:C_vars){
    venus_C[,i] <- rbeta(N, cult.alpha[step,i], cult.beta[step,i])} # cultural traits, each distribution moving in the direction of positive phenotypic outcome
  venus_P <- matrix(rep(0,N), N, P_vars)
  venus <- data.frame(cbind(venus_G,venus_E,venus_C,venus_P))
  
  # assigning names
  names(venus)[(GEC_vars+1):(GEC_vars+P_vars)] <- pheno.names
  venus[,1:(G_vars+E_vars)] <- round(venus[,1:(G_vars+E_vars)], 2)
  for(i in 1:G_vars) { names(venus)[i] <- paste0("gene_",as.character(i)) }
  for(i in 1:E_vars) { names(venus)[i+G_vars] <- paste0("eco_",as.character(i)) }
  for(i in 1:C_vars) { names(venus)[i+G_vars+E_vars] <- paste0("cult_",as.character(i)) }
  names(rankings)[3:dim(rankings)[2]] <- names(venus)
  
  
#@@@@@@ INTERACTION EFFECTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # each gene has N=(length(eco)+length(cult)) interactions;
  # each eco factor has N=length(cult) interactions (eco X gene is redundant with above); 10 eco, 100 interactions
  interactions.wg.rowsums <- matrix(NA,N,GEC_vars)
  for (agent in 1:N){
    interactions.wg <- matrix(0, GEC_vars, GEC_vars )
    for (i in 1:GEC_vars){ 
      for (j in 1:GEC_vars){ # {genes X eco} and {genes X culture}
        interactions.wg[i,j] <- venus[agent,i] * venus[agent,j] * interactions[i,j] }}
    interactions.wg.rowsums[agent,] <- rowSums(interactions.wg, na.rm=T)
    if (agent %% 5000 == 0 ) {print(agent)}
    }
  maineffects.wg <- sweep(venus[,1:GEC_vars], MARGIN=2, maineffects, '*')
  phenotype <-  rowSums(maineffects.wg) + rowSums(interactions.wg.rowsums) 
  
  if (step==1){  # setting the scale of the phenotype during the first iteration, and continuing to use it to show Flynn effect 
    Pscale.sd <- sd(phenotype)
    Pscale.mean <- mean(phenotype/Pscale.sd *15) # taking mean only after rescaling variability
  }
  phenotype.IQformat <- (phenotype/Pscale.sd *15) - Pscale.mean + 100
  venus[,GEC_vars+1] <- phenotype.IQformat
  
  venus.st <- data.frame(matrix(NA,dim(venus)[1],dim(venus)[2]))
  venus.st[,1:5] <- scale(venus[,1:5])
  venus.st[,2:dim(venus)[2]] <- scale(venus[,2:dim(venus)[2]]) 
  colnames(venus.st) <- colnames(venus)
  
#@@@@@@ GENERATE REGRESSION MODELS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    
  lm_G <- lm( eval(parse(text=modelstring_G)), data=venus.st)
  lm_E <- lm( eval(parse(text=modelstring_E)), data=venus.st)
  lm_C <- lm( eval(parse(text=modelstring_C)), data=venus.st)
  lm_EC <- lm( eval(parse(text=modelstring_EC)), data=venus.st)
  lm_GEC <- lm( eval(parse(text=modelstring_GEC)), data=venus.st)
  lm_I <- lm( eval(parse(text=modelstring_I)), data=venus.st)
  lm_full <- lm( eval(parse(text=modelstring_full)), data=venus.st)
  
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
  
  
  for (i in 1:C_vars){ # entering properties of each cultural distribution 
    beta.a <- round(cult.alpha[step,i], 3)
    beta.b <- round(cult.beta[step,i], 3)
    cult.mean <- round(beta.a/(beta.a +beta.b), 3)
    cult.var <- round(beta.a *beta.b/((beta.a +beta.b)^2 *(beta.a +beta.b +1)), 6)
    cultparams[(i*4-3):(i*4), 2+step] <- c(beta.a, beta.b, cult.mean, cult.var)
  }
   
  results[step,] <- c(step, cult.mode, interaction.ratio, interaction.strength, init_alpha_beta,   
                   effect_G, effect_E, effect_C, effect_EC, effect_GEC,
                   effect_I, effect_full,  h2, round(mean(phenotype.IQformat),1),
                   round(var(phenotype.IQformat),1), round(skewness(phenotype.IQformat),3), 
                    round(kurtosis(phenotype.IQformat),3), N, G_vars, E_vars, C_vars)
  
  sortedagents <- cbind(sort(phenotype.IQformat, decreasing=T, index.return=T)$x[c(1:ranking_length,(N-ranking_length+1):N)], 
                  sort(phenotype.IQformat, decreasing=T, index.return=T)$ix[c(1:50,(N-49):N)])
  rankings[which(rankings$cult.steps==step),(3:dim(rankings)[2])] <- round(venus[sortedagents[,2],],3)

  # # reformatting full model parameters onto a matrix (matrix diagonal is main effects)  
  # qmat <- matrix(NA , GEC_vars, GEC_vars)
  # qoo <- round(model_full$coef[2:(dim(model_full$coef)[1]),1], 3)
  # diag(qmat) <- qoo[1:GEC_vars]
  # GxE_ind <- GEC_vars+G_vars*(E_vars+C_vars)
  # qq <- matrix(qoo[(GEC_vars+1):GxE_ind],G_vars,E_vars+C_vars)
  # qr <- t(qq)
  # qmat[ind_G, ind_E] <- qr[seq(1, (E_vars+C_vars) ,by=2),]
  # qmat[ind_G, ind_C] <- qr[seq(1, (E_vars+C_vars) ,by=2)+1,]
  # qmat[ind_E, ind_C] <- t(matrix(qoo[(GxE_ind+1):length(qoo)],C_vars,E_vars))
  
  # paste model parameters into cumulative matrices
  modelparam_list_full[step, 2:(GEC_vars+1)] <- t(round(model_full$coef,3)[2:(GEC_vars+1)])
  modelparam_list_GEC[step, 2:(GEC_vars+1)] <- t(round(model_GEC$coef,3)[2:(GEC_vars+1)])
  ####modelparam_matrix[which(modelparam_matrix[,1]==step),2:(GEC_vars+1)] <- qmat
  
  print(paste0("iteration: ",step, " out of ", cult.steps))

  
} # MAIN LOOP END

#@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# prepare both ground-truth effects and model estimates for output 
maineffects.out <- data.frame(t(round(data.frame(maineffects),3))) 
colnames(maineffects.out) <- names(venus)[1:(dim(venus)[2]-1)]
alleffects.beta <- round(data.frame(interactions),3) 
diag(alleffects.beta) <- round(maineffects,3)
colnames(alleffects.beta) <- names(venus)[1:(dim(venus)[2]-1)]
alleffects.beta$names <- names(venus)[1:(dim(venus)[2]-1)]
alleffects.beta$vars_influence <- round(maineffects + colSums(interactions.double), 3)

colnames(modelparam_list_full) <- c("cult.step", names(venus)[1:(dim(venus)[2]-1)])
colnames(modelparam_list_GEC) <- c("cult.step", names(venus)[1:(dim(venus)[2]-1)])
colnames(modelparam_matrix) <- c("cult.step", names(venus)[1:(dim(venus)[2]-1)])

time <- as.POSIXlt(Sys.time())
timeB <- paste0("-",time$hour,"h", time$min,"m")

# output files
write.csv(results, paste0(filestring, timeB, "-summary.csv"))
write.csv(rankings, paste0(filestring, timeB,"-rankings.csv"))
###write.csv(maineffects.out, paste0(filestring, timeB, "-groundtruth-maineffects.csv"))
write.csv(alleffects.beta, paste0(filestring, timeB, "-trueeffects.csv"))
write.csv(modelparam_list_full, paste0(filestring, timeB, "-maineffects-fullmodel.csv"))
write.csv(modelparam_list_GEC, paste0(filestring, timeB, "-maineffects-GECmodel.csv"))
####write.csv(modelparam_matrix, paste0(filestring, timeB, "-alleffects-fullmodel.csv"))
write.csv(cultparams, paste0(filestring, timeB, "-cultparams.csv"))








