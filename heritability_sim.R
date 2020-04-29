

# =*=*=*=*=*=*=*=*=*=*=*
library(tidyverse)
library(ggplot2)

coh <- 1 ####
N <- 50000
meanscore <- 10
#var.gene <- 5
var.eco <- 1
G_cols <- 3
E_cols <- 3
C_cols <- 3
P_cols <- 1
culturalvars <- 8
pheno.names <- c("smart")
cult.coherence <- seq(-6,6,by=0.5)  #  #c(-6,-4,-2,0,2,4,6)
filestring <- "0429-N50000"

results.variables <- c("cult.coherence","entropy","gene.effect","eco.effect","cult.effect",
                       "EandC.effect","GandEandC.effect","interaction.effect","full.effect",
                    "phenotype.var","h2","cult_freq-1","cult_freq-2","cult_freq-3","cult_freq-4",
                    "cult_freq-5","cult_freq-6","cult_freq-7","cult_freq-8","N","meanscore",
                    "var.gene")
results <- data.frame(matrix(NA,length(cult.coherence),length(results.variables))) # row number = cult.coherence levels * number of pheno traits 
names(results) <- results.variables

rankings <- data.frame(matrix(NA,20*length(cult.coherence),G_cols+E_cols+C_cols+P_cols+2)) 
names(rankings)[1:2] <- c("cult.coherence","rank")
rankings$cult.coherence <- rep(cult.coherence, each=20)
rankings$rank <- rep(c(1:10,(N-9):N),length(cult.coherence))

x <- 1 # initializing counter for 'results' matrix
GEC_cols <- G_cols + E_cols + C_cols
ind_G <- 1:G_cols
ind_E <- (G_cols+1):(G_cols+E_cols)
ind_C <- (G_cols+E_cols+1):GEC_cols
ind_GC <- c(1:G_cols, (G_cols+E_cols+1):GEC_cols)
ind_GE <- 1:(G_cols+E_cols)
ind_EC <- (G_cols+1):GEC_cols

mainarray <- seq(1:G_cols)^-1.5
maineffects.beta <- rep(0,GEC_cols) #betas
maineffects.beta[ind_G] <- sample(mainarray/mean(mainarray)) *(meanscore/4)   #main effects of genes.  maineffects.beta[ind_E] <- sample(mainarray/mean(mainarray)) #main effects of eco-factors. 
maineffects.beta[ind_E] <- sample(mainarray/mean(mainarray)) ##main effects of eco factors
maineffects.beta[ind_C] <- rep(1,C_cols) ##main effects of cultural traits
interactions.beta <- matrix(0, GEC_cols, GEC_cols )
interactions.beta[ind_G, ind_EC] <- matrix(rnorm(G_cols*(E_cols+C_cols),0,0.1),G_cols,E_cols+C_cols) 
interactions.beta[ind_E, ind_C] <- matrix(rnorm(E_cols*C_cols,0,0.1),E_cols,C_cols) # populating upper triangle

modelparam_list_full <- data.frame(matrix(NA,length(cult.coherence),GEC_cols+1)) 
modelparam_list_full[,1] <- cult.coherence
modelparam_list_GEC <- data.frame(matrix(NA,length(cult.coherence),GEC_cols+1)) 
modelparam_list_GEC[,1] <- cult.coherence 
modelparam_matrix <- matrix(NA , GEC_cols*length(cult.coherence), GEC_cols+1)
modelparam_matrix[,1] <- rep(cult.coherence, each=GEC_cols)

#@@@@@@ MAIN LOOP START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
for (coh in cult.coherence){
 
  # initial scores for each genetic, ecological, and cultural factor
  venus_G <- matrix(0, N, G_cols)
  for(k in 1:G_cols){
    a1 <- sample(1:9,1); a2 <- (10-a1) # All genes have 2 alleles, integer ratio between them are randomly selected
    venus_G[,k] <- sample(c(rep(0,a1/10*N), rep(1,a2/10*N)))
#    a1 <- sample(1:8,1); a2 <- sample(1:(10-a1-1),1); a3 <- (10-a1-a2) # All genes have 3 alleles, integer ratio between them are randomly selected
#    venus_G[,k] <- sample(c(rep(1,a1/10*N), rep(2,a2/10*N), rep(3,a3/10*N)))
      }   
  venus_E <- matrix(0, N, E_cols)
  venus_E <- apply(venus_E, 2, function(x) x + rnorm(N, meanscore, var.eco))
  venus_CP <- matrix(rep(0,N), N, C_cols+P_cols)
  venus <- data.frame(cbind(venus_G,venus_E,venus_CP))

  # assigning names
  names(venus)[(GEC_cols+1):(GEC_cols+P_cols)] <- pheno.names
  venus[,1:(G_cols+E_cols)] <- round(venus[,1:(G_cols+E_cols)], 2)
  for(i in 1:G_cols) { names(venus)[i] <- paste0("gene_",as.character(i)) }
  for(i in 1:E_cols) { names(venus)[i+G_cols] <- paste0("eco_",as.character(i)) }
  for(i in 1:C_cols) { names(venus)[i+G_cols+E_cols] <- paste0("cult_",as.character(i)) }
  names(rankings)[3:dim(rankings)[2]] <- names(venus)
  
#@@@@@@ CULTURAL TRAITS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # setting distribution of cultural traits
  cult_freq <- seq(1:culturalvars) ^ coh
  cult_freq <- cult_freq/sum(cult_freq)
  cult_freq <- round(cult_freq*N)
  
  # if round() pushes it above or below total number of agents, this fixes it by adding or subtracting from highest freq variant
  if(sum(cult_freq)< N){
    cult_freq[which.max(cult_freq)] <- cult_freq[which.max(cult_freq)] + (N-sum(cult_freq))
  } else if(sum(cult_freq)> N){   
    cult_freq[which.max(cult_freq)] <- cult_freq[which.max(cult_freq)] - (sum(cult_freq)-N)}  
  
  # generate array of cultural variants to match cult_freq
  cult_array <- rep(NA,N)
  for(i in 1:culturalvars){
    if(cult_freq[i]>0){
      unpop <- which(is.na(cult_array)) # subset just the slots that are unpopulated
      cult_array[unpop][1:cult_freq[i]] <- rep(i, cult_freq[i]) }}  
  
  for(i in ind_C){ venus[,i] <- sample(cult_array) } # randomized versions of the cultural array into  
  
#@@@@@@ INTERACTION EFFECTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # each gene has N=(length(eco)+length(cult)) interactions;
  # each eco factor has N=length(cult) interactions (eco X gene is redundant with above); 10 eco, 100 interactions
  interactions.wg.rowsums <- matrix(NA,N,GEC_cols)
  for (agent in 1:N){
    interactions.wg <- matrix(0, GEC_cols, GEC_cols )
    for (i in ind_G){ 
      for (j in ind_EC){ # {genes X eco} and {genes X culture}
    interactions.wg[i,j] <- venus[agent,i] * venus[agent,j] * interactions.beta[i,j] }}
    for (i in ind_E){ 
      for (j in ind_C){ # {eco X culture} 
    interactions.wg[i,j] <- venus[agent,i] * venus[agent,j] * interactions.beta[i,j] }}
#    interactions.wg <- (interactions.wg + t(interactions.wg))/2 # fill in lower triangle and HALVE (because each interaction will be counted twice ([gene1,eco4] and [eco4,gene1]) when doing rowSums. Could have also just summed the upper or lower triangle)
    interactions.wg.rowsums[agent,] <- rowSums(interactions.wg, na.rm=T)
    if (agent %% 10000 == 0 ) {print(agent)}
    }
  maineffects.wg <- sweep(venus[,1:GEC_cols], MARGIN=2, maineffects.beta, '*')
    #### ◊◊◊◊◊◊ how should I scale magnitudes of the interactions and the main effects?? Need to???
  
  phenotype <-  rowSums(maineffects.wg) + rowSums(interactions.wg.rowsums) # WHAT HAPPENS WITH WEIGHTS THAT DEVIATE FROM 1:1??
  venus[,GEC_cols+1] <- phenotype

#@@@@@@ GENERATE MODELS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # generate strings to enter into lm()
  string_G <- character()
  string_E <- character()
  string_C <- character()
  for (i in 1:G_cols){ string_G <- str_c(string_G, paste0("gene_", as.character(i),"+ ")) }
  for (i in 1:E_cols){ string_E <- str_c(string_E, paste0("eco_", as.character(i),"+ ")) }
  for (i in 1:C_cols){ string_C <- str_c(string_C, paste0("cult_", as.character(i),"+ ")) }
  string_G <- substr(string_G, 1, nchar(string_G)-2) # shave off the last "+ " in string
  string_E <- substr(string_E, 1, nchar(string_E)-2) 
  string_C <- substr(string_C, 1, nchar(string_C)-2) 
  
  string_I <- character() #string of all interactions
  for (i in 1:G_cols){ 
    for (j in 1:E_cols){string_I <- str_c(string_I, paste0("gene_",as.character(i),"*eco_",as.character(j),"+ "))}
    for (k in 1:C_cols){string_I <- str_c(string_I, paste0("gene_",as.character(i),"*cult_",as.character(k),"+ "))}}
  for (j in 1:E_cols){
    for (k in 1:C_cols){string_I <- str_c(string_I, paste0("eco_",as.character(j),"*cult_",as.character(k),"+ "))}}
  string_I <- substr(string_I, 1, nchar(string_I)-2) # shave off the last "+ " in string
  
  # generate full string for evaluating model
  modelstring_G <- paste0("smart ~ ", string_G) 
  modelstring_E <- paste0("smart ~ ", string_E)
  modelstring_C <- paste0("smart ~ ", string_C)
  modelstring_EC <- paste0("smart ~ ", string_E,"+",string_C)
  modelstring_GEC <- paste0("smart ~ ", string_G,"+",string_E,"+",string_C)
  modelstring_I <- paste0("smart ~ ", string_I)
  modelstring_full <- paste0("smart ~ ", string_G,"+",string_E,"+",string_C,"+",string_I)
    
  lm_G <- lm( eval(parse(text=modelstring_G)), data=venus)
  lm_E <- lm( eval(parse(text=modelstring_E)), data=venus)
  lm_C <- lm( eval(parse(text=modelstring_C)), data=venus)
  lm_EC <- lm( eval(parse(text=modelstring_EC)), data=venus)
  lm_GEC <- lm( eval(parse(text=modelstring_GEC)), data=venus)
  lm_I <- lm( eval(parse(text=modelstring_I)), data=venus)
  lm_full <- lm( eval(parse(text=modelstring_full)), data=venus)
  
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
  info <- rep(NA,length(cult_freq))      
  for (i in 1:length(cult_freq)){ 
    info[i] <- (cult_freq[i]/N) * log2(cult_freq[i]/N)}
  entropy <- round(-sum(info, na.rm=T), 3)
  
  results[x,] <- c(coh, entropy, effect_G, effect_E, effect_C, effect_EC, effect_GEC,
                       effect_I, effect_full, (effect_G + effect_E + effect_C) , h2,
                       cult_freq[1],cult_freq[2],cult_freq[3],cult_freq[4],cult_freq[5],cult_freq[6],
                       cult_freq[7],cult_freq[8], N, meanscore, var.eco)
  
  sortedagents <- cbind(sort(phenotype, decreasing=T, index.return=T)$x[c(1:10,(N-9):N)], 
                  sort(phenotype, decreasing=T, index.return=T)$ix[c(1:10,(N-9):N)])
  rankings[which(rankings$cult.coherence==coh),(3:dim(rankings)[2])] <- round(venus[sortedagents[,2],],3)
  
  
  # reformatting model parameters onto a matrix (matrix diagonal is main effects)  
  qmat <- matrix(NA , GEC_cols, GEC_cols)
  qoo <- round(model_full$coef[2:(dim(model_full$coef)[1]),1], 3)
  diag(qmat) <- qoo[1:GEC_cols]
  GxE_ind <- GEC_cols+G_cols*(E_cols+C_cols)
  qq <- matrix(qoo[(GEC_cols+1):GxE_ind],G_cols,E_cols+C_cols)
  qr <- t(qq)
  qmat[ind_G, ind_E] <- qr[seq(1, (E_cols+C_cols) ,by=2),]
  qmat[ind_G, ind_C] <- qr[seq(1, (E_cols+C_cols) ,by=2)+1,]
  qmat[ind_E, ind_C] <- t(matrix(qoo[(GxE_ind+1):length(qoo)],C_cols,E_cols))
  
  # paste model parameters into cumulative matrices
  modelparam_list_full[coh,2:(GEC_cols+1)] <- data.frame(t(round(model_full$coef,3)[2:(GEC_cols+1)]))
  modelparam_list_GEC[coh,2:(GEC_cols+1)] <- data.frame(t(round(model_GEC$coef,3)[2:(GEC_cols+1)]))
  modelparam_matrix[which(modelparam_matrix[,1]==coh),2:(GEC_cols+1)] <- qmat
  
  print(paste0("iteration: ",x, " out of ",length(cult.coherence)))
  x <- x+1
} # MAIN LOOP END

#@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# prepare both ground-truth effects and model estimates for output 
maineffects.beta <- data.frame(t(round(data.frame(maineffects.beta),3))) 
colnames(maineffects.beta) <- names(venus)[1:(dim(venus)[2]-1)]
alleffects.beta <- round(data.frame(interactions.beta),3) 
diag(alleffects.beta) <- as.numeric(maineffects.beta)
colnames(alleffects.beta) <- names(venus)[1:(dim(venus)[2]-1)]

colnames(modelparam_list_full) <- c("cult.coherence", names(venus)[1:(dim(venus)[2]-1)])
colnames(modelparam_list_GEC) <- c("cult.coherence", names(venus)[1:(dim(venus)[2]-1)])
colnames(modelparam_matrix) <- c("cult.coherence", names(venus)[1:(dim(venus)[2]-1)])

# output files
write.csv(results, paste0(filestring, "-summary.csv"))
write.csv(rankings, paste0(filestring,"-rankings.csv"))
write.csv(maineffects.beta, paste0(filestring, "-groundtruth-maineffects.csv"))
write.csv(alleffects.beta, paste0(filestring, "-groundtruth-alleffects.csv"))
write.csv(modelparam_list_full, paste0(filestring, "-modelparam-list-fullmodel.csv"))
write.csv(modelparam_list_GEC, paste0(filestring, "-modelparam-list-GECmodel.csv"))
write.csv(modelparam_matrix, paste0(filestring, "-modelparam-matrix-fullmodel.csv"))


