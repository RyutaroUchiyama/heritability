

# =*=*=*=*=*=*=*=*=*=*=*
library(tidyverse)
library(ggplot2)

coh <- 2 ####
N <- 20000
meanscore <- 40
#var.gene <- 5
var.eco <- 5
G_cols <- 3
E_cols <- 3
C_cols <- 3
P_cols <- 1
culturalvars <- 8
pheno.names <- c("smart")
cult.coherence <- c(-6:10,-0.1,0.1,0.2,0.4)

results.variables <- c("cult.coherence","entropy","gene.effect","eco.effect","cult.effect",
                       "EandC.effect","GandEandC.effect","interaction.effect","full.effect",
                    "phenotype.var","h2","cult_freq-1","cult_freq-2","cult_freq-3","cult_freq-4",
                    "cult_freq-5","cult_freq-6","cult_freq-7","cult_freq-8","N","meanscore",
                    "var.gene")
results <- data.frame(matrix(NA,length(cult.coherence),length(results.variables))) # row number = cult.coherence levels * number of pheno traits 
names(results) <- results.variables

x <- 1 # initializing counter for 'results' matrix
GEC_cols <- G_cols + E_cols + C_cols
ind_G <- 1:G_cols
ind_E <- (G_cols+1):(G_cols+E_cols)
ind_C <- (G_cols+E_cols+1):GEC_cols
ind_GC <- c(1:G_cols, (G_cols+E_cols+1):GEC_cols)
ind_GE <- 1:(G_cols+E_cols)
ind_EC <- (G_cols+1):GEC_cols

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
  names(venus)[(GEC_cols+1):(GEC_cols+P_cols)] <- pheno.names
  
  # assigning names
  venus[,1:(G_cols+E_cols)] <- round(venus[,1:(G_cols+E_cols)], 2)
  for(i in 1:G_cols) { names(venus)[i] <- paste0("gene_",as.character(i)) }
  for(i in 1:E_cols) { names(venus)[i+G_cols] <- paste0("eco_",as.character(i)) }
  for(i in 1:C_cols) { names(venus)[i+G_cols+E_cols] <- paste0("cult_",as.character(i)) }
  
  
  #€€€€€€€€€€€€€€€€€
  
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
  ####++++ 
  mainarray <- seq(1:G_cols)^-1.5
  maineffects.beta <- rep(0,GEC_cols) #betas
  maineffects.beta[ind_G] <- sample(mainarray/mean(mainarray))*20 #main effects of genes. Multiplying by 20 so that expected distribution mean comes out to 10.
  maineffects.beta[ind_E] <- sample(mainarray/mean(mainarray))*(10/meanscore) #main effects of eco-factors. Scaling down by (10/meanscore) so that expected distribution mean comes out to 10 (coordinating magnitude with genes)
  #maineffects.beta[ind_C] <- sample(mainarray/mean(mainarray)) ##main effects of cultural traits

  interactions.beta <- matrix(0, GEC_cols, GEC_cols )
  interactions.beta[ind_G, ind_E] <- matrix(rnorm(G_cols*(E_cols),0, 0.05),G_cols,E_cols) 
###  interactions.beta[ind_E, ind_C] <- matrix(rnorm(E_cols*C_cols,0, 0.05),E_cols,C_cols) # populating upper triangle
  for (i in 1:length(ind_E)){
    interactions.beta[ind_E[i], ind_C[i]] <- -maineffects.beta[ind_E[i]]/(culturalvars) } 

    
  # each gene has N=(length(eco)+length(cult)) interactions;
  # each eco factor has N=length(cult) interactions (eco X gene is redundant with above); 10 eco, 100 interactions
  interactions.wg.rowsums <- matrix(NA,N,GEC_cols)
  for (agent in 1:N){
    interactions.wg <- matrix(0, GEC_cols, GEC_cols )
    for (i in ind_G){ 
      for (j in ind_E){ # {genes X eco}
    interactions.wg[i,j] <- venus[agent,i] * venus[agent,j] * interactions.beta[i,j] }}
    for (i in ind_E){ 
      for (j in ind_C){ # {eco X culture} 
    interactions.wg[i,j] <- venus[agent,i] * venus[agent,j] * interactions.beta[i,j] }}
#    interactions.wg <- (interactions.wg + t(interactions.wg))/2 # fill in lower triangle and HALVE (because each interaction will be counted twice ([gene1,eco4] and [eco4,gene1]) when doing rowSums. Could have also just summed the upper or lower triangle)
    interactions.wg.rowsums[agent,] <- rowSums(interactions.wg, na.rm=T)
    if (agent %% 1000 == 0 ) {print(agent)}
    }

  maineffects.wg <- sweep(venus[,1:GEC_cols], MARGIN=2, maineffects.beta, '*')
    #### ◊◊◊◊◊◊ how should I scale magnitudes of the interactions and the main effects?? Need to???
  
  
  # standardize
  ###venus.st <- venus
  ###venus.st[,1:G_cols] <- apply(venus.st[,1:G_cols], 2, scale) 
  ###venus.st[,(G_cols+1):(GEC_cols+P_cols)] <- apply(venus.st[,(G_cols+1):(GEC_cols+P_cols)], 2, scale)
  # ^ for some reason, need to separate scale function into two chunks in order to prevent venus.st from collapsing from a data frame to a matrix
  
  phenotype <-  rowSums(maineffects.wg) + rowSums(interactions.wg.rowsums) # WHAT HAPPENS WITH WEIGHTS THAT DEVIATE FROM 1:1??
  venus[,GEC_cols+1] <- phenotype
  venus[,ind_E] <- scale(venus[,ind_E])
  
  #venus[,c(ind_E, GEC_cols+1)] <- scale(venus[,c(ind_E, GEC_cols+1)])
  
  
  # ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # generate strings to enter into regression
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
  modelstring_G <- paste0("summary(lm(smart ~ ", string_G, ", data= venus))" ) 
  modelstring_E <- paste0("summary(lm(smart ~ ", string_E, ", data= venus))" )
  modelstring_C <- paste0("summary(lm(smart ~ ", string_C, ", data= venus))" )
  modelstring_EC <- paste0("summary(lm(smart ~ ", string_E,"+",string_C, ", data= venus))" )
  modelstring_GEC <- paste0("summary(lm(smart ~ ", string_G,"+",string_E,"+",string_C, ", data= venus))" )
  modelstring_I <- paste0("summary(lm(smart ~ ", string_I, ", data= venus))" )
  modelstring_full <- paste0("summary(lm(smart ~ ", string_G,"+",string_E,"+",string_C,"+",string_I, ", data= venus))" )
  
  model_G <- eval(parse(text= modelstring_G)) # eval()
  model_E <- eval(parse(text= modelstring_E))
  model_C <- eval(parse(text= modelstring_C))
  model_EC <- eval(parse(text= modelstring_EC))
  model_GEC <- eval(parse(text= modelstring_GEC))
  model_I <- eval(parse(text= modelstring_I))
  model_full <- eval(parse(text= modelstring_full))
  
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
  
  
  
  
  
  print(x)
  x <- x+1
}


write.csv(results,"output-0427-02.csv")
