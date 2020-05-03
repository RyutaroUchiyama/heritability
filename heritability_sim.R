


rm(list = ls())
# =*=*=*=*=*=*=*=*=*=*=*
library(tidyverse)
library(ggplot2)
library(moments)


N <- 5000
G_cols <- 6
E_cols <- 3
C_cols <- 3
P_cols <- 1
ranking_topK <- 50  
pheno.names <- c("smart")
gen_scale <- 1 # Genetic main effects get multiplied by (this value)
cult.mode <- "sel_dir"  # { "sel_pos", "sel_neg","sel_stab", "sel_ent", "sel_div" } 
cult.steps <- 8
filestring <- "0502-fullinteractionmatrix"

results.variables <- c("step","cult.mode", "beta.paramA", "beta.paramB", "cult.mean", "cult.var", 
                       "gene.effect","eco.effect","cult.effect",
                       "EandC.effect","GandEandC.effect","interaction.effect","full.effect","h2", 
                       "P.IQmean", "P.IQvar","P.IQskewness","P.IQkurtosis","N",
                       "variables_G","variables_E","variables_C")

results <- data.frame(matrix(NA,1:cult.steps,length(results.variables))) 
names(results) <- results.variables

rankings <- data.frame(matrix(NA, ranking_topK *2 *cult.steps, G_cols+E_cols+C_cols+P_cols+2)) 
names(rankings)[1:2] <- c("cult.steps","rank")
rankings$cult.steps <- rep(1:cult.steps, each=ranking_topK*2)
rankings$rank <- rep(c(1:ranking_topK,(N-ranking_topK+1):N), cult.steps)

step <- 1 
GEC_cols <- G_cols + E_cols + C_cols
ind_G <- 1:G_cols
ind_E <- (G_cols+1):(G_cols+E_cols)
ind_C <- (G_cols+E_cols+1):GEC_cols
ind_GC <- c(1:G_cols, (G_cols+E_cols+1):GEC_cols)
ind_GE <- 1:(G_cols+E_cols)
ind_EC <- (G_cols+1):GEC_cols

# setting up parameters of the generative model
###mainarray <- seq(1:G_cols)^-1.5
maineffects.beta <- rep(0,GEC_cols) #betas
maineffects.beta[ind_G] <- rnorm(G_cols) *gen_scale
maineffects.beta[ind_E] <- rnorm(E_cols)  
maineffects.beta[ind_C] <- rnorm(C_cols)

interactions.beta <- matrix(0, GEC_cols, GEC_cols )
interactions.beta[upper.tri(interactions.beta)] <- rnorm(sum(upper.tri(interactions.beta)))
interactions.beta <- interactions.beta/(GEC_cols-1) # reducing magnitude of interactions according to number of interactions (balancing 1:1 with main effects)
interactions.beta.double <- (interactions.beta + t(interactions.beta))

cult.sign <- sign( maineffects.beta[ind_C] + colSums(interactions.beta.double[,ind_C])) # which direction do the cult distribuitions need to move
cult.sign <- cult.sign*-1/2+1.5 # for indexing. 1 is for positive slopes, 2 is for negative (d*pheno/d*culturallevel)
# sel_dir <- cbind(2^seq(3,6,length.out=cult.steps),2^seq(3,0,length.out=cult.steps)) #positive directional
# sel_stab <- cbind(2^seq(3,6,length.out=cult.steps),2^seq(3,6,length.out=cult.steps)) #stabilizing
# sel_ent <-  cbind(2^seq(3,0,length.out=cult.steps),2^seq(3,0,length.out=cult.steps)) #entropic
# sel_div <- cbind(2^seq(3,-2,length.out=cult.steps),2^seq(3,-2,length.out=cult.steps)) #diversifying

sel_dir <- c(98,1); sel_stab <- c(1250,1250); sel_ent <- (2,2); 
sel <- eval(parse(text=cult.mode))
 
cult.alpha <- matrix(NA, cult.steps, C_cols)
cult.beta <- matrix(NA, cult.steps, C_cols)
for (i in C_cols){

    cult.alpha[,i] <-  exp(seq(log(12), log(sel_dir),length.out=cult.steps))
    cult.beta[,i] <-   exp(seq(log(12), log(sel_dir),length.out=cult.steps))

    
        cult.alpha[,i] <-  exp(seq(12, 1,length.out=cult.steps))
    cult.beta[,i] <-  exp(seq(12, 98,length.out=cult.steps)) }
  
  }

sel_dir <-  cbind(2^seq(3, 3+3*cult.sign[i] ,length.out=cult.steps),
                  2^seq(3, 3+3*-cult.sign[i] ,length.out=cult.steps))


# when alpha=beta=12= exp(2.4849), variance is at 0.01 (initial distribution)
# when alpha=beta=1250 = exp(7.1305) , variance is at 0.0001
# when alpha=beta= 2 = exp(0.6931) , variance is at 0.05
# when alpha=98=exp(4.5849) & beta=1, variance is at 0.0001 
# when alpha=beta=0.75 , variance is at 0.1

ll <- cbind(2^seq(3,6,length.out=cult.steps),2^seq(3,0,length.out=cult.steps)) 

k <- 8
h <- rbeta(N, ll[k,1], ll[k,2])
hist(h, breaks=70)


[,1]     [,2]
[1,]  8.00000 8.000000
[2,] 10.76720 5.943977
[3,] 14.49158 4.416358
[4,] 19.50422 3.281341
[5,] 26.25073 2.438027
[6,] 35.33086 1.811447
[7,] 47.55182 1.345900
[8,] 64.00000 1.000000

h <- rbeta(N, 2^12 , 2^12 )


step <- 1
beta.a <- round(cultparams[step,1], 3); beta.b <- round(cultparams[step,2], 8)

beta.a <- 2^8; beta.b <- 2^8
cult.mean <- round(beta.a/(beta.a +beta.b), 3)
cult.var <- round(beta.a *beta.b/((beta.a +beta.b)^2 *(beta.a +beta.b +1)), 8)
print(paste(cult.mean, cult.var))
0.015
0.5


# xxvar <- function(beta.a, beta.b) {round(beta.a *beta.b/((beta.a +beta.b)^2 *(beta.a +beta.b +1)), 8) }
# xxvar <- function(beta.a) {round(beta.a *beta.a/((beta.a +beta.a)^2 *(beta.a +beta.a +1)), 8) - 0.5  }
# xxvar(200)
# h <- rbeta(N, 90, 1  )
# hist(h, breaks=200, xlim=c(0,1))
# uniroot(xxvar, lower=0.7, upper=0.8)$root
# when alpha=beta=12= exp(2.4849), variance is at 0.01
# when alpha=beta=1250 = exp(7.1305) , variance is at 0.0001
# when alpha=beta= 2 = exp(0.6931) , variance is at 0.05
# when alpha=98=exp(4.5849) & beta=1, variance is at 0.0001
# when alpha=beta=0.75 , variance is at 0.1



- 3e-05


cultparams <- eval(parse(text=cult.mode))

modelparam_list_full <- data.frame(matrix(NA,cult.steps,GEC_cols+1)) 
modelparam_list_full[,1] <- 1:cult.steps
modelparam_list_GEC <- data.frame(matrix(NA,cult.steps,GEC_cols+1)) 
modelparam_list_GEC[,1] <- 1:cult.steps 
modelparam_matrix <- matrix(NA , GEC_cols*cult.steps, GEC_cols+1)
modelparam_matrix[,1] <- rep(1:cult.steps, each=GEC_cols)

#@@@@@@ MAIN LOOP START @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
for (step in 1:cult.steps){
 
  # initial scores for each genetic, ecological, and cultural factor
  venus_G <- matrix(0, N, G_cols)
  for(k in 1:G_cols){
#    a1 <- sample(1:9,1); a2 <- (10-a1) # All genes have 2 alleles, integer ratio between them are randomly selected
#    venus_G[,k] <- sample(c(rep(0,a1/10*N), rep(1,a2/10*N)))
     venus_G[,k] <- sample(c(rep(0,0.5*N), rep(1,0.5*N)))
      }   
  venus_E <- matrix(0, N, E_cols)
  venus_E <- apply(venus_E, 2, function(x) x + rbeta(N, cultparams[1,1], cultparams[1,2])) # assuming eco distribuition is same as initial dist for culture
  venus_C <- matrix(0, N, C_cols)
  venus_C <- apply(venus_C, 2, function(x) x + rbeta(N, cultparams[step,1], cultparams[step,2])) # cultural traits
  venus_P <- matrix(rep(0,N), N, P_cols)
  venus <- data.frame(cbind(venus_G,venus_E,venus_C,venus_P))
  
  # assigning names
  names(venus)[(GEC_cols+1):(GEC_cols+P_cols)] <- pheno.names
  venus[,1:(G_cols+E_cols)] <- round(venus[,1:(G_cols+E_cols)], 2)
  for(i in 1:G_cols) { names(venus)[i] <- paste0("gene_",as.character(i)) }
  for(i in 1:E_cols) { names(venus)[i+G_cols] <- paste0("eco_",as.character(i)) }
  for(i in 1:C_cols) { names(venus)[i+G_cols+E_cols] <- paste0("cult_",as.character(i)) }
  names(rankings)[3:dim(rankings)[2]] <- names(venus)
  
  
#@@@@@@ INTERACTION EFFECTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # each gene has N=(length(eco)+length(cult)) interactions;
  # each eco factor has N=length(cult) interactions (eco X gene is redundant with above); 10 eco, 100 interactions
  interactions.wg.rowsums <- matrix(NA,N,GEC_cols)
  for (agent in 1:N){
    interactions.wg <- matrix(0, GEC_cols, GEC_cols )
    for (i in 1:GEC_cols){ 
      for (j in 1:GEC_cols){ # {genes X eco} and {genes X culture}
        interactions.wg[i,j] <- venus[agent,i] * venus[agent,j] * interactions.beta[i,j] }}
    interactions.wg.rowsums[agent,] <- rowSums(interactions.wg, na.rm=T)
    if (agent %% 5000 == 0 ) {print(agent)}
    }
  maineffects.wg <- sweep(venus[,1:GEC_cols], MARGIN=2, maineffects.beta, '*')
  phenotype <-  rowSums(maineffects.wg) + rowSums(interactions.wg.rowsums) 
  
  if (step==1){ 
    Pscale.sd <- sd(phenotype)
    Pscale.mean <- mean(phenotype/Pscale.sd *15) # taking mean only after rescaling variability
  }
  phenotype.IQformat <- (phenotype/Pscale.sd *15) - Pscale.mean + 100
  venus[,GEC_cols+1] <- phenotype.IQformat
  
  venus.st <- data.frame(matrix(NA,dim(venus)[1],dim(venus)[2]))
  venus.st[,1:5] <- scale(venus[,1:5])
  venus.st[,2:dim(venus)[2]] <- scale(venus[,2:dim(venus)[2]]) 
  colnames(venus.st) <- colnames(venus)
  
#@@@@@@ GENERATE REGRESSION MODELS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
  
  beta.a <- round(cultparams[step,1], 3); beta.b <- round(cultparams[step,2], 3)
  cult.mean <- round(beta.a/(beta.a +beta.b), 3)
  cult.var <- round(beta.a *beta.b/((beta.a +beta.b)^2 *(beta.a +beta.b +1)), 3)

  #  info <- rep(NA,length(cult_freq))      
#  for (i in 1:length(cult_freq)){ 
#    info[i] <- (cult_freq[i]/N) * log2(cult_freq[i]/N)}
#  entropy <- round(-sum(info, na.rm=T), 3)
  
  results[step,] <- c(step, cult.mode, beta.a, beta.b, cult.mean, cult.var,  
                   effect_G, effect_E, effect_C, effect_EC, effect_GEC,
                   effect_I, effect_full,  h2, round(mean(phenotype.IQformat),1),
                   round(var(phenotype.IQformat),1), round(skewness(phenotype.IQformat),3), 
                    round(kurtosis(phenotype.IQformat),3), N, G_cols, E_cols, C_cols)
  
  sortedagents <- cbind(sort(phenotype.IQformat, decreasing=T, index.return=T)$x[c(1:ranking_topK,(N-ranking_topK+1):N)], 
                  sort(phenotype.IQformat, decreasing=T, index.return=T)$ix[c(1:50,(N-49):N)])
  rankings[which(rankings$cult.steps==step),(3:dim(rankings)[2])] <- round(venus[sortedagents[,2],],3)

  # reformatting full model parameters onto a matrix (matrix diagonal is main effects)  
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
  modelparam_list_full[step, 2:(GEC_cols+1)] <- t(round(model_full$coef,3)[2:(GEC_cols+1)])
  modelparam_list_GEC[step, 2:(GEC_cols+1)] <- t(round(model_GEC$coef,3)[2:(GEC_cols+1)])
  modelparam_matrix[which(modelparam_matrix[,1]==step),2:(GEC_cols+1)] <- qmat

  print(paste0("iteration: ",step, " out of ", cult.steps))

  
} # MAIN LOOP END

#@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# prepare both ground-truth effects and model estimates for output 
maineffects.beta.out <- data.frame(t(round(data.frame(maineffects.beta),3))) 
colnames(maineffects.beta.out) <- names(venus)[1:(dim(venus)[2]-1)]
alleffects.beta <- round(data.frame(interactions.beta),3) 
diag(alleffects.beta) <- round(maineffects.beta,3)
colnames(alleffects.beta) <- names(venus)[1:(dim(venus)[2]-1)]
alleffects.beta$cult_sums <- round(rowSums(alleffects.beta[,ind_C]),3)

colnames(modelparam_list_full) <- c("cult.step", names(venus)[1:(dim(venus)[2]-1)])
colnames(modelparam_list_GEC) <- c("cult.step", names(venus)[1:(dim(venus)[2]-1)])
colnames(modelparam_matrix) <- c("cult.step", names(venus)[1:(dim(venus)[2]-1)])

time <- as.POSIXlt(Sys.time())
timeB <- paste0("-",time$hour,"h", time$min,"m")

# output files
write.csv(results, paste0(filestring, timeB, "-summary.csv"))
write.csv(rankings, paste0(filestring, timeB,"-rankings.csv"))
###write.csv(maineffects.beta.out, paste0(filestring, timeB, "-groundtruth-maineffects.csv"))
write.csv(alleffects.beta, paste0(filestring, timeB, "-groundtruth-alleffects.csv"))
write.csv(modelparam_list_full, paste0(filestring, timeB, "-modelparam-list-fullmodel.csv"))
write.csv(modelparam_list_GEC, paste0(filestring, timeB, "-modelparam-list-GECmodel.csv"))
write.csv(modelparam_matrix, paste0(filestring, timeB, "-modelparam-matrix-fullmodel.csv"))







