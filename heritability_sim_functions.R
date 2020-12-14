



get_interactions <- function(GEC_vars, init_alpha_beta_sd, interaction.ratio, interaction.strength){
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
  interactions <- interactions /((GEC_vars-1) *interaction.ratio) * interaction.strength # first controlling for number of interactions culled by interaction.ratio, then multiplying by interaction.strength. Result: [ sum(abs(interaction terms)) /interaction.strength] is expected to be equal to the absolute magnitude of the main effect   
  interactions
}

get_cult <- function(cult.mode, maineffects, interactions.double, ind_C, C_vars, cult.steps, init_alpha_beta, alpha_beta){
  sel_95tight <- c(800,42); sel_95mid <- c(80,4.2); sel_95loose <- c(8,0.42);
  sel <- eval(parse(text=cult.mode))
  cult.sign <- sign( maineffects[ind_C] + colSums(interactions.double[,ind_C])) # which direction do the cult distribuitions need to move
  cult <- matrix(NA, cult.steps, C_vars)
  if (alpha_beta == "alpha") {
    for (i in 1:C_vars){
      if (cult.sign[i]== 1){
        cult[,i] <-  exp(seq(log(init_alpha_beta), log(sel[1]),length.out=cult.steps))}
      if (cult.sign[i]== -1){
        cult[,i] <-  exp(seq(log(init_alpha_beta), log(sel[2]),length.out=cult.steps))}
    }
  } else if (alpha_beta == "beta") {
    for (i in 1:C_vars){
      if (cult.sign[i]== 1){
        cult[,i] <-   exp(seq(log(init_alpha_beta), log(sel[2]),length.out=cult.steps))}
      if (cult.sign[i]== -1){
        cult[,i] <-   exp(seq(log(init_alpha_beta), log(sel[1]),length.out=cult.steps))}
    }
  }   
  cult
}

get_pop <- function(N,G_vars,E_vars,C_vars,init_alpha_beta,step,phenotype.name){ 
  # initial scores for each genetic, ecological, and cultural factor
  pop_G <- matrix(0, N, G_vars)
  pop_G <- apply(pop_G, 2, function(x) x + rbeta(N, init_alpha_beta, init_alpha_beta)) # eco distribuition is same as initial dist for culture
  pop_E <- matrix(0, N, E_vars)
  pop_E <- apply(pop_E, 2, function(x) x + rbeta(N, init_alpha_beta, init_alpha_beta)) # eco distribuition is same as initial dist for culture
  pop_C <- matrix(0, N, C_vars)
  for (i in 1:C_vars){
    pop_C[,i] <- rbeta(N, cult.alpha[step,i], cult.beta[step,i])} # cultural traits, each distribution moving in the direction of positive phenotypic outcome
  pop_P <- matrix(rep(0,N), N, 1)
  pop <- data.frame(cbind(pop_G,pop_E,pop_C,pop_P))
  # assigning names
  names(pop)[(GEC_vars+1):(GEC_vars+1)] <- phenotype.name
  pop[,1:(G_vars+E_vars)] <- round(pop[,1:(G_vars+E_vars)], 2)
  for(i in 1:G_vars) { names(pop)[i] <- paste0("gene_",as.character(i)) }
  for(i in 1:E_vars) { names(pop)[i+G_vars] <- paste0("eco_",as.character(i)) }
  for(i in 1:C_vars) { names(pop)[i+G_vars+E_vars] <- paste0("cult_",as.character(i)) }
  pop
}



get_interactions.wg.rowsums <- function(N,GEC_vars,interactions ){
  interactions.wg.rowsums <- matrix(NA,N,GEC_vars)
  for (agent in 1:N){
    interactions.wg <- matrix(0, GEC_vars, GEC_vars )
    for (i in 1:GEC_vars){ 
      for (j in 1:GEC_vars){ # {genes X eco} and {genes X culture}
        interactions.wg[i,j] <- pop[agent,i] * pop[agent,j] * interactions[i,j] }}
    interactions.wg.rowsums[agent,] <- rowSums(interactions.wg, na.rm=T)
    if (agent %% 5000 == 0 ) {print(paste0("agent: ",agent))}
  }
  interactions.wg.rowsums
}


get_pop.st <- function(pop){
  # standardized population scores for use in regression models
  pop.st <- data.frame(matrix(NA,dim(pop)[1],dim(pop)[2]))
  pop.st[,1:5] <- scale(pop[,1:5])
  pop.st[,2:dim(pop)[2]] <- scale(pop[,2:dim(pop)[2]]) 
  colnames(pop.st) <- colnames(pop)
  pop.st
}


get_cultparams <- function(cult.alpha, cult.beta, step){
  for (i in 1:C_vars){ # entering properties of each cultural distribution 
    beta.a <- round(cult.alpha[step,i], 3)
    beta.b <- round(cult.beta[step,i], 3)
    cult.mean <- round(beta.a/(beta.a +beta.b), 3)
    cult.var <- round(beta.a *beta.b/((beta.a +beta.b)^2 *(beta.a +beta.b +1)), 6)
    cultparams[(i*4-3):(i*4), 2+step] <- c(beta.a, beta.b, cult.mean, cult.var)
  }
  cultparams 
}
