simulate <-function(cluster1_size, cluster2_size, cluster3_size, min_exp, max_exp, Delta, nb_spec_genes,nb_ns_genes, beta_param){
  #### Simulation ####
  
  ## Building the clusters (spec. genes) ##
  
  #Generate matrix of mean gene expression
  cluster1_spec_mu_matrix = matrix(runif(nb_spec_genes,min_exp,max_exp), byrow = T, cluster1_size , nb_spec_genes )
  cluster2_spec_mu_matrix = matrix(runif(nb_spec_genes,min_exp,max_exp), byrow = T, cluster2_size , nb_spec_genes )
  #Generate sigma vector for specific genes
  sigma_vect = (abs(cluster1_spec_mu_matrix[1,] - cluster2_spec_mu_matrix[1,])/Delta)
  #Generate matrix of variation
  #cluster1_spec_sd_matrix = cbind(sapply(1:nb_spec_genes, function(x) {matrix(rnorm(cluster1_size,0,sigma_vect[x]),cluster1_size,1) }))
  #cluster2_spec_sd_matrix = cbind(sapply(1:nb_spec_genes, function(x) {matrix(rnorm(cluster2_size,0,sigma_vect[x]),cluster2_size,1) }))
  #Merge matrix of median expression and variance
  #cluster1_spec_matrix = cluster1_spec_mu_matrix + cluster1_spec_sd_matrix
  #cluster2_spec_matrix = cluster2_spec_mu_matrix + cluster2_spec_sd_matrix
  
  cluster1_spec_matrix = matrix( rnorm(length(cluster1_spec_mu_matrix) , cluster1_spec_mu_matrix , rep( sigma_vect , each=cluster1_size) ) , cluster1_size, nb_spec_genes )
  cluster2_spec_matrix = matrix( rnorm(length(cluster2_spec_mu_matrix) , cluster2_spec_mu_matrix , rep( sigma_vect , each=cluster2_size) ) , cluster2_size, nb_spec_genes )
  
  ## Building the Intermediates (spec. genes) ##
  
  #Generate vector of mixture coeficient (shiftf) per sample : #mixture coeficient values
  cluster3_shiftf_list = rbeta(1 , beta_param,beta_param , n = cluster3_size)
  #Generate matrix of mean gene expression
  #cluster3_spec_mu_matrix = rbind(sapply(1:nb_spec_genes, function(x) { sapply(1:cluster3_size, function(y) {
  #  mu1 = cluster1_spec_mu_matrix[1,x]; #mean expression of gene x in cluster1
  #  mu2 = cluster2_spec_mu_matrix[1,x]; #mean expression of gene x in cluster2
  #  shift = abs(mu1-mu2)*cluster3_shiftf_list[y]; # mixture coeficient
  #  if (mu1 > mu2) {mu3 = mu1 - shift;};
  #  if (mu1 <= mu2) {mu3 = mu1 + shift;};
  #  mu3; #mean expression of intermediate y
  #} ) } ))
  cluster3_spec_mu_matrix = matrix(1-cluster3_shiftf_list,ncol=1)%*%cluster1_spec_mu_matrix[1,] + matrix(cluster3_shiftf_list,ncol=1)%*%cluster2_spec_mu_matrix[1,]
  #Generate matrix of variation (interm. spec. genes)
  #cluster3_spec_sd_matrix = matrix(rnorm(cluster3_size,0,sigma_vect[x]), byrow = T, cluster3_size, nb_spec_genes )
  #Merge matrix of median expression and variation (interm. spec. genes)
  #cluster3_spec_matrix = cluster3_spec_mu_matrix + cluster3_spec_sd_matrix
  
  cluster3_spec_matrix = matrix( rnorm(length(cluster3_spec_mu_matrix) , cluster3_spec_mu_matrix , rep( sigma_vect , each=cluster3_size) ) , cluster3_size, nb_spec_genes )
  
  
  ##
  
  #total number of samples (number of rows in matrix)
  totalS = cluster1_size + cluster2_size + cluster3_size
  
  ## Adding non specific genes ##
  
  ns_matrix = cbind(sapply(1:nb_ns_genes, function(x) {sigma =runif(1,min_sigma,max_sigma) ; matrix(runif(1,min_exp,max_exp),totalS,1) + matrix(rnorm(totalS,0,sigma),totalS,1)} ))
  sample_spec_matrix = rbind( cluster1_spec_matrix, cluster2_spec_matrix, cluster3_spec_matrix)
  full_matrix = cbind(sample_spec_matrix,ns_matrix)
  full_matrix[full_matrix<0]<-0
  admixt3 <- full_matrix

  returnList <- list( admixt3 = admixt3, cluster3_shiftf_list = cluster3_shiftf_list)
  
  #### end simulation
return(returnList)
}