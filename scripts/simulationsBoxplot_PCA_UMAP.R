###############################################################################
### Script for simulation of high-dimensional (transcriptomic) data         ###
### and assessment of dimension reduction methods for the detection of      ###
### intermediate phenotypes                                                 ###
###                                                                         ###
### Generates PCA and UMAP plots, and boxplot of simulation R squarres      ###
###																			                                    ###
###					by Jerome Poizat										                            ###
###############################################################################

require(ade4)
require(umap)

nb_sim = 5 #number of simulation 
testparem = "Spread" #tested parameter (used as text for plots)
param_vector = c(1,2,5,10) #values of the parameter to test (for boxplot)

#### Constant simulation parameters ####

nb_spec_genes = 500 #number of specific genes
nb_ns_genes = 500 #number of non specific genes
cluster3_prop = 0.33 #proportion of intermediates
data_size = 500
cluster1_size = ceiling((data_size - data_size*cluster3_prop)/2)
cluster2_size = ceiling((data_size - data_size*cluster3_prop)/2)
cluster3_size = ceiling(data_size*cluster3_prop)
if (cluster3_size < 2) {cluster3_size = 2}
min_exp = 0 #minimum gene expression median value
max_exp = 100 #minimum gene expression median value
min_sigma = 1
max_sigma = 12 #possible sigmas for non spec. genes, maybe use uniform distribution instead of sampling ?
Delta = 2
beta_param = 1 #parameters a & b of the beta distribution of the shift factor

# Dimensionality reduction #
custom.config = umap.defaults
custom.config$n_neighbors = 15
custom.config$min_dist = 0.1
#custom.config$spread = 10
custom.config$metric = "euclidean"

#Loop through values of parameter
boxplot_matrix <- sapply(param_vector, simplify = 'array', function(param_value) { 
  
  #### Variable simulation parameter ####
  #nb_spec_genes = param_value 
  #nb_ns_genes = param_value
  #cluster1_size = param_value
  #cluster2_size = param_value
  #cluster3_size = param_value
  #min_exp = param_value 
  #max_exp = param_value
  #sigmas = param_value
  #beta_param = param_value 
  
  #Dimensionality reduction #
  #custom.config$n_neighbors = param_value
  #custom.config$min_dist = param_value
  custom.config$spread = param_value
  
  #### ####
  
  #Build color matrix for simulation plot
  color1_matrix = sapply(1:cluster1_size, function(x) {c(0.9)})  #cluster 1 hue color
  color2_matrix = sapply(1:cluster2_size, function(x) {c(0.55)}) #cluster 2 hue color
  color3m1_matrix = sapply(1:round(cluster3_size/2), function(x) {c(0.1)})  #intermediate hue color
  color3m2_matrix = sapply(1:( cluster3_size - round(cluster3_size/2) ), function(x) {c(0.2)})  #intermediate hue color
  color_matrix = c(color1_matrix,color2_matrix,color3m1_matrix,color3m2_matrix)
  
  #Loop through simulations
  Rsquare_matrix = sapply(1:nb_sim, function(x) {

    #### Simulation ####
    source("D:/Jerome/master bio-info/Docs_stage_IARC/final_work/results/Scripts_and_Plots/simulation2b_script.R")
    simulation <- simulate(cluster1_size, cluster2_size, cluster3_size, min_exp, max_exp, Delta, nb_spec_genes,nb_ns_genes, beta_param)
    admixt3 <- simulation$admixt3
    cluster3_shiftf_list <- simulation$cluster3_shiftf_list
    #### end simulation

    #### Dimensionality reduction ####
    admixt3.pca  = dudi.pca(admixt3,scannf = F,nf = 2,center = T,scale = F)
    admixt3.pca$li[,1]= admixt3.pca$li[,1]*(admixt3.pca$eig[1]/admixt3.pca$eig[2]) #adjust axe scale for variance explained
    admixt3.umap = umap(admixt3, custom.config)

    ####


    #### Vector projections ####
    ## PCA projections ##
    pca_centroid_cluster1 = colMeans(admixt3.pca$li[1:cluster1_size,])
    pca_centroid_cluster2 = colMeans(admixt3.pca$li[cluster1_size+1:cluster2_size,])
    pca_centroid_vector = (pca_centroid_cluster2 - pca_centroid_cluster1)
    pca_centroid_vector_lenght = sqrt((pca_centroid_vector[1]^2)+(pca_centroid_vector[2]^2))
    #PCA normalized projection lenght
    pca_norm_proj_list = rbind(sapply(cluster1_size+cluster2_size+1:cluster3_size, function(x) {
      pca_sample_prevector <- admixt3.pca$li[x,];
      pca_sample_vector = pca_sample_prevector - pca_centroid_cluster1;
      pca_sample_projection = ((sum(pca_sample_vector*pca_centroid_vector))/(pca_centroid_vector_lenght^2))*pca_centroid_vector;
      pca_sample_projection_lenght = sqrt((pca_sample_projection[1]^2)+(pca_sample_projection[2]^2));
      pca_norm_proj = pca_sample_projection_lenght/pca_centroid_vector_lenght;
      #sign correction of norm_proj
      if ( (sign(pca_sample_projection[1]) != sign(pca_centroid_vector[1])) || (sign(pca_sample_projection[2])!=sign(pca_centroid_vector[2])) ) {
        pca_norm_proj <- -pca_norm_proj;
      };
      pca_norm_proj;
    } ))
    #PCA projections (for plot)
    pca_proj_list = rbind(sapply(cluster1_size+cluster2_size+1:cluster3_size, function(x) {
      pca_sample_prevector <- admixt3.pca$li[x,];
      pca_sample_vector = pca_sample_prevector - pca_centroid_cluster1;
      pca_sample_projection = ((sum(pca_sample_vector*pca_centroid_vector))/(pca_centroid_vector_lenght^2))*pca_centroid_vector;
      pca_sample_projection;
    } ))

    
    ## UMAP projections ##
    umap_centroid_cluster1 = colMeans(admixt3.umap$layout[1:cluster1_size,])
    umap_centroid_cluster2 = colMeans(admixt3.umap$layout[cluster1_size+1:cluster2_size,])
    umap_centroid_vector = (umap_centroid_cluster2 - umap_centroid_cluster1)
    umap_centroid_vector_lenght = sqrt((umap_centroid_vector[1]^2)+(umap_centroid_vector[2]^2))
    #UMAP normalized projection lenght
    umap_norm_proj_list = rbind(sapply(cluster1_size+cluster2_size+1:cluster3_size, function(x) {
      umap_sample_prevector <- admixt3.umap$layout[x,];
      umap_sample_vector = (umap_sample_prevector - umap_centroid_cluster1);
      umap_sample_projection = ((sum(umap_sample_vector*umap_centroid_vector))/(umap_centroid_vector_lenght^2))*umap_centroid_vector;
      umap_sample_projection_lenght = sqrt((umap_sample_projection[1]^2)+(umap_sample_projection[2]^2));
      umap_norm_proj = (umap_sample_projection_lenght/umap_centroid_vector_lenght);
      #sign correction of norm_proj
      if ( (sign(umap_sample_projection[1]) != sign(umap_centroid_vector[1])) || (sign(umap_sample_projection[2])!=sign(umap_centroid_vector[2])) ) {
        umap_norm_proj <- -umap_norm_proj;
      };
      umap_norm_proj;
    } ))
    #UMAP projections (for plot)
    umap_proj_list = rbind(sapply(cluster1_size+cluster2_size+1:cluster3_size, function(x) {
      umap_sample_prevector <- admixt3.umap$layout[x,];
      umap_sample_vector = umap_sample_prevector - umap_centroid_cluster1;
      umap_sample_projection = ((sum(umap_sample_vector*umap_centroid_vector))/(umap_centroid_vector_lenght^2))*umap_centroid_vector;
      umap_sample_projection;
    } ))


    #### end vector projections

    #### Coefficent of determination : R squared ####
    # PCA #
    pca_norm_proj_yob = pca_norm_proj_list;
    pca_norm_proj_yth = cluster3_shiftf_list;
    pca_norm_proj_yme = mean(pca_norm_proj_list);
    pca_SS1 = sum((pca_norm_proj_yob - pca_norm_proj_yth)^2); #residual sum of squares
    pca_SS2 = sum((pca_norm_proj_yob - pca_norm_proj_yme)^2); #total sum of squares
    pca_Rsquare = 1 - (pca_SS1/pca_SS2); #coefficient of determination
    #Adjust Rsquare
    pca_Rsquare = 1-(1-pca_Rsquare)*((cluster3_size-1)/(cluster3_size-1-1));
    
    # UMAP #
    umap_norm_proj_yob = umap_norm_proj_list;
    umap_norm_proj_yth = cluster3_shiftf_list;
    umap_norm_proj_yme = mean(umap_norm_proj_list);
    umap_SS1 = sum((umap_norm_proj_yob - umap_norm_proj_yth)^2);
    umap_SS2 = sum((umap_norm_proj_yob - umap_norm_proj_yme)^2);
    umap_Rsquare = 1 - (umap_SS1/umap_SS2);
    #Adjust Rsquare
    umap_Rsquare = 1-(1-umap_Rsquare)*((cluster3_size-1)/(cluster3_size-1-1));
    
    Rsquare_tuple = cbind(pca_Rsquare , umap_Rsquare);
    Rsquare_tuple[Rsquare_tuple<0]<-0;
    #### end R squared
    
    
    #### simulation Plots ####
    par(mfrow=c(2,3),las=1);
    
    legend_str= paste0(
      #'number of simulations : ',nb_sim,
      '\n SIMULATION PARAMETERS \n',
      '\n cluster1 size : ',cluster1_size,
      '\n cluster2 size : ',cluster2_size,
      '\n number of intermediates : ',cluster3_size,
      '\n specific genes : ',nb_spec_genes,
      '\n non specific genes : ',nb_ns_genes,
      '\n max expression : ',max_exp,
      '\n min expression : ',min_exp,
      '\n Delta : ',Delta,
      '\n beta dist. param. : ',beta_param,
      '\n \n UMAP PARAMETERS \n',
      '\n n_neighbors : ',custom.config$n_neighbors,
      '\n min_dist : ',custom.config$min_dist,
      '\n spread : ',custom.config$spread
    )
    
    plot(admixt3[,1:2] , col =  hsv(color_matrix,1,1),main="2D of 2 spec.gene expr.", xlab="Gene expr. dim. 1", ylab="Gene expr. dim. 2" ,pch=16);
    ##PCA plot
    plot(admixt3.pca$li  , col =  hsv(color_matrix,1,1),main=paste0("PCA "), xlab="PC1", ylab="PC2",pch=16 );
    segments(pca_centroid_cluster1[1],pca_centroid_cluster1[2],pca_centroid_cluster2[1],pca_centroid_cluster2[2],col="orange",lwd=2);
    if (cluster3_size <= 100){ #plot projection lines
      p=sapply(1:cluster3_size, function(sample) {segments(admixt3.pca$li[cluster1_size+cluster2_size+sample,1],admixt3.pca$li[cluster1_size+cluster2_size+sample,2],pca_centroid_cluster1[1]+pca_proj_list[1,sample],pca_centroid_cluster1[2]+pca_proj_list[2,sample],col="goldenrod1",lty=3)})
    }
    plot(cluster3_shiftf_list,pca_norm_proj_list, asp=1, xlim = c(0, 1), ylim = c(0, 1), xlab="Expected projection", ylab="Observed projection", main=paste0("PCA projections, adj.R squared=",round(Rsquare_tuple[1],digits=4)));
    abline(0,1,col="red",lwd=2);
    
    #Empty plot to draw text on
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') ;
    text(x = 0.5, y = 0.5,legend_str)
    
    #UMAP plot
    plot(admixt3.umap$layout , col =  hsv(color_matrix,1,1), asp=1,main="UMAP", xlab="UMAP dim 1", ylab="UMAP dim 2" ,pch=16);
    segments(umap_centroid_cluster1[1],umap_centroid_cluster1[2],umap_centroid_cluster2[1],umap_centroid_cluster2[2],col="orange",lwd=2);
    if (cluster3_size <= 100){ #plot projection lines
      p=sapply(1:cluster3_size, function(sample) {segments(admixt3.umap$layout[cluster1_size+cluster2_size+sample,1],admixt3.umap$layout[cluster1_size+cluster2_size+sample,2],umap_centroid_cluster1[1]+umap_proj_list[1,sample],umap_centroid_cluster1[2]+umap_proj_list[2,sample],col="goldenrod1",lty=3)})
    }
    plot(cluster3_shiftf_list,umap_norm_proj_list, xlim = c(0, 1), ylim = c(0, 1), xlab="Expected projection", ylab="Observed projection", main=paste0("UMAP projections, adj.R squared=",round(Rsquare_tuple[2],digits=4)));
    abline(0,1,col="red",lwd=2);
    
    #### end simulation Plots
    #Sys.sleep(5); #wait 5 seconds before next simulation
        
    #Return tuple of Rsquare
    Rsquare_tuple;
  } )

  Rsquare_matrix_pca = cbind(Rsquare_matrix[1,])
  Rsquare_matrix_umap = cbind(Rsquare_matrix[2,])#umap
  Rsquare_matrix = cbind(Rsquare_matrix_pca , Rsquare_matrix_umap)

} )

#### R squared Boxplot ####
boxplot_matrix
legend_str= paste0(
  'number of simulations : ',nb_sim,
  '\n SIMULATION PARAMETERS \n',
  '\n cluster1 size : ',cluster1_size,
  '\n cluster2 size : ',cluster2_size,
  '\n number of intermediates : ',cluster3_size,
  '\n specific genes : ',nb_spec_genes,
  '\n non specific genes : ',nb_ns_genes,
  '\n max expression : ',max_exp,
  '\n min expression : ',min_exp,
  '\n Delta : ',Delta,
  '\n beta dist. param. : ',beta_param,
  '\n \n UMAP PARAMETERS \n',
  '\n n_neighbors : ',custom.config$n_neighbors,
  '\n min_dist : ',custom.config$min_dist
  )

legend_str
par(mfrow=c(1,3),las=1)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')#Empty plot to draw text on
text(x = 0.5, y = 0.5,legend_str)
boxplot(boxplot_matrix[,1,],ylim = c(0, 1),main="PCA projections adj.R squared",names=param_vector,xlab = testparem ,ylab = "Adjusted R squared",col = "orange",border = "brown") #PCA
boxplot(boxplot_matrix[,2,],ylim = c(0, 1),main="UMAP projections adj.R squared",names=param_vector,xlab = testparem,ylab = "Adjusted R squared",col = "orange",border = "brown") #UMAP

