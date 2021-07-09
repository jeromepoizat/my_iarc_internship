###############################################################################
### Script for simulation of high-dimensional (transcriptomic) data         ###
### and assessment of dimension reduction methods for the detection of      ###
### intermediate phenotypes                                                 ###
###                                                                         ###
### Generates a heatmap with with R² values for UMAP method         		    ###
###																			                                    ###
###					by Jerome Poizat										                            ###
###############################################################################

#dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
#install.packages("ade4") 
#install.packages("umap") 
#install.packages("gplots") 
#install.packages("optparse") 


library("optparse") 
# get options
option_list = list(
  make_option(c("-n", "--name"), type="character", default="HeatMap_DataSize_x_Delta", help="output file name [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-s", "--sim"), type="integer", default=2, help="Number of simulations [default= %default]", metavar="integer"),
  make_option(c("-a", "--var1"), type="character", default="param_vector = c(50, 100, 200, 500, 1000)", help="set values of first variable (DataSize) [default= %default]", metavar="character"),  
  make_option(c("-b", "--var2"), type="character", default="param_vector2 = c(0.5, 1, 2, 5)", help="set values of second variable (Delta: distance between clusters) [default= %default]", metavar="character"),   
  make_option(c("-x", "--spread"), type="integer", default=1, help="value of first constant (UMAP Spread) [default= %default]", metavar="integer"),
  make_option(c("-y", "--propint"), type="integer", default=10, help="value of second constant (% proportion of intermediates) [default= %default]", metavar="integer"),
  make_option(c("-z", "--const3"), type="integer", default=0, help="value of third constant if applicable [default= %default]", metavar="integer"),
  make_option(c("-g", "--sgenes"), type="integer", default=2000, help="Number of specific genes [default= %default]", metavar="integer"),
  make_option(c("-e", "--nsgenes"), type="integer", default=4000, help="Number of non specific genes [default= %default]", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


library(ade4) 
library(umap) 
library(gplots) 


nb_sim = opt$sim #number of simulation 
outputname = paste0(opt$name,"_Spread_",opt$spread,"_PropInterm_",opt$propint)
eval(parse(text=opt$var1))  #vector of values of the first parameter (data size) (param_vector)
eval(parse(text=opt$var2)) #vector of values of the second parameter (Delta: distance between clusters) (param_vector2)

#### Constant simulation parameters ####
nb_spec_genes = opt$sgenes #number of specific genes
nb_ns_genes = opt$nsgenes #number of non specific genes
cluster3_prop = (opt$propint/100) #proportion of intermediates
min_exp = 0 #minimum gene expression median value
max_exp = 100 #minimum gene expression median value
min_sigma = 1
max_sigma = 12#vector of possible sigmas values for non spec. genes, maybe use distribution instead of random sampling
beta_param = 1 #parameters a & b of the beta distribution of the mixture coeficient

#### Dimensionality reduction parameters ####
# UMAP Custom configs
custom.config = umap.defaults
custom.config$n_neighbors = 15
custom.config$spread = opt$spread
custom.config$min_dist = 0.1
custom.config$metric = "euclidean"

#### ####

current_sim <<- 0
total_sim <- (nb_sim*length(param_vector2)*length(param_vector))

#### heatmap_matrix ####
print("Building heatmap matrix")
#Loop through values of parameter2"
heatmap_matrix <- sapply(param_vector2, simplify = 'array', function(param_value2) { 
  
  Delta = param_value2
  
  #Loop through values of parameter1"
  boxplot_matrix <- sapply(param_vector, simplify = 'array', function(param_value) { 
    
    #### Variable simulation parameter ####
    data_size = param_value
    cluster1_size = (data_size - data_size*cluster3_prop)%/%2
    cluster2_size = (data_size - data_size*cluster3_prop)%/%2
    cluster3_size = (data_size - (cluster1_size+cluster2_size))
    if (cluster3_size < 2) {cluster3_size = 2}
    
    print(paste0("simulation : ",current_sim,"/",total_sim))
    
    #### Rsquare_matrix ####
    Rsquare_matrix = sapply(1:nb_sim, function(x) {
      
      current_sim <<- (current_sim + 1)
      #### Simulation ####
      
      source("//inti/Gcs/MESOMICS/Internship/Jerome/results/Scripts_and_Plots/simulation_script.R")
      simulation <- simulate(cluster1_size, cluster2_size, cluster3_size, min_exp, max_exp, Delta, nb_spec_genes,nb_ns_genes, beta_param)
      admixt3 <- simulation$admixt3
      cluster3_shiftf_vect <- simulation$cluster3_shiftf_list
      
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
      
      #### end vector projections
      
      
      
      #### Coefficent of determination : R² ####
      # PCA #
      pca_norm_proj_yob = pca_norm_proj_list;
      pca_norm_proj_yth = cluster3_shiftf_vect;
      pca_norm_proj_yme = mean(pca_norm_proj_list);
      pca_SS1 = sum((pca_norm_proj_yob - pca_norm_proj_yth)^2); #residual ++------------sum of squares
      pca_SS2 = sum((pca_norm_proj_yob - pca_norm_proj_yme)^2); #total sum of squares
      pca_Rsquare = 1 - (pca_SS1/pca_SS2); #coefficient of determination
      #Adjust Rsquare
      #pca_Rsquare = 1-(1-pca_Rsquare)*((cluster3_size-1)/(cluster3_size-1-1));
      
      #UMAP
      umap_norm_proj_yob = umap_norm_proj_list;
      umap_norm_proj_yth = cluster3_shiftf_vect;
      umap_norm_proj_yme = mean(umap_norm_proj_list);
      umap_SS1 = sum((umap_norm_proj_yob - umap_norm_proj_yth)^2);
      umap_SS2 = sum((umap_norm_proj_yob - umap_norm_proj_yme)^2);
      umap_Rsquare = 1 - (umap_SS1/umap_SS2);
      #Adjust Rsquare
      #umap_Rsquare = 1-(1-umap_Rsquare)*((cluster3_size-1)/(cluster3_size-1-1));
      
      Rsquare_tuple = cbind(pca_Rsquare , umap_Rsquare);
      Rsquare_tuple[Rsquare_tuple<0]<-0;
      #### end R²
      
      Rsquare_tuple;
    } )
    
    Rsquare_matrix_pca = cbind(Rsquare_matrix[1,])
    Rsquare_matrix_umap = cbind(Rsquare_matrix[2,])#umap
    Rsquare_matrix = cbind(Rsquare_matrix_pca , Rsquare_matrix_umap)
    
  } )
  
  
  pca_means = sapply(1:length(boxplot_matrix[1,1,]), function(x) {mean(boxplot_matrix[,1,][,x])})
  pca_25 = sapply(1:length(boxplot_matrix[1,1,]), function(x) {quantile(boxplot_matrix[,1,][,x],0.25)})
  pca_75 = sapply(1:length(boxplot_matrix[1,1,]), function(x) {quantile(boxplot_matrix[,1,][,x],0.75)})
  pca_var_diff = pca_75 - pca_25 # interquartile
  
  umap_means = sapply(1:length(boxplot_matrix[1,1,]), function(x) {mean(boxplot_matrix[,2,][,x])})
  umap_25 = sapply(1:length(boxplot_matrix[1,1,]), function(x) {quantile(boxplot_matrix[,2,][,x],0.25)})
  umap_75 = sapply(1:length(boxplot_matrix[1,1,]), function(x) {quantile(boxplot_matrix[,2,][,x],0.75)})
  umap_var_diff = umap_75 - umap_25 # interquartile
  
  
  meansNdiffs = cbind(pca_means, umap_means, pca_var_diff, umap_var_diff);
  meansNdiffs; #for heatmap_matrix
})

save(heatmap_matrix, file = paste0(outputname,".RData"))

#### plot heatmaps ####
#plot heatmaps umap R squarre means
png(filename=paste0(outputname,"_umap_Rsquarred.png"),width = 800, height = 600, pointsize=16)
par(mfrow=c(1,1),las=1)
heatmap.2(heatmap_matrix[,2,], cellnote = round(heatmap_matrix[,2,],digits=2),scale = "none", col = colorRampPalette(c("white","white","yellow","red"))(60), 
          trace = "none", density.info = "none",Rowv=NA, Colv=NA, notecol="black",margins =c(12,12),
          labRow = param_vector,labCol = param_vector2,
          main = "Projection's R squared",
          breaks = seq(0,1,1/60),
          
          revC =TRUE
)

dev.off()

#plot heatmap umap R squarre inter-quartile
png(filename=paste0(outputname,"_umap_interQ.png"),width = 800, height = 600, pointsize=16)
par(mfrow=c(1,1),las=1)
heatmap.2(heatmap_matrix[,4,], cellnote = round(heatmap_matrix[,4,],digits=2),scale = "none", col = colorRampPalette(c("white","cadetblue1","cyan3","darkcyan","dodgerblue4","dodgerblue4","dodgerblue4","navyblue","navyblue"))(60), 
          trace = "none", density.info = "none",Rowv=NA, Colv=NA, notecol="black",margins =c(12,12),
          labRow = param_vector,labCol = param_vector2,
          main = "Projection's R squared Q1-Q3",
          breaks = seq(0,1,1/60),
          
          revC =TRUE
)

dev.off()

#plot heatmaps pca R squarre means
png(filename=paste0(outputname,"_pca_Rsquarred.png"),width = 800, height = 600, pointsize=16)
par(mfrow=c(1,1),las=1)
heatmap.2(heatmap_matrix[,1,], cellnote = round(heatmap_matrix[,1,],digits=2),scale = "none", col = colorRampPalette(c("white","white","yellow","red"))(60), 
          trace = "none", density.info = "none",Rowv=NA, Colv=NA, notecol="black",margins =c(12,12),
          labRow = param_vector,labCol = param_vector2,
          main = "Projection's R squared",
          breaks = seq(0,1,1/60),
          
          revC =TRUE
)

dev.off()

#plot heatmap pca R squarre inter-quartile
png(filename=paste0(outputname,"_pca_interQ.png"),width = 800, height = 600, pointsize=16)
par(mfrow=c(1,1),las=1)
heatmap.2(heatmap_matrix[,3,], cellnote = round(heatmap_matrix[,3,],digits=2),scale = "none", col = colorRampPalette(c("white","cadetblue1","cyan3","darkcyan","dodgerblue4","dodgerblue4","dodgerblue4","navyblue","navyblue"))(60), 
          trace = "none", density.info = "none",Rowv=NA, Colv=NA, notecol="black",margins =c(12,12),
          labRow = param_vector,labCol = param_vector2,
          main = "Projection's R squared Q1-Q3",
          breaks = seq(0,1,1/60),
          
          revC =TRUE
)

dev.off()


print("Heatmaps done")