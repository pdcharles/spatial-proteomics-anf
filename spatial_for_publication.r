# # # Header # # #

#Code to generate Figures 3E and 6A in Davis et al, "Deep topographic proteomics of a human brain tumour"
#Author: Phil Charles
#Date: 15.02.22

# !! Main analysis starts on line 224, everything prior is setup !!

# # # Setup # # # 

#Load libraries
library(data.table)
library(ANF)
library(Spectrum)
library(magick)
library(tools)

#Define a generic image writing function that's compatible with notebooks and markdown processing
img <- function(filename,code,width=7,height=7,units='in',res=144,env=.GlobalEnv,...) {
 is_jupyter <- 'IRkernel' %in% loadedNamespaces() 
 is_knitr <- !is.null(getOption("knitr.in.progress"))
 is_png <- file_ext(filePath)=='png'
 is_jpg <- file_ext(filePath)=='jpg' || file_ext(filePath)=='jpeg'
 dir.create('out',showWarnings=F)
 if (is_jupyter || is_knitr || is_png) {
  fn <- paste0('out/',if (is_png) filename else sub('\\.[^.]+$','.png',filename))
  png(fn,width=width,height=height,units=units,res=res,...)
  eval(substitute(code),envir=env)
  dev.off()
  if (is_jupyter) IRdisplay::display_png(file=fn,width=width,height=height,units=units,res=res)
  else if (is_knitr) knitr::include_graphics(fn,dpi=NA)
 }
 else if (is_jpg) {
  jpg(paste0('out/',filename),width=width,height=height,units=units,res=res,...)
  eval(substitute(code),envir=env)
  dev.off()
 }
 else {
  if (units=='px') {
   width <- width/72
   height <- height/72
  }
  if (!grepl('\\..+$',filename)) {
   filename <- paste0(filename,'.pdf')
  }
  else {
   filename <- sub('\\.[^.]+$','.pdf',filename)
  }
  cairo_pdf(paste0('out/',filename),width=width,height=height,...)
  eval(substitute(code),envir=env)
  dev.off()
 }
}

#Define an object class, 'slide' to hold slide data and relative distance data
setClass('slide', slots=list(layout='matrix', expr='matrix', vars='data.frame', meta='data.frame', cells='list', cellDims='list'))
setGeneric('plot_similarity',function(obj,distance,cols=NULL,clusterNumbers=NULL,annotate='cluster',borderCol=NULL,borderInset=0,borderWidth=4,labelCex,images=F,filePath=NULL,...) warning('No plot_similarity method exists for objs of class ',paste0(class(obj),collapse=',')))
setMethod('plot_similarity','slide',function(obj,distance,cols=NULL,clusterNumbers=NULL,annotate='cluster',borderCol=NULL,borderInset=0,borderWidth=4,labelCex,images=F,filePath=NULL,...) {
 if (images && (missing(filePath) || is.null(filePath))) stop('Must supply filePath if images=TRUE',call.=F)
 if (images && file_ext(filePath)=='png') stop('PNG Format does not seem to work; use JPEG',call.=F)
 if (missing(distance)) distance <- 0
 if (typeof(distance)=='double') { #a distance measure
  if (length(distance) == 1) {
   distance <- matrix(distance,nrow=length(obj@layout),ncol=length(obj@layout),dimnames=list(c(obj@layout),c(obj@layout)))
  }
  else if (!(length(dim(distance)) > 1 && all(dim(distance) <= length(obj@layout)) && all(unlist(dimnames(distance)) %in% slide@layout) )) {
   stop('Distance parameter array does not match slide layout entries',call.=F)
  }
 }
 else if (typeof(distance)=='logical' && length(distance) == length(slide@layout)) { #an include / exclude bool vector
  distance <- matrix(0,nrow=length(slide@layout),ncol=length(slide@layout),dimnames=list(colnames(slide@expr),colnames(slide@expr)))[distance,distance]
 }
 else {
  stop('Unknown distance parameter format',call.=F)
 }
 x_distance <- sapply(2:ncol(obj@layout),function(j) {
  sapply(1:nrow(obj@layout),function(i) {
   if (obj@layout[i,j-1] %in% rownames(distance) && obj@layout[i,j] %in% colnames(distance)) distance[obj@layout[i,j-1],obj@layout[i,j]]
   else 0
  })
 })
 y_distance <- sapply(1:ncol(obj@layout),function(j) {
  sapply(2:nrow(obj@layout),function(i) {
   if (obj@layout[i-1,j] %in% rownames(distance) && obj@layout[i,j] %in% colnames(distance)) distance[obj@layout[i-1,j],obj@layout[i,j]]
   else 0
  })
 })
 mid_x <- floor(ncol(obj@layout)/2)
 mid_y <- floor(nrow(obj@layout)/2)
 len_x <- max(apply(x_distance,1,sum,na.rm=T))+ncol(obj@layout)+8
 len_y <- max(apply(y_distance,2,sum,na.rm=T))+nrow(obj@layout)+8
 if (images) {
  fig <- image_blank(width=ceiling(len_x)*obj@cellDims$width
                    ,height=ceiling(len_y)*obj@cellDims$height
                    ,color='white'
                    )
  if (missing(labelCex)) labelCex <- 3
 }
 else {
  if (names(dev.cur()) %in% c('quartz','X11','windows','null device')) dev.new(width=7.5,height=7.5/4*3)
  par(mar=c(4,4,4,4))
  plot(NA,xlim=c(-len_x/2,len_x/2),ylim=c(-len_y/2,len_y/2),bty='n',xlab='',ylab='',xaxt='n',yaxt='n',...)
  if (missing(labelCex)) labelCex <- 0.5
 }

 borders_and_numbers <- list()

 for (j in 1:ncol(obj@layout)) {
  for (i in 1:nrow(obj@layout)) {
   if (j < mid_x) x <- -(sum(x_distance[i,j:(mid_x-1)]) + (mid_x-j))
   else if (j > mid_x) x <- sum(x_distance[i,mid_x:(j-1)]) + (j-mid_x)
   else x <- 0
   if (i < mid_y) y <- sum(y_distance[i:(mid_y-1),j]) + (mid_y-i)
   else if (i > mid_y) y <- -(sum(y_distance[mid_y:(i-1),j]) + (i-mid_y))
   else y <- 0
   if (obj@layout[i,j] %in% rownames(distance)) {
    cell_clus_number <- NULL
    border_col_override <- NULL
    shading_col <- NULL
    bn <- list(layoutCode=obj@layout[i,j])
    if (!is.null(clusterNumbers) && obj@layout[i,j] %in% names(clusterNumbers)) cell_clus_number <- clusterNumbers[[ obj@layout[i,j] ]]
    if (!is.null(cols) && obj@layout[i,j] %in% names(cols)) {
     if(is.list(cols[[ obj@layout[i,j] ]])) {
      shading_col <- cols[[ obj@layout[i,j] ]]$cell
      border_col_override <- cols[[ obj@layout[i,j] ]][c('bottom','top','left','right')]
     }
     else shading_col <- cols[[ obj@layout[i,j] ]]
    }
    if (!is.null(cell_clus_number)) {
     bn$number <- list(label=cell_clus_number,x=x,y=y)
     if (!is.null(borderCol)) {
      if (i < nrow(obj@layout) && (obj@layout[i+1,j] %in% names(clusterNumbers)) && cell_clus_number != clusterNumbers[[ obj@layout[i+1,j] ]]) {
       bn$bottom <- list(x=c(x+borderInset,x+1-borderInset)-0.5,y=c(y+borderInset,y+borderInset)-0.5,col=borderCol)
      }
      if (i > 1 && (obj@layout[i-1,j] %in% names(clusterNumbers)) && cell_clus_number != clusterNumbers[[ obj@layout[i-1,j] ]]) {
       bn$top <- list(x=c(x+borderInset,x+1-borderInset)-0.5,y=c(y+1-borderInset,y+1-borderInset)-0.5,col=borderCol)
      }
      if (j > 1 && (obj@layout[i,j-1] %in% names(clusterNumbers)) && cell_clus_number != clusterNumbers[[ obj@layout[i,j-1] ]]) {
       bn$left <-  list(x=c(x+borderInset,x+borderInset)-0.5,y=c(y+borderInset,y+1-borderInset)-0.5,col=borderCol)
      }
      if (j < ncol(obj@layout) && (obj@layout[i,j+1] %in% names(clusterNumbers)) && cell_clus_number != clusterNumbers[[ obj@layout[i,j+1] ]]) {
       bn$right <- list(x=c(x+1-borderInset,x+1-borderInset)-0.5,y=c(y+borderInset,y+1-borderInset)-0.5,col=borderCol)
      }
     }
     borders_and_numbers <- append(borders_and_numbers,list(bn))
    }
    else if (!is.null(shading_col) && !is.null(borderCol)) {
     if (i < nrow(obj@layout) && (obj@layout[i+1,j] %in% names(clusterNumbers)) && shading_col != cols[[ obj@layout[i+1,j] ]]) {
      bn$bottom <- list(x=c(x+borderInset,x+1-borderInset)-0.5,y=c(y+borderInset,y+borderInset)-0.5,col=borderCol)
     }
     if (i > 1 && (obj@layout[i-1,j] %in% names(clusterNumbers)) && shading_col != cols[[ obj@layout[i-1,j] ]]) {
      bn$top <- list(x=c(x+borderInset,x+1-borderInset)-0.5,y=c(y+1-borderInset,y+1-borderInset)-0.5,col=borderCol)
     }
     if (j > 1 && (obj@layout[i,j-1] %in% names(clusterNumbers)) && shading_col != cols[[ obj@layout[i,j-1] ]]) {
      bn$left <-  list(x=c(x+borderInset,x+borderInset)-0.5,y=c(y+borderInset,y+1-borderInset)-0.5,col=borderCol)
     }
     if (j < ncol(obj@layout) && (obj@layout[i,j+1] %in% names(clusterNumbers)) && shading_col != cols[[ obj@layout[i,j+1] ]]) {
      bn$right <- list(x=c(x+1-borderInset,x+1-borderInset)-0.5,y=c(y+borderInset,y+1-borderInset)-0.5,col=borderCol)
     }
    }
    if (!is.null(border_col_override)) {
     for (name in names(bn)) {
      if (name %in% names(border_col_override)) bn[[name]]$col <- border_col_override[[name]]
     }
    }
    borders_and_numbers <- append(borders_and_numbers,list(bn))
    if (images) {
     cell <- obj@cells[[ obj@layout[i,j] ]]
     if (!is.null(shading_col)) cell <- image_colorize(cell, 50, shading_col)
     fig <- image_composite(fig
                           ,cell
                           ,offset = paste0(if(x>=0) '+' else '-',abs(x)*obj@cellDims$width,if(y>=0) '-' else '+',abs(y)*obj@cellDims$height)
                           ,gravity = 'center'
                           )
    }
    else {
     rect(x-0.5,y-0.5,x+0.5,y+0.5,col=if(is.null(shading_col)) 'white' else shading_col)
    }
   }
  }
 }
 if (images) {
  fig <- image_draw(fig,bg='white',xlim=c(-len_x/2,len_x/2),ylim=c(-len_y/2,len_y/2),...)
 }
 if (length(borders_and_numbers) > 0) {
  for (bn in borders_and_numbers) {
   if (!is.null(annotate)) {
    if (annotate=='cluster' && !is.null(bn$number)) {
     text(bn$number$x,bn$number$y,bn$number$label,cex=labelCex,font=2,col='black')
    }
    else if (annotate=='layout') {
     text(x,y,bn$layoutCode,cex=0.5)
    }
   }
   for (border in c('bottom','top','left','right')) {
    if (!is.null(bn[[border]])) {
     lines(bn[[border]]$x,bn[[border]]$y,col=bn[[border]]$col,lwd=borderWidth)
    }
   }
  }
 }
 if (images) {
  dev.off()
  image_write(fig,path=filePath,format=file_ext(filePath))
 }
})

#Define a helper function to identify the optimum number of clusters for spectral clustering, using maximal eigengap heuristic
get_clusters <- function(affinity_matrix) {
 dv <- 1/sqrt(rowSums(affinity_matrix))
 L <- dv * affinity_matrix %*% diag(dv) #symmetrically normalised laplacian
 evals <- suppressWarnings(as.numeric(eigen(L)$values))
 n_clusters <- which.max(abs(diff(evals)[-1])[1:10])+1
 cat(paste0('Number of clusters: ',n_clusters,'\n'))
 n_clusters
}

#Define a helper function to return cluster wise colors
clust_cols <- RColorBrewer::brewer.pal(12,'Set3')
get_cols <- function(clustering) {
 sapply(names(clustering),function(s) clust_cols[as.numeric(clustering[[s]])])
}

# # # Main Analysis # # #

dir.create('out/', showWarnings=F)

#Load the quantitation data into a slide object
slide <- (function() {
 d <- data.table::fread('sources/ATRT-384-LFQ.csv')

 plate1 <- matrix(rev(1:96), nrow=8, ncol=12, byrow=T)
 plate2 <- matrix(rev(97:192), nrow=8, ncol=12, byrow=T)
 plate3 <- matrix(rev(193:288), nrow=8, ncol=12, byrow=T)
 plate4 <- matrix(rev(289:384), nrow=8, ncol=12, byrow=T)
 lyt <- rbind(cbind(plate1,plate2),cbind(plate3,plate4))
 lyt <- apply(lyt,2,as.character)

 expr <- as.matrix(d[,1:384])
 rownames(expr) <- d$Gene.names
 
 meta <- data.frame(type=c(rep('Tissue',84),rep('Haemorrhage',2),rep('Tissue',10) #plate1
                          ,rep('Tissue',52),rep('Haemorrhage',3),rep('Tissue',5),'Empty',rep('Tissue',2),rep('Haemorrhage',5),rep('Tissue',4),'Empty','Tissue',rep('Haemorrhage',9),'Tissue','Empty',rep('Haemorrhage',11) #plate2
                          ,rep(c(rep('Tissue',11),'Empty'),7),rep('Tissue',12) #plate3
                          ,'Empty',rep('Tissue',11),'Empty',rep('Tissue',11),rep('Tissue',72)#plate4 
                          )
                   ,plate=c(rep(1,length(plate1))
                           ,rep(2,length(plate2))
                           ,rep(3,length(plate3))
                           ,rep(4,length(plate4))
                           )
                   ,sample=colnames(expr)
                   )

 rownames(meta) <- colnames(expr) <- as.character(1:384)

 cat('Reading image\n')
 img <- image_read('sources/SampledRegion.jpg')
 img_left_mar <- 10
 img_right_mar <- 10
 img_top_mar <- 7
 img_bottom_mar <- 7

 img <- image_crop(img
                  ,geometry_area(width = image_info(img)$width-img_left_mar-img_right_mar
                                ,height = image_info(img)$height-img_top_mar-img_bottom_mar
                                ,x_off = img_left_mar
                                ,y_off = img_top_mar
                                )
                  )
 img <- image_scale(img
                   ,geometry_area(width = image_info(img)$width/10
                                 ,height = image_info(img)$height/10
                                 )
                   )

 img_block_x <- (image_info(img)$width-img_left_mar-img_right_mar)/dim(lyt)[2]
 img_block_y <- (image_info(img)$height-img_top_mar-img_bottom_mar)/dim(lyt)[1]

 imgCells <- list()
 cat('Splitting image\n')
 for (s in sort(lyt)) {
  coords <- which(lyt==s,arr.ind=T)
  imgCells[[s]] <- image_crop(img,paste0(img_block_x,'x',img_block_y,'+',img_left_mar+img_block_x*(coords[2]-1),'+',img_top_mar+img_block_y*(coords[1]-1)))
 }
 new('slide',layout=lyt,expr=expr,vars=data.frame(d[,385:388]),meta=meta,cells=imgCells,cellDims=list(width=img_block_x-1,height=img_block_y))
})()


#Generate a distance matrix for each LCM-sampled location relative to other locations, using euclidean distances
set.seed(0)
layout_dist <- matrix(NA
                     ,nrow=length(slide@layout)
                     ,ncol=length(slide@layout)
                     ,dimnames=list(colnames(slide@expr),colnames(slide@expr))
                     )

for (j in 1:ncol(slide@layout)) {
 for (i in 1:nrow(slide@layout)) {
  for (b in 1:ncol(slide@layout)) {
   for (a in 1:nrow(slide@layout)) {
    layout_dist[ slide@layout[i,j], slide@layout[a,b] ] <- sqrt((i-a)^2 + (j-b)^2)
   }
  }
 }
}

#Filter the distance matrix to remove 'empty' (no material) locations to prevent data skew. Empty locations are not considered for clustering 
non_empty <- slide@meta[colnames(slide@expr),'type']!='Empty'
layout_dist_non_empty <- layout_dist[non_empty,non_empty]

#Generate a symmetric affinity matrix based on pairwise distances
layout_dist_affinity <- affinity_matrix(layout_dist_non_empty/max(layout_dist_non_empty), 12)

#Basic normalisation on MS protein-wise quantiation values - log2 transform, set zero values as missing and median centre
exprs_processed <- log2(slide@expr)[,non_empty]
exprs_processed[is.infinite(exprs_processed)] <- NA
exprs_processed <- sweep(exprs_processed,2,apply(exprs_processed,2,median,na.rm=T))

#Calculate protein-wise variance across locations
variances <- apply(exprs_processed,1,var,na.rm=T)

#Use protein-wise variance to identify informative (high variance, top 25th percentile) protein features
features <- !is.na(variances) & variances > quantile(variances,0.75,na.rm=T)

#Generate a second symmetric affinity matrix based on pairwise spearman corrleation across informative features
expr_features_affinity <- affinity_matrix(1-cor(slide@expr[features,non_empty],method='spearman'), 12)

#Fuse both affinity matrices using ANF
fusion_features = ANF(list(expr_features_affinity,layout_dist_affinity),K=20)

#Perform spectral clustering on the ANF fusion matrix using optimum cluster number
cat('Naive clustering\n')
clust_features_spectral_anf = spectral_clustering(fusion_features, get_clusters(fusion_features))
names(clust_features_spectral_anf) <- colnames(slide@expr[,non_empty])

#Generate the slide plot using the spectral clustering result
plot_similarity(slide,non_empty,cols=get_cols(clust_features_spectral_anf),clusterNumbers=clust_features_spectral_anf,annotate=NULL,borderCol='black',images=T,filePath='out/figforpaper - Naive.jpg')
plot_similarity(slide,non_empty,cols=get_cols(clust_features_spectral_anf),clusterNumbers=clust_features_spectral_anf,borderCol='black',images=T,filePath='out/figforpaper - Naive annotated.jpg')
write.table(data.frame(Sample=names(clust_features_spectral_anf),Cluster=clust_features_spectral_anf),file='out/cluster_assignment - Naive.tsv',sep='\t',row.names=F,col.names=T)

#Load known matrisome ids
ids_of_interest <- sapply(scan('sources/Matrisome.csv',what=list('character','character'),skip=1,quiet=T,sep=','),function(ids) ids[ids!=''])
names(ids_of_interest) <- scan('sources/Matrisome.csv',what=list('character','character'),nlines=1,quiet=T,sep=',')

# # # Core Matrisome Component Analysis # # # 

#Select quantified proteins matching core matrisome ids
id_set <- 'Core Matrisome'
rows_of_interest <- sapply(strsplit(rownames(slide@expr),';'),function(id) any(id %in% ids_of_interest[[id_set]]))

#Further filter the distance matrix to remove locations where no matrisome proteins change
cols_of_interest <- non_empty & apply(slide@expr[rows_of_interest,],2,sd) > 0
layout_dist_interest <- layout_dist[cols_of_interest,cols_of_interest]

#Generate a new symmetric affinity matrix based on pairwise distances restricted just to the filtered distance matrix
layout_dist_interest_affinity <- affinity_matrix(layout_dist_interest/max(layout_dist_interest), 12)

#Generate a new second symmetric affinity matrix based on pairwise spearman corrleation across core matrisome features
expr_interest_affinity <- affinity_matrix(1-cor(slide@expr[rows_of_interest,cols_of_interest],method='spearman'), 12)

#Fuse both matrisome affinity matrices using ANF
fusion_features_interest = ANF(list(expr_interest_affinity,layout_dist_interest_affinity),K=20)

#Perform spectral clustering on the matrisome ANF fusion matrix using optimum cluster number
cat(paste0(id_set,' clustering\n'))
clust_features_spectral_anf_interest <- spectral_clustering(fusion_features_interest, get_clusters(fusion_features_interest))
names(clust_features_spectral_anf_interest) <- colnames(slide@expr)[cols_of_interest]

#Generate the slide plot using the matrisome spectral clustering result
plot_similarity(slide,cols_of_interest,cols=get_cols(clust_features_spectral_anf_interest),clusterNumbers=clust_features_spectral_anf_interest,annotate=NULL,borderCol='black',images=T,filePath=paste0('out/figforpaper - ',id_set,'.jpg'))
plot_similarity(slide,cols_of_interest,cols=get_cols(clust_features_spectral_anf_interest),clusterNumbers=clust_features_spectral_anf_interest,borderCol='black',images=T,filePath=paste0('out/figforpaper - ',id_set,' annotated.jpg'))
write.table(data.frame(Sample=names(clust_features_spectral_anf_interest),Cluster=clust_features_spectral_anf_interest),file=paste0('out/cluster_assignment - ',id_set,'.tsv'),sep='\t',row.names=F,col.names=T)
