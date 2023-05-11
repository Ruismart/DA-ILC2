library(CountClust)

#### ref to https://bitbucket.org/jerry00/mouse_small_intestine_immune_cell_atlas/src/master/script/topicModeling_PP.R 

# Function for plotting cell specific values on tSNE, used in the above topic modeling function 
# obj = Seurat object 
# features = gene name or meta data column label to plot (can input several as a vector c(a,b,c))
# feature.type = either "meta" or "gene", indiciating whether you are plotting genes or meta data 
# ncols = number of columns you want in  your final plot, defaults to making the grid as square as possible  
# pt.size = size of dots on your tSNE
# title = title of your legend 
# lower = lower bound of your legend
# upper = upper bound of your legend 
# na.color = what colors for na values should be set to (useful when dealing with outliers)
plot.tsne.feature <- function(obj=countData, features="nGene", feature.type="meta", 
                              ncols=ceiling(sqrt(length(features))), pt.size=1, title="val", 
                              lower=NULL, upper=NULL, na.color="gray") {
  
  library(reshape2)
  library(plotly)
  library(tidyr)
  
  # Set color gradient
  colors<-c("#191970","#121285","#0C0C9A","#0707B0","#0101C5","#0014CF","#0033D3","#0053D8","#0072DD","#0092E1","#00B2E6",
            "#00D1EB","#23E8CD","#7AF17B","#D2FA29","#FFEB00","#FFC300","#FF9B00","#FF8400","#FF7800","#FF6B00","#FF5F00","#FF5300",
            "#FF4700","#F73B00","#EF2E00","#E62300","#DD1700","#D50B00","#CD0000")
  
  # Error out if user doesn't specify correct feature type 
  if (feature.type!="meta" & feature.type!="gene"){stop("feature type must be 'meta' or 'gene'")}
  
  # Collect feature info
  if (feature.type=="meta"){ 
    feature.info <- as.matrix(obj@meta.data[,features]) 
    if (length(features)==1){
      colnames(feature.info) <- features
    } 
    
    # Build data frame of feature info and tsne coordinates 
    #tmp.df <- data.frame(feature.info, obj@dr$tsne@cell.embeddings)
    tmp.df <- data.frame(feature.info, obj@reductions[['tsne']]@cell.embeddings)
    plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
    
  } else if (feature.type=="gene"){ 
    feature.info <- as.matrix(obj@assays[['RNA']]@data[features,]) 
    if (length(features)==1){
      colnames(feature.info) <- features
      #feature.info <- t(feature.info)
    }   
    
    # Build data frame of feature info and tsne coordinates 
    #tmp.df <- data.frame(t(feature.info), obj@dr$tsne@cell.embeddings)
    tmp.df <- data.frame(feature.info, obj@reductions[['tsne']]@cell.embeddings)
    plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
  }
  
  # Set scales for color bar
  if (!length(lower) & !length(upper)){
    lower = min(plot.df$val)
    upper = max(plot.df$val)
  }
  
  # Plot the data on tSNE
  p <- ggplot(plot.df, aes(x=tSNE_1, y=tSNE_2)) + geom_point(aes(color=val), alpha=0.8, shape=16, size=pt.size) + 
    theme(aspect.ratio = 1) + scale_color_gradientn(colors=colors, limits=c(lower, upper), na.value=na.color)
  p <- p + theme(aspect.ratio=1, text = element_text(size=10), axis.text=element_text(size=6), 
                 strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")), 
                 panel.background = element_blank(), plot.background=element_blank(), 
                 axis.line = element_line(size=0.5)) + 
    labs(color = title) +
    facet_wrap( ~ name, ncol=ncols) + 
    theme(strip.text = element_text(size=10))
  return(p)
}

plot.umap.feature <- function(obj=countData, features="nGene", feature.type="meta", 
                              ncols=ceiling(sqrt(length(features))), pt.size=1, title="val", 
                              lower=NULL, upper=NULL, na.color="gray", friendly=FALSE) {
  
  library(reshape2)
  library(plotly)
  library(tidyr)
  
  # Set color gradient
  colors<-c("#191970","#121285","#0C0C9A","#0707B0","#0101C5","#0014CF","#0033D3","#0053D8","#0072DD","#0092E1","#00B2E6",
            "#00D1EB","#23E8CD","#7AF17B","#D2FA29","#FFEB00","#FFC300","#FF9B00","#FF8400","#FF7800","#FF6B00","#FF5F00","#FF5300",
            "#FF4700","#F73B00","#EF2E00","#E62300","#DD1700","#D50B00","#CD0000")
  
  # Error out if user doesn't specify correct feature type 
  if (feature.type!="meta" & feature.type!="gene"){stop("feature type must be 'meta' or 'gene'")}
  
  # Collect feature info
  if (feature.type=="meta"){ 
    feature.info <- as.matrix(obj@meta.data[,features]) 
    if (length(features)==1){
      colnames(feature.info) <- features
    } 
    
    # Build data frame of feature info and umap coordinates 
    #tmp.df <- data.frame(feature.info, obj@dr$umap@cell.embeddings)
    tmp.df <- data.frame(feature.info, obj@reductions[['umap']]@cell.embeddings)
    plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
    
  } else if (feature.type=="gene"){ 
    feature.info <- as.matrix(obj@assays[['RNA']]@data[features,]) 
    if (length(features)==1){
      colnames(feature.info) <- features
      #feature.info <- t(feature.info)
    }   
    
    # Build data frame of feature info and umap coordinates 
    #tmp.df <- data.frame(t(feature.info), obj@dr$umap@cell.embeddings)
    tmp.df <- data.frame(feature.info, obj@reductions[['umap']]@cell.embeddings)
    plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
  }
  
  # Set scales for color bar
  if (!length(lower) & !length(upper)){
    lower = min(plot.df$val)
    upper = max(plot.df$val)
  }
  
  # Plot the data on UMAP
  p <- ggplot(plot.df, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=val), alpha=0.8, shape=16, size=pt.size) + 
    theme(aspect.ratio = 1) + scale_color_gradientn(colors=colors, limits=c(lower, upper), na.value=na.color)
  p <- p + theme(aspect.ratio=1, text = element_text(size=10), axis.text=element_text(size=6), 
                 strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")), 
                 panel.background = element_blank(), plot.background=element_blank(), 
                 axis.line = element_line(size=0.5)) + 
    labs(color = title) +
    facet_wrap( ~ name, ncol=ncols) + 
    theme(strip.text = element_text(size=10), strip.background = element_rect(color = "white", fill = "white"))
  if(friendly==FALSE){
      return(p)   
  }else if(friendly==TRUE){
      return(Seurat::AugmentPlot(plot = p)) 
  }

}

