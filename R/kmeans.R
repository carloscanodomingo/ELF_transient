#library(easyGgplot2)
library(latex2exp)
library(R.matlab)
library(ggplot2)
library(qqplotr)
#library(qualityTools)
library(goft)
library(fitdistrplus)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

matlabFile <- readMat("r_am_sr.mat")




varNames    <- names(matlabFile$data[,,1])
datList     <- matlabFile$data
datList     <- lapply(datList, unlist, use.names=FALSE)
data        <- as.data.frame(datList)
names(data) <- varNames


res.pca <- prcomp(df, scale = TRUE, rank. = 2)
res.pca
data_clustering <- function(data)
{
  #data$start.time
  data_subset <- subset( data, select = c(peak.value,PSD.in.band,FWHMT,RiseTime,PSD.rest, T1) )
  
  # REmove unused values
  df = data_subset;
  
  
  df <- na.omit(df);
  
  df <- scale(df)
  
  qburst_bool = (data$peak.value > 2e6)
  
  
  k2 <- kmeans(df, centers =2, nstart = 25, iter.max = 20)
  p <- fviz_cluster(k2, data = df, choose.vars = c("peak.value", "PSD.in.band"),
               palette = "data1", 
               geom = "point",
               ellipse.type = "convex", 
               ggtheme = theme_bw()
  )
  p2 <- fviz_cluster(k2, data = df,
                    palette = "data1", 
                    geom = "point",
                    ellipse.type = "convex", 
                    ggtheme = theme_bw()
  )
  k2$cluster[k2$cluster == 2] = "Normal";
  k2$cluster[k2$cluster == 1] = "Qburst"
  fviz_cluster(object = k2, data = df,geom = "point",shape = 3, show.clust.cent = TRUE, labelsize =0,
               ellipse.level = 0.90, choose.vars = c("peak.value", "PSD.in.band"), star.plot = TRUE, repel = TRUE) +
   labs( x="Normalized PSD [6 - 10]Hz AU", y = "Normalized cPeak Value (pT)") + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"))
    theme_bw() +
    theme(legend.position = "none", legend.text = c("Qburst","Normal"))
  
  
  data <- p$data # this is all you need
  data_2<- p2$data 
  
  # calculate the convex hull using chull(), for each cluster
  hull_data <-  data %>%
    group_by(cluster) %>%
    slice(chull(x, y))
  # plot: you can now customize this by using ggplot sintax
  ggplot(data, aes(x, y)) + geom_point() +
  geom_polygon(data = hull_data, alpha = 0.5, aes(fill=cluster, linetype=cluster)) + 
    labs(x=expression("Normalized Peak Value"), y=expression("Normalized Band Power")) + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      legend.position = c(0.60, .95),
      legend.direction = "horizontal",
      legend.justification = c("top"),
      legend.title = element_blank(),#element_text(face = "bold",hjust = 0.5),
      legend.box.just = "left",
      legend.box.background = element_rect(color = 'black', size=2),
      legend.text = element_text(size = 10, colour = "black", face="bold"),
      legend.margin = margin(6, 6, 6, 6),
      axis.text.x = element_text(angle=0, size = 12, color = "black"),
      axis.text.y = element_text(angle=0, size = 12, color = "black"),
    )
  path = paste("cluster/", "cluster_dim.pdf", sep = "")
  
  ggsave(path, width = 12, height = 8, units = "cm")
  
  
  # calculate the convex hull using chull(), for each cluster
  hull_data_2 <-  data_2 %>%
    group_by(cluster) %>%
    slice(chull(x, y))
  # plot: you can now customize this by using ggplot sintax
  ggplot(data_2, aes(x, y)) + geom_point() +
    geom_polygon(data = hull_data_2, alpha = 0.5, aes(fill=cluster, linetype=cluster)) + 
    labs(x=bquote(.(p2$labels$x)), y=bquote(.(p2$labels$y))) + theme(
      plot.title = element_text(color="black", size=14, hjust = 0.5),
      axis.title.y = element_text(color="black", size=14),
      axis.title.x = element_text(color="black", size=14),
      legend.position = c(0.60, .95),
      legend.direction = "horizontal",
      legend.justification = c("top"),
      legend.title = element_blank(),#element_text(face = "bold",hjust = 0.5),
      legend.box.just = "left",
      legend.box.background = element_rect(color = 'black', size=2),
      legend.text = element_text(size = 10, colour = "black", face="bold"),
      legend.margin = margin(6, 6, 6, 6),
      axis.text.x = element_text(angle=0, size = 12, color = "black"),
      axis.text.y = element_text(angle=0, size = 12, color = "black"),
    )
  path = paste("cluster/", "cluster_adim.pdf", sep = "")
  
  ggsave(path, width = 12, height = 8, units = "cm")
  
  
  #sil <- silhouette(k2$cluster, dist(data), ordered = FALSE)
  #row.names(sil) <- row.names(data) # Needed to use label option
  #fviz_silhouette(sil, label = TRUE)
  val = fviz_nbclust(data, kmeans, method = "wss", k.max = 8)
  result <- cbind(val$data, val$data$y/max(val$data$y))
  return(k2$cluster)
}

