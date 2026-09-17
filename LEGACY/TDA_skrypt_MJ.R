library(TDA)
library(TDAmapper)
library(igraph)
library(dplyr)
library(kableExtra)

set.seed(123)


# Przygotowanie danych
data.prep <- function(path, typ = "original")
{
  dane <- as.data.frame(read.csv2(path))
  dane.num <- dane[, -c(1,4,5,6)]
  dane.num[] <- lapply(dane.num, as.numeric)
  dane.scale <- scale(dane.num)
  if(typ == "original") return(dane)
  if(typ == "numeric") return(dane.num)
  if(typ == "scale") return(dane.scale)
}



# Homologia trwała - metoda landmarków
landmarks <- function(data.scale,
                      lands,
                      iter.max = 50,
                      nstart,
                      maxdimension = 1,
                      library = "GUDHI")
{
  clusters <- kmeans(data.scale, centers = lands, iter.max = iter.max, nstart = nstart)
  landmarks <- clusters$centers
  
  odleglosci <- dist(landmarks)
  
  maxscale <- max(odleglosci)
  
  homology.landmarks <- ripsDiag(X = landmarks, maxdimension = maxdimension, maxscale = maxscale, library = library)
  
  plot(homology.landmarks[["diagram"]], main = "Diagram persystencji")
}

# Metoda landscape
landscapes <- function(homolody.landmarks, 
                       len,
                       dimension,
                       xlim)
{
  max.death <- max(homology.landmarks[["diagram"]][,3])
  lands1 <- TDA::landscape(homology.landmarks[["diagram"]], dimension = dimension, KK = 1, tseq = seq(0, max.death, length = len))
  
  lands2 <- TDA::landscape(homology.landmarks[["diagram"]], dimension = dimension, KK = 2, tseq = seq(0, max.death, length = len))
  
  plot(seq(0, max.death, length = len),lands1, type= "l", main = "Landscapes", col = "maroon", ylab = "lands", xlim = xlim)
  lines(seq(0, max.death, length = len),lands2, type= "l", col = "green")
  legend("topright",                   
         legend = c("Lands1", "Lands2"), 
         col = c("maroon", "green"),  
         lty = 1,                    
         lwd = 2)
}

# Algorytm mapper
Mapper <- function(data.num,
                   num.intervals,
                   percent.overlap,
                   num.bins,
                   table = F)
{
  data.scale <- scale(data.num)
  pca.wynik <- prcomp(data.scale, rank. = 1)
  soczewka <- pca.wynik$x[,1]
  
  wynik.mapper <- mapper1D(distance_matrix = dist(data.scale),
                           filter_values = soczewka,
                           num_intervals = num.intervals,
                           percent_overlap = percent.overlap,
                           num_bins_when_clustering = num.bins)
  
  graf.mapper <- graph.adjacency(wynik.mapper$adjacency, mode = "undirected")
  
  srednie.edycje <- sapply(wynik.mapper$points_in_vertex, function(ind.points) {mean(data.num$edits[ind.points], na.rm = T)})
  
  paleta <- colorRampPalette(c('blue', 'yellow', 'red'))(100)
  min.ed <- min(srednie.edycje, na.rm = T)
  max.ed <- max(srednie.edycje, na.rm  = T)
  
  scale.edits <- round((srednie.edycje - min.ed) / (max.ed - min.ed) * 99) + 1
  
  col.vertices <- paleta[scale.edits]
  
  plot(graf.mapper,
       layout = layout_with_fr(graf.mapper),
       vertex.color = col.vertices,
       vertex.size = log(sapply(wynik.mapper$points_in_vertex, length) + 1) * 2.5,
       vertex.label = NA,
       edge.color = "gray50")
  
  if (table == T)
  {
    prog.odciecia <- quantile(srednie.edycje, probs = 0.8, na.rm = TRUE)
    
    top.wezly <- which(srednie.edycje >= prog.odciecia)
    
    ind.red <- unlist(wynik.mapper$points_in_vertex[top.wezly])
    ind.red <- unique(ind.red)
    
    
    mean.red <- colMeans(dane.num[ind.red, ], na.rm = T)
    mean.rest <- colMeans(dane.num[-ind.red, ], na.rm = T)
    
    diff <- mean.red - mean.rest
    imp.features <- sort(diff, decreasing = T)
    
    kbl(head(imp.features, 10), col.names = c("Zmienna", "Różnica względem tła"), caption = "Tabela zmiennych wykazujących największy wpływ na występowanie edycji oraz różnica średnich tych zmiennych w stosunku do tła biologicznego") %>% kable_styling(latex_options = "H")
    
  }
}















