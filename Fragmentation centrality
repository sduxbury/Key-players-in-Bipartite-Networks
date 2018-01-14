
##====================Fragmentation centrality function in a bipartite network

#=========Intramodal fragmentation is the fragmentation score for the row x row and column by column projections

#========Crossmodal fragmentation is the fragmentation score for first mode actors on the second mode projection and vice versa

##=======Total fragmentation is the average of these two fragmentation scores




#===================================================================================================
#                                            WARNING                                              
#===================================================================================================
##     Function is not optimized for efficiency--computation can take a long time in large networks
#===================================================================================================






##affiliation.mat is a binary affiliation matrix with first mode actors on the rows and second mode actors on the columns

fragmentation.centrality.bipartite<-function(affiliation.mat){
require(statnet)
require(igraph)
require(tnet)

mafia<-as.matrix(affiliation.mat)
fakenet<-graph_from_incidence_matrix(mafia)

fakenet<-fakenet

###intramodal fragmentation--takes ~1 minute
nact<-nrow(mafia)+ncol(mafia)
mafiafrag1<-vector(length=nact)
for (i in 1:nrow(mafia)){
  
  projnet<-fakenet
  fragnet<-delete_vertices(projnet, i) ##delete vertex for focal actor, fragnet is subgraph without i
  
  fragnet<-get.edgelist(fragnet)
  fragnet<-as.tnet(fragnet, type="binary two-mode tnet")
  projfirstmode<-projecting_tm(fragnet, method="Newman") ##project with Newman weights
  distfirstproj<-distance_w(projfirstmode, directed=NULL, gconly=FALSE) #shortest path algorithm
  
  distnet<-(sum(1/distfirstproj, na.rm=TRUE)) ##sum 1/distances, all distances are counted twice, so no need to multiply by two (as in equation)
  distnet<-((nrow(mafia)^2)-nrow(mafia))/distnet ##numerator is total number of potential edges in projection, ensures increasing values indicate increasing vulnerability
  
  mafiafrag1[i]<-distnet}

for (i in nrow(mafia):nact){
  if(i==nrow(mafia)){i<-i+1}
  projnet<-fakenet
  fragnet<-delete_vertices(projnet, i)
  
  fragnet<-t(get.incidence(fragnet))###transpose to project second mode
  projsecmode<-projecting_tm(fragnet, method="Newman")
  distsecproj<-distance_w(projsecmode, directed=NULL, gconly=FALSE)
  
  #projsecmode$w
  
  distnet2<-(sum(1/distsecproj, na.rm=TRUE))
  distnet2<-((ncol(mafia)^2)-ncol(mafia))/distnet2
  
  mafiafrag1[i]<-distnet2
}

intramodalmafiafrag<-mafiafrag1



###cross modal fragmentation

crossmafia<-vector(length=nact)
for (i in nrow(mafia):nact){
  if(i==nrow(mafia)){i<-i+1}
  
  projnet<-fakenet
  fragnet<-delete_vertices(projnet, i) ##delete focal vertex
  
  fragnet<-get.edgelist(fragnet)
  projfirstmode<-projecting_tm(fragnet, method="Newman") ##project first mode
  distfirstproj<-distance_w(projfirstmode, directed=NULL, gconly=FALSE)
  
  
  distnet<-(sum(1/distfirstproj, na.rm=TRUE))
  distnet<-((nrow(mafia)^2)-nrow(mafia))/distnet ###numerator is size of first mode network
  
  crossmafia[i]<-distnet
  if(i==nact){break}
  }

for (i in 1:nrow(mafia)){##first mode
  
  projnet<-fakenet
  fragnet<-delete_vertices(projnet, i)
  
  fragnet<-t(get.incidence(fragnet))##take transpose to project second mode
  projsecmode<-projecting_tm(fragnet, method="Newman")
  distsecproj<-distance_w(projsecmode, directed=NULL, gconly=FALSE)
  
  #projsecmode$w
  
  distnet2<-(sum(1/distsecproj, na.rm=TRUE))
  distnet2<-((ncol(mafia)^2)-ncol(mafia))/distnet2 ##constrain by total potential edges in second mode network
  
  crossmafia[i]<-distnet2
}

crossnodalmafia<-crossmafia

##total fragmentation
totalmafia<-(crossmafia+mafiafrag1)/2

fragmentation.output<-list(intramodalmafiafrag,crossnodalmafia,totalmafia)
names(fragmentation.output)<-c("Intramodal Fragmentation Centrality","Crossmodal Fragmentation Centrality","Overall Fragmentation Centrality")


fragmentation.output
}





###Example--not run
library(igraph)
fakenet<-sample_bipartite(50, 15, type="gnp", runif(1, min=0.1, max=0.4), directed=FALSE)

my.affiliation.matrix<-as.matrix(as_incidence_matrix(fakenet))




fragmentation.centrality.bipartite(my.affiliation.matrix)
