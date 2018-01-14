


#========================Cohesion centrality scores for bipartite networks

##======================intramodal cohesion is the influence of an actor in the actor x actor projection or an event in the event x event projection
##======================Crossmodal cohesion is the influence of an actor in hte event x event projection, or vice versa
##======================Total cohesion is the average of intramodal and cross modal cohesion






#===================================================================================================
#                                            WARNING                                              
#===================================================================================================
##     Function is not optimized for efficiency--computation can take a long time in large networks
#===================================================================================================








##affiliation.mat is an affiliation matrix with the first mode in the rows and the second mode in the columns
cohesion.centrality.bipartite<-function(affiliation.mat){
require(statnet)
require(igraph)
require(tnet)
  
  
davis<-affiliation.mat
projnet<-as.edgelist.sna(davis)


projfirstmode<-projecting_tm(projnet, method="Newman") #projects first mode with distances proposed by Newman
distfirstproj<-distance_w(projfirstmode, directed=NULL, gconly=FALSE) #weights distances with shortest path algorithm

nact<-nrow(davis)+ncol(davis)

daviscoh<-vector(length=nact)
###calculate intramodal cohesion
for (i in 1:nrow(davis)){
  
  distnet<-(sum(1/distfirstproj[i,], na.rm=TRUE)) #sum 1/path lengths for i
  distnet<-distnet/(nrow(davis)) #divide by # of actors in projection
  
  daviscoh[i]<-distnet}


projnet<-t(davis)
projsecmode<-projecting_tm(projnet, method="Newman")
distsecproj<-distance_w(projsecmode, directed=NULL, gconly=FALSE)


for (i in nrow(davis):nact){
  if(i==nrow(davis)){i<-i+1}
  
  distnet2<-(sum(1/distsecproj[i-nrow(davis),], na.rm=TRUE))
  distnet2<-distnet2/(ncol(davis))
  
  daviscoh[i]<-distnet2
}


daviscohesioncentrality<-daviscoh


###cross modal cohesion
crossdav<-vector(length=nact)
davisgraph<-graph_from_incidence_matrix(davis)

a<-matrix(0, nr=1, nc=nrow(davis))

for(i in 1:nrow(davis)){
  b<-neighbors(davisgraph, i)##index neighbors of i in two-mode network
  b<-as.vector(b)
  b<-as.matrix(b)
  
  dimen<-nrow(b)
  c<-matrix(0, nrow=1, ncol=dimen)
  
  for(j in 1:dimen){
    c[1,j]<-daviscoh[b[j,]] ##retrieve neighbors' intramodal cohesion scores
  }
  
  sumscore<-(sum(c, na.rm=TRUE))/ncol(davis) ###sum neighbors'cohesion, divide by number of actors in the alternate mode
  crossdav[i]<-sumscore
}

for(i in nrow(davis):nact){
  if(i==nrow(davis)){i<-i+1}
  b<-neighbors(davisgraph, i)
  b<-as.vector(b)
  b<-as.matrix(b)
  
  dimen<-nrow(b)
  c<-matrix(0, nrow=1, ncol=dimen)
  
  for(j in 1:dimen){
    c[1,j]<-daviscoh[b[j,]]
  }
  
  sumscore<-(sum(c, na.rm=TRUE))/nrow(davis)
  crossdav[i]<-sumscore
}

crosscohesiondavis<-crossdav

totdav<-(crossdav+daviscoh)/2 ###total cohesion

totalcohesiondavis<-totdav


cohesion.centrality.output<-list(daviscohesioncentrality,crosscohesiondavis,totalcohesiondavis)
names(cohesion.centrality.output)<-c("Intramodal cohesion centrality","Crossmodal cohesion centrality","Total cohesion centrality")

cohesion.centrality.output
}



#####Example--not run
library(igraph)
fakenet<-sample_bipartite(50, 15, type="gnp", runif(1, min=0.1, max=0.4), directed=FALSE)

my.affiliation.matrix<-as.matrix(as_incidence_matrix(fakenet))


cohesion.centrality.bipartite(my.affiliation.matrix)
