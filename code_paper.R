# Extracting information from gene coexpression networks of Rhizobium leguminosarum
# Javier Pardo-Diaz et al 2021
# jdiaz@stats.ox.ac.uk

# LOAD REQUIERED PACKAGES

require(igraph)
require(parallel)
require(latex2exp)
require(psych)
require(plot.matrix)
library(zoo)
library(preprocessCore)

load("colour_paletter.RData")


# CODE TO OBTAIN A MAPPING MATRIX INDICATING THE COMMUNITY TO WHICH EACH GENE BELONGS AT EACH OF THE STUDIED PARTITIONS.

# INPUT: igraph network (undirected, weighted). For reference on how to construct a network using signed distance correlation: https://github.com/javier-pardodiaz/sdcorGCN

## Previously, run the code in comm_det.ipynb which takes as input the network in GML format

### To export the network in GML format

write_graph(network,file,format =  "gml")

## Now, obtain the mapping matrix "mapping_comms*

info_comms=read.csv2("Comms_info.txt")
number_partitions=nrow(info_comms)

load_partition=function(res){
  communities<-read.csv2(paste("partition_",res,".txt",sep=""), header=FALSE)
  v=rep(NA,vcount(network))
  for (j in c(1:length(communities))){
    v[as.vector(strsplit(as.character(communities[[j]]),",")[[1]],mode="numeric")+1]=j
  }
  return(v)
}

load_community=function(j,community,network){
  nodes=V(WSwS)$name[as.vector(strsplit(as.character(community[[1]]),",")[[1]],mode="numeric")+1]
  return(nodes)
}

mapping_comms=sapply(c(0:(number_partitions-1)),load_partition)
colnames(mapping_comms_2)=c(1:number_partitions)
rownames(mapping_comms_2)=V(WSwS)$name


# FUNCTIONS TO COMPUTE S SCORES (those involving a seed list)

CD_score_function_seedlist=function(network,seed_list,map){
  # Arguments: network, seed_list, and a map containing the list of partitions
  i=V(network)$name
  return(sapply(i, function(i) get_comm_det_score(i,map,seed_list,network)))
} 

get_comm_det_score=function(gene,mapping_comms,seed_list,network){
  j=c(1:ncol(mapping_comms))
  return(sum(sapply(j, function(j) get_comm_det_score_part_adj(gene,network,mapping_comms[,j],seed_list)),na.rm=TRUE))
}

get_comm_det_score_part_adj=function(gene,network,vector,seed_list){
  return((sum(vector[seed_list[seed_list!=gene]]==vector[gene]))/((sum(vector==vector[gene])-1*(gene%in%seed_list))*length(seed_list)))
}

PR_score_function_seedlist=function(network,seed_list){
  PR=page_rank(network)$vector
  a=sapply(seed_list, function(seed_list) page_rank_seed_list(network,seed_list))
  a[sapply(colnames(a), `==`, rownames(a))] <- NA
  a=apply(a,1,mean,na.rm=TRUE)
  return(a/PR)
}

page_rank_seed_list=function(network,seed_list){
  v=rep(0,vcount(network))
  v[V(network)$name%in%seed_list]=1/length(seed_list)
  return(page_rank(network,personalized=v)$vector)
}

EW_score_function_seedlist=function(network,seed_list){
  adj_matrix_NA_diag=network[V(network)$name,V(network)$name]
  diag(adj_matrix_NA_diag)=NA
  if(length(seed_list)>1){
    a=apply(adj_matrix_NA_diag[V(network)$name,seed_list],1,mean,na.rm=TRUE)
  }else{
    a=adj_matrix_NA_diag[V(network)$name,seed_list]
  }
  return(a)
}


# FUNCTIONS TO COMPUTE Z SCORES (those involiving a numeric vector such as gene expression)


CD_score_function_vector=function(network,seed_list,map,names){
  names(seed_list)=names
  i=V(network)$name
  return(sapply(i, function(i) get_comm_det_score_cont(i,map,seed_list,network)))
}


get_comm_det_score_cont=function(gene,mapping_comms,seed_list,network){
  j=c(1:ncol(mapping_comms))
  return(sum(sapply(j, function(j) get_comm_det_score_part_adj_continuos(gene,network,mapping_comms[,j],seed_list))))
}


get_comm_det_score_part_adj_continuos=function(gene,network,vector,seed_list){
  return((sum(seed_list[names(vector[vector==vector[gene]])],na.rm = TRUE))*1000/(sum(vector==vector[gene])*sum(seed_list,na.rm = TRUE)))
}

PR_score_function_vector=function(v,network){
  return(page_rank(network,personalized = v[V(network)$name])$vector)
}

EW_score_function_vector=function(network,seed_list){
  
  a=apply(WSwS[V(WSwS)$name,V(WSwS)$name],1,multiply_edgeweight_expression,seed_list)
  
  return(a)
}

multiply_edgeweight_expression=function(weights,expression){
  return(sum(weights*expression[names(weights)])/(length(weights)-1))
}


# FUNCTION TO SELECT AND PLOT THE GENES IN THE TOP M FOR AT LEAST ONE OF THE RANKINGS AND T FOR ALL OF THEM, AND THOSE IN THE TOP S FOR ANY OF THEM


selection_top=function(CD,PR,WE,T=100,S=10,network,M=50,CS=NULL){
  inter=T
  union=S
  n=M
  n2=M
  CD=CD[names(CD)%in%CS==FALSE]
  PR=PR[names(PR)%in%CS==FALSE]
  WE=WE[names(WE)%in%CS==FALSE]
  CD_top_inter=names(sort(CD,decreasing = TRUE))[1:inter]
  PR_top_inter=names(sort(PR,decreasing = TRUE))[1:inter]
  WE_top_inter=names(sort(WE,decreasing = TRUE))[1:inter]
  CD_top_union=names(sort(CD,decreasing = TRUE))[1:union]
  PR_top_union=names(sort(PR,decreasing = TRUE))[1:union]
  WE_top_union=names(sort(WE,decreasing = TRUE))[1:union]
  genes=sort(unique(c(intersect(intersect(CD_top_inter,WE_top_inter),PR_top_inter),CD_top_union,PR_top_union,WE_top_union)))
  genes=genes[(genes%in%names(sort(CD,decreasing = TRUE))[1:n2])|(genes%in%names(sort(PR,decreasing = TRUE))[1:n2])|(genes%in%names(sort(WE,decreasing = TRUE))[1:n2])]
  sg=induced_subgraph(network,genes)
  V(sg)$color="white"
  V(sg)$color[V(sg)$name%in%names(sort(CD,decreasing = TRUE))[1:n]]=palette_ch4["CD"]
  V(sg)$color[V(sg)$name%in%names(sort(PR,decreasing = TRUE))[1:n]]=palette_ch4["PR"]
  V(sg)$color[V(sg)$name%in%names(sort(WE,decreasing = TRUE))[1:n]]=palette_ch4["Weight"]
  V(sg)$color[V(sg)$name%in%names(sort(CD,decreasing = TRUE))[1:n]&V(sg)$name%in%names(sort(PR,decreasing = TRUE))[1:n]]=palette_ch4["CD+PR"]
  V(sg)$color[V(sg)$name%in%names(sort(CD,decreasing = TRUE))[1:n]&V(sg)$name%in%names(sort(WE,decreasing = TRUE))[1:n]]=palette_ch4["CD+W"]
  V(sg)$color[V(sg)$name%in%names(sort(PR,decreasing = TRUE))[1:n]&V(sg)$name%in%names(sort(WE,decreasing = TRUE))[1:n]]=palette_ch4["PR+W"]
  V(sg)$color[V(sg)$name%in%names(sort(CD,decreasing = TRUE))[1:n]&V(sg)$name%in%names(sort(PR,decreasing = TRUE))[1:n]&V(sg)$name%in%names(sort(WE,decreasing = TRUE))[1:n]]=palette_ch4["All2"]
  
  plot(sg)
  order_CD=match(genes,names(sort(CD,decreasing = TRUE)))
  order_PR=match(genes,names(sort(PR,decreasing = TRUE)))
  order_WE=match(genes,names(sort(WE,decreasing = TRUE)))
  
  return(list("genes"=genes,"ranks_CD"=order_CD,"ranks_PR"=order_PR,"ranks_WE"=order_WE))
}


# FUNCTIONS TO COMPARE CONDITIONS


selection_top_conditions=function(cond,T=20,S=5){
  CD_top_inter=names(cond$CD)[1:inter]
  PR_top_inter=names(cond$PR)[1:inter]
  WE_top_inter=names(cond$WE)[1:inter]
  CD_top_union=names(cond$CD)[1:union]
  PR_top_union=names(cond$PR)[1:union]
  WE_top_union=names(cond$WE)[1:union]
  genes=sort(unique(c(intersect(intersect(CD_top_inter,WE_top_inter),PR_top_inter),CD_top_union,PR_top_union,WE_top_union)))
  return(genes)
}



top_genes_condition_vs_condition=function(CD,PR,WE,Exp,c1,c2,c3=NULL,c4=NULL,m=NULL,n=NULL,sorted=TRUE,rm_low=0.2){
  if (is.null(n)==TRUE){
    n=nrow(CD)
  }
  if (is.null(m)==TRUE){
    m=1
  }
  if(is.null(c3)==TRUE){
    FC_Ordered_genes_c1_c2=(Exp[,c1]/Exp[,c2])[m:n]
    CD_Ordered_genes_c1_c2=(CD[,c1]/CD[,c2])[m:n]
    PR_Ordered_genes_c1_c2=(PR[,c1]/PR[,c2])[m:n]
    WE_Ordered_genes_c1_c2=(WE[,c1]/WE[,c2])[m:n]
  }else{
    FC_Ordered_genes_c1_c2=(Exp[,c1]/Exp[,c2]*Exp[,c3]/Exp[,c4])[m:n]
    CD_Ordered_genes_c1_c2=(CD[,c1]/CD[,c2]*CD[,c3]/CD[,c4])[m:n]
    PR_Ordered_genes_c1_c2=(PR[,c1]/PR[,c2]*PR[,c3]/PR[,c4])[m:n]
    WE_Ordered_genes_c1_c2=(WE[,c1]/WE[,c2]*WE[,c3]/WE[,c4])[m:n]
  }
  thr=quantile(Exp[,c1],rm_low)
  FC_Ordered_genes_c1_c2[Exp[,c1]<thr]=NA
  CD_Ordered_genes_c1_c2[Exp[,c1]<thr]=NA
  PR_Ordered_genes_c1_c2[Exp[,c1]<thr]=NA
  WE_Ordered_genes_c1_c2[Exp[,c1]<thr]=NA
  
  if (sorted==TRUE){
    FC_Ordered_genes_c1_c2=(sort(FC_Ordered_genes_c1_c2,decreasing = TRUE))[m:n]
    CD_Ordered_genes_c1_c2=(sort(CD_Ordered_genes_c1_c2,decreasing = TRUE))[m:n]
    PR_Ordered_genes_c1_c2=(sort(PR_Ordered_genes_c1_c2,decreasing = TRUE))[m:n]
    WE_Ordered_genes_c1_c2=(sort(WE_Ordered_genes_c1_c2,decreasing = TRUE))[m:n]
  }
  return(list("FC"=FC_Ordered_genes_c1_c2,"CD"=CD_Ordered_genes_c1_c2,"PR"=PR_Ordered_genes_c1_c2,"WE"=WE_Ordered_genes_c1_c2))
}

plot_intersections_selection=function(conditions,U,n2=U,cnames=NULL,selected_genes=NULL){
  genes=selected_genes
  
  genes=unique(genes)
  genes=sort(genes)
  m=matrix(nrow=length(genes),ncol=length(conditions))
  rownames(m)=genes
  colnames(m)=cnames
  for (i in c(1:length(conditions))){
    v=rep(8,length(genes))
    CD=genes%in%names(conditions[[i]]$CD)[1:n2]
    WE=genes%in%names(conditions[[i]]$WE)[1:n2]
    PR=genes%in%names(conditions[[i]]$PR)[1:n2]
    v[CD]=1
    v[WE]=2
    v[PR]=3
    v[CD&PR]=4
    v[CD&WE]=5
    v[WE&PR]=6
    v[CD&WE&PR]=7
    m[,i]=v
  }
  par(mfrow=c(1,3))
  #plot(induced_subgraph(WSwS,genes),vertex.size=5)
  par(mar=c(11, 8, 1, 4.1))   # adapt margins
  plot(m[c(1:ceiling(length(genes)/3)),], breaks=c(1:9),col=c(palette_ch4[c("CD", "Weight",'PR',"CD+PR", "CD+W","PR+W","All1")],"white"),las=2,key=NULL,ylab="",xlab="",main="")
  plot(m[c((ceiling(length(genes)/3)+1):(2*ceiling(length(genes)/3))),], breaks=c(1:9),col=c(palette_ch4[c("CD", "Weight",'PR',"CD+PR", "CD+W","PR+W","All1")],"white"),las=2,key=NULL,ylab="",xlab="",main="")
  
  plot(m[c(((2*ceiling(length(genes)/3))+1):length(genes)),], breaks=c(1:9),col=c(palette_ch4[c("CD", "Weight",'PR',"CD+PR", "CD+W","PR+W","All1")],"white"),las=2,key=NULL,ylab="",xlab="",main="")
  
  #plot_gene_pattern(genes)
  
}


# ANALYSIS OF THE PERFORMANCE OF THE SEED LIST SCORES

control_seed=c("RL0106", "RL1773", "RL1797", "RL1769" ,"RL1796", "RL1787", "RL4549", "RL1783",
               "RL1554" ,"RL1778", "RL2221", "RL0374" ,"RL4130", "RL1780", "RL2624", "RL1791",
               "RL1555" ,"RL1770", "RL1788", "RL1672" ,"RL1762", "RL1764", "RL1761", "RL1784",
               "RL1793" ,"RL1781", "RL1799", "RL1790" ,"RL4552", "RL1777", "RL0268", "RL4676",
               "RL1779" ,"RL1776", "RL1785", "RL3471" ,"RL4677", "RL1782", "RL1774", "RL1792",
               "RL4017" ,"RL4624", "RL1731", "RL0267", "RL1775" ,"RL1786", "RL1789", "RL1765",
               "RL1552"
)

# We use 25 random seed lists of size 5

random_control_seeds=list()
for(i in c(1:25)){
  set.seed(i)
  random_control_seeds[[i]]=sample(control_seed,5)
}


## STD NETWORK

load("RL3841_std_mapping_comms.RData")
load("RL3841_std_WSwS.RData")

### Obtain random rank

random_ranks_std=list()
for(i in c(1:25)){
  set.seed(i)
  random_ranks_std[[i]]=sample(V(WSwS)$name,vcount(WSwS),replace = FALSE)
}

### Obtain scores

cl1= makeCluster(5, "FORK")
CD_scores_ribosomal_iterations_5_std<<-parSapply(cl1,random_control_seeds,function(random_control_seeds) CD_score_function_seedlist(network=WSwS,seed_list = random_control_seeds, map=mapping_comms))
stopCluster(cl1)

PR=page_rank(WSwS)$vector
PR_scores_ribosomal_iterations_5_std=sapply(random_control_seeds,function(random_control_seeds) PR_score_function_seedlist(network=WSwS,seed_list = random_control_seeds,PR))

WE_scores_ribosomal_iterations_5_std=sapply(random_control_seeds,function(random_control_seeds) EW_score_function_seedlist(network=WSwS,seed_list = random_control_seeds))

### Plot ROC

plot(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[1]]],names(sort(PR_scores_ribosomal_iterations_5_std[,1],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[1]]])/(length(control_seed[!control_seed%in% random_control_seeds[[1]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["PR"], xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",main="Recovery of ribosomal genes \n Complete Network")

for (i in c(1:25)){
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(PR_scores_ribosomal_iterations_5_std[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["PR"])
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(CD_scores_ribosomal_iterations_5_std[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["CD"])
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(WE_scores_ribosomal_iterations_5_std[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["Weight"])
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],random_ranks_std[[i]]),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["All2"],lty=3)
  
}
legend(0.65,0.5,col=c(palette_ch4["CD"],palette_ch4["PR"],palette_ch4["Weight"],palette_ch4["All2"]),legend=c("CommDet","PageRank","EdgeWeight","RandomOrder"),lwd=2,lty=c(1,1,1,3))

### Compute AUROC and boxplot

AUC_CD_std=c()
AUC_PR_std=c()
AUC_WE_std=c()
AUC_RN_std=c()
for (i in c(1:25)){
  AUC_PR_std=c(AUC_PR_std, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(PR_scores_ribosomal_iterations_5_std[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  AUC_CD_std=c(AUC_CD_std, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(CD_scores_ribosomal_iterations_5_std[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  AUC_WE_std=c(AUC_WE_std, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(WE_scores_ribosomal_iterations_5_std[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  AUC_RN_std=c(AUC_RN_std, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],random_ranks_std[[i]]),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  
}

par(mfrow=c(2,2))
boxplot(list(AUC_CD_std,AUC_PR_std,AUC_WE_std,AUC_RN_std),names=c("CommDet","PageRank","EdgeWeight","RandomOrder"),col=c(palette_ch4["CD"],palette_ch4["PR"],palette_ch4["Weight"],palette_ch4["All2"]),ylab="Area under the curve",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
boxplot(list(AUC_CD_std,AUC_PR_std,AUC_WE_std),names=c("CommDet","PageRank","EdgeWeight"),col=c(palette_ch4["CD"],palette_ch4["PR"],palette_ch4["Weight"]),main="",ylab="Area under the curve",cex.axis=1.5,cex.lab=1.5)

## CPL NETWORK

load("RL3841_cpl_mapping_comms.RData")
load("RL3841_cpl_WSwS.RData")

### Obtain random rank

random_ranks_cpl=list()
for(i in c(1:25)){
  set.seed(i)
  random_ranks_cpl[[i]]=sample(V(WSwS)$name,vcount(WSwS),replace = FALSE)
}

### Obtain scores

cl1= makeCluster(5, "FORK")
CD_scores_ribosomal_iterations_5_cpl<<-parSapply(cl1,random_control_seeds,function(random_control_seeds) CD_score_function_seedlist(network=WSwS,seed_list = random_control_seeds, map=mapping_comms))
stopCluster(cl1)

PR=page_rank(WSwS)$vector
PR_scores_ribosomal_iterations_5_cpl=sapply(random_control_seeds,function(random_control_seeds) PR_score_function_seedlist(network=WSwS,seed_list = random_control_seeds,PR))

WE_scores_ribosomal_iterations_5_cpl=sapply(random_control_seeds,function(random_control_seeds) EW_score_function_seedlist(network=WSwS,seed_list = random_control_seeds))

### Plot ROC

plot(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[1]]],names(sort(PR_scores_ribosomal_iterations_5_cpl[,1],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[1]]])/(length(control_seed[!control_seed%in% random_control_seeds[[1]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["PR"], xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",main="Recovery of ribosomal genes \n Complete Network")

for (i in c(1:25)){
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(PR_scores_ribosomal_iterations_5_cpl[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["PR"])
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(CD_scores_ribosomal_iterations_5_cpl[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["CD"])
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(WE_scores_ribosomal_iterations_5_cpl[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["Weight"])
  points(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],random_ranks_cpl[[i]]),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1),c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),type="s",xlim=c(0,1),ylim=c(0,1),col=palette_ch4["All2"],lty=3)
}
legend(0.65,0.5,col=c(palette_ch4["CD"],palette_ch4["PR"],palette_ch4["Weight"],palette_ch4["All2"]),legend=c("CommDet","PageRank","EdgeWeight","RandomOrder"),lwd=2,lty=c(1,1,1,3))

### Compute AUROC and boxplot

AUC_CD_cpl=c()
AUC_PR_cpl=c()
AUC_WE_cpl=c()
AUC_RN_cpl=c()

for (i in c(1:25)){
  AUC_PR_cpl=c(AUC_PR_cpl, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(PR_scores_ribosomal_iterations_5_cpl[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  AUC_CD_cpl=c(AUC_CD_cpl, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(CD_scores_ribosomal_iterations_5_cpl[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  AUC_WE_cpl=c(AUC_WE_cpl, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],names(sort(WE_scores_ribosomal_iterations_5_cpl[,i],decreasing = T))),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  AUC_RN_cpl=c(AUC_RN_cpl, sum(diff(c(0,sort(match(control_seed[!control_seed%in% random_control_seeds[[i]]],random_ranks_cpl[[i]]),decreasing = F)/(vcount(WSwS)-(length(control_seed)-length(random_control_seeds[[1]]))),1))*rollmean(c(0:length(control_seed[!control_seed%in% random_control_seeds[[i]]])/(length(control_seed[!control_seed%in% random_control_seeds[[i]]])),1),2)))
  
}

boxplot(list(AUC_CD_cpl,AUC_PR_cpl,AUC_WE_cpl,AUC_RN_cpl),names=c("CommDet","PageRank","EdgeWeight","RandomOrder"),col=c(palette_ch4["CD"],palette_ch4["PR"],palette_ch4["Weight"],palette_ch4["All2"]),ylab="Area under the curve",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
boxplot(list(AUC_CD_cpl,AUC_PR_cpl,AUC_WE_cpl),names=c("CommDet","PageRank","EdgeWeight"),col=c(palette_ch4["CD"],palette_ch4["PR"],palette_ch4["Weight"]),main="",ylab="Area under the curve",cex.axis=1.5,cex.lab=1.5)


# STUDY OF RIBOSOMAL GENES FOR THE STD AND CPL NETWORKS

control_seed=c("RL0106", "RL1773", "RL1797", "RL1769" ,"RL1796", "RL1787", "RL4549", "RL1783",
               "RL1554" ,"RL1778", "RL2221", "RL0374" ,"RL4130", "RL1780", "RL2624", "RL1791",
               "RL1555" ,"RL1770", "RL1788", "RL1672" ,"RL1762", "RL1764", "RL1761", "RL1784",
               "RL1793" ,"RL1781", "RL1799", "RL1790" ,"RL4552", "RL1777", "RL0268", "RL4676",
               "RL1779" ,"RL1776", "RL1785", "RL3471" ,"RL4677", "RL1782", "RL1774", "RL1792",
               "RL4017" ,"RL4624", "RL1731", "RL0267", "RL1775" ,"RL1786", "RL1789", "RL1765",
               "RL1552"
)

## STD NETWORK

load("RL3841_std_mapping_comms.RData")
load("RL3841_std_WSwS.RData")

scores_CD_control_std=CD_score_function_seedlist(network=WSwS,seed_list=control_seed,map=mapping_comms)
scores_PR_control_std=PR_score_function_seedlist(network = WSwS,seed_list = control_seed)
scores_WE_control_std=EW_score_function_seedlist(network = WSwS,seed_list = control_seed)

ribosome_selection_top_std=selection_top(CD=scores_CD_control_std,PR=scores_PR_control_std,WE=scores_WE_control_std,network=WSwS,CS=control_seed, T=100, S=10, M=50)
print(ribosome_selection_top_std$genes)
print(ribosome_selection_top_std$ranks_CD)
print(ribosome_selection_top_std$ranks_PR)
print(ribosome_selection_top_std$ranks_WE)

## CPL NETWORK

load("RL3841_cpl_mapping_comms.RData")
load("RL3841_cpl_WSwS.RData")

scores_CD_control_cpl=CD_score_function_seedlist(network=WSwS,seed_list=control_seed,map=mapping_comms)
scores_PR_control_cpl=PR_score_function_seedlist(network = WSwS,seed_list = control_seed)
scores_WE_control_cpl=EW_score_function_seedlist(network = WSwS,seed_list = control_seed)

ribosome_selection_top_cpl=selection_top(CD=scores_CD_control_cpl,PR=scores_PR_control_cpl,WE=scores_WE_control_cpl,network=WSwS,CS=control_seed, T=100, S=10, M=50)
print(ribosome_selection_top_cpl$genes)
print(ribosome_selection_top_cpl$ranks_CD)
print(ribosome_selection_top_cpl$ranks_PR)
print(ribosome_selection_top_cpl$ranks_WE)

# STUDY OF DIFFERENT GROWTH CONDITIONS

## We take as an input a matrix (expression_matrix_conditions) that has in its rows all the genes in the network and in the columns the values for all the microarrays studied (both channels)

load("expression_matrix_conditions.RData")

expression_matrix_conditions_qnorm=normalize.quantiles(expression_matrix_conditions)
rownames(expression_matrix_conditions_qnorm)=rownames(expression_matrix_conditions)
colnames(expression_matrix_conditions_qnorm)=colnames(expression_matrix_conditions)

## Compute scores for the standard network

load("RL3841_std_mapping_comms.RData")
load("RL3841_std_WSwS.RData")

cl1= makeCluster(4, "FORK")
CD_scores_conditions_std=parApply(cl1,expression_matrix_conditions_qnorm, 2, CD_score_function_vector,network=WSwS,map=mapping_comms,names=rownames(expression_matrix_conditions))
stopCluster(cl1)

PR_scores_conditions_std=apply(expression_matrix_conditions_qnorm, 2, PR_score_function_vector,network=WSwS)

WE_scores_conditions_std=apply(expression_matrix_conditions_qnorm, 2, EW_score_function_vector,network=WSwS)

## Refer the scores of conditions to the controls

Pyr_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=1,c2=2,rm_low = 0.2)
Succ_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=3,c2=4,rm_low = 0.2)
Acetoacetate_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=5,c2=6,rm_low = 0.2)
Protocatechuate_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=7,c2=8,c3=1,c4=2,rm_low = 0.2)
Acetate_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=13,c2=14,rm_low = 0.2)
Formate_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=24,c2=26,c3=1,c4=2,rm_low = 0.2)
Hydroxybenzoate_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=39,c2=40,c3=1,c4=2,rm_low = 0.2)
Ara_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=7,c2=8,c3=43,c4=44,rm_low = 0.2)
Inositol_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=11,c2=12,rm_low = 0.2)
Gal_vs_Glu_std=top_genes_condition_vs_condition(CD=CD_scores_conditions_std,PR=PR_scores_conditions_std,WE=WE_scores_conditions_std,Exp=expression_matrix_conditions_qnorm, c1=9,c2=10,rm_low = 0.2)

## Selection of the genes in each of the conditions (using the parameters described in the paper)

selected_genes_Pyr_vs_Glu_std=selection_top_conditions(Pyr_vs_Glu_std,T=20,S=5)
selected_genes_Succ_vs_Glu_std=selection_top_conditions(Succ_vs_Glu_std,T=20,S=5)
selected_genes_Acetoacetate_vs_Glu_std=selection_top_conditions(Acetoacetate_vs_Glu_std,T=20,S=5)
selected_genes_Protocatechuater_vs_Glu_std=selection_top_conditions(Protocatechuate_vs_Glu_std,T=20,S=5)
selected_genes_Acetate_vs_Glu_std=selection_top_conditions(Acetate_vs_Glu_std,T=20,S=5)
selected_genes_Formate_vs_Glu_std=selection_top_conditions(Formate_vs_Glu_std,T=20,S=5)
selected_genes_Hydroxybenzoate_vs_Glu_std=selection_top_conditions(Hydroxybenzoate_vs_Glu_std,T=20,S=5)
selected_genes_Ara_vs_Glu_std=selection_top_conditions(Ara_vs_Glu_std,T=20,S=5)
selected_genes_Inositol_vs_Glu_std=selection_top_conditions(Inositol_vs_Glu_std,T=20,S=5)
selected_genes_Gal_vs_Glu_std=selection_top_conditions(Gal_vs_Glu_std,T=20,S=5)

selected_genes_all_std=unique(c(selected_genes_Gal_vs_Glu_std,selected_genes_Inositol_vs_Glu_std,selected_genes_Pyr_vs_Glu_std,selected_genes_Succ_vs_Glu_std,selected_genes_Acetoacetate_vs_Glu_std,selected_genes_Protocatechuater_vs_Glu_std,selected_genes_Acetate_vs_Glu_std,selected_genes_Formate_vs_Glu_std,selected_genes_Hydroxybenzoate_vs_Glu_std,selected_genes_Ara_vs_Glu_std))

## Plot of genes across the conditions

plot_intersections_selection(list(Pyr_vs_Glu_std,Succ_vs_Glu_std,Acetate_vs_Glu_std,Acetoacetate_vs_Glu_std,Protocatechuate_vs_Glu_std,Formate_vs_Glu_std,Hydroxybenzoate_vs_Glu_std,Ara_vs_Glu_std,Inositol_vs_Glu_std,Gal_vs_Glu_std),U=50,cnames=c("Pyr_vs_Glu","Succ_vs_Glu","Acetate_vs_Glu","Acetoacetate_vs_Glu","Protocatechuate_vs_Glu","Formate_vs_Glu","Hydroxybenzoate_vs_Glu","Ara_vs_Glu","Inositol_vs_Glu","Gal_vs_Glu"),selected_genes=selected_genes_all_std)


## Compute scores for the complete network

load("RL3841_cpl_mapping_comms.RData")
load("RL3841_cpl_WSwS.RData")

cl1= makeCluster(4, "FORK")
CD_scores_conditions_cpl=parApply(cl1,expression_matrix_conditions_qnorm, 2, CD_score_function_vector,network=WSwS,map=mapping_comms,names=rownames(expression_matrix_conditions))
stopCluster(cl1)

PR_scores_conditions_cpl=apply(expression_matrix_conditions_qnorm, 2, PR_score_function_vector,network=WSwS)

WE_scores_conditions_cpl=apply(expression_matrix_conditions_qnorm, 2, EW_score_function_vector,network=WSwS)

## Refer the scores of conditions to the controls

Pyr_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=1,c2=2,rm_low = 0.2)
Succ_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=3,c2=4,rm_low = 0.2)
Acetoacetate_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=5,c2=6,rm_low = 0.2)
Protocatechuate_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=7,c2=8,c3=1,c4=2,rm_low = 0.2)
Acetate_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=13,c2=14,rm_low = 0.2)
Formate_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=24,c2=26,c3=1,c4=2,rm_low = 0.2)
Hydroxybenzoate_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=39,c2=40,c3=1,c4=2,rm_low = 0.2)
Ara_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=7,c2=8,c3=43,c4=44,rm_low = 0.2)
Inositol_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=11,c2=12,rm_low = 0.2)
Gal_vs_Glu_cpl=top_genes_condition_vs_condition(CD=CD_scores_conditions_cpl,PR=PR_scores_conditions_cpl,WE=WE_scores_conditions_cpl,Exp=expression_matrix_conditions_qnorm, c1=9,c2=10,rm_low = 0.2)

## Selection of the genes in each of the conditions (using the parameters described in the paper)

selected_genes_Pyr_vs_Glu_cpl=selection_top_conditions(Pyr_vs_Glu_cpl,T=20,S=5)
selected_genes_Succ_vs_Glu_cpl=selection_top_conditions(Succ_vs_Glu_cpl,T=20,S=5)
selected_genes_Acetoacetate_vs_Glu_cpl=selection_top_conditions(Acetoacetate_vs_Glu_cpl,T=20,S=5)
selected_genes_Protocatechuater_vs_Glu_cpl=selection_top_conditions(Protocatechuate_vs_Glu_cpl,T=20,S=5)
selected_genes_Acetate_vs_Glu_cpl=selection_top_conditions(Acetate_vs_Glu_cpl,T=20,S=5)
selected_genes_Formate_vs_Glu_cpl=selection_top_conditions(Formate_vs_Glu_cpl,T=20,S=5)
selected_genes_Hydroxybenzoate_vs_Glu_cpl=selection_top_conditions(Hydroxybenzoate_vs_Glu_cpl,T=20,S=5)
selected_genes_Ara_vs_Glu_cpl=selection_top_conditions(Ara_vs_Glu_cpl,T=20,S=5)
selected_genes_Inositol_vs_Glu_cpl=selection_top_conditions(Inositol_vs_Glu_cpl,T=20,S=5)
selected_genes_Gal_vs_Glu_cpl=selection_top_conditions(Gal_vs_Glu_cpl,T=20,S=5)

selected_genes_all_cpl=unique(c(selected_genes_Gal_vs_Glu_cpl,selected_genes_Inositol_vs_Glu_cpl,selected_genes_Pyr_vs_Glu_cpl,selected_genes_Succ_vs_Glu_cpl,selected_genes_Acetoacetate_vs_Glu_cpl,selected_genes_Protocatechuater_vs_Glu_cpl,selected_genes_Acetate_vs_Glu_cpl,selected_genes_Formate_vs_Glu_cpl,selected_genes_Hydroxybenzoate_vs_Glu_cpl,selected_genes_Ara_vs_Glu_cpl))

## Plot of genes across the conditions

plot_intersections_selection(list(Pyr_vs_Glu_cpl,Succ_vs_Glu_cpl,Acetate_vs_Glu_cpl,Acetoacetate_vs_Glu_cpl,Protocatechuate_vs_Glu_cpl,Formate_vs_Glu_cpl,Hydroxybenzoate_vs_Glu_cpl,Ara_vs_Glu_cpl,Inositol_vs_Glu_cpl,Gal_vs_Glu_cpl),U=50,cnames=c("Pyr_vs_Glu","Succ_vs_Glu","Acetate_vs_Glu","Acetoacetate_vs_Glu","Protocatechuate_vs_Glu","Formate_vs_Glu","Hydroxybenzoate_vs_Glu","Ara_vs_Glu","Inositol_vs_Glu","Gal_vs_Glu"),selected_genes=selected_genes_all_cpl)




