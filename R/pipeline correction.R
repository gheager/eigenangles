library(magrittr)
library(stringr)
library(purrr)
library(ggplot2)
library(igraph)
library(gtools)

#batch effect correction packages
library(sva) #for ComBat
library(RUVSeq) #for RUVs
library(batchelor) #for mnnCorrect

load_experiments<-function(directory, item.SimpleList='rnaseq'){
  directory %>% dir %>% map(~get(load(paste0(directory,'/',.x))) %>% (if(class(.)=='SimpleList') as_mapper(~.[[item.SimpleList]]) else identity)) %>% 
    set_names(directory %>% dir %>% str_split('-') %>% map(~.x[2:3] %>% paste(collapse='')) %>% unlist)
}

download_experiments_from_ExpressionAtlas<-function(..., destdir=getwd() %>% paste('experiments',sep='/')){
  if(destdir %>% dir.exists){
    stop('Attempted to create directory `',desdir,'`` but this address already exists. Rename it or modify `destdir` argument in `download_experiments`.')
  }else{
    destdir %>% dir.creates
  }
  for(experiment in list(...)){
    paste0('https://www.ebi.ac.uk/gxa/experiments-content/',experiment,'/static/',experiment,'-atlasExperimentSummary.Rdata') %>% download.file(destdir)
  }
  desdir %>% load_experiments
}

remove_isolated_experiments<-function(experiments, biological.group){
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[biological.group]]) %>% unlist(use.names=FALSE)
  groups <- group %>% split(batch)
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  intersections %>% graph_from_adjacency_matrix %>% plot
  experiments[rownames(intersections)[rowSums(intersections!=0)<=1]]<-NULL
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[biological.group]]) %>% unlist(use.names=FALSE)
  groups <- group %>% split(batch)
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  intersections %>% graph_from_adjacency_matrix %>% plot
  return(experiments)
}

merge_experiments<-function(experiments, log, filter.unexpressed.genes){
  UseMethod("merge_experiments")
}

merge_experiments.SummarizedExperiment <- function(experiments, log=TRUE, filter.unexpressed.genes=TRUE){
  genes<-experiments %>% map(rownames)
  shared.genes<-genes %>% purrr::reduce(intersect)
  unshared.genes<-genes %>% map(setdiff %>% partial(y=shared.genes))
  data<-experiments %>% map(~.x@assays$data$counts %>% extract(rownames(.x)%in%shared.genes,)) %>% purrr::reduce(cbind)
  warning(unshared.genes %>% unlist %>% unique %>% length,' genes have been removed as they are not shared across the batches.')
  batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
  if(filter.unexpressed.genes){
    unexpressed.genes <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)==0) %>% purrr::reduce(`&`)
    data%<>%extract(!unexpressed.genes,)
    message(sum(unexpressed.genes),' genes have been removed as they were unexpressed across the samples of a batch.')
  }
  if(log) data%<>%log1p
  return(SummarizedExperiment(
    assays=if(log) list(log_counts=data) else list(counts=data),
    colData=experiments %>% map(colData) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
    metadata=list(batch=batch)
  ))
}

merge_experiments.ExpressionSet <- function(experiments, log=FALSE, filter.unexpressed.genes=TRUE){
  genes<-experiments %>% map(rownames)
  shared.genes<-genes %>% purrr::reduce(intersect)
  unshared.genes<-genes %>% map(setdiff %>% partial(y=shared.genes))
  data<-experiments %>% map(~.x@assayData$exprs %>% extract(rownames(.x)%in%shared.genes,)) %>% purrr::reduce(cbind)
  warning(unshared.genes %>% unlist %>% unique %>% length,' genes have been removed as they are not shared across the batches.')
  batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
  if(filter.unexpressed.genes){
    unexpressed.genes <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)==0) %>% purrr::reduce(`&`)
    data%<>%extract(!unexpressed.genes,)
    message(sum(unexpressed.genes),' genes have been removed as they were unexpressed across the samples of a batch.')
  }
  if(log) data%<>%log1p
  return(ExpressionSet(
    assayData=if(log) list(log_exprs=data) else list(exprs=data),
    phenoData=experiments %>% map(phenoData) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
    experimentData=list(batch=batch)
  ))
}

correct_batch_effect<-function(experiment, model, method=c('ComBat','RUV','MNN'), k){
  UseMethod("correct_batch_effect")
}

correct_batch_effect.SummarizedExperiment<-function(experiment, model, method=c('ComBat','RUV','MNN'), k){
  log<-experiment@assays$data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
  model.data<-model.frame(model, experiment@colData[all.vars(model)])
  if(length(k)==1|method=='ComBat'){
    return(SummarizedExperiment(
      assays = list(switch(
        method,
        ComBat = ComBat(experiment@assays$data[[1]], experiment$batch, mod=model.matrix(model, data=model.data)),
        RUV = RUVs(experiment@assays$data[[1]], cIdx=seq_len(nrow(experiment@assays$data[[1]])), k=k, 
                   scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=log)$normalizedCounts,
        MNN = mnnCorrect(experiment@assays$data[[1]], batch=experiment$batch, k=k)@assays$data$corrected
      )) %>% set_names(if(log) 'corrected_log_counts' else 'corrected_counts'),
      colData = experiment@colData,
      metadata = experiment@metadata
    ))
  }else{
    return(k %>% map(correct_batch_effect %>% partial(experiment=experiment, model=model, method=method)))
  }
}

correct_batch_effect.ExpressionSet<-function(experiment, model, method=c('ComBat','RUV','MNN'), k){
  log<-experiment@assayData$exprs %>% names %>% switch(log_exprs=TRUE, exprs=FALSE)
  model.data<-model.frame(model, experiment@phenoData[all.vars(model)])
  if(length(k)==1|method=='ComBat'){
    return(ExpressionSet(
      assayData = list(switch(
        method,
        ComBat = ComBat(experiment@assayData$exprs[[1]], experiment$batch, mod=model.matrix(model, data=model.data)),
        RUV = RUVs(experiment@assayData$exprs[[1]], cIdx=seq_len(nrow(experiment@assayData$exprs[[1]])), k=k, 
                   scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=log)$normalizedCounts,
        MNN = mnnCorrect(experiment@assayData$exprs[[1]], batch=experiment$batch, k=k)@assayData$exprs$corrected
      )) %>% set_names(if(log) 'corrected_log_exprs' else 'corrected_exprs'),
      phenoData = experiment@phenoData,
      experimentData = experiment@experimentData
    ))
  }else{
    return(k %>% map(correct_batch_effect %>% partial(experiment=experiment, model=model, method=method)))
  }
}