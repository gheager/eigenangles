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

load_experiments<-function(directory, names=dir(directory), item.SimpleList='rnaseq'){
  directory %>% dir %>% 
    map(~get(load(paste0(directory,'/',.x)))) %>% 
    map(~if(class(.x)=='SimpleList') .x[[item.SimpleList]] else .x) %>% 
    set_names(names)
}

download_experiments_from_ExpressionAtlas<-function(..., destdir=getwd() %>% paste('experiments',sep='/')){
  if(!(destdir %>% dir.exists) ){
    destdir %>% dir.create
  }
  for(experiment in list(...)){
    paste0('https://wwwdev.ebi.ac.uk/gxa/experiments-content/',experiment,'/static/',experiment,'-atlasExperimentSummary.Rdata') %>% download.file(destfile = paste0(destdir,"/",experiment,".Rdata"))
  }
  destdir %>% load_experiments
}

# download_experiments_from_ExpressionAtlas<-function(..., destdir=getwd() %>% paste('experiments',sep='/')){
#   if(destdir %>% dir.exists){
#     stop('Attempted to create directory `',desdir,'`` but this address already exists. Rename it or modify `destdir` argument in `download_experiments`.')
#   }else{
#     destdir %>% dir.creates
#   }
#   for(experiment in list(...)){
#     paste0('https://www.ebi.ac.uk/gxa/experiments-content/',experiment,'/static/',experiment,'-atlasExperimentSummary.Rdata') %>% download.file(destdir)
#   }
#   desdir %>% load_experiments
# }

remove_isolated_experiments<-function(experiments, biological.group){
  warning('The following batches were removed as they do not have a ',biological.group,' column:\n',
          names(experiments)[experiments %>% map(~is.null(.x[[biological.group]])) %>% unlist])
  experiments[experiments %>% map(~is.null(.x[[biological.group]])) %>% unlist]<-NULL
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
  message('The following batches were removed as they do not share any common point on ',biological.group,' column:\n',
          rownames(intersections)[rowSums(intersections!=0)<=1])
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

access_data<-function(experiment) UseMethod("access_data")
access_data.SummarizedExperiment<-function(experiment){
  return(experiment@assays$data$counts)
}
access_data.ExpressionSet<-function(experiment){
  return(experiment@assayData$exprs)
}

access_pheno<-function(experiment) UseMethod("access_pheno")
access_pheno.SummarizedExperiment<-function(experiment){
  return(experiment@colData)
}
access_pheno.ExpressionSet<-function(experiment){
  return(experiment@phenoData@data)
}

access_meta<-function(experiment) UseMethod("access_meta")
access_meta.SummarizedExperiment<-function(experiment){
  return(experiment@metadata)
}
access_meta.ExpressionSet<-function(experiment){
  return(experiment@protocolData@data)
}

merge_experiments <- function(experiments, filter.unexpressed.genes=TRUE, log, force=FALSE){
  if(experiments %>% map(class) %>% unlist %>% unique %>% length %>% is_greater_than(1) & !force) stop("The experiments must have the same class. Their classes are :\n", experiments %>% map(class) %>% unlist)
  genes<-experiments %>% map(rownames)
  shared.genes<-genes %>% purrr::reduce(intersect)
  unshared.genes<-genes %>% map(setdiff %>% partial(y=shared.genes))
  data<-experiments %>% map(~.x %>% access_data %>% extract(shared.genes,)) %>% purrr::reduce(cbind)
  warning(unshared.genes %>% unlist %>% unique %>% length,' genes have been removed as they are not shared across the batches.')
  batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
  if(filter.unexpressed.genes){
    unexpressed.genes <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)==0) %>% purrr::reduce(`&`)
    data%<>%extract(!unexpressed.genes,)
    message(sum(unexpressed.genes),' genes have been removed as they were unexpressed across the samples of a batch.')
  }
  if(missing(log)){
    log<-switch(
      class(experiments[[1]]),
      SummarizedExperiment=, RangedSummarizedExperiment=TRUE,
      ExpressionSet=FALSE
    )
  }
  if(log) data %<>% log1p
  return(switch(
    class(experiments[[1]]),
    SummarizedExperiment=,RangedSummarizedExperiment=
      SummarizedExperiment(
        assays=if(log) list(log_counts=data) else list(counts=data),
        colData=experiments %>% map(colData) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
        metadata=list(batch=batch)
      ),
    ExpressionSet=
      ExpressionSet(
        assayData=data,
        phenoData=experiments %>% map(~.x@phenoData@data) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch) %>% AnnotatedDataFrame
      )
  ))
}

# merge_experiments.ExpressionSet <- function(experiments, log=FALSE, filter.unexpressed.genes=TRUE){
#   genes<-experiments %>% map(rownames)
#   shared.genes<-genes %>% purrr::reduce(intersect)
#   unshared.genes<-genes %>% map(setdiff %>% partial(y=shared.genes))
#   data<-experiments %>% map(~.x@assayData$exprs %>% extract(shared.genes,)) %>% purrr::reduce(cbind)
#   warning(unshared.genes %>% unlist %>% unique %>% length,' genes have been removed as they are not shared across the batches.')
#   batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
#   if(filter.unexpressed.genes){
#     unexpressed.genes <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)==0) %>% purrr::reduce(`&`)
#     data%<>%extract(!unexpressed.genes,)
#     message(sum(unexpressed.genes),' genes have been removed as they were unexpressed across the samples of a batch.')
#   }
#   if(log) data%<>%log1p
#   return(ExpressionSet(
#     assayData=if(log) list(log_exprs=data) else list(exprs=data),
#     phenoData=experiments %>% map(phenoData) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
#     experimentData=list(batch=batch)
#   ))
# }

correct_batch_effect<-function(experiment, model, method=c('ComBat','RUV','MNN'), k, batch="batch"){
  UseMethod("correct_batch_effect")
}

correct_batch_effect.SummarizedExperiment<-function(experiment, model, method, k=NULL, batch="batch"){
  log<-experiment@assays$data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
  data <- experiment@assays$data[[1]]
  model.data<-model.frame(model, experiment@colData[all.vars(model)])
  if(length(k)==1|method=='ComBat'){
    return(SummarizedExperiment(
      assays = list(switch(
        method,
        ComBat = ComBat(data, experiment[[batch]], mod=model.matrix(model, data=model.data)),
        RUV = RUVs(data, cIdx=seq_len(nrow(data)), k=k,
                   scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=TRUE)$normalizedCounts,
        MNN = mnnCorrect(data, batch=experiment[[batch]], k=k)@assays$data$corrected
      )) %>% set_names(if(log) 'corrected_log_counts' else 'corrected_counts'),
      colData = experiment@colData,
      metadata = experiment@metadata
    ))
  }else{
    return(k %>% map(~correct_batch_effect.SummarizedExperiment(experiment=experiment, model=model, method=method, k=.x, batch=batch)))
  }
}

correct_batch_effect.ExpressionSet<-function(experiment, model, method, k=NULL, batch="batch"){
  #log<-experiment@assayData$exprs %>% names %>% switch(log_exprs=TRUE, exprs=FALSE)
  data <- experiment@assayData$exprs
  model.data<-model.frame(model, experiment@phenoData@data[all.vars(model)])
  if(length(k)==1|method=='ComBat'){
    return(ExpressionSet(
      assayData = switch(
        method,
        ComBat = ComBat(data, experiment[[batch]], mod=model.matrix(model, data=model.data)),
        RUV = RUVs(data, cIdx=seq_len(nrow(data)), k=k,
                   scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=TRUE)$normalizedCounts,
        MNN = mnnCorrect(data, batch=experiment[[batch]], k=k)@assays$data$corrected
      ),
      phenoData = experiment@phenoData
    ))
  }else{
    return(k %>% map(~correct_batch_effect.ExpressionSet(experiment=experiment, model=model, method=method, k=.x, batch=batch)))
  }
}

# correct_batch_effect<-function(experiment, model, method=c('ComBat','RUV','MNN'), k=NULL){
#   #log<-experiment@assays$data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
#   model.data<-model.frame(model, experiment %>% access_pheno %>% extract(all.vars(model)))
#   for(kk in k){
#     
#   }
#   if(length(k)==1|method=='ComBat'){
#     return(ExpressionSet(
#       assayData = switch(
#         method,
#         ComBat = ComBat(experiment %>% access_data, experiment$batch, mod=model.matrix(model, data=model.data)),
#         RUV = RUVs(experiment %>% access_data, cIdx=seq_len(nrow(experiment %>% access_data)), k=k, 
#                    scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=log)$normalizedCounts,
#         MNN = mnnCorrect(experiment %>% access_data, batch=experiment$batch, k=k)@assays$data$corrected
#       ),# %>% set_names(if(log) 'corrected_log_counts' else 'corrected_counts'),
#       phenoData = experiment %>% access_pheno %>% AnnotatedDataFrame,
#       protocolData = experiment %>% access_meta %>% AnnotatedDataFrame
#     ))
#   }else{
#     return(k %>% map(~correct_batch_effect(experiment=experiment, model=model, method=method, k=.x)))
#   }
# }
# 
