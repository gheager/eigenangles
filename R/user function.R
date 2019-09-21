# correct.batch.effect<-function(data,batch,
#                                method=c('none','combat','ruv','mnn'),
#                                model=NULL,model.data,log=TRUE,k=1){
#   for(u in model.data %>% seq_along){
#     eval(parse(text=paste0(
#       names(model.data)[u],'<-model.data[[u]]'
#     )))
#   }
#   batch%<>%factor
#   if(log) data%<>%log1p
#   if(method=='none') return(data)
#   else if(method=='bmc'){
#     if(is.null(model)) return(pamr.batchadjust(list(x=data,batchlabels=batch))$x)
#   }
#   else if(method=='combat'){
#     if(is.null(model)) return(ComBat(data,batch))
#     else return(ComBat(data,batch,mod=model.matrix(model,data=model.data)))
#   }
#   else if(method=='ruv'){
#     if(is.null(model)) stop('Biological model is needed to run RUVs')
#     else return(RUVs(data,cIdx=seq_len(nrow(data)),k=k,
#                      scIdx=makeGroups(expand.grid(model.data) %>% apply(1,paste)),isLog=log)$normalizedCounts)
#   }
#   else if(method=='mnn'){
#     return(mnnCorrect(data,batch=batch,k=k) %>% assays %>% use_series(corrected))
#   }
# }
# 
# integrate.experiments<-function(...,list=NULL,model=NULL,method=c('none','combat','ruv','mnn'),log=TRUE,k=1){
#   if(is.null(list)) experiments<-list(...) else experiments<-list; if(is.null(names(experiments))) names(experiments)<-paste0('batch',experiments %>% seq_along)
#   genes<-experiments %>% map(rownames)
#   common.genes<-genes %>% purrr::reduce(intersect)
#   filtered.genes<-genes %>% map(setdiff %>% partial(y=common.genes))
#   data<-NULL;batch<-NULL;vars<-list()
#   for(i in experiments %>% seq_along){
#     experiments[[i]]->exp
#     data%<>%cbind(exp %>% assays %>% use_series(counts) %>% extract(rownames(exp)%in%common.genes,))
#     batch%<>%c(names(experiments)[[i]] %>% rep(dim(exp)[2]))
#     for(v in all.vars(model)) vars[[v]]%<>%c(exp[[v]])
#   }; batch%<>%factor
#   filter <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
#   data%<>%extract(filter,)
#   data%<>%correct.batch.effect(batch,method,model,log,model.data=if(is.null(model)) NULL else model.frame(model,vars),k=k)
#   return(SummarizedExperiment(
#     assays=if(log) list(corrected_log_counts=data) else list(corrected_counts=data),
#     colData=gtools::smartbind(list=experiments %>% map(colData)) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
#     metadata=list(batch=batch)
#   ))
# }

benchmark_eigenangles<-function(directory, group){
  directory %>% dir %>% 
    map(~get(load(paste0(directory,'/',.x)))) %>% set_names(dir(directory)) %>% 
    do.call(apply_eigenangles %>% partial(group=group), .)
}

apply_eigenangles<-function(...,group){
  list(...) %>% map(do_eigenangles %>% partial(group=group)) %>% 
    imap(~mutate(.x,algorithm_=.y)) %>% 
    purrr::reduce(rbind.fill) %>% as_tibble %>% structure(class=c('eigenangles',class(.)))
}

do_eigenangles<-function(experiment,group){
  if(!is.list(experiment)){
    return(eigenangles(
      data=experiment %>% assays %>%extract2(1),
      batch=experiment %>% metadata %>% use_series(batch),
      group=experiment[[group]]
    ))
  }else{
    return(experiment %>% 
             map(do_eigenangles %>% partial(group=group)) %>% 
             imap(~mutate(.x,k_=.y)) %>% purrr::reduce(rbind))
  }
}

extract_dim<-function(tbl,dim){
  tbl$angles_ %<>% map(~.x[dim %>% min(length(.x))]) %<>% unlist
  return(tbl)
}

# do_eigenangles<-function(experiment,group,...){
#   eigenangles(
#     data=experiment %>% assays %>% use_series(corrected),
#     batch=experiment %>% metadata %>% use_series(batch),
#     group=experiment[[group]],
#     ...
#   )
# }
# 
do_gPCA<-function(experiment,...){
  experiment %>% assays %>% use_series(corrected) %>% gPCA(batch=experiment$batch,...)
}
