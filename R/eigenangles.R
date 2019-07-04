angledet<-function(U,V){
  U%<>%as.matrix; V%<>%as.matrix
  cbind(U,V) %>% qr %>% qr.R -> new.coord
  cbind(
    new.coord[,1:ncol(U)] %>% qr %>% qr.Q,
    new.coord[,ncol(U)+(1:ncol(V))] %>% qr %>% qr.Q
  ) %>% det %>% abs %>% asin/pi
}

angles<-function(U,V){
  a<-NULL
  for(i in seq_len(min(ncol(U),ncol(V)))){
    a%<>%c(angledet(U[,1:i],V[,1:i]))
  }
  return(a)
}

eigenangles<-function(data,batch,group,scale=FALSE,center=TRUE,average=TRUE,verbose=TRUE){
  data%<>%t%<>%na.omit%<>%t
  batch%<>%factor; batch %>% levels -> batches
  ###########################
  angles_batch_vs_all<-list()
  for(b in batches){
    data.all<-NULL; data.batch<-NULL
    for(g in group[batch==b] %>% unique){
      data.all%<>%cbind(data[,group==g] %>% as.matrix %>% (if(average) rowMeans else identity))
      data.batch%<>%cbind(data[,group==g & batch==b] %>% as.matrix %>% (if(average) rowMeans else identity))
    }
    filter<-if(scale) rowVars(data.all)!=0 & rowVars(data.batch)!=0 else TRUE
    angles_batch_vs_all[[b]]<-angles(
      data.all[filter,] %>% t %>% prcomp(center=center,scale=scale) %>% use_series(rotation),
      data.batch[filter,] %>% t %>% prcomp(center=center,scale=scale) %>% use_series(rotation)
    )
    if(verbose) b %>% paste('\n') %>% cat
  }
  ##########################
  angles_inter_batch<-list()
  for(i in 1:(nlevels(batch)-1)){
    for(j in (i+1):nlevels(batch)){
      data.i<-NULL; data.j<-NULL
      for(g in intersect(group[batch==batches[i]],group[batch==batches[j]])){
        data.i%<>%cbind(data[,group==g & batch==batches[i]] %>% as.matrix %>% (if(average) rowMeans else identity))
        data.j%<>%cbind(data[,group==g & batch==batches[j]] %>% as.matrix %>% (if(average) rowMeans else identity))
      }
      if(length(intersect(group[batch==batches[i]],group[batch==batches[j]]))>ifelse(scale,1,0)){
        filter<-if(scale) rowVars(data.i)!=0 & rowVars(data.j)!=0 else TRUE
        angles_inter_batch[[batches[i]]][[batches[j]]]<-angles(
          data.i[filter,] %>% t %>% prcomp(center=center,scale=scale) %>% use_series(rotation),
          data.j[filter,] %>% t %>% prcomp(center=center,scale=scale) %>% use_series(rotation)
        )
        angles_inter_batch[[batches[j]]][[batches[i]]]<-angles_inter_batch[[batches[i]]][[batches[j]]]
      }
      if(verbose) batches[i] %>% paste('vs',batches[j],'\n') %>% cat
    }
  }
  for(b in batches) angles_inter_batch[[b]][[b]]<-0
  return(list(
    batch_vs_all=angles_batch_vs_all %>% structure(class='eigenangles.batch_vs_all'),
    inter_batch=angles_inter_batch %>% structure(class='eigenangles.inter_batch')
  ) %>% structure(class='eigenangles'))
}

plot.eigenangles.batch_vs_all<-function(...,angles=NULL){
  if(is.null(angles)) list(...) %>% transpose -> angles
  angles.df<-data.frame(
    cosinus=angles %>% unlist %>% cospi,
    sinus=angles %>% unlist %>% sinpi,
    method=angles %>% map(imap %>% partial(...=,~rep(..2,each=length(..1)))) %>% unlist %>% factor,
    batch=angles %>% imap(~rep(..2,each=length(unlist(..1)))) %>% unlist
  )
  ggplot(angles.df,aes(x=cosinus,y=sinus,colour=method))+
    xlim(c(0,1.25))+ylim(c(0,1))+
    geom_segment(xend=0,yend=0)+
    geom_text(label=angles %>% unlist %>% round(3) %>% paste('$\\pi$') %>% TeX,parse=TRUE,nudge_x=.1)+
    geom_arc(aes(x0=0,y0=0,r=1,start=0,end=pi/2),colour='black',inherit.aes = FALSE)+coord_fixed()+
    facet_wrap(~batch)
}

plot.eigenangles.inter_batch<-function(...,angles=NULL){
  if(is.null(angles)) list(...) -> angles
  angles.df<-data.frame(
    cosinus=angles %>% unlist %>% cospi,
    sinus=angles %>% unlist %>% sinpi,
    method=angles %>% imap(~rep(.y,each=length(unlist(.x)))) %>% unlist %>% factor,
    batch1=angles %>% map(imap %>% partial(...=,~rep(.y,each=length(unlist(.x))))) %>% unlist,
    batch2=angles %>% map(map %>% partial(...=,imap %>% partial(...=,~rep(.y,each=length(.x))))) %>% unlist
  )
  ggplot(angles.df,aes(x=cosinus,y=sinus,colour=method))+
    xlim(c(0,1.25))+ylim(c(0,1))+
    geom_segment(xend=0,yend=0)+
    geom_text(label=angles %>% unlist %>% round(3) %>% paste('$\\pi$') %>% TeX,parse=TRUE,nudge_x=.1)+
    geom_arc(aes(x0=0,y0=0,r=1,start=0,end=pi/2),colour='black',inherit.aes = FALSE)+coord_fixed()+
    facet_grid(batch2~batch1)
}

plot.eigenangles<-function(...){
  plot(plot.eigenangles.batch_vs_all(angles=list(...) %>% map(~.$batch_vs_all) %>% transpose))
  plot(plot.eigenangles.inter_batch(angles=list(...) %>% map(~.$inter_batch)))
}