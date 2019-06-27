viz_eigenangles_batch_vs_all<-function(...){
  list(...) %>% map(~.$batch_vs_all) %>% transpose -> angles
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

viz_eigenangles_inter_batch<-function(...){
  list(...) %>% map(~.$inter_batch) -> angles
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
