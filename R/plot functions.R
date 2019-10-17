plot.eigenangles<-function(tbl,component=1, colour='blue'){
  ggplot(tbl %>% extract_component(component))+
    geom_hline(aes(yintercept = integration_angles), colour=colour)+
    geom_hline(aes(yintercept = -conservation_angles), colour=colour)+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y',start=pi,direction=-1)+ylim(c(-1,1))+xlim(c(0,1))+
    annotate(label='integration',y=2/3,x=1/2,geom='text',size=2)+
    annotate(label='conservation',y=-2/3,x=1/2,geom='text',size=2)+
    facet_wrap(~batch_)
}

plot.eigenangles.parametric<-function(tbl,component=1, colour='blue'){
  ggplot(tbl %>% extract_component(component))+
    geom_point(aes(x=k, y=integration_angles), colour=colour)+
    geom_point(aes(x=k, y=-conservation_angles), colour=colour)+
    geom_hline(aes(yintercept = integration_angles*is.na(k)), colour=colour)+
    geom_hline(aes(yintercept = -conservation_angles*is.na(k)), colour=colour)+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y',start=pi,direction=-1)+ylim(c(-1,1))+xlim(c(0,max(tbl$k,na.rm=TRUE)))+
    annotate(label='integration',y=2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text',size=2)+
    annotate(label='conservation',y=-2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text',size=2)+
    facet_wrap(~batch_)
}

plot.eigenangles.benchmark<-function(tbl,component=1){
  ggplot(tbl %>% extract_component(component))+
    geom_point(aes(x=k, y=integration_angles, colour=algorithm))+
    geom_point(aes(x=k, y=-conservation_angles, colour=algorithm))+
    geom_hline(aes(yintercept = integration_angles*is.na(k), colour=algorithm))+
    geom_hline(aes(yintercept = -conservation_angles*is.na(k), colour=algorithm))+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y',start=pi,direction=-1)+ylim(c(-1,1))+xlim(c(0,max(tbl$k,na.rm=TRUE)))+
    annotate(label='integration',y=2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text',size=2)+
    annotate(label='conservation',y=-2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text',size=2)+
    facet_wrap(~batch_)
}
#
# plot.eigenangles<-function(tbl,ref='none',dim=1){
#   ggplot(tbl %>% extract_dim(dim))+
#     aes(x=k_, y=angles_, colour=algorithm_)+
#     geom_point()+
#     geom_hline(aes(yintercept = angles_*is.na(k_), colour=algorithm_))+
#     geom_hline(yintercept=0, colour='black')+
#     coord_polar(theta='y')+ylim(c(0,2))+
#     facet_wrap(~batch_)
# }

tanmean<-function(tbl) UseMethod("tanmean")

tanmean.eigenangles<-function(tbl, component=1){
  tbl %>% extract_component(component) %>% 
    group_by(k) %>% 
    dplyr::summarise(
      integration_angles=atan(mean(tanpi(integration_angles)))/pi,
      conservation_angles=atan(mean(tanpi(conservation_angles)))/pi
    ) %>% structure(class=c('eigenangles.tanmean',class(.)))
}

tanmean.eigenangles.benchmark<-function(tbl, component=1){
  tbl %>% extract_component(component) %>% 
    group_by(algorithm,k) %>% 
    dplyr::summarise(
      integration_angles=atan(mean(tanpi(integration_angles)))/pi,
      conservation_angles=atan(mean(tanpi(conservation_angles)))/pi
    ) %>% structure(class=c('eigenangles.benchmark.tanmean',class(.)))
}

plot.eigenangles.tanmean<-function(tbl, colour='blue'){
  ggplot(tbl)+
    geom_point(aes(x=k, y=integration_angles), colour=colour)+
    geom_point(aes(x=k, y=-conservation_angles), colour=colour)+
    geom_hline(aes(yintercept = integration_angles*is.na(k)), colour=colour)+
    geom_hline(aes(yintercept = -conservation_angles*is.na(k)), colour=colour)+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y',start=pi,direction=-1)+ylim(c(-1,1))+xlim(c(0,max(tbl$k,na.rm=TRUE)))+
    annotate(label='integration',y=2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text')+
    annotate(label='conservation',y=-2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text')
}

plot.eigenangles.benchmark.tanmean<-function(tbl){
  ggplot(tbl)+
    geom_point(aes(x=k, y=integration_angles, colour=algorithm))+
    geom_point(aes(x=k, y=-conservation_angles, colour=algorithm))+
    geom_hline(aes(yintercept = integration_angles*is.na(k), colour=algorithm))+
    geom_hline(aes(yintercept = -conservation_angles*is.na(k), colour=algorithm))+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y',start=pi,direction=-1)+ylim(c(-1,1))+xlim(c(0,max(tbl$k,na.rm=TRUE)))+
    annotate(label='integration',y=2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text')+
    annotate(label='conservation',y=-2/3,x=max(tbl$k,na.rm=TRUE)/2,geom='text')
}

extract_component<-function(tbl,component){
  tbl$integration_angles %<>% map(~.x[component %>% min(length(.x))]) %<>% unlist
  tbl$conservation_angles %<>% map(~.x[component %>% min(length(.x))]) %<>% unlist
  return(tbl)
}

#
# plot.eigenangles.tanmean<-function(tbl,ref='none',dim=1){
#   ggplot(tbl)+
#     aes(x=k_, y=angles_, colour=algorithm_)+
#     geom_point()+
#     geom_hline(aes(yintercept = angles_*is.na(k_), colour=algorithm_))+
#     geom_hline(yintercept=0, colour='black')+
#     coord_polar(theta='y')+ylim(c(0,2))+xlim(c(min(tbl$k_,na.rm=TRUE)-1,max(tbl$k_,na.rm=TRUE)))
# }

# viz_meantan<-function(tbl,ref='none',dim=1){
#   ggplot(tbl %>% extract_dim(dim) %>% 
#            group_by(algorithm_,k_) %>% 
#            dplyr::summarise(angles_=atan(mean(tanpi(angles_)))/pi))+
#     aes(x=k_, y=angles_, colour=algorithm_)+
#     geom_point()+
#     geom_hline(aes(yintercept = angles_*is.na(k_), colour=algorithm_))+
#     geom_hline(yintercept=0, colour='black')+
#     coord_polar(theta='y')+ylim(c(0,2))+xlim(c(min(tbl$k_,na.rm=TRUE)-1,max(tbl$k_,na.rm=TRUE)))
# }
# 
# viz_meansincos<-function(tbl,ref='none',dim=1){
#   ggplot(tbl %>% extract_dim(dim) %>% 
#            group_by(algorithm_,k_) %>% 
#            dplyr::summarise(angles_=atan(mean(sinpi(angles_))/mean(cospi(angles_)))/pi))+
#     aes(x=k_, y=angles_, colour=algorithm_)+
#     geom_point()+
#     geom_hline(aes(yintercept = angles_*is.na(k_), colour=algorithm_))+
#     geom_hline(yintercept=0, colour='black')+
#     coord_polar(theta='y')+ylim(c(0,2))+xlim(c(min(tbl$k_,na.rm=TRUE)-1,max(tbl$k_,na.rm=TRUE)))
# }
#
# viz_angles_corrections<-function(ref, ...){
#   list(...) -> corrections
#   corrections %>% map(~class(.x[[1]])) %>% equals("list") -> parametric
#   ggplot(data.frame(
#     k = corrections[parametric] %>% map(seq_along) %>% unlist,#à changer
#     correction = corrections[parametric] %>% imap(~.y %>% rep(length(.x))) %>% unlist,#à changer
#     angle = corrections[parametric] %>% map(~.x %>% map(~.x %>% map(~.x[[1]]))) %>% unlist,
#     batch = corrections[parametric] %>% map(~.x %>% map(names)) %>% unlist
#   ))
# }