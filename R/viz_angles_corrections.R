library(magrittr)
library(purrr)
library(ggplot2)

plot.eigenangles<-function(tbl,ref='none',dim=1){
  ggplot(tbl %>% extract_dim(dim))+
    aes(x=k_, y=angles_, colour=algorithm_)+
    geom_point()+
    geom_hline(aes(yintercept = angles_*is.na(k_), colour=algorithm_))+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y')+ylim(c(0,2))+
    facet_wrap(~batch_)
}

tanmean<-function(tbl, dim=1){
  tbl %>% extract_dim(dim) %>% 
    group_by(algorithm_,k_) %>% 
    dplyr::summarise(angles_=atan(mean(tanpi(angles_)))/pi) %>% 
    structure(class=c('eigenangles.tanmean',class(.)))
}

plot.eigenangles.tanmean<-function(tbl,ref='none',dim=1){
  ggplot(tbl)+
    aes(x=k_, y=angles_, colour=algorithm_)+
    geom_point()+
    geom_hline(aes(yintercept = angles_*is.na(k_), colour=algorithm_))+
    geom_hline(yintercept=0, colour='black')+
    coord_polar(theta='y')+ylim(c(0,2))+xlim(c(min(tbl$k_,na.rm=TRUE)-1,max(tbl$k_,na.rm=TRUE)))
}

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