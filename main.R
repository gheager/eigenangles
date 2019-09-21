library(magrittr)
library(purrr)
library(eigenangles)

#Mouse
experiments <- load_experiments('Mouse Expression Atlas/batch data/')
experiments %<>% remove_isolated_experiments('organism_part')
experiments %<>% merge_experiments
experiments %>% save(file='Mouse Expression Atlas/integrated data/uncorrected.Rdata')

experiments %>% correct_batch_effect(~organism_part, 'ComBat') %>% 
  save(file='Mouse Expression Atlas/integrated data/ComBat_corrected.Rdata')
1:20 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'RUV', k=.x)) %>% 
  save(file='Mouse Expression Atlas/integrated data/RUV_corrected_k=1:20.Rdata')
1:20 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'MNN', k=.x)) %>% 
  save(file='Mouse Expression Atlas/integrated data/MNN_corrected_k=1:20.Rdata')
experiments %>% correct_batch_effect(~organism_part, 'BMC') %>% 
  save(file='Mouse Expression Atlas/integrated data/BMC_corrected.Rdata')

apply_eigenangles(
  uncorrected = get(load('Mouse Expression Atlas/integrated data/uncorrected.Rdata')),
  ComBat = get(load('Mouse Expression Atlas/integrated data/ComBat_corrected.Rdata')),
  RUVs = get(load('Mouse Expression Atlas/integrated data/RUVs_corrected_k=1:20.Rdata')),
  MNN = get(load('Mouse Expression Atlas/integrated data/MNN_corrected_k=1:20.Rdata')),
  BMC = get(load('Mouse Expression Atlas/integrated data/BMC_corrected.Rdata')),
  group = 'organism_part'
) -> ea

ea %>% plot
ea %>% tanmean %>% plot

#Human
experiments <- load_experiments('Human Expression Atlas/batch data/')
experiments %<>% remove_isolated_experiments('organism_part')
experiments %<>% merge_experiments
experiments %>% save(file='Human Expression Atlas/integrated data/uncorrected.Rdata')

experiments<-get(load('Human Expression Atlas/integrated data/uncorrected.Rdata'))
experiments %>% correct_batch_effect(~organism_part, 'ComBat') %>% 
  save(file='Human Expression Atlas/integrated data/ComBat_corrected.Rdata')
1:5 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'RUV', k=.x)) %>% 
  save(file='Human Expression Atlas/integrated data/RUVs_corrected_k=1:5.Rdata')
1:5 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'MNN', k=.x)) %>% 
  save(file='Human Expression Atlas/integrated data/MNN_corrected_k=1:5.Rdata')
experiments %>% correct_batch_effect(~organism_part, 'BMC') %>% 
  save(file='Human Expression Atlas/integrated data/BMC_corrected.Rdata')

apply_eigenangles(
  uncorrected = get(load('Human Expression Atlas/integrated data/uncorrected.Rdata')),
  ComBat = get(load('Human Expression Atlas/integrated data/ComBat_corrected.Rdata')),
  RUVs = get(load('Human Expression Atlas/integrated data/RUVs_corrected_k=1:5.Rdata')),
  MNN = get(load('Human Expression Atlas/integrated data/MNN_corrected_k=1:5.Rdata')),
  BMC = get(load('Human Expression Atlas/integrated data/BMC_corrected.Rdata')),
  group = 'organism_part'
) -> ea
ea %>% save(file='Human Expression Atlas/integrated data/eigenangles.Rdata')
ea %>% plot
ea %>% tanmean %>% plot
