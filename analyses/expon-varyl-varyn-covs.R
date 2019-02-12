library(tidyverse)
devtools::load_all()

# this full dataset is too large to summarize on my laptop, 
# so this script does the necessary work and produces some data files for figures

load('../simdata/expfit-varyl-varyn-covs.Rdata')

# Create the N=1000 subset of data, which is used for main analyses.
# We do this in local so as not to change the original dataframe 
# (and annoyingly you can't pipe to save()).
local({
   expfit_varyl_res <- expfit_varyl_res %>% filter(N == 1000) 
   save(expfit_varyl_res, file='../simdata/expfit-varyl-covs.Rdata')
})

vark <- expfit_varyl_res %>% 
          unnest(stats, genic_va, allelic_cov) %>% 
          {if (all(.$gen == .$gen1)) return(.) else stop()}  %>%
          group_by(L, N, Va, alpha, genlen, r, gen) %>% 
          summarize(kvar=mean(kvar), 
                    zvar=mean(zvar), 
                    zbar=mean(zbar), 
                    genic_va=mean(genic_va, na.rm=TRUE),
                    ac_var=mean(ac_var),
                    ac_cov=mean(ac_cov)) %>% 
          ungroup()



covs <- expfit_varyl_res %>% 
         unnest(covs) %>% 
         # see note in temp_cov() for why these times are incremented
         mutate(t0=t0+1L, t1=t1+1L) %>%  
         group_by(L, N, Va, genlen, r, t0, t1) %>% 
         summarize(cov=mean(cov))

rec_params <- c(0.5, 0.1, 0.01, 1.5)
va_params <- c(0.01, 0.02, 0.05)


pd5 <- covs %>% filter((t1 == 5 & t1 != t0) | (t0 == 5 & t1 != t0)) %>% 
          ungroup() %>% mutate(L=factor(L), gen=ifelse(t1 > 5, t1, t0))

pred_all <- vark %>% ungroup() %>% 
          mutate(N=as.integer(as.character(N)), 
                        pred_zvar=pmap_dbl(list(genlen, gen-5, zvar, N), 
                                      analytic_cov),
                        pred_genic=pmap_dbl(list(genlen, gen-5, genic_va, N), 
                                      analytic_cov))

pd5f_varyn <- pd5 %>% filter(Va %in% va_params, genlen %in% rec_params)
predf_varyn <- pred_all %>% filter(Va %in% va_params, genlen %in% rec_params) %>% 
             gather(type, pred, pred_zvar, pred_genic) %>%
             mutate(L=factor(L), cov=pred)

# test with N
#ggplot() +  geom_point(data=pd5f %>% filter(L == 500), aes(t1, cov, color=as.factor(N)), 
#                       alpha=0.5, size=0.8) +
#  geom_hline(yintercept=0, color='red') + 
#  geom_line(data=predf %>%  filter(L==500, type=='pred_zvar'), aes(gen, pred, color=as.factor(N)))  + 
#  facet_grid(Va ~ genlen, scales='free_y')  

# tables and Rda files for plots:
write_tsv(pd5f_varyn, '../data/pd5f-varyn.tsv')
write_tsv(predf_varyn, '../data/predf-varyn.tsv')
devtools::use_data(pd5f_varyn, overwrite=TRUE)
devtools::use_data(predf_varyn, overwrite=TRUE)

