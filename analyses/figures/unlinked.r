library(devtools)
library(cowplot)
library(tidyverse)
library(viridis)
library(scales)
load_all()

inverse_haldane <- function(r) log(1-2*r) * -0.5

linked <- function(t=2, N=1e3) integrate(function(r) ok71(4*N*r)*(1-r)^t, 0, 0.5)$value
unlinked <- function(M, t=2, N=1e3) ok71(4*N*0.5) * 0.5^t * (M)


d1e3 <- crossing(M=seq(0, 40, length.out=500), t=seq(0, 7, length.out=500)) %>% 
  mutate(f=map2_dbl(M, t, function(M, t) unlinked(M, t) - linked(t))) %>% 
  mutate(g=map2_dbl(M, t, function(M, t) unlinked(M, t) / linked(t))) %>% arrange(M, t) 
 
# OLD NAME: we use 1e5 here now
d1e4 <- crossing(M=seq(0, 40, length.out=500), t=seq(0, 7, length.out=500)) %>% 
  mutate(f=map2_dbl(M, t, function(M, t) unlinked(M, t, N=1e5) - linked(t, N=1e5))) %>% 
  mutate(g=map2_dbl(M, t, function(M, t) unlinked(M, t, N=1e5) / linked(t, N=1e5))) %>% 
  arrange(M, t) 


#xy <- d %>% filter(abs(f) < 0.00001)
#d %>% ggplot(aes(M, t, fill=f))  + geom_tile()  +  
#  #scale_fill_gradient2(high='dark blue', mid='white', low='dark red', space='Lab') + 
#  scale_fill_gradientn(name='unlinked - linked', colours = c("blue","white","red"), 
#                         values = rescale(c(-.003, 0, .02)),
#                         guide = "colorbar", limits=c(-.003, .02)) +
#  geom_line(data=xy, aes(M, t), linetype='dashed')  + 
#  geom_point(data=xy[1, ], aes(M, t)) + theme_classic()  +
#  geom_hline(data=tibble(x=0:7), aes(yintercept=x), linetype='dashed', size=.1, alpha=0.4) + 
#  annotate('text', x=30, y=1.5, label='unlinked > linked') + 
#  annotate('text', x=20, y=4.5, label='linked > unlinked') + 
#  xlab('size of genome, excluding focal chromosome (Morgans)') + 
#  ylab('time (generations)') 

xy1e3 <- d1e3 %>% filter(abs(g-1) < 0.001)
xy1e4 <- d1e4 %>% filter(abs(g-1) < 0.001)

p1e3 <- d1e3 %>% ggplot(aes(M, t, fill=log10(g)))  + geom_tile()  +  
  #scale_fill_gradient2(high='dark blue', mid='white', low='dark red', space='Lab') + 
  scale_fill_gradientn(name='unlinked / linked', colours = c("blue","white","red"), 
                         values = rescale(c(-4, 0, 1.1)),
                         guide = "colorbar",
                         limits=c(-4, 1.1))  +
  geom_line(data=xy1e3, aes(M, t), linetype='dashed')  + 
  geom_point(data=xy1e3[1, ], aes(M, t)) + theme_classic()  +
  geom_hline(data=tibble(x=0:7), aes(yintercept=x), linetype='dashed', size=.1, alpha=0.4) + 
  annotate('text', x=30, y=1.5, label='unlinked > linked') + 
  annotate('text', x=20, y=4.5, label='linked > unlinked') + 
  xlab('') + 
  ylab('generations elapsed (|s-t|)') + theme(legend.position = "none")

p1e4 <- d1e4 %>% ggplot(aes(M, t, fill=log10(g)))  + geom_tile()  +  
  #scale_fill_gradient2(high='dark blue', mid='white', low='dark red', space='Lab') + 
  scale_fill_gradientn(name='unlinked / linked', colours = c("blue","white","red"), 
                         values = rescale(c(-4, 0, 1.1)),
                         guide = "colorbar",
                         limits=c(-4, 1.1))  +
  geom_line(data=xy1e4, aes(M, t), linetype='dashed')  + 
  geom_point(data=xy1e4[1, ], aes(M, t)) + theme_classic()  +
  geom_hline(data=tibble(x=0:7), aes(yintercept=x), linetype='dashed', size=.1, alpha=0.4) + 
  annotate('text', x=30, y=1.5, label='unlinked > linked') + 
  annotate('text', x=20, y=4.5, label='linked > unlinked') + 
  xlab('') + 
  ylab('')+ 
  theme(legend.position = "none")


p <- plot_grid(p1e3, p1e4, #labels = c('N = 1,000', 'N = 100,000'), 
          label_size = 12, label_x=10)


ggsave('unlinked_linked_contributions.png', p, width=6, height=.4*6)
