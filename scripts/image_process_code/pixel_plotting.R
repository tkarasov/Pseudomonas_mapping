#This script takes the pixel info 
library(dplyr)
library(tidyverse)
library(ggtree)
library(phytools)
library(heritability)
library(lme4)
'%!in%' <- function(x,y)!('%in%'(x,y))
#now combine pixel information with phylogeny
phy = read.tree("~/Dropbox/Pseudomonas_1524/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")

#drop irrelevant branches
pixel = read.table("~/work_remote/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/image_analysis/cell_pixel_strain.txt", header = T, sep = " ")

pixel_plot = pixel %>% group_by(actual_strain) %>% summarise(mean = mean(pixel_count, na.rm=T), se = sd(pixel_count, na.rm=T) / length(pixel_count))
tips=phy$tip.label
not_in = tips[tips %!in% pixel_plot$actual_strain]
red_phy = drop.tip(phy, not_in)

#rownames(pixel_plot) = pixel_plot$actual_strain
pixel_order = pixel_plot[red_phy$tip.label,]


pixel_keep = pixel[pixel$actual_strain %in% red_phy$tip.label,]
pixel_plot_keep = pixel_keep %>% group_by(actual_strain) %>% summarise(mean = mean(pixel_count, na.rm=T), se = sd(pixel_count, na.rm=T) / length(pixel_count))

#now plot the death stuff
p <- ggtree(red_phy)
p2 <- facet_plot(p, panel = 'dot', data = pixel_plot_keep, geom = geom_segment,aes(x=0, xend=mean, y=y, yend=y), size=.5, color='cyan4')


#Estimate heritability
model1 <- lmer(pixel_count ~ (1|actual_strain), data = pixel, REML = T)
vc =VarCorr(model1)
residual_var = attr(vc, 'sc')^2
intercept_var = attr(vc$actual_strain,'stddev')[1]^2
H = intercept_var/(residual_var + intercept_var)

#or not by hand
herit = repeatability(pixel$pixel_count, pixel$actual_strain)






