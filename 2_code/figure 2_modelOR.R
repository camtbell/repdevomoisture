##### -------- Analysis for Reptile Moisture Meta-analysis  -------- ######
#####                   Trait: Hatching Survival                     ######
#####     Analysis lead by Cameron Bell                              ######
#####     Reviewed with Daniel Noble edits; 24 Jan 2023 ****         ######
##### -------------------------------------------------------------- ######

#clear
rm(list=ls())

#load packages
pacman::p_load(metaAidR, metafor, corrplot, MCMCglmm, ape, phytools, stats, rotl, readr, metafor, MuMIn, clubSandwich, dplyr, cowplot, ggplot2, ggpubr, dplyr, hrbrthemes, viridis, RColorBrewer, patchwork)

#grab the files
survival <- read_csv("3_trait data/survival_data.csv")
sex <- read.csv("3_trait data/sex_data.csv")

#escalc 
SurvivalOR <- escalc(measure="OR", ai=alive.2, bi=dead.2, ci =alive.1, di =dead.1 , data=survival)
SexOR <- escalc(measure="OR", ai=male.2, bi=female.2, ci =male, di =female , data=sex)

###Survival Plots-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Survival plots (S1 and S2)
S1<-SurvivalOR %>%
  ggplot(aes(x=T, y=yi, size=N, color=order)) +
  labs(x="Temperature (°C)", y=" Effect Size (lnOR)", title= "Hatching Success") +
  labs(shape="Sample Size", color="Order") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  geom_hline(aes(fill=yi),yintercept =0, linetype=2)+
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 10), name="Sample size") +
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=10)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.position = "none") +
  geom_smooth(method=lm, lwd=0.5, colour="black")

S2<-SurvivalOR %>%
  filter(waterpot_diff < 4000) %>%
  ggplot(aes(x=waterpot_diff, y=yi, size=N, color=order)) +
  labs(x="Moisture Difference (kPa)", y=" Effect size (lnOR)", title= "") +
  labs(shape="Sample Size", color="Order") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  geom_hline(aes(fill=yi),yintercept =0, linetype=2)+
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 10), name="Sample size") +
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=10)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 22)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.position = "none") +
  geom_smooth(method=lm, lwd=0.5, colour="black")

###Sex Plots-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Sex Plots (Sex1 and Sex2)
Sex1 <-SexOR %>%
  ggplot(aes(x=T, y=yi, size=N, color=order)) +
  labs(x="Temperature (°C)", y=" Effect size (lnOR)", title= "Sex Ratio") +
  labs(shape="Sample Size", color="Order") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  geom_hline(aes(fill=yi),yintercept =0, linetype=2)+
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 10), name="Sample size") +
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=10)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.position = "none") 

Sex2 <-SexOR %>%
  filter(waterpot_diff < 4000) %>%
  ggplot(aes(x=waterpot_diff, y=yi, size=N, color=order)) +
  labs(x="Moisture Difference (kPa)", y=" Effect size (lnOR)", title= "") +
  labs(shape="Sample Size", color="Order") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  geom_hline(aes(fill=yi),yintercept =0, linetype=2)+
  geom_point(alpha=0.5) +
  scale_size(range = c(.1, 10), name="Sample size") +
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=10)) +
  theme(axis.title = element_text(size = 12)) +
  theme(plot.title = element_text(size = 22)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.position = "none")

###Patchwork Codes-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
patchwork <- S1 + Sex1 + S2 + Sex2 + theme(legend.position = "right") + plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
patchwork

