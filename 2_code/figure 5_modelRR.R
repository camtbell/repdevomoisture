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
mass <- read_csv('3_trait data/mass_data.csv')
length <- read.csv('3_trait data/length_data.csv')
ID <- read.csv('3_trait data/incubation_data.csv')

#escalc 
MassRR <- escalc(measure = "ROM", n1i = N.2, n2i =N, m1i =mean.2, m2i =mean, sd1i =sd.2, sd2i =sd, data=mass)
LengthRR<-escalc(measure = "ROM", n1i = N.2, n2i =N, m1i =mean.2, m2i =mean, sd1i =sd.2, sd2i =sd, data = length)
IDRR<-escalc(measure = "ROM", n1i = N.2, n2i =N, m1i =mean.2, m2i =mean, sd1i =sd.2, sd2i =sd, data = ID)
###Mass Plots-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Mass plots (Mass1 and Mass2)
Mass1<-MassRR %>%
  ggplot(aes(x=T, y=yi, size=N, color=order)) +
  labs(x="Temperature (°C)", y=" Effect Size (lnRR)", title= "Mass") +
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

Mass2<-MassRR %>%
  filter(waterpot_diff < 4000) %>%
  ggplot(aes(x=waterpot_diff, y=yi, size=N, color=order)) +
  labs(x="Moisture Difference (kPa)", y=" Effect size (lnRR)", title= "") +
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


###Length Plots-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Length Plots (Length1 and Length2)
Length1 <-LengthRR %>%
  ggplot(aes(x=T, y=yi, size=N, color=order)) +
  labs(x="Temperature (°C)", y=" Effect size (lnRR)", title= "Length") +
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

Length2 <-LengthRR %>%
  filter(waterpot_diff < 4000) %>%
  ggplot(aes(x=waterpot_diff, y=yi, size=N, color=order)) +
  labs(x="Moisture Difference (kPa)", y=" Effect size (lnRR)", title= "") +
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

###Incubation Duration Plots-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Incubation Duration Plots (ID1 and ID2)
ID1 <-IDRR %>%
  ggplot(aes(x=T, y=yi, size=N, color=order)) +
  labs(x="Temperature (°C)", y=" Effect size (lnRR)", title= "Incubation Duration") +
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

ID2 <-IDRR %>%
  filter(waterpot_diff < 4000) %>%
  ggplot(aes(x=waterpot_diff, y=yi, size=N, color=order)) +
  labs(x="Moisture Difference (kPa)", y=" Effect size (lnRR)", title= "") +
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

###Patchwork Plot Code--------------------------------------------------------------------------------------------------------------------------------------------------------------------

patchwork <- Mass1 + Length1 + ID1 + Mass2 + Length2 + ID2 + theme(legend.position = "right") + plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
patchwork

