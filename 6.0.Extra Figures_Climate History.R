### FIGURE 2: Historical Conditions- Climate Normals and Fire History 
rm(list=ls())
load(file='/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/Annual_Climate_Summary_ENP.RDATA' )

# Climate Normals ..........................................................................................
summary.tot$Year <- summary.tot$Year %>% as.numeric()

plot.temp <- ggplot(data =summary.tot) + 
  geom_line(aes(x = Year, y = TMAX.mean), linetype = "dashed", size = 2, alpha = 0.8) + 
  geom_line(aes(x = Year, y = TMIN.mean), linetype = "dotted", size = 2, alpha = 0.8) +
  #geom_point( aes(x=Year, y = TMAX.mean), col="#e76155") + 
  #geom_line( aes(x=Year, y = TMAX.mean), col="#e76155", size=2, alpha=.5) + 
  #geom_point( aes(x=Year, y = TMIN.mean), col="#376795") + 
  #geom_line( aes(x=Year, y = TMIN.mean), col="#376795", size=2, alpha=.5) +
  labs(x="", y="Temperature (C)", tag="B") + 
  ylim(15, 40) +
  theme_bw() + 
  theme(text = element_text(size = 25), 
        legend.position = "top") +
  annotate("text", x = 2010, y = 32, label = "Maximum", size = 8, fontface ="plain" ) +
  annotate("text", x = 2010, y = 21, label = "Minimum", size = 8,  fontface = "plain") 

plot.prcp <- ggplot(data =summary.tot) + 
  geom_line( aes(x=Year, y = PRCP), size=3, alpha=.5) +
  geom_point( aes(x=Year, y = PRCP)) + 
  labs(x="", y="Precipitation (mm)", tag="C") +
  ylim(1200, 2050) +
  theme_bw() + 
  theme(text = element_text(size = 25))

plot.pdsi <- ggplot(data =summary.tot) + 
  geom_line( aes(x=Year, y = PDSI.mean), size=3, alpha=.5) +
  geom_point( aes(x=Year, y = PDSI.mean)) + 
  labs(x="", y="PDSI", tag="D") +
  ylim(-4, 4) +
  theme_bw() + 
  theme(text = element_text(size = 25))

# combine 
plot_climate <- ggarrange(plot.temp,
                          plot.prcp,
                          plot.pdsi,
                          nrow=3, ncol=1,
                          #labels=c('B','C', 'D'),
                          font.label = list(size = 30))
