# get data from bm_PlotResponseCurves
response_curves<- bm_PlotResponseCurves(bm.out = myBiomodModelOut5,
                                        models.chosen = get_built_models(myBiomodModelOut5)[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39)],
                                        fixed.var = 'median')


response_curves_data <- response_curves$plot$data

library("viridis")
library("foreach")
palette <- viridis_pal(option = "viridis")(20)

# plot
response_bat = response_curves_data %>%
  filter(expl.name == "bat") %>% 
  ggplot( aes(x = expl.val, y = pred.val, 
              color = pred.name)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette) +
  scale_y_continuous(position="left", n.breaks = 10, expand = c(0, 0), limits = c(0,1)) +
  
  labs(y = "Predição", x = "Batimetria (m)") +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "grey",
                                linewidth = 0.8, linetype = "solid"),
    axis.line.y = element_line(colour = "grey",
                               linewidth = 0.8, linetype = "solid"),
    axis.ticks.x = element_line(colour = "grey",
                                linewidth = 0.8, linetype = "solid"),
    axis.line.x = element_line(colour = "grey",
                               linewidth = 0.8, linetype = "solid"),
    axis.text.x = element_text(size = 13,  color = "#284b80" ),
    axis.text.y = element_text(size = 15,  color = "grey" ),
    axis.title.y = element_text(size = 14,  color = "#284b80" ),
    axis.title.x = element_text(size = 14,  color = "#284b80" ),
    legend.position="none"
  )

response_bat
ggsave("pacs_figs/reponse_bat.png", width = 10, height = 5, dpi = 300)

response_dist_inv = response_curves_data %>%
  filter(expl.name == "dist_inv") %>% 
  ggplot( aes(x = expl.val/1000, y = pred.val, 
              color = pred.name)) +
  geom_line(size = 1) +
  scale_color_manual(values = palette) +
  scale_y_continuous(position="left", n.breaks = 10, expand = c(0, 0), limits = c(0,1)) +
  scale_x_continuous(position="bottom", n.breaks = 20, expand = c(0, 0)) +
  
  labs(y = "Predição", x = "Distancia dos focos RN e Engenho (Km)") +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "grey",
                                linewidth = 0.8, linetype = "solid"),
    axis.line.y = element_line(colour = "grey",
                               linewidth = 0.8, linetype = "solid"),
    axis.ticks.x = element_line(colour = "grey",
                                linewidth = 0.8, linetype = "solid"),
    axis.line.x = element_line(colour = "grey",
                               linewidth = 0.8, linetype = "solid"),
    axis.text.x = element_text(size = 13,  color = "#284b80" ),
    axis.text.y = element_text(size = 15,  color = "grey" ),
    axis.title.y = element_text(size = 14,  color = "#284b80" ),
    axis.title.x = element_text(size = 14,  color = "#284b80" ),
    legend.position="none"
  )  

response_dist_inv
ggsave("pacs_figs/response_dist_inv.png", width = 10, height = 5, dpi = 300)
