# Figure one - combined

library(ggpubr)
library(patchwork)
library(here)
library(yaml)

cfg <- yaml::read_yaml(here("config.yaml"))


wel_contr <- readRDS(file.path(cfg$paths$output_figures, "sad_uni_model", "uni_surv_well controlled asthma_B_R520ABNRp.rds"))

poor_contr <- readRDS(file.path(cfg$paths$output_figures, "sad_uni_model", "uni_surv_poorly controlled asthma_B_R520ABNRp.rds"))

# 1. Prepare Plot A components
col_a <- wrap_elements(panel = (wel_contr$plot / wel_contr$table + plot_layout(heights = c(3, 1)))) + 
  labs(title = "Figure 1A - Well-controlled asthma (ACQ6 < .75)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )
  
# 2. Prepare Plot B components
col_b <- wrap_elements(panel = (poor_contr$plot / poor_contr$table + plot_layout(heights = c(3, 1)))) + 
  labs(title = "Figure 1B - Poorly-controlled asthma (ACQ6 > 1.5)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )
# 3. Combine and collect the legend
combined_fig <- (col_a | col_b) + 
  plot_layout(guides = "collect") &   # This creates the shared legend
  theme(legend.position = "bottom")

# 5. Display and save
png(file.path(cfg$paths$output_figures, "sad_uni_model", "uni_surv_figure1.png"), 
    width = 1200, 
    height = 600)
print(combined_fig)
dev.off()

