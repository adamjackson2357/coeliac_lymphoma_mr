# PCs

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(ggplot2)
})

config <- read_yaml('../../configs/main.yml')
case_covars_output <- config$case_covars_output

case_covars <- readRDS(case_covars_output)
case_covars$outcome <- ifelse(case_covars$outcome == 1, "NHL", "Control")
table(case_covars$outcome)

pc1_pc2 <- case_covars %>% 
  ggplot() +
  geom_point(data=subset(case_covars, outcome == "Control"), aes(x=PC.1, y=PC.2, color=outcome), size = 0.5) +
  geom_point(data=subset(case_covars, outcome == "NHL"), aes(x=PC.1, y=PC.2, color=outcome), size = 0.5) +
  scale_color_manual(values=c("Control" = "grey", "NHL" = "dodgerblue4")) +
  theme_minimal()

pc1_pc3 <- case_covars %>% 
  ggplot() +
  geom_point(data=subset(case_covars, outcome == "Control"), aes(x=PC.1, y=PC.3, color=outcome), size = 0.5) +
  geom_point(data=subset(case_covars, outcome == "NHL"), aes(x=PC.1, y=PC.3, color=outcome), size = 0.5) +
  scale_color_manual(values=c("Control" = "grey", "NHL" = "dodgerblue4")) +
  theme_minimal()

pc2_pc3 <- case_covars %>% 
  ggplot() +
  geom_point(data=subset(case_covars, outcome == "Control"), aes(x=PC.2, y=PC.3, color=outcome), size = 0.5) +
  geom_point(data=subset(case_covars, outcome == "NHL"), aes(x=PC.2, y=PC.3, color=outcome), size = 0.5) +
  scale_color_manual(values=c("Control" = "grey", "NHL" = "dodgerblue4")) +
  theme_minimal()

pc1_pc4 <- case_covars %>% 
  ggplot() +
  geom_point(data=subset(case_covars, outcome == "Control"), aes(x=PC.1, y=PC.4, color=outcome), size = 0.5) +
  geom_point(data=subset(case_covars, outcome == "NHL"), aes(x=PC.1, y=PC.4, color=outcome), size = 0.5) +
  scale_color_manual(values=c("Control" = "grey", "NHL" = "dodgerblue4")) +
  theme_minimal()
print(pc1_pc4)

pcs <- ggarrange(
  pc1_pc2, pc1_pc3, pc2_pc3, ncol=3, nrow=1,
  common.legend = TRUE, legend = "bottom"
)
ggsave(plot = pcs, filename="../../figures/pcs.png")
print(pcs)