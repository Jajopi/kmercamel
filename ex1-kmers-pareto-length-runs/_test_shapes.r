# install.packages(c("ggplot2","ggsci")) # run if needed

library(ggplot2)
library(ggsci)

# shapes available in ggplot
shapes <- 0:25

# NEJM palette colors
cols <- pal_nejm("default")(8)

# dataframe
df <- data.frame(
  shape = shapes,
  x = shapes %% 6,
  y = shapes %/% 6
)

# plot shapes
ggplot(df, aes(x, y)) +
  geom_point(
    aes(shape = factor(shape), color = factor(shape)),
    size = 6,
    stroke = 1.2
  ) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = rep(cols, length.out = length(shapes))) +
  geom_text(aes(label = shape), vjust = 2.2, size = 3) +
  labs(
    title = "ggplot2 Point Shapes (0–25) using NEJM palette",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

  ggsave("_test_shapes.png",
       height = 10,
       width = 10,
       unit = "cm")
