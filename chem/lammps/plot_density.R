read_csv('data.csv') %>%
ggplot() +
aes(Step,Density) +
geom_line() +
theme_tom() +
scale_x_continuous(labels=scales::comma) +
ggsave('plot.png', dpi=300)
