library(shiny)
library(tidyverse)
library(ggiraph)

# https://rviews.rstudio.com/2019/10/09/building-interactive-world-maps-in-shiny/
world_data <- ggplot2::map_data('world')
happy_select <- read_csv("happy_select.csv")

# combine
combo_df <- happy_select %>% 
  dplyr::inner_join(
    world_data,
    join_by(country_name == region)
  )
