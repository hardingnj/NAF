library(tidyverse)
library(rhdf5)

# Tidies coach data from .h5 file for plotting
# Rows are the races and columns the dates
tidy_data <- function(data, r_name, c_name) {
  data <- data %>% 
    as_data_frame() %>% 
    set_colnames(c_name) %>% 
    mutate(race = r_name) %>% 
    gather(date, ranking, -race) %>% 
    filter(!is.nan(ranking)) %>% 
    mutate(date = lubridate::ymd(date))
}

# Given a rankings.h5 file and a coach name, extracts the coach data
extract_coach_data <- function(coach, file) {
  accessor = paste0("coaches/", coach)
  
  data <- h5read(file, accessor) %>% 
    map(tidy_data, h5read(file, "race_ids"), h5read(file, "date")) %>% 
    cbind() %>% 
    bind_cols() %>% 
    select(race, date, mu = "ranking", phi = "ranking1") %>% 
    group_by(race) %>% 
    arrange(date) %>% 
    mutate(last_date = lag(date))
  
  H5close()
  
  data
}

# Produces the underlying data for the density plots.
# For each date/race combination, calculates the normal density of the specified number of bins from mu +/- 3 phi 
density_data <- function(d, bins = 100) {
  dens_data <- data_frame(
    date      = rep(d$date, each = bins), 
    last_date = rep(d$last_date, each = bins), 
    race      = rep(d$race, each = bins), 
    mu        = rep(d$mu, each=bins), 
    phi       = rep(d$phi, each = bins), 
    y         = map2(d$mu, d$phi, ~seq(.x-.y*3, .x+.y*3, length.out = bins)) %>% unlist # $bins points evenly spaced mu +/- 3 phi  
  )
  
  dens_data %>% 
    group_by(date,race) %>% 
    mutate(dens = dnorm(y, mean=mu, sd=phi)) # normal density at each bin
}

top10_coaches <-c("Jimjimany", "AndyDavo", "Pipey", "straume", "Phoenix11", "Jeff", "delevus", "Karaak", "delevus", "Harti")

# Given a coach name, make a plot of historical ratings
plot_coach <- function(coach) {
  extract_coach_data(coach, "../output/rankings.h5") %>% 
    ggplot(aes(x = date, fill = race)) +
    geom_rect( # Hidden skill probability distribution
      data = function(d) density_data(d, bins = 40), 
      aes(xmin = last_date, xmax = date, ymin = y, ymax = lead(y), fill = race, alpha = dens, group=race), 
      colour=NA
    ) +
    geom_point(aes(y = mu, colour = race), alpha = 0.2, size = 0.5) + # Points to show mu
    geom_point(aes(y = mu-(2.5*phi)), colour = "grey40", size = 0.5) + # Points to show rating
    geom_smooth(aes(y = mu-(2.5*phi)), colour = "black", se = F, span = 0.2) + # Smoothed ranking
    labs(title = "Glicko Ratings", subtitle = coach, x = NULL, y = "Rating") +
    facet_wrap(~race) +
    viridis::scale_color_viridis(discrete = T) + 
    viridis::scale_fill_viridis(discrete = T) +
    theme_linedraw() +
    theme(legend.position= "none")
}

plots <- map(top10_coaches, plot_coach)
