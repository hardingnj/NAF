library(tidyverse)
library(rhdf5)

d = H5Fopen("~/Downloads/rankings.h5")

d$"coaches/Pipey" %>% 
  map(
    ~{as.data.frame(.) %>% 
        `colnames<-`(d$date) %>% 
        `rownames<-`(d$race_ids) %>% 
        mutate(race = rownames(.)) %>% 
        gather(date, ranking, -race) %>% 
        filter(!is.nan(ranking)) %>% 
        mutate(date = lubridate::ymd(date))
    } 
  ) %>%  
  bind_cols() %>% 
  select(date, race, mu = "ranking", phi="ranking1") %>% 
  group_by(race) %>% 
  arrange(date) %>% 
  mutate(next_date = lead(date, default = NA)) %>% 
  ggplot(aes(date, mu-2.5*phi, fill = race)) + 
  geom_rect(data = function(d) {e = 70; data_frame(date = rep(d$date, each = e+1), next_date = rep(d$next_date, each = e+1), race = rep(d$race, each = e+1), mu = rep(d$mu, each=e+1), phi = rep(d$phi, each = e+1), y =map2(d$mu, d$phi, ~seq(.x-.y*3, .x+.y*3, length.out = e+1)) %>% unlist) %>% group_by(date,race)%>% mutate(f = dnorm(y, mean=mu, sd=phi), width = y-lead(y))}, aes(xmin= date, xmax = next_date, ymin = y-width/2, ymax =y+(width/2) , fill=race, alpha = f, group=race), height = 3, colour=NA) + 
  geom_point(aes(y=mu, colour = race), alpha = 0.4, size = 0.5) + 
  facet_wrap(~race) + 
  geom_point(span = 0.2, colour="grey70", se=F, size = 0.5)+
  geom_smooth(span = 0.2, colour= "black", se = F) + 
  theme_nufflytics() + 
  xlab(NULL) + 
  ylab("mu") + 
  ggtitle("Rating changes","Pipey") + 
  viridis::scale_color_viridis(discrete = T) + 
  viridis::scale_fill_viridis(discrete = T) +
  theme(legend.position= "none")
