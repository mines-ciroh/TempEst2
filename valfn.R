library(tidyverse)
library(hydroGOF)
library(sf)
library(ggspatial)

states <- map_data("state")
eco.raw <- if(!exists("eco.raw")) {
  st_read("../Ecoregions/NA_CEC_Eco_Level1.shp", quiet=TRUE)
} else eco.raw

plot.eco <- function(vir.random=FALSE) {
  ecos <- c("Northern Forests", "Eastern Temperate Forests", "Tropical Wet Forests", "Northwestern Forested Mountains", "Great Plains", "North American Deserts", "Mediterranean California", "Marine West Coast Forest", "Temperate Sierras", "Southern Semiarid Highlands")
  ecoreg <- eco.raw %>% st_transform("WGS84") %>%
    rename(ecoregion = NA_L1NAME) %>%
    mutate(ecoregion = str_to_title(ecoregion)) %>% filter(ecoregion %in% ecos)
  # Order: ETF, GP, MWCF, MC, NAD, NF, NWFM, SSH, TS, TWF
  plt <- ggplot(ecoreg) +
    geom_sf(aes(fill=if (vir.random) fct_shuffle(as.factor(ecoregion)) else ecoregion), color=NA) +
    geom_polygon(aes(long, lat, group=group), data=states, fill=NA,
                 color=if (vir.random) "grey" else "black", size=1) +
    scale_y_continuous(limits=c(
      min(states$lat),
      max(states$lat)
    )) +
    scale_x_continuous(limits=c(
      min(states$long),
      max(states$long)
    )) +
    (if (vir.random) {
      scale_fill_viridis_d()
    } else {
      scale_fill_discrete(
        # Order: ETF, GP, MWCF, MC, NAD, NF, NWFM, SSH, TS, TWF
        type = c(
          "#AA22AA", "#999955", "#BBBBDD", "orange", "#FF6666", "#D5D5BC", "skyblue",
          "lightpink", "darkgrey", "lightblue"
        ),
        guide = guide_legend(nrow=3, title.position = "top")
      )
    })
    # annotation_north_arrow(which_north = "grid", location="br") +
    # annotation_scale(location="bl", height=unit(1, "cm"), text_cex=2)
  
  plt
}

theme_std <- function(...) {
  theme_bw() +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(
        fill="white",
        color="white"
      ),
      ...
    )
}
theme_vert <- function(...) {
  theme_std(
    axis.text.x = element_text(
      angle=90, vjust=0.5, hjust=1
    ),
    ...
  )
}

twrap <- scale_x_discrete(labels = function(x) str_wrap(str_to_title(x), width=15))
