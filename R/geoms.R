### ggplot2 theme element for no space between panels.
theme_no_space <- function(...){
  tryCatch({
    theme(panel.spacing=grid::unit(0, "lines"), ...)
  }, error=function(e){
    theme(panel.margin=grid::unit(0, "lines"), ...)
  })
}

geom_tallrect <- function
### ggplot2 geom with xmin and xmax aesthetics that covers the entire
### y range, useful for clickSelects background elements.
(mapping = NULL,
 data = NULL,
 stat = "identity", position = "identity",
 ...,
 na.rm = FALSE,
 show.legend = NA,
 inherit.aes = TRUE) {
  ggplot2::layer(
    geom = GeomTallRect,
    data = data,
    mapping = mapping,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

### ggproto object for geom_tallrect
GeomTallRect <- ggplot2::ggproto(
  "GeomTallRect", ggplot2::Geom,
  default_aes = ggplot2::aes(colour = "black",
                    fill = "grey35", 
                    size = 1, 
                    linetype = 1,
                    alpha = 0.5),
  
  required_aes = c("xmin", "xmax"),
  
  draw_panel = function(self, data, 
                        panel_scales, coord) {
    coords <- coord$transform(data, panel_scales)
    ymax <- grid::unit(1, "npc")
    ymin <- grid::unit(0, "npc")
    grid::rectGrob(
      coords$xmin, ymax,
      width = coords$xmax - coords$xmin,
      height = ymax - ymin,
      default.units = "native",
      just = c("left", "top"),
      gp = grid::gpar(
        col = coords$colour,
        fill = scales::alpha(coords$fill, 
                             coords$alpha), 
        lwd = coords$size * .pt,
        lty = coords$linetype,
        lineend = "butt"
      )
    )
  },
  
  draw_key=ggplot2::draw_key_polygon
  
)

