

# variables

line_col_palette <- c("#ff00ff", "#ff2400", "#66ff66", "#00cc99", "#009900", "#558000")
names(line_col_palette) <- exp_line_factor
breaks_step <- 50


# define functions

multiply_nested_list <- function(nested_list, multiplier) {
  lapply(nested_list, function(inner) {
    lapply(inner, function(x) x * multiplier)
  })
}


# functions
custom_plot_theme <- function(application = "publication", scale = 1) {
    size_set <- list(
    publication = list(
        title = 0,
        axis_title = 40,
        axis_text = 25,
        strip_text = 40,
        legend_title = 40,
        legend_text = 25
    ),
    presentation = list(
        title = 50,
        axis_title = 40,
        axis_text = 30,
        strip_text = 40,
        legend_title = 40,
        legend_text = 30
    ))

    size_set <- multiply_nested_list(size_set, scale)
    
    
    theme_bw() +
    theme(
        plot.title = element_text(size = size_set[[application]]["title"], hjust = 0.5),
        axis.title = element_text(size = size_set[[application]]["axis_title"]),
        axis.text.x = element_text(size = size_set[[application]]["axis_text"], color = "black", angle = 90),
        axis.text.y = element_text(size = size_set[[application]]["axis_text"], color = "black"),
        strip.text.x = element_text(size = 50, angle=0),
        legend.title = element_text(size = size_set[[application]]["legend_title"]),
        legend.text = element_text(size = size_set[[application]]["legend_text"])
    )
}



custom_plot_save <- function(fig, name, width = 14, height = 14) {
    ## 1
    date <- str_split(Sys.time(), pattern = " ")[[1]][1]
    fig_dir1 <- paste0("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/publication/2e-vs-btk/provisional/", date, "/")
    ifelse(!dir.exists(file.path(fig_dir1)), dir.create(file.path(fig_dir1)), FALSE)
    pdf(file = paste0(fig_dir1, name), height = height, width = width)
    print(fig)
    dev.off()
    ## 2
    fig_dir2 <- "/Users/alimos313/Documents/studies/phd/research/manuscripts/2e-vs-btk-manuscript/figures/"
    pdf(file = paste0(fig_dir2, name), height = height, width = width)
    print(fig)
    dev.off()
}