

# variables
line_col_palette <- c("#ff00ff", "#ff2400", "#6600cc", "#0000ff")
names(line_col_palette) <- c("MT-2_1", "MT-2_2", "MT-4_1", "MT-4_2")
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
        axis_title = 45,
        axis_text = 35,
        strip_text = 45,
        legend_title = 45,
        legend_text = 35,
        top_gap = 20,
        right_gap = 20
    ),
    presentation = list(
        title = 50,
        axis_title = 40,
        axis_text = 30,
        strip_text = 40,
        legend_title = 40,
        legend_text = 30,
        top_gap = 20,
        right_gap = 20
    ))

    size_set <- multiply_nested_list(size_set, scale)
    
    
    theme_bw() +
    theme(
        plot.title = element_text(size = size_set[[application]]["title"], hjust = 0.5),
        axis.title.x = element_text(
            size = size_set[[application]]["axis_title"],
            margin = margin(t = size_set[[application]]["top_gap"])  # space above x-axis title
        ),
        axis.title.y = element_text(
            size = size_set[[application]]["axis_title"],
            margin = margin(r = size_set[[application]]["right_gap"])  # space to the right of y-axis title
        ),
        axis.text.x = element_text(size = size_set[[application]]["axis_text"], color = "black", angle = 90),
        axis.text.y = element_text(size = size_set[[application]]["axis_text"], color = "black"),
        strip.text.x = element_text(size = 50, angle=0),
        strip.text.y = element_text(size = 50, angle=270),
        legend.title = element_text(size = size_set[[application]]["legend_title"]),
        legend.text = element_text(size = size_set[[application]]["legend_text"])
    )
}



custom_plot_save <- function(fig, name, width = 14, height = 14) {
    ## 1
    date <- str_split(Sys.time(), pattern = " ")[[1]][1]
    fig_dir1 <- paste0("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/publication/expiii/provisional/", date, "/")
    ifelse(!dir.exists(file.path(fig_dir1)), dir.create(file.path(fig_dir1)), FALSE)
    pdf(file = paste0(fig_dir1, name), height = height, width = width)
    print(fig)
    dev.off()
    ## 2
    fig_dir2 <- "/Users/alimos313/Documents/studies/phd/research/manuscripts/decelerating-adaptation/figures/"
    pdf(file = paste0(fig_dir2, name), height = height, width = width)
    print(fig)
    dev.off()
}
