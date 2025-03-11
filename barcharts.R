# Clear workspace
rm(list=ls())

# Load required packages
library(lattice)

# Load helping functions
source(file.path('utilities', 'output_tools.R'))

# Choose settings
SAVE_TO_PDF <- TRUE
LABELS <- TRUE
STD_WIDTH <- 5
STD_HEIGHT <- 4

label_string <- if (LABELS) {
    'labels'
} else {
    'no_labels'
}

# Load data
top_cited_since_1990 <- read.delim(file.path('data', 'top_cited_since_1990.tsv'))
top_cited_since_2013 <- read.delim(file.path('data', 'top_cited_since_2013.tsv'))
fitting_tools <- read.delim(file.path('data', 'curve_fitting_tools.tsv'))

# Get paper and tool counts
since_1990_paper_count <- with(top_cited_since_1990, {
  tapply(equations_paper, equations_paper, length)
})

since_2013_paper_count <- with(top_cited_since_2013, {
  tapply(equations_paper, equations_paper, length)
})

curve_fitting_tool_count <- with(fitting_tools, {
  tapply(limiting_factor_determination, limiting_factor_determination, length)
})

# Get citation counts
since_1990_citation_count <- with(top_cited_since_1990, {
  tapply(times_cited, equations_paper, sum)
})

since_2013_citation_count <- with(top_cited_since_2013, {
  tapply(times_cited, equations_paper, sum)
})

# Form data frames
paper_count_fractions <- rbind(
  data.frame(
    model = names(since_1990_paper_count),
    frac = as.numeric(since_1990_paper_count) / nrow(top_cited_since_1990),
    count = as.numeric(since_1990_paper_count),
    type = 'since_1990'
  ),
  data.frame(
    model = names(since_2013_paper_count),
    frac = as.numeric(since_2013_paper_count) / nrow(top_cited_since_2013),
    count = as.numeric(since_2013_paper_count),
    type = 'since_2013'
  )
)

citation_count_fractions <- rbind(
  data.frame(
    model = names(since_1990_citation_count),
    frac = as.numeric(since_1990_citation_count) / sum(top_cited_since_1990$times_cited),
    type = 'since_1990'
  ),
  data.frame(
    model = names(since_2013_citation_count),
    frac = as.numeric(since_2013_citation_count) / sum(top_cited_since_2013$times_cited),
    type = 'since_2013'
  )
)

tool_count_totals <- data.frame(
  model = names(curve_fitting_tool_count),
  count = as.numeric(curve_fitting_tool_count),
  type = 'tools'
)

col_to_keep <- c('citation', 'tpu_c_lim')

# Get info about C restrictions on TPU
tpu_c_lim_info <- rbind(
    top_cited_since_1990[, col_to_keep],
    top_cited_since_2013[, col_to_keep],
    fitting_tools[, col_to_keep]
)

tpu_c_lim_info <- unique(tpu_c_lim_info)

tpu_c_lim_info <- tpu_c_lim_info[order(tpu_c_lim_info$tpu_c_lim, tpu_c_lim_info$citation), ]

tpu_c_lim_info <- tpu_c_lim_info[tpu_c_lim_info$citation != 'Tu (2013)', ]

# Set column order for plotting
paper_count_fractions$model <- factor(
  paper_count_fractions$model,
  levels = c('original FvCB', 'other', 'min-A')
)

citation_count_fractions$model <- factor(
  citation_count_fractions$model,
  levels = c('original FvCB', 'other', 'min-A')
)

tool_count_totals$model <- factor(
  tool_count_totals$model,
  levels = c(
    'original FvCB',
    'original FvCB+FTT',
    'manual',
    'exhaustive',
    'min-A',
    'min-A+FRL+FTT'
  )
)

# Define plots
to_plot <- list(

  list(
    plot_name = paste0('since_1990_paper_counts_', label_string),
    plot_obj = barchart(
      count ~ type,
      groups = model,
      data = paper_count_fractions[paper_count_fractions$type == 'since_1990', ],
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, nrow(top_cited_since_1990) * 1.1),
      ylab = 'Number of papers'
    )
  ),

  list(
    plot_name = paste0('since_1990_paper_count_fractions_', label_string),
    plot_obj = barchart(
      frac ~ type,
      groups = model,
      data = paper_count_fractions[paper_count_fractions$type == 'since_1990', ],
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, 1.1),
      ylab = 'Fraction of papers'
    )
  ),

  list(
    plot_name = paste0('since_1990_citation_count_fractions_', label_string),
    plot_obj = barchart(
      frac ~ type,
      groups = model,
      data = paper_count_fractions[paper_count_fractions$type == 'since_1990', ],
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, 1.1),
      ylab = 'Fraction of citations'
    )
  ),

  list(
    plot_name = paste0('since_2013_paper_counts_', label_string),
    plot_obj = barchart(
      count ~ type,
      groups = model,
      data = paper_count_fractions[paper_count_fractions$type == 'since_2013', ],
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, nrow(top_cited_since_2013) * 1.1),
      ylab = 'Number of papers'
    )
  ),

  list(
    plot_name = paste0('since_2013_paper_count_fractions_', label_string),
    plot_obj = barchart(
      frac ~ type,
      groups = model,
      data = paper_count_fractions[paper_count_fractions$type == 'since_2013', ],
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, 1.1),
      ylab = 'Fraction of papers'
    )
  ),

  list(
    plot_name = paste0('since_2013_citation_count_fractions_', label_string),
    plot_obj = barchart(
      frac ~ type,
      groups = model,
      data = paper_count_fractions[paper_count_fractions$type == 'since_2013', ],
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, 1.1),
      ylab = 'Fraction of citations'
    )
  ),

  list(
    plot_name = paste0('tool_count_totals_', label_string),
    plot_obj = barchart(
      count ~ type,
      groups = model,
      data = tool_count_totals,
      stack = TRUE,
      auto = LABELS,
      ylim = c(0, nrow(fitting_tools) * 1.1),
      ylab = "Number of fitting tools"
    )
  )
)

# Make plots
pdf_print_list(to_plot, SAVE_TO_PDF, STD_WIDTH, STD_HEIGHT)

# Save citation counts
write.csv(as.data.frame(since_1990_paper_count),   file = file.path(OUTPUT_DIR, TABLE_DIR, 'since_1990_paper_count.csv'))
write.csv(as.data.frame(since_2013_paper_count),   file = file.path(OUTPUT_DIR, TABLE_DIR, 'since_2013_paper_count.csv'))
write.csv(as.data.frame(curve_fitting_tool_count), file = file.path(OUTPUT_DIR, TABLE_DIR, 'curve_fitting_tool_count.csv'))
write.csv(as.data.frame(tpu_c_lim_info),           file = file.path(OUTPUT_DIR, TABLE_DIR, 'tpu_c_lim.csv'), row.names = FALSE)
