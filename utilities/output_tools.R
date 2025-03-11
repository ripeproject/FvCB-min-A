library(PhotoGEA) # for pdf_print

OUTPUT_DIR <- 'output'
TABLE_DIR <- 'tables'
FIGURE_DIR <- 'figures'
RDATA_DIR <- 'rdata'

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

if (!dir.exists(file.path(OUTPUT_DIR, TABLE_DIR))) {
  dir.create(file.path(OUTPUT_DIR, TABLE_DIR))
}

if (!dir.exists(file.path(OUTPUT_DIR, FIGURE_DIR))) {
  dir.create(file.path(OUTPUT_DIR, FIGURE_DIR))
}

if (!dir.exists(file.path(OUTPUT_DIR, RDATA_DIR))) {
  dir.create(file.path(OUTPUT_DIR, RDATA_DIR))
}

pdf_print_list <- function(list_to_print, save_as_pdf, width = 3, height = 4) {
  invisible(lapply(seq_along(list_to_print), function(i) {
    plot_obj <- list_to_print[[i]]$plot_obj
    plot_name <- list_to_print[[i]]$plot_name

    pdf_print(
      plot_obj,
      save_to_pdf = save_as_pdf,
      width = width,
      height = height,
      file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(plot_name, '.pdf'))
    )
  }))
}
