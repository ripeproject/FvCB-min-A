# Clear workspace
rm(list=ls())

# Load required packages
library(lattice)
library(PhotoGEA) # for multi_curve_colors

# Load helping functions
source(file.path('utilities', 'model_definitions.R'))
source(file.path('utilities', 'output_tools.R'))

# Choose settings
STD_WIDTH <- 6
STD_HEIGHT <- 6
SAVE_TO_PDF <- TRUE

# Make a function that runs both versions of the model across a range of C
# values and then plots the results
fvcb_compare <- function(
  Vcmax, # micromol / m^2 / s
  Kc,    # microbar
  O,     # mbar
  Ko,    # mbar
  J,     # micromol / m^2 / s
  Gstar, # microbar
  Tp,    # micromol / m^2 / s
  alpha, # dimensionless
  Rd,    # micromol / m^2 / s
  plot_type,
  C_seq = seq(1, 800, by = 2),
  clim = c(0, max(C_seq)),
  alim = c(-20, 50),
  use_FRL = FALSE,
  TPU_threshold = -1
)
{
  min_W_res <- min_W(
    C_seq,
    Vcmax,
    Kc,
    O,
    Ko,
    J,
    Gstar,
    Tp,
    alpha,
    Rd
  )

  min_A_res <- min_A(
    C_seq,
    Vcmax,
    Kc,
    O,
    Ko,
    J,
    Gstar,
    Tp,
    alpha,
    Rd,
    use_FRL,
    TPU_threshold
  )

  if (plot_type == 1) {
    xyplot(
      An + Vc + Wc + Wj + Wp ~ C_seq,
      data = min_W_res,
      type = 'l',
      lty = c(1, 1, 2, 2, 2),
      lwd = c(3, 3, 2, 2, 2),
      par.settings = list(
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
      ),
      grid = TRUE,
      xlim = clim,
      ylim = alim,
      auto.key = list(
        lines = TRUE,
        points = FALSE
      ),
      xlab = 'C (microbar)',
      ylab = 'Rates (micromol / m^2 / s)',
      main = 'min-W'
    )
  } else if (plot_type == 2) {
    xyplot(
      An + Ac + Aj + Ap ~ C_seq,
      data = min_A_res,
      type = 'l',
      lty = c(1, 2, 2, 2),
      lwd = c(3, 2, 2, 2),
      par.settings = list(
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
      ),
      grid = TRUE,
      xlim = clim,
      ylim = alim,
      auto.key = list(
        lines = TRUE,
        points = FALSE
      ),
      xlab = 'C (microbar)',
      ylab = 'Rates (micromol / m^2 / s)',
      main = 'min-A'
    )
  } else if (plot_type == 3) {
    xyplot(
      min_W_res$An + min_A_res$An ~ C_seq,
      type = 'l',
      lty = c(1, 2),
      lwd = c(3, 2),
      par.settings = list(
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
      ),
      grid = TRUE,
      xlim = clim,
      ylim = alim,
      xlab = 'C (microbar)',
      ylab = 'An (micromol / m^2 / s)',
      main = 'min-W vs. min-A',
      auto.key = list(
        text = c('min-W', 'min-A'),
        lines = TRUE,
        points = FALSE
      )
    )
  }
}

# Define plots
to_plot <- list(
  list(plot_name = "bright_1", plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 1}))),
  list(plot_name = "bright_2", plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 2}))),
  list(plot_name = "bright_3", plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 3}))),
  list(plot_name = "alpha_1",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 1; alpha = 0.5; Tp = 8; Vcmax = 120}))),
  list(plot_name = "alpha_2",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 2; alpha = 0.5; Tp = 8; Vcmax = 120}))),
  list(plot_name = "alpha_3",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 3; alpha = 0.5; Tp = 8; Vcmax = 120}))),
  list(plot_name = "alpha_4",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 3; alpha = 0.5; Tp = 8; Vcmax = 120; TPU_threshold = 400}))),
  list(plot_name = "dim_1",    plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 1; J = 10; C_seq = seq(1, 200, by = 0.5); alim = c(-10, 10)}))),
  list(plot_name = "dim_2",    plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 2; J = 10; C_seq = seq(1, 200, by = 0.5); alim = c(-10, 10)}))),
  list(plot_name = "dim_3",    plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 3; J = 10; C_seq = seq(1, 200, by = 0.5); alim = c(-10, 10)}))),
  list(plot_name = "weird_1",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 1; Kc = 50; Gstar = 150; Vcmax = 30; J = 180}))),
  list(plot_name = "weird_2",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 2; Kc = 50; Gstar = 150; Vcmax = 30; J = 180}))),
  list(plot_name = "weird_3",  plot_obj = do.call(fvcb_compare, within(defaults, {plot_type = 3; Kc = 50; Gstar = 150; Vcmax = 30; J = 180})))
)

# Make plots
pdf_print_list(to_plot, SAVE_TO_PDF, STD_WIDTH, STD_HEIGHT)
