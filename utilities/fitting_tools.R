# Note: These functions assume that `model_definitions.R` will also be sourced
# before any of them are called.

# Choose min-A type
MIN_A_WITH_FRL <- FALSE

# Don't fit TPU, alpha, or RL, and use default GammaStar
FIXED_ALPHA <- 0
FIXED_GSTAR <- defaults$Gstar
FIXED_RD    <- defaults$Rd
FIXED_TP    <- 1000

FIT_OPTIONS <- list(
  alpha_old = FIXED_ALPHA,
  Gamma_star = FIXED_GSTAR,
  RL_at_25 = FIXED_RD,
  Tp_at_25 = FIXED_TP
)

# Returns a function that calculates error values using the original FvCB
fvcb_error <- function(
  C_seq,  # microbar
  An_seq, # micromol / m^2 / s
  Kc,     # microbar
  O,      # mbar
  Ko,     # mbar
  sd_A    # micromol / m^2 / s
)
{
  # Make an exdf input table
  inputs <- input_exdf(
    C_seq,
    Kc,
    O,
    Ko,
    PRESSURE,
    'Ci',
    return_exdf = TRUE
  )

  # Add the assimilation rates
  inputs <- set_variable(
    inputs,
    'A',
    'micromol m^(-2) s^(-1)',
    value = An_seq
  )

  # Get the error function, which returns the negative of the logarithm of the
  # likelihood
  error_function_c3_aci(
    inputs,
    fit_options = FIT_OPTIONS,
    sd_A = sd_A,
    atp_use = Wj_C_coef,
    nadph_use = Wj_G_coef
  )
}

# Returns a function that calculates error values using the min-A variant
min_A_error <- function(
  C_seq,    # microbar
  An_seq,   # micromol / m^2 / s
  Kc,       # microbar
  O,        # mbar
  Ko,       # mbar
  force_Ac, # set to TRUE for min-A+Ac, FALSE for min-A
  sd_A      # micromol / m^2 / s
)
{
  # Make an exdf input table
  inputs <- input_exdf(
    C_seq,
    Kc,
    O,
    Ko,
    PRESSURE,
    'Ci',
    return_exdf = TRUE
  )

  # Add the assimilation rates
  inputs <- set_variable(
    inputs,
    'A',
    'micromol m^(-2) s^(-1)',
    value = An_seq
  )

  # Get the error function, which returns the negative of the logarithm of the
  # likelihood
  error_function_c3_aci(
    inputs,
    fit_options = FIT_OPTIONS,
    sd_A = sd_A,
    atp_use = Wj_C_coef,
    nadph_use = Wj_G_coef,
    use_min_A = TRUE,
    use_FRL = force_Ac
  )
}

# Fits an A-C curve using the FvCB and psudeo-original FvCBs
fit_model <- function(
  C_seq,          # microbar (a vector to be fitted)
  An_seq,         # microbar (a vector to be fitted)
  Kc,             # microbar
  O,              # mbar
  Ko,             # mbar
  force_Ac        # boolean
)
{
  # Make an exdf input table
  inputs <- input_exdf(
    C_seq,
    Kc,
    O,
    Ko,
    PRESSURE,
    'Ci',
    return_exdf = TRUE
  )

  # Add the assimilation rates
  inputs <- set_variable(
    inputs,
    'A',
    'micromol m^(-2) s^(-1)',
    value = An_seq
  )

  # Fit the curve using the original FvCB
  fvcb_fit <- fit_c3_aci(
    inputs,
    atp_use = Wj_C_coef,
    nadph_use = Wj_G_coef,
    fit_options = FIT_OPTIONS,
    calculate_confidence_intervals = TRUE,
    remove_unreliable_param = 0
  )

  # Fit the curve using the min-A variant
  min_A_fit <- fit_c3_aci(
    inputs,
    atp_use = Wj_C_coef,
    nadph_use = Wj_G_coef,
    fit_options = FIT_OPTIONS,
    calculate_confidence_intervals = TRUE,
    remove_unreliable_param = 0,
    use_min_A = TRUE,
    use_FRL = force_Ac
  )

  # Return important results
  list(
    fvcb_fit = list(
      success = fvcb_fit$parameters[1, 'convergence'],
      feval = fvcb_fit$parameters[1, 'feval'],
      value = fvcb_fit$parameters[1, 'optimum_val'],
      J = fvcb_fit$parameters[1, 'J_at_25'],
      J_lower = fvcb_fit$parameters[1, 'J_at_25_lower'],
      J_upper = fvcb_fit$parameters[1, 'J_at_25_upper'],
      Vcmax = fvcb_fit$parameters[1, 'Vcmax_at_25'],
      Vcmax_lower = fvcb_fit$parameters[1, 'Vcmax_at_25_lower'],
      Vcmax_upper = fvcb_fit$parameters[1, 'Vcmax_at_25_upper'],
      RMSE = fvcb_fit$parameters[1, 'RMSE']
    ),
    min_A_fit = list(
      success = min_A_fit$parameters[1, 'convergence'],
      feval = min_A_fit$parameters[1, 'feval'],
      value = min_A_fit$parameters[1, 'optimum_val'],
      J = min_A_fit$parameters[1, 'J_at_25'],
      J_lower = min_A_fit$parameters[1, 'J_at_25_lower'],
      J_upper = min_A_fit$parameters[1, 'J_at_25_upper'],
      Vcmax = min_A_fit$parameters[1, 'Vcmax_at_25'],
      Vcmax_lower = min_A_fit$parameters[1, 'Vcmax_at_25_lower'],
      Vcmax_upper = min_A_fit$parameters[1, 'Vcmax_at_25_upper'],
      RMSE = min_A_fit$parameters[1, 'RMSE']
    )
  )
}

fit_and_plot <- function(
  C_seq,          # microbar (a vector to be fitted)
  An_seq,         # microbar (a vector to be fitted)
  Kc,             # microbar
  O,              # mbar
  Ko,             # mbar
  force_Ac,       # boolean
  bn,             # basename for output plot files
  main,           # title for plot
  xlim,
  ylim,
  make_plots      # boolean
)
{
  # Fit the A-Ci curve with both models
  fitres <- fit_model(
    C_seq,
    An_seq,
    Kc,
    O,
    Ko,
    force_Ac
  )

  # Calculate fitted An values
  an_seq_fvcb_fit <- fvcb_model(
    C_seq,
    fitres$fvcb_fit$Vcmax,
    Kc,
    O,
    Ko,
    fitres$fvcb_fit$J,
    FIXED_GSTAR,
    FIXED_TP,
    FIXED_ALPHA,
    FIXED_RD
  )

  an_seq_min_A_fit <- min_A(
    C_seq,
    fitres$min_A_fit$Vcmax,
    Kc,
    O,
    Ko,
    fitres$min_A_fit$J,
    FIXED_GSTAR,
    FIXED_TP,
    FIXED_ALPHA,
    FIXED_RD,
    force_Ac,
    TPU_threshold = -1
  )

  # Set the axis labels
  xlab <- 'Cc [ micromol / mol]'
  ylab <- 'Rates [micromol / m^2 / s]'

  # Specify line plots
  low_j_p1 <- xyplot(
    An_seq + an_seq_fvcb_fit$An + an_seq_min_A_fit$An ~ C_seq,
    type = 'b',
    grid = TRUE,
    auto.key = list(lines = TRUE, points = TRUE),
    ylim = ylim,
    xlim = xlim,
    ylab = ylab,
    xlab = xlab,
    main = main,
    par.settings = list(
      superpose.line = list(col = multi_curve_line_colors()),
      superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
    )
  )

  low_j_p2 <- xyplot(
    An_seq + an_seq_fvcb_fit$An + an_seq_fvcb_fit$Wc + an_seq_fvcb_fit$Wj ~ C_seq,
    type = 'b',
    pch = 16,
    grid = TRUE,
    auto.key = list(lines = TRUE, points = TRUE),
    ylim = ylim,
    xlim = xlim,
    ylab = ylab,
    xlab = xlab,
    main = main,
    par.settings = list(
      superpose.line = list(col = multi_curve_line_colors()),
      superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
    )
  )

  low_j_p3 <- xyplot(
    An_seq + an_seq_min_A_fit$An + an_seq_min_A_fit$Ac + an_seq_min_A_fit$Aj ~ C_seq,
    type = 'b',
    pch = 16,
    grid = TRUE,
    auto.key = list(lines = TRUE, points = TRUE),
    ylim = ylim,
    xlim = xlim,
    ylab = ylab,
    xlab = xlab,
    main = main,
    par.settings = list(
      superpose.line = list(col = multi_curve_line_colors()),
      superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
    )
  )

  # Make line plots
  to_plot <- list(
    list(plot_name = paste(bn, 'fit_p1'), plot_obj = low_j_p1),
    list(plot_name = paste(bn, 'fit_p2'), plot_obj = low_j_p2),
    list(plot_name = paste(bn, 'fit_p3'), plot_obj = low_j_p3)
  )

  if (make_plots) {
    pdf_print_list(to_plot, SAVE_TO_PDF, STD_WIDTH, STD_HEIGHT)
  }

  # Return the fit result
  fitres
}
