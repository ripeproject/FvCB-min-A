library(PhotoGEA) # for calculate_c3_assimilation

# Decide whether to use Wj coefficients from Bernacchi et al. (2003) or von
# Caemmerer (2000). This will also determine which temperature response
# parameters to use for curve fitting.
JEQ_FROM_BERNACCHI <- FALSE

Wj_C_coef <- if (JEQ_FROM_BERNACCHI) {
  4.5
} else {
  4
}

Wj_G_coef <- if (JEQ_FROM_BERNACCHI) {
  10.5
} else {
  8
}

TEMPERATURE_PARAM <- if (JEQ_FROM_BERNACCHI) {
  c3_temperature_param_bernacchi
} else {
  c3_temperature_param_sharkey
}

PRESSURE <- 1 # bar

LIMIT_TOL <- 1e-10

CHECK_UNITS <- FALSE # Set to FALSE for a small speed boost

# Helping function to make an exdf object suitable for
# PhotoGEA:calculate_c3_assimilation
input_exdf <- function(
  C,        # microbar
  Kc,       # microbar
  O,        # mbar
  Ko,       # mbar
  pressure, # bar
  c_column_name = 'Cc',
  return_exdf = CHECK_UNITS
)
{
  # Make an exdf input table
  input_table <- exdf(
    data.frame(
      gmc_norm = rep_len(1, length(C)),
      J_norm = 1,
      Kc = Kc / pressure,
      Ko = Ko / pressure,
      oxygen = 0.1 * O / pressure,
      RL_norm = 1,
      total_pressure = pressure,
      Tp_norm = 1,
      Vcmax_norm = 1,
      Ca = 420
    ),
    units = data.frame(
      gmc_norm = 'normalized to gmc at 25 degrees C',
      J_norm = 'normalized to J at 25 degrees C',
      Kc = 'micromol mol^(-1)',
      Ko = 'mmol mol^(-1)',
      oxygen = 'percent',
      RL_norm = 'normalized to RL at 25 degrees C',
      total_pressure = 'bar',
      Tp_norm = 'normalized to Tp at 25 degrees C',
      Vcmax_norm = 'normalized to Vcmax at 25 degrees C',
      Ca = 'micromol mol^(-1)',
      stringsAsFactors = FALSE
    )
  )

  # If we aren't checking units, there's not much point in keeping track of
  # them, so we can just use data frames instead of exdf
  if (!return_exdf) {
    input_table <- input_table$main_data
  }

  set_variable(
    input_table,
    c_column_name,
    'micromol mol^(-1)',
    value = C / pressure
  )
}

# Use PhotoGEA to run the original FvCB
fvcb_model_exdf <- function(
  inputs,
  Vcmax, # micromol / m^2 / s
  J,     # micromol / m^2 / s
  Gstar, # microbar
  Tp,    # micromol / m^2 / s
  alpha, # dimensionless
  Rd     # micromol / m^2 / s
)
{
  # Calculate assimilation rates
  outputs <- calculate_c3_assimilation(
    inputs,
    0,                # alpha_g
    alpha,            # alpha_old
    0,                # alpha_s
    0,                # alpha_t
    Gstar / PRESSURE, # Gamma_star
    J,                # J_at_25
    Rd,               # RL_at_25
    Tp,               # Tp_at_25
    Vcmax,            # Vcmax_at_25
    atp_use = Wj_C_coef,
    nadph_use = Wj_G_coef,
    perform_checks = CHECK_UNITS
  )

  # Calculate limiting processes
  outputs <- identify_c3_limiting_processes(
    outputs,
    a_column_name = 'An',
    tol = LIMIT_TOL
  )

  # Return results as a list
  list(
    Wc = outputs[, 'Wc'],
    Wj = outputs[, 'Wj'],
    Wp = outputs[, 'Wp'],
    Vc = outputs[, 'Vc'],
    An = outputs[, 'An'],
    limit = outputs[, 'limiting_process']
  )
}

fvcb_model <- function(
  C,     # microbar
  Vcmax, # micromol / m^2 / s
  Kc,    # microbar
  O,     # mbar
  Ko,    # mbar
  J,     # micromol / m^2 / s
  Gstar, # microbar
  Tp,    # micromol / m^2 / s
  alpha, # dimensionless
  Rd     # micromol / m^2 / s
)
{
  # Make an exdf input table
  inputs <- input_exdf(
    C,
    Kc,
    O,
    Ko,
    PRESSURE
  )

  # Calculate assimilation rates
  fvcb_model_exdf(
    inputs,
    Vcmax, # micromol / m^2 / s
    J,     # micromol / m^2 / s
    Gstar, # microbar
    Tp,    # micromol / m^2 / s
    alpha, # dimensionless
    Rd     # micromol / m^2 / s
  )
}

# Use PhotoGEA to run the min-A variant
min_A_exdf <- function(
  inputs,
  Vcmax,        # micromol / m^2 / s
  J,            # micromol / m^2 / s
  Gstar,        # microbar
  Tp,           # micromol / m^2 / s
  alpha,        # dimensionless
  Rd,           # micromol / m^2 / s
  use_FRL,      # TRUE or FALSE
  TPU_threshold # microbar
)
{
  # Calculate assimilation rates
  outputs <- calculate_c3_assimilation(
    inputs,
    0,                # alpha_g
    alpha,            # alpha_old
    0,                # alpha_s
    0,                # alpha_t
    Gstar / PRESSURE, # Gamma_star
    J,                # J_at_25
    Rd,               # RL_at_25
    Tp,               # Tp_at_25
    Vcmax,            # Vcmax_at_25
    atp_use = Wj_C_coef,
    nadph_use = Wj_G_coef,
    use_min_A = TRUE,
    use_FRL = use_FRL,
    TPU_threshold = TPU_threshold / PRESSURE,
    perform_checks = CHECK_UNITS
  )

  # Calculate limiting processes
  outputs <- identify_c3_limiting_processes(
    outputs,
    a_column_name = 'An',
    tol = LIMIT_TOL
  )

  # Return results as a list
  list(
    Ac = outputs[, 'Ac'],
    Aj = outputs[, 'Aj'],
    Ap = outputs[, 'Ap'],
    An = outputs[, 'An'],
    limit = outputs[, 'limiting_process']
  )
}

min_A <- function(
  C,            # microbar
  Vcmax,        # micromol / m^2 / s
  Kc,           # microbar
  O,            # mbar
  Ko,           # mbar
  J,            # micromol / m^2 / s
  Gstar,        # microbar
  Tp,           # micromol / m^2 / s
  alpha,        # dimensionless
  Rd,           # micromol / m^2 / s
  use_FRL,      # TRUE or FALSE
  TPU_threshold # microbar
)
{
  # Make an exdf input table
  inputs <- input_exdf(
    C,
    Kc,
    O,
    Ko,
    PRESSURE
  )

  # Calculate assimilation rates
  min_A_exdf(
    inputs,
    Vcmax,        # micromol / m^2 / s
    J,            # micromol / m^2 / s
    Gstar,        # microbar
    Tp,           # micromol / m^2 / s
    alpha,        # dimensionless
    Rd,           # micromol / m^2 / s
    use_FRL,      # TRUE or FALSE
    TPU_threshold # microbar
  )
}

# Default settings for most figures
defaults <- list(
  Vcmax = 100,  # micromol / m^2 / s
  Kc = 259,     # microbar
  O = 210,      # mbar
  Ko = 179,     # mbar
  J = 170,      # micromol / m^2 / s
  Gstar = 38.6, # microbar
  Tp = 11.8,    # micromol / m^2 / s
  alpha = 0,    # dimensionless
  Rd = 1        # micromol / m^2 / s
)

# Define settings for supplemental figure S2
supplement <- list(
  Kc = 50,
  Gstar = 150
)
