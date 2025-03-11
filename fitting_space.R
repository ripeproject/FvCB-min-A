# Clear workspace
rm(list=ls())

# Load required packages
library(lattice)

# Load helping functions
source(file.path('utilities', 'model_definitions.R'))
source(file.path('utilities', 'fitting_tools.R'))
source(file.path('utilities', 'error_space_tools.R'))
source(file.path('utilities', 'output_tools.R'))

# Choose settings
STD_WIDTH <- 6
STD_HEIGHT <- 6
SAVE_TO_PDF <- TRUE

NPTS <- 501

# LOW LIGHT A-Ci

# Set the sequence of C values to use for simulating measured A-Ci curves.
c_seq <- c(5, 10, 20, 30, 50, 70, 90, 110, 130, 175, 200, 300, 500) # microbar

# Simulate An values using a low J
an_seq <- fvcb_model(
  c_seq,
  defaults$Vcmax,
  defaults$Kc,
  defaults$O,
  defaults$Ko,
  35,
  FIXED_GSTAR,
  FIXED_TP,
  FIXED_ALPHA,
  FIXED_RD
)$An

# Add some random noise to the simulated data (up to a specified percentage of
# the true An)
set.seed(1234)
noise_frac <- 0.2
an_seq <- an_seq * (1 + runif(length(an_seq), -noise_frac, noise_frac))

# Fit the A-Ci curve with both models and print the results
res <- fit_and_plot(
  c_seq,
  an_seq,
  defaults$Kc,
  defaults$O,
  defaults$Ko,
  MIN_A_WITH_FRL,
  'simulated_low_j',
  'simulated_low_j',
  c(-20, 520),
  c(-5, 7),
  TRUE
)

# Choose Vcmax and J sequences
v_seq <- seq(10, 100, length.out = NPTS)
j_seq <- seq(25, 45, length.out = NPTS)

# Explore error space
error_res <- explore_error_space(
  c_seq,       # microbar
  an_seq,      # micromol / m^2 / s
  defaults$Kc, # microbar
  defaults$O,  # mbar
  defaults$Ko, # mbar
  MIN_A_WITH_FRL,
  v_seq,
  j_seq,
  res
)

# Plot the images
jlab <- 'J [ micromol / m^2 / s ]'
vcmaxlab <- 'Vcmax [micromol / m^2 / s]'

plot_error_space_image(make_colormap_image(), 'error_colormap')
plot_error_space_image(error_res$fvcb_image,  'simulated_aci_curve fvcb_error',                           xlab = vcmaxlab, ylab = jlab)
plot_error_space_image(error_res$min_A_image, 'simulated_aci_curve min_A_error',                          xlab = vcmaxlab, ylab = jlab)
plot_error_space_image(error_res$fvcb_image,  'simulated_aci_curve fvcb_error_contour',  boundary = TRUE, xlab = vcmaxlab, ylab = jlab)
plot_error_space_image(error_res$min_A_image, 'simulated_aci_curve min_A_error_contour', boundary = TRUE, xlab = vcmaxlab, ylab = jlab)

# Extract and reorganize the fitting results
fit_res <- rbind(
  within(as.data.frame(res$fvcb_fit),  {type = 'original FvCB'}),
  within(as.data.frame(res$min_A_fit), {type = 'min_A'})
)

# Save the fit results
write.csv(fit_res, file.path(OUTPUT_DIR, TABLE_DIR, 'simulated_aci_curve_fitting_best_estimates.csv'), row.names = FALSE)
