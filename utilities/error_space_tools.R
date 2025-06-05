# Note: These functions assume that `model_definitions.R` and `output_tools.R`
# will also be sourced before any of them are called.

library(RColorBrewer)

MAX_ERROR_TO_PLOT <- 1
DEFAULT_ERROR_STEP_SIZE <- 0.05
DEFAULT_ERROR_THRESH <- 0.147

explore_error_space <- function(
  pcc_values, # microbar
  an_values,  # micromol / m^2 / s
  Kc,         # microbar
  O,          # mbar
  Ko,         # mbar
  force_Ac,
  v_seq,
  j_seq,
  fitres # output from fit_model
)
{
  # Get ready to explore error space
  z_min_W <- matrix(nrow = length(v_seq), ncol = length(j_seq))
  z_min_A <- z_min_W

  # Create error functions.
  min_W_error_func <- min_W_error(
    pcc_values, # microbar
    an_values,  # micromol / m^2 / s
    Kc,         # microbar
    O,          # mbar
    Ko,         # mbar
    fitres$min_W_fit$RMSE # micromol / m^2 / s
  )

  min_A_error_func <- min_A_error(
    pcc_values, # microbar
    an_values,  # micromol / m^2 / s
    Kc,         # microbar
    O,          # mbar
    Ko,         # mbar
    force_Ac,   # set to TRUE for min-A+Ac, FALSE for min-A
    fitres$min_A_fit$RMSE # micromol / m^2 / s
  )

  # Fill in the matrix with likelihood values
  for (i in seq_along(v_seq)) {
    for (j in seq_along(j_seq)) {
      param <- c(j_seq[j], v_seq[i])

      z_min_W[i,j]  <- exp(-min_W_error_func(param))
      z_min_A[i,j] <- exp(-min_A_error_func(param))
    }
  }

  # Get relative likelihood values
  z_min_W  <- z_min_W / exp(-fitres$min_W_fit$value)
  z_min_A <- z_min_A / exp(-fitres$min_A_fit$value)

  # Make the images
  min_W_image  <- list(x = v_seq, y = j_seq, z = z_min_W)
  min_A_image <- list(x = v_seq, y = j_seq, z = z_min_A)

  # Return the images in a list
  list(
    min_W_image = min_W_image,
    min_A_image = min_A_image
  )
}

# Make a colormap function
make_colormap_image <- function(step_size = DEFAULT_ERROR_STEP_SIZE) {
  x <- 0
  y <- seq(-step_size, MAX_ERROR_TO_PLOT, by = step_size) + step_size / 2
  z <- matrix(ncol = length(y), nrow = 1)
  z[1,] <- y
  return(list(x=x, y=y, z=z))
}

# Helping function for plotting (Vcmax, J)-space images, either as PDFs or to
# the current R session
plot_error_space_image <- function(
  F_image,
  name,
  step_size = DEFAULT_ERROR_STEP_SIZE,
  boundary = FALSE,
  boundary_threshold = DEFAULT_ERROR_THRESH,
  ...
)
{
  # Breaks specifies the intervals used for setting colors based on z values.
  # Intervals are closed on the right and open on the left except for the lowest
  # interval which is closed at both ends. There must be one more breakpoint
  # than the number of colors.
  breaks <- if (!boundary) {
    c(-100, seq(0, MAX_ERROR_TO_PLOT, by = step_size), 100)
  } else {
    c(-100, boundary_threshold, 100)
  }

  # Specify the colors to use
  cols <- if (!boundary) {
    c(
      '#ffffff',                                                    # white for z < 0
      colorRampPalette(brewer.pal(9, "Blues"))(length(breaks) - 3), # blue for intermediate values
      '#000000'                                                     # black for z > MAX_ERROR_TO_PLOT
    )
  } else {
    c('#ffffff', '#000000')
  }

  if (SAVE_TO_PDF) {
    pdf(
      file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(name, '.pdf')),
      width = STD_WIDTH,
      height = STD_HEIGHT,
      useDingbats = FALSE
    )
  } else {
    dev.new(width = STD_WIDTH, height = STD_HEIGHT)
  }

  image(
    F_image,
    breaks = breaks,
    col = cols,
    useRaster = TRUE,
    ...
  )

  if (SAVE_TO_PDF) {
    dev.off()
  }
}
