# Clear workspace
rm(list=ls())

# Load required packages
library(lattice)

# Load helping functions
source(file.path('utilities', 'model_definitions.R'))
source(file.path('utilities', 'output_tools.R'))

# Choose settings
STD_WIDTH <- 6
STD_HEIGHT <- 6
SAVE_TO_PDF <- TRUE

MAKE_NEW_CALCULATIONS <- TRUE

NPTS_VJ     <- 601
C_STEP_SIZE <- 0.25

if (MAKE_NEW_CALCULATIONS){
  # Set the sequence of C values
  C_seq <- seq(1, 2000, by = C_STEP_SIZE) # microbar

  # Define a function that determines the sequence of limiting factors as C
  # changes according to C_seq. The output is determined as follows:
  #  0: Co-limitation across all C, or any situation not covered below
  #  1: Ac
  #  2: Ac -> Aj
  #  3: Ac -> Ap
  #  4: Ac -> Aj -> Ap
  #  5: Aj
  #  6: Aj -> Ac
  #  7: Aj -> Ap
  #  8: Aj -> Ac -> Ap
  #  9: Aj -> Ac -> Aj       (only possible in min-A)
  # 10: Aj -> Ac -> Aj -> Ap (only possible in min-A)
  # 11: Ac -> Aj -> Ac       (only possible in min-A)
  # 12: Ac -> Aj -> Ac -> Ap (only possible in min-A)
  sequence_id <- function(limits) {
    limit_seq <- limits[1]

    for (i in seq(2, length(limits))) {
      if (limit_seq[length(limit_seq)] != limits[i]) {
        limit_seq <- append(limit_seq, limits[i])
      }
    }

    if (identical(limit_seq, 'Ac') || identical(limit_seq, c('Ac', 'co-limited (Ac and Aj)')) || identical(limit_seq, c('Ac', 'co-limited (Ac and Aj)', 'Ac'))) {
      1
    } else if (identical(limit_seq, c('Ac', 'Aj'))) {
      2
    } else if (identical(limit_seq, c('Ac', 'Ap'))) {
      3
    } else if (identical(limit_seq, c('Ac', 'Aj', 'Ap'))) {
      4
    } else if (identical(limit_seq, 'Aj') || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj'))) {
      5
    } else if (identical(limit_seq, c('Aj', 'Ac')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj', 'Ac')) || identical(limit_seq, c('Aj', 'Ac', 'co-limited (Ac and Aj)', 'Ac')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Ac')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj', 'co-limited (Ac and Aj)', 'Ac'))) {
      6
    } else if (identical(limit_seq, c('Aj', 'Ap')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj', 'Ap')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj', 'co-limited (Aj and Ap)', 'Ap'))) {
      7
    } else if (identical(limit_seq, c('Aj', 'Ac', 'Ap')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj', 'Ac', 'Ap')) || identical(limit_seq, c('Aj', 'Ac', 'co-limited (Ac and Aj)', 'Ac', 'Ap')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Ac', 'Ap')) || identical(limit_seq, c('Aj', 'co-limited (Ac and Aj)', 'Aj', 'co-limited (Ac and Aj)', 'Ac', 'Ap'))) {
      8
    } else if (identical(limit_seq, c('Aj', 'Ac', 'Aj'))) {
      9
    } else if (identical(limit_seq, c('Aj', 'Ac', 'Aj', 'Ap'))) {
      10
    } else if (identical(limit_seq, c('Ac', 'Aj', 'Ac'))) {
      11
    } else if (identical(limit_seq, c('Ac', 'Aj', 'Ac', 'Ap'))) {
      12
    } else if (identical(limit_seq, 'co-limited (Ac and Aj)')) {
      0
    } else {
      print(limit_seq)
      0
    }
  }

  # Set up inputs for FvCB calculations
  default_inputs <- input_exdf(
    C_seq,       # microbar
    defaults$Kc, # microbar
    defaults$O,  # mbar
    defaults$Ko, # mbar
    PRESSURE,    # bar
    c_column_name = 'Cc'
  )

  supplement_inputs <- input_exdf(
    C_seq,         # microbar
    supplement$Kc, # microbar
    defaults$O,    # mbar
    defaults$Ko,   # mbar
    PRESSURE,      # bar
    c_column_name = 'Cc'
  )

  # Define functions for extracting sequences of limiting factors
  min_W_limits <- function(
      input_exdf,
      Vcmax,
      J,
      Gstar,
      Tp,
      alpha,
      Rd
  )
  {
      full_res <- min_W_exdf(
          input_exdf,
          Vcmax,
          J,
          Gstar,
          Tp,
          alpha,
          Rd
      )

      sequence_id(full_res$limit)
  }

  min_A_limits <- function(
      input_exdf,
      Vcmax,
      J,
      Gstar,
      Tp,
      alpha,
      Rd,
      use_FRL,
      TPU_threshold
  )
  {
      full_res <- min_A_exdf(
          input_exdf,
          Vcmax,
          J,
          Gstar,
          Tp,
          alpha,
          Rd,
          use_FRL,
          TPU_threshold
      )

      sequence_id(full_res$limit)
  }

  # Set up the input values and output matrices
  v_seq <- seq(0, 300, length.out = NPTS_VJ)
  j_seq <- seq(0, 300, length.out = NPTS_VJ)

  z_cc_normal     <- matrix(nrow = length(v_seq), ncol = length(j_seq))
  z_cc_supplement <- z_cc_normal
  z_cc_min_A      <- z_cc_normal
  z_cc_min_A_FRL  <- z_cc_normal

  # Fill in the matrices
  for (i in seq_along(v_seq)) {
    print(i / NPTS_VJ * 100)
    for (j in seq_along(j_seq)) {
      z_cc_normal[i,j]     <- min_W_limits(default_inputs,    v_seq[i], j_seq[j], defaults$Gstar,   defaults$Tp, defaults$alpha, defaults$Rd)
      z_cc_supplement[i,j] <- min_W_limits(supplement_inputs, v_seq[i], j_seq[j], supplement$Gstar, defaults$Tp, defaults$alpha, defaults$Rd)
      z_cc_min_A[i,j]      <- min_A_limits(default_inputs,   v_seq[i], j_seq[j], defaults$Gstar,   defaults$Tp, defaults$alpha, defaults$Rd, FALSE, -1)
      z_cc_min_A_FRL[i,j]  <- min_A_limits(default_inputs,   v_seq[i], j_seq[j], defaults$Gstar,   defaults$Tp, defaults$alpha, defaults$Rd, TRUE,  -1)
    }
  }

  # Make the images
  cc_image_normal     <- list(x = v_seq, y = j_seq, z = z_cc_normal)
  cc_image_supplement <- list(x = v_seq, y = j_seq, z = z_cc_supplement)
  cc_image_min_A      <- list(x = v_seq, y = j_seq, z = z_cc_min_A)
  cc_image_min_A_FRL  <- list(x = v_seq, y = j_seq, z = z_cc_min_A_FRL)

  # Save results
  save(
      cc_image_normal, cc_image_supplement,
      cc_image_min_A, cc_image_min_A_FRL,
      file = file.path(OUTPUT_DIR, RDATA_DIR, 'transition_space_cjp.RData')
  )
} else {
  load(file.path(OUTPUT_DIR, RDATA_DIR, 'transition_space_cjp.RData'))
}

# Make a custom image plotting function
plot_image <- function(F_image, image_name, ...) {
    # Specify the colors to use for each z value
    cols <- c(
        '#000000', # case 0:  black
        '#0066cc', # case 1:  dark blue
        '#e6f2ff', # case 2:  light blue
        '#b78f5e', # case 3:  tan
        '#ddd1c3', # case 4:  light tan
        '#cc3300', # case 5:  dark red
        '#ffece6', # case 6:  light red
        '#493a2d', # case 7:  brown
        '#5eb5b5', # case 8:  teal
        '#88c987', # case 9:  light green
        '#406b3f', # case 10: dark green
        '#ba87c9', # case 11: light purple
        '#4A3650'  # case 12: dark purple
    )

    # Breaks specifies the intervals used for setting colors based on z values.
    # Intervals are closed on the right and open on the left except for the
    # lowest interval which is closed at both ends. There must be one more
    # breakpoint than the number of colors.
    breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5)

    if (SAVE_TO_PDF) {
        pdf(
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(image_name, '.pdf')),
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

# Plot the images
plot_image(
    cc_image_normal,
    'transition_space_cjp_normal_image',
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)

plot_image(
    cc_image_supplement,
    'transition_space_cjp_supplement_image',
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)

plot_image(
    cc_image_min_A,
    'transition_space_cjp_min_A_image',
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)


plot_image(
    cc_image_min_A_FRL,
    'transition_space_cjp_min_A_FRL_image',
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)
