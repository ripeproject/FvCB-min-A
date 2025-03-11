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

CONSIDER_PRACTICALITY <- FALSE

# Define a function that determines the progression of limiting factors as C
# increases. The return values have the following meanings:
# 0: None of the other situations apply (because the two processes are equally
#    limiting at C = 0)
# 1: Carboxylation is always Rubisco-limited (because Cc is negative)
# 2: Carboxylation is always Rubisco-limited (because Cc is too high)
# 3: Carboxylation changes from Rubisco-limited to RuBP-regeneration-limited
# 4: Carboxylation is always RuBP-regeneration-limited (because Cc is negative)
# 5: Carboxylation is always RuBP-regeneration-limited (because Cc is too high)
# 6: Carboxylation changes from RuBP-regeneration-limited to Rubisco-limited
crossover <- function(
    Vcmax, # micromol / m^2 / s
    J,     # micromol / m^2 / s
    Kc,    # microbar
    O,     # mbar
    Ko,    # mbar
    Gstar, # microbar
    return_vals = FALSE
)
{
    cc_threshold <- 2000 # microbar

    jv <- J / (Wj_C_coef * Vcmax) # dimensionless

    Cc <- if (Vcmax == 0) {
        -1
    } else {
        (Kc * (1.0 + O / Ko) * jv - (Wj_G_coef / Wj_C_coef) * Gstar) / (1.0 - jv) # microbar
    }

    Acc0 <- -Gstar * Vcmax / (Kc * (1 + O / Ko)) # micromol / m^2 / s

    Ajc0 <- -J / Wj_G_coef # micromol / m^2 / s

    if (return_vals) {
      return(data.frame(
        Cc = Cc,
        Acc0 = Acc0,
        Ajc0 = Ajc0
      ))
    }

    if (Acc0 > Ajc0) {
        # In this case, carboxylation is Rubisco-limited at C = 0
        if (Cc < 0) {
            # In this case, no crossover occurs
            1
        } else if (Cc > cc_threshold) {
            # In this case, no crossover occurs at a reasonable value of C
            2
        } else {
            # In this case, a crossover occurs to RuBP-regeneration-limited carboxylation
            3
        }
    } else if (Acc0 < Ajc0) {
        # In this case, carboxylation is RuBP-regeneration-limited at C = 0
        if (Cc < 0) {
            # In this case, no crossover occurs
            4
        } else if (Cc > cc_threshold) {
            # In this case, no crossover occurs at a resonable value of C
            5
        } else {
            # In this case, a crossover occurs to Rubisco-limited carboxylation.
            # Note: this should never happen.
            6
        }
    } else {
        # In this case, the factors are equally limiting at C = 0. Note: this is
        # unlikely to ever happen.
        0
    }
}

# Calculate crossover values using default parameters
default_crossover_info <- crossover(
    defaults$Vcmax,
    defaults$J,
    defaults$Kc,
    defaults$O,
    defaults$Ko,
    defaults$Gstar,
    return_vals = TRUE
)

write.csv(default_crossover_info, file = file.path(OUTPUT_DIR, TABLE_DIR, 'default_crossover_info.csv'), row.names = FALSE)

# Set up the input values and the output matrices
npts <- 1001

v_seq <- seq(0, 300, length.out = npts)
j_seq <- seq(0, 300, length.out = npts)

z_cc_normal <- matrix(nrow = length(j_seq), ncol = length(v_seq))
z_cc_supplement <- z_cc_normal

# Fill in the matrices
for (i in seq_along(v_seq)) {
  for (j in seq_along(j_seq)) {
    z_cc_normal[i,j] <- crossover(v_seq[i], j_seq[j], defaults$Kc, defaults$O, defaults$Ko, defaults$Gstar)
    z_cc_supplement[i,j] <- crossover(v_seq[i], j_seq[j], supplement$Kc, defaults$O, defaults$Ko, supplement$Gstar)
  }
}

# Make the images
cc_image_normal <- list(x = v_seq, y = j_seq, z = z_cc_normal)
cc_image_supplement <- list(x = v_seq, y = j_seq, z = z_cc_supplement)

# Make a custom image plotting function. We cannot use `pdf_print_list` here
# because the `image` function does not return a plotting object. Here we
# (optionally) use different colors to distinguish regions where a transition
# occurs at a very high value of Cc (cases 2 and 5) and where it occurs at a
# more reasonable value of Cc (cases 3 and 6). If consider_practicality is
# FALSE, different colors will be used for cases 2 and 3 (and 5 and 6) to
# indicate whether a transition is completely forbidden or whether it is simply
# unlikely to be observed in practice.
plot_image <- function(F_image, image_name, consider_practicality, ...) {
    # Specify the colors to use for each z value
    cols <- c(
        '#808080', # case 0: gray
        '#0066cc', # case 1: dark blue
        if (consider_practicality) {'#99ccff'} else {'#e6f2ff'}, # case 2: middle blue or light blue
        '#e6f2ff', # case 3: light blue
        '#cc3300', # case 4: dark red
        if (consider_practicality) {'#ffb399'} else {'#ffece6'}, # case 5: middle red or light red
        '#ffece6'  # case 6: light red
    )

    # Breaks specifies the intervals used for setting colors based on z values.
    # Intervals are closed on the right and open on the left except for the
    # lowest interval which is closed at both ends. There must be one more
    # breakpoint than the number of colors.
    breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5)

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
    'transition_space_cj_normal_image',
    CONSIDER_PRACTICALITY,
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)

plot_image(
    cc_image_supplement,
    'transition_space_cj_supplement_image',
    CONSIDER_PRACTICALITY,
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)

# Define functions for Jn and Jd
Jn <- function(Vcmax, Gstar, Kc, O, Ko) {
    Wj_G_coef * Gstar * Vcmax / (Kc * (1 + O / Ko))
}

Jd <- function(Vcmax) {
    Wj_C_coef * Vcmax
}

Jd_seq      <- Jd(v_seq)
Jn_seq      <- Jn(v_seq, defaults$Gstar, defaults$Kc, defaults$O, defaults$Ko)
Jn_seq_supp <- Jn(v_seq, supplement$Gstar, supplement$Kc, defaults$O, defaults$Ko)

# Define plots of Jn and Jd
plot_obj_normal <- xyplot(
    Jd_seq + Jn_seq ~ v_seq,
    type = 'l',
    auto = TRUE,
    xlim = c(0, 100),
    ylim = c(0, 100),
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)

plot_obj_supplement <- xyplot(
    Jd_seq + Jn_seq_supp ~ v_seq,
    type = 'l',
    auto = TRUE,
    xlim = c(0, 100),
    ylim = c(0, 100),
    xlab = 'Vcmax (micromol / m^2 / s)',
    ylab = 'J (micromol / m^2 / s)'
)

to_plot <- list(
  list(plot_name = 'transition_space_cj_normal_curves', plot_obj = plot_obj_normal),
  list(plot_name = 'transition_space_cj_supplement_curves', plot_obj = plot_obj_supplement)
)

# Make plots
pdf_print_list(to_plot, SAVE_TO_PDF, STD_WIDTH, STD_HEIGHT)
