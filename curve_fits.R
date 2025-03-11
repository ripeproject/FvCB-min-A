###               ###
### PRELIMINARIES ###
###               ###

# Load libraries
library(PhotoGEA)
library(lattice)
library(latticeExtra)
library(RColorBrewer)

# Clear workspace
rm(list=ls())

# Load some helping functions
source(file.path('utilities', 'model_definitions.R'))
source(file.path('utilities', 'output_tools.R'))

# Load data
load(file.path('data', 'tobacco_aci_raw_data.RData')) # load licor_data; created by `load_tobacco_data.R`

# Choose some settings
SAVE_TO_PDF <- TRUE
MAKE_NEW_CALCULATIONS <- TRUE

LIMIT_BY_GSTAR <- FALSE
CI_LOW_THRESHOLD <- 45

TPU_THRESHOLD <- 400

ITERMAX <- 200

PLOT_OPERATING_POINT <- FALSE
PLOT_AD <- TRUE

BASE_NAME <- 'tobacco_aci'

A_LIM          <- c(-5, 50)
A_LIM_ZOOM     <- c(-5, 5)
CI_LIM         <- c(-100, 1500)
CI_LIM_ZOOM    <- c(0, 100)
J_XLIM         <- c(0, 225)
J_YLIM_ABS     <- c(-25, 5)
J_YLIM_CI      <- c(0, 60)
J_YLIM_REL     <- c(-15, 15)
RL_XLIM        <- c(0, 3)
RL_YLIM_ABS    <- c(-2.5, 1)
RL_YLIM_CI     <- c(0, 2)
RL_YLIM_REL    <- c(-150, 150)
RMSE_XLIM      <- c(0, 1.8)
RMSE_YLIM_ABS  <- c(-0.25, 2.25)
RMSE_YLIM_REL  <- c(-20, 1700)
VCMAX_XLIM     <- c(0, 120)
VCMAX_YLIM_ABS <- c(-15, 10)
VCMAX_YLIM_CI  <- c(0, 15)
VCMAX_YLIM_REL <- c(-30, 30)

# Define a helping function that fits a curve with and without considering Ad,
# returning the option with a smaller RMSE
fit_c3_aci_Ad <- function(exdf_obj, ...) {
    res_without_Ad <- fit_c3_aci(
        exdf_obj,
        ...,
        consider_depletion = FALSE
    )

    res_without_Ad$parameters[, 'fit_type'] <- 'without Ad'

    RMSE_without_Ad <- res_without_Ad$parameters[1, 'RMSE']

    res_with_Ad <- fit_c3_aci(
        exdf_obj,
        ...,
        consider_depletion = TRUE
    )

    res_with_Ad$parameters[, 'fit_type'] <- 'with Ad'

    if (res_with_Ad$parameters[1, 'RMSE'] < RMSE_without_Ad) {
        res_with_Ad
    } else {
        res_without_Ad
    }
}

if (MAKE_NEW_CALCULATIONS) {
    ###                                ###
    ### FIT A-CI CURVES USING PHOTOGEA ###
    ###                                ###

    # Calculate total pressure in the Licor chamber
    licor_data <- calculate_total_pressure(licor_data)

    # Calculate temperature-dependent values of C3 photosynthetic parameters
    licor_data <- calculate_temperature_response(licor_data, TEMPERATURE_PARAM)

    # Get a version excluding points with C < Gstar
    licor_data[, 'Ci_low'] <- if (LIMIT_BY_GSTAR) {
        licor_data[, 'Ci'] < licor_data[, 'Gamma_star']
    } else {
        licor_data[, 'Ci'] < CI_LOW_THRESHOLD
    }

    licor_data_subset <- remove_points(
        licor_data,
        list(Ci_low = TRUE),
        method = 'exclude'
    )

    # Set a seed number because the deoptiom algorithm uses randomness
    set.seed(123)

    # Fit the A-Ci curves (original FvCB, all Ci)
    c3_aci_results <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci_Ad,                    # The function to apply to each chunk of `licor_data`
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(ITERMAX),
        atp_use = Wj_C_coef,
        nadph_use = Wj_G_coef
    ))

    print('done fitting set 1')

    # Fit the A-Ci curves (original FvCB, Ci above threshold)
    c3_aci_results_no_low <- consolidate(by(
        licor_data_subset,                       # The `exdf` object containing the curves
        licor_data_subset[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci,                              # The function to apply to each chunk of `licor_data`
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(ITERMAX),
        atp_use = Wj_C_coef,
        nadph_use = Wj_G_coef
    ))

    print('done fitting set 2')

    # Fit the A-Ci curves (min-A variant, all Ci)
    c3_aci_results_min_A <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(ITERMAX),
        atp_use = Wj_C_coef,
        nadph_use = Wj_G_coef,
        use_min_A = TRUE,
        TPU_threshold = TPU_THRESHOLD
    ))

    print('done fitting set 3')

    # Fit the A-Ci curves (min-A variant, Ci above threshold)
    c3_aci_results_no_low_min_A <- consolidate(by(
        licor_data_subset,                       # The `exdf` object containing the curves
        licor_data_subset[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci,                              # The function to apply to each chunk of `licor_data`
        Ca_atmospheric = 420,
        optim_fun = optimizer_deoptim(ITERMAX),
        atp_use = Wj_C_coef,
        nadph_use = Wj_G_coef,
        use_min_A = TRUE,
        TPU_threshold = TPU_THRESHOLD
    ))

    print('done fitting set 4')

    # Save the results
    write.csv.exdf(
        c3_aci_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters.csv'))
    )

    write.csv.exdf(
        c3_aci_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits.csv'))
    )

    write.csv.exdf(
        c3_aci_results_no_low$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_no_low.csv'))
    )

    write.csv.exdf(
        c3_aci_results_no_low$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_no_low.csv'))
    )

    write.csv.exdf(
        c3_aci_results_min_A$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_min_A.csv'))
    )

    write.csv.exdf(
        c3_aci_results_min_A$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_min_A.csv'))
    )

    write.csv.exdf(
        c3_aci_results_no_low_min_A$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_no_low_min_A.csv'))
    )

    write.csv.exdf(
        c3_aci_results_no_low_min_A$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_no_low_min_A.csv'))
    )

    save(
        c3_aci_results, c3_aci_results_no_low,
        c3_aci_results_min_A, c3_aci_results_no_low_min_A,
        file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData'))
    )
} else {
    load(file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData')))
}

# Plot the FvCB fits (including limiting rates and operating point)
pdf_print(
  plot_c3_aci_fit(
    c3_aci_results,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM,
    ylim = A_LIM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_fits.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM_ZOOM,
    ylim = A_LIM_ZOOM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_fits_zoom.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results_no_low,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM,
    ylim = A_LIM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_no_low_fits.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results_no_low,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM_ZOOM,
    ylim = A_LIM_ZOOM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_no_low_fits_zoom.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results_min_A,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM,
    ylim = A_LIM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_min_A_fits.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results_min_A,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM_ZOOM,
    ylim = A_LIM_ZOOM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_min_A_fits_zoom.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results_no_low_min_A,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM,
    ylim = A_LIM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_min_A_no_low_fits.pdf'))
)

pdf_print(
  plot_c3_aci_fit(
    c3_aci_results_no_low_min_A,
    'curve_identifier',
    'Ci',
    plot_operating_point = PLOT_OPERATING_POINT,
    plot_Ad = PLOT_AD,
    xlim = CI_LIM_ZOOM,
    ylim = A_LIM_ZOOM
  ),
  width = 15,
  height = 15,
  save_to_pdf = SAVE_TO_PDF,
  file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_min_A_no_low_fits_zoom.pdf'))
)


# Plot RMSE and parameter comparisons for each curve
plot_diff <- function(pname, xlim, ylim_rel, ylim_abs, ylim_ci, ordering = NULL, x_curve = FALSE, alternating = 0) {
    val_no_low     <- c3_aci_results_no_low$parameters[, pname]
    val_fvcb_all   <- c3_aci_results$parameters[, pname]
    val_min_A_all <- c3_aci_results_min_A$parameters[, pname]

    abs_diff_fvcb   <- val_fvcb_all - val_no_low
    abs_diff_min_A <- val_min_A_all - val_no_low

    rel_diff_fvcb   <- 100 * abs_diff_fvcb / val_no_low
    rel_diff_min_A <- 100 * abs_diff_min_A / val_no_low

    dataf <- data.frame(
        curve_identifier = c3_aci_results$parameters[, 'curve_identifier'],
        abs_diff_fvcb = abs_diff_fvcb,
        abs_diff_min_A = abs_diff_min_A,
        fit_type = c3_aci_results$parameters[, 'fit_type'],
        rel_diff_fvcb = rel_diff_fvcb,
        rel_diff_min_A = rel_diff_min_A,
        val_no_low = val_no_low
    )

    if (pname != 'RMSE') {
        pname_upper <- paste0(pname, '_upper')
        pname_lower <- paste0(pname, '_lower')

        upper_no_low <- c3_aci_results_no_low$parameters[, pname_upper]
        lower_no_low <- c3_aci_results_no_low$parameters[, pname_lower]
        upper_fvcb   <- c3_aci_results$parameters[, pname_upper]
        lower_fvcb   <- c3_aci_results$parameters[, pname_lower]
        upper_min_A <- c3_aci_results_min_A$parameters[, pname_upper]
        lower_min_A <- c3_aci_results_min_A$parameters[, pname_lower]

        ci_width_no_low <- upper_no_low - lower_no_low
        ci_width_fvcb   <- upper_fvcb - lower_fvcb
        ci_width_min_A <- upper_min_A - lower_min_A

        ci_width_no_low[is.infinite(ci_width_no_low)] <- NA
        ci_width_fvcb[is.infinite(ci_width_fvcb)]     <- NA
        ci_width_min_A[is.infinite(ci_width_min_A)] <- NA

        dataf <- cbind(dataf, data.frame(
            ci_width_no_low = ci_width_no_low,
            ci_width_fvcb = ci_width_fvcb,
            ci_width_min_A = ci_width_min_A
        ))
    }

    ordering <- if (is.null(ordering)) {
        order(dataf$fit_type, dataf$abs_diff_fvcb)
    } else {
        ordering
    }

    dataf <- dataf[ordering, ]

    if (x_curve) {
        pdf_print(
            xyplot(
                val_fvcb_all + val_min_A_all + val_no_low ~ seq_along(val_no_low) | fit_type,
                data = dataf,
                scales = list(alternating = alternating),
                type = 'p',
                auto.key = list(text = c('original FvCB', 'min-A variant', 'Ci above threshold'), space = 'top'),
                par.settings = list(superpose.symbol = list(pch = c(16, 1))),
                grid = TRUE,
                xlim = c(0, length(val_no_low) + 1),
                #ylim = ylim_abs,
                xlab = 'Curve',
                ylab = pname
            ),
            width = 12,
            save_to_pdf = SAVE_TO_PDF,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_comp_', pname, '_', alternating, '_raw.pdf'))
        )
    }

    pdf_print(
        xyplot(
            if (x_curve) {
                rel_diff_fvcb + rel_diff_min_A ~ seq_along(val_no_low) | fit_type
            } else {
                rel_diff_fvcb + rel_diff_min_A ~ val_no_low
            },
            data = dataf,
            scales = list(alternating = alternating),
            type = 'p',
            auto.key = list(text = c('original FvCB', 'min-A variant'), space = 'top'),
            par.settings = list(superpose.symbol = list(pch = c(16, 1))),
            grid = !x_curve,
            xlim = if (x_curve) {
                c(0, length(val_no_low) + 1)
            } else {
                xlim
            },
            ylim = ylim_rel,
            xlab = if (x_curve) {
                'Curve'
            } else {
                paste(pname, ' (Ci above threshold)')
            },
            ylab = paste('Relative change in', pname, 'when including all Ci'),
            panel = function(...) {
                if (x_curve) {
                    for (i in seq_along(rel_diff_fvcb)) {
                        panel.lines(
                            c(dataf$rel_diff_fvcb[i], dataf$rel_diff_min_A[i]) ~ c(i, i),
                            col = 'gray30'
                        )
                    }
                }

                panel.xyplot(...)
            }
        ),
        width = if (x_curve) {
            12
        } else {
            7
        },
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_comp_', pname, '_', alternating, '_rel.pdf'))
    )

    pdf_print(
        xyplot(
            if (x_curve) {
                abs_diff_fvcb + abs_diff_min_A ~ seq_along(val_no_low) | fit_type
            } else {
                abs_diff_fvcb + abs_diff_min_A ~ val_no_low
            },
            data = dataf,
            scales = list(alternating = alternating),
            type = 'p',
            auto.key = list(text = c('original FvCB', 'min-A variant'), space = 'top'),
            par.settings = list(superpose.symbol = list(pch = c(16, 1))),
            grid = !x_curve,
            xlim = if (x_curve) {
                c(0, length(val_no_low) + 1)
            } else {
                xlim
            },
            ylim = ylim_abs,
            xlab = if (x_curve) {
                'Curve'
            } else {
                paste(pname, ' (Ci above threshold)')
            },
            ylab = paste('Change in', pname, 'when including all Ci'),
            panel = function(...) {
                if (x_curve) {
                    for (i in seq_along(rel_diff_fvcb)) {
                        panel.lines(
                            c(dataf$abs_diff_fvcb[i], dataf$abs_diff_min_A[i]) ~ c(i, i),
                            col = 'gray30'
                        )
                    }
                }

                panel.xyplot(...)
            }
        ),
        width = if (x_curve) {
            12
        } else {
            7
        },
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_comp_', pname, '_', alternating, '_abs.pdf'))
    )

    if (x_curve && pname != 'RMSE') {
        pdf_print(
            xyplot(
                ci_width_fvcb + ci_width_min_A + ci_width_no_low ~ seq_along(val_no_low) | fit_type,
                data = dataf,
                scales = list(alternating = alternating),
                type = 'p',
                auto.key = list(text = c('FvCB (all Ci)', 'min-A (all Ci)', 'Ci above threshold'), space = 'top'),
                par.settings = list(superpose.symbol = list(pch = c(16, 1, 1))),
                xlim = c(0, length(val_no_low) + 1),
                ylim = ylim_ci,
                xlab = 'Curve',
                ylab = paste(pname, 'confidence interval width'),
                panel = function(...) {
                    for (i in seq_along(rel_diff_fvcb)) {
                        ci_widths <- c(
                          dataf$ci_width_fvcb[i],
                          dataf$ci_width_min_A[i],
                          dataf$ci_width_no_low[i]
                        )

                        minwidth <- if (all(is.na(ci_widths))) {
                          NA
                        } else {
                          min(ci_widths, na.rm = TRUE)
                        }

                        maxwidth <- if (all(is.na(ci_widths))) {
                          NA
                        } else {
                          max(ci_widths, na.rm = TRUE)
                        }

                        panel.lines(
                            c(minwidth, maxwidth) ~ c(i, i),
                            col = 'gray30'
                        )
                    }

                    panel.xyplot(...)
                }
            ),
            width = 12,
            save_to_pdf = SAVE_TO_PDF,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_comp_', pname, '_', alternating, '_ci_width.pdf'))
        )
    }

    invisible(ordering)
}

ord <- plot_diff('RMSE', RMSE_XLIM, RMSE_YLIM_REL, RMSE_YLIM_ABS, NULL, x_curve = TRUE, alternating = 0)

plot_diff('J_at_25',     J_XLIM,     J_YLIM_REL,     J_YLIM_ABS,     J_YLIM_CI,     ordering = ord, x_curve = TRUE, alternating = 0)
plot_diff('J_at_25',     J_XLIM,     J_YLIM_REL,     J_YLIM_ABS,     J_YLIM_CI,     ordering = ord, x_curve = TRUE, alternating = 1)
plot_diff('RL_at_25',    RL_XLIM,    RL_YLIM_REL,    RL_YLIM_ABS,    RL_YLIM_CI,    ordering = ord, x_curve = TRUE, alternating = 0)
plot_diff('RL_at_25',    RL_XLIM,    RL_YLIM_REL,    RL_YLIM_ABS,    RL_YLIM_CI,    ordering = ord, x_curve = TRUE, alternating = 1)
plot_diff('RMSE',        RMSE_XLIM,  RMSE_YLIM_REL,  RMSE_YLIM_ABS,  NULL,          ordering = ord, x_curve = TRUE, alternating = 1)
plot_diff('Vcmax_at_25', VCMAX_XLIM, VCMAX_YLIM_REL, VCMAX_YLIM_ABS, VCMAX_YLIM_CI, ordering = ord, x_curve = TRUE, alternating = 0)
plot_diff('Vcmax_at_25', VCMAX_XLIM, VCMAX_YLIM_REL, VCMAX_YLIM_ABS, VCMAX_YLIM_CI, ordering = ord, x_curve = TRUE, alternating = 1)

# Summarize this information in a table
fit_result_summary <- data.frame(
    ord, # just a placeholder, will be replaced later
    as.character(c3_aci_results$parameters[, 'curve_identifier']),
    c3_aci_results$parameters[, 'fit_type'],

    c3_aci_results$parameters[, 'RMSE'],
    c3_aci_results$parameters[, 'Vcmax_at_25'],
    c3_aci_results$parameters[, 'J_at_25'],
    c3_aci_results$parameters[, 'Tp_at_25'],
    c3_aci_results$parameters[, 'alpha_old'],
    c3_aci_results$parameters[, 'RL_at_25'],

    c3_aci_results_no_low$parameters[, 'RMSE'],
    c3_aci_results_no_low$parameters[, 'Vcmax_at_25'],
    c3_aci_results_no_low$parameters[, 'J_at_25'],
    c3_aci_results_no_low$parameters[, 'Tp_at_25'],
    c3_aci_results_no_low$parameters[, 'alpha_old'],
    c3_aci_results_no_low$parameters[, 'RL_at_25'],

    c3_aci_results_min_A$parameters[, 'RMSE'],
    c3_aci_results_min_A$parameters[, 'Vcmax_at_25'],
    c3_aci_results_min_A$parameters[, 'J_at_25'],
    c3_aci_results_min_A$parameters[, 'Tp_at_25'],
    c3_aci_results_min_A$parameters[, 'alpha_old'],
    c3_aci_results_min_A$parameters[, 'RL_at_25'],

    c3_aci_results_no_low_min_A$parameters[, 'RMSE'],
    c3_aci_results_no_low_min_A$parameters[, 'Vcmax_at_25'],
    c3_aci_results_no_low_min_A$parameters[, 'J_at_25'],
    c3_aci_results_no_low_min_A$parameters[, 'Tp_at_25'],
    c3_aci_results_no_low_min_A$parameters[, 'alpha_old'],
    c3_aci_results_no_low_min_A$parameters[, 'RL_at_25']
)

# We supply column names this way so R doesn't replace spaces with periods
colnames(fit_result_summary) <- c(
    'Curve number',
    'Curve identifier',
    'Curve type',

    'RMSE (original FvCB, all Ci)',
    'Vcmax_at_25 (original FvCB, all Ci)',
    'J_at_25 (original FvCB, all Ci)',
    'Tp_at_25 (original FvCB, all Ci)',
    'alpha_old (original FvCB, all Ci)',
    'RL_at_25 (original FvCB, all Ci)',

    'RMSE (original FvCB, Ci above threshold)',
    'Vcmax_at_25 (original FvCB, Ci above threshold)',
    'J_at_25 (original FvCB, Ci above threshold)',
    'Tp_at_25 (original FvCB, Ci above threshold)',
    'alpha_old (original FvCB, Ci above threshold)',
    'RL_at_25 (original FvCB, Ci above threshold)',

    'RMSE (min-A variant, all Ci)',
    'Vcmax_at_25 (min-A variant, all Ci)',
    'J_at_25 (min-A variant, all Ci)',
    'Tp_at_25 (min-A variant, all Ci)',
    'alpha_old (min-A variant, all Ci)',
    'RL_at_25 (min-A variant, all Ci)',

    'RMSE (min-A variant, Ci above threshold)',
    'Vcmax_at_25 (min-A variant, Ci above threshold)',
    'J_at_25 (min-A variant, Ci above threshold)',
    'Tp_at_25 (min-A variant, Ci above threshold)',
    'alpha_old (min-A variant, Ci above threshold)',
    'RL_at_25 (min-A variant, Ci above threshold)'
)

fit_result_summary <- fit_result_summary[ord, ]

fit_result_summary[['Curve number']] <- seq_len(nrow(fit_result_summary))

row.names(fit_result_summary) <- NULL

write.csv(fit_result_summary, file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_fit_result_summary.csv')), row.names = FALSE)

# Print mean RMSE values
print(mean(c3_aci_results$parameters[, 'RMSE']))
print(mean(c3_aci_results_no_low$parameters[, 'RMSE']))
print(mean(c3_aci_results_min_A$parameters[, 'RMSE']))
print(mean(c3_aci_results_no_low_min_A$parameters[, 'RMSE']))

# Plot RMSE histograms
rmse_df <- rbind(
    within(c3_aci_results$parameters$main_data,              {fit_type = 'Original FvCB (all Ci)'}),
    within(c3_aci_results_no_low$parameters$main_data,       {fit_type = 'Original FvCB (Ci above threshold)'}),
    within(c3_aci_results_min_A$parameters$main_data,        {fit_type = 'Min-A+FTT (all Ci)'}),
    within(c3_aci_results_no_low_min_A$parameters$main_data, {fit_type = 'Min-A+FTT (Ci above threshold)'})
)

pdf_print(
    histogram(
        ~ RMSE | fit_type,
        data = rmse_df,
        breaks = seq(min(RMSE_YLIM_ABS), max(RMSE_YLIM_ABS), by = 0.1),
        layout = c(1, 4),
        xlim = RMSE_YLIM_ABS,
        ylim = c(0, 14),
        type = 'count'
    ),
    width = 4.5,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_comp_RMSE.pdf'))
)
