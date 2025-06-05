# Check to make sure we have the correct version of PhotoGEA
library(PhotoGEA)

# Check versions of key R packages
print('running: check_versions.R')
source('check_versions.R')

# Make sure the output directories exist
source(file.path('utilities', 'output_tools.R'))

# Download the tobacco A-Ci curve data
print('running: download_tobacco_data.R')
source('download_tobacco_data.R')

# Make barcharts used in Figure 1
print('running: barcharts.R')
source('barcharts.R')

# Make simulated A-Ci curves used in Figures 2, 3, and 4
print('running: model_comparison.R')
source('model_comparison.R')

# Make C-J limitation space maps and curves (Figure 4)
print('running: transition_space_cj.R')
source('transition_space_cj.R')

# Make C-J-P limitation space maps (Figure 4)
print('running: transition_space_cjp.R')
source('transition_space_cjp.R')

# Make components of Figure 5
print('running: fitting_space.R')
source('fitting_space.R')

# Make components of Figure 6
print('running: curve_fits.R')
source('curve_fits.R')
