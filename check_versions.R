# Specify some required packages and their versions
PKGS <- list(
  lattice      = NA,
  latticeExtra = NA,
  PhotoGEA     = '1.2.0',
  RColorBrewer = NA
)

# Helping function for trying to load a package. Returns an error message if
# the package could not be loaded, and an empty string otherwise.
try_loading_package <- function(name) {
  if (require(name, character.only = TRUE, quietly = TRUE)) {
    ''
  } else {
    paste0(
      'The `',
      name,
      '` package could not be loaded; check `README.md` for ',
      'installation instructions.'
    )
  }
}

# Helper function for checking a package version. If the actual version is less
# than the expected one, it returns an error message string; otherwise, it
# returns an empty string. If the actual version is greater than the expected
# one, it sends a warning message.
check_package_version <- function(name, expected, actual) {
  if (is.na(expected)) {
    # No need to check the version for this package
    ''
  } else {
    comp_result <- compareVersion(
      as.character(packageVersion(name)),
      expected
    )

    if (comp_result > 0) {
      msg <- paste0(
        'Version `',
        expected,
        '` of the `',
        name,
        '` package is expected, but version `',
        actual, '` is installed. This could potentially cause issues, ',
        'but it may also be okay.'
      )
      warning(msg)
    }

    if (comp_result < 0) {
      paste0(
        'Version `',
        expected,
        '` of the `',
        name,
        '` package is expected, but an insufficient version (`',
        actual, '`) is installed.'
      )
    } else {
      ''
    }
  }
}

# Load and check packages
msgs <- sapply(seq_along(PKGS), function(i) {
  pkg_name         <- names(PKGS)[i]
  expected_version <- PKGS[[i]]

  msg <- try_loading_package(pkg_name)

  if (msg == '') {
    msg <- check_package_version(
      pkg_name,
      expected_version,
      packageVersion(pkg_name)
    )
  }

  msg
})

# Send messages if any errors occurred
pkg_problem <- msgs != ''

if (any(pkg_problem)) {
  msg <- paste0(
    'Issues were encountered with some required packages:\n  ',
    paste(msgs[pkg_problem], collapse = '\n  ')
  )
  stop(msg)
}
