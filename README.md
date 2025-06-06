## FvCB-min-A
This repository includes R scripts and input data that were used to produce
the figures of the manuscript titled "Widely Used Variants of the
Farquhar-von-Caemmerer-Berry Model Can Cause Errors in Parameter Estimation,"
currently available on [bioRxiv](https://doi.org/10.1101/2025.03.11.642611).

These scripts have been tested using the following installations:
- Windows:
  - R version 4.3.2 (2023-10-31 ucrt)
  - Windows 10 Enterprise version 22H2
- Windows:
  - R version 4.4.2 (2024-10-30 ucrt)
  - Windows 11 Enterprise version 10.0.26100
- Linux:
  - R version 4.3.3 (2024-02-29)
  - Ubuntu 24.04.2 LTS (GNU/Linux 6.8.0-49-generic x86_64)

All outputs were generated using these scripts and stored in the
`output_archive` directory.

Figures and tables in the manuscript were produced from the files in
`output_archive/figures` and `output_archive/tables`. For figures, Adobe
Illustrator was used to perform the final combinations, recoloring, annotations,
etc.

### Reproducing the Outputs

#### Requirements
- The [R environment](https://cran.r-project.org/)
- The `lattice` R package
- The `latticeExtra` R package
- The `RColorBrewer` R package
- The `PhotoGEA` R package (version `1.2.0`)

#### Installing the R packages

Most of the required packages can be installed as usual by calling
`install.packages` from within R:

```R
install.packages('lattice')
install.packages('latticeExtra')
install.packages('RColorBrewer')
```

Version `1.2.0` of the
[PhotoGEA R package](https://eloch216.github.io/PhotoGEA/) can be installed
from within R by calling:

```R
remotes::install_github('eloch216/PhotoGEA', ref = 'v1.2.0')
```

Note that this command required the `remotes` package, which can be installed
via `install.packages('remotes')`.

#### Steps
- Start a fresh R session and set the working directory to this one
- Run the main script: `source('run_all.R')`
- All outputs (will be generated in a new `output` directory, whose content
  should be identical to the `output_archive` directory provided with the
  repository
- **Note:** The `transition_space_cjp.R` script may take a long time to run. It
  prints progress updates to the R terminal; the script will finish when this
  number reaches 100.

#### Additional Comments
- The main script (```run_all.R```) simply calls several other scripts. As an
  alternative to sourcing that script, individual analysis scripts can be
  sourced instead to produce a subset of all outputs.
- Many settings can be changed in these scripts, as well as in the settings and
  helping functions defined in the `utilities` directory.

## License
This repository is licensed under the CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/)
