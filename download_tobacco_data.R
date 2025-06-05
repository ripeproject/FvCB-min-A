TOBACCO_ACI_URL <-
  'https://github.com/ripeproject/PhotoGEA-paper/raw/refs/heads/main/output_archive/rdata/tobacco_aci_raw_data.RData'

TOBACCO_FPATH <- file.path(OUTPUT_DIR, RDATA_DIR, 'tobacco_aci_raw_data.RData')

if (!file.exists(TOBACCO_FPATH)) {
  success <- tryCatch(
    download.file(TOBACCO_ACI_URL, TOBACCO_FPATH, 'auto'),
    error = function(e) {
      print(e)
      return(1)
    }
  )

  if (success != 0) {
    msg <- paste(
      'Could not automatically donwnload tobacco A-Ci data.',
      'There are three ways to proceed:\n',
      '(1) Try again\n',
      '(2) Manually download the file from', TOBACCO_ACI_URL,
      'and save it to ', TOBACCO_FPATH, '\n',
      '(3) Copy the archived version of the file from the `output_archive`',
      'directory to', TOBACCO_FPATH
    )
    stop(msg)
  }
}
