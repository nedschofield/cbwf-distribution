#' author: Ned Ryan-Schofield
#' date: 5/02/2026

#' Function for downloading files from the NCI thredds database (or other data repositories probably)
#' Catalogue url to browse the server is https://thredds.nci.org.au/thredds/catalog/catalog.html
#' To download, have to use https://thredds.nci.org.au/thredds/fileServer/... instead of catalog
#' 
#' @param download_urls character vector of download links starting with format above
#' @param out_dir character string of output file directory location. Function creates it if it doesn't exist

nci_thredds_download <- function(
    download_urls,
    out_dir,
    overwrite = FALSE,
    timeout = 180,
    quiet = FALSE
    ) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = timeout)
  
  out_paths <- file.path(out_dir, basename(download_urls))
  
  for (i in seq_along(download_urls)) {
    if (!overwrite && file.exists(out_paths[i])) {
      if (!quiet) message("Skipping: ", basename(out_paths[i]))
      next
    }
    
    if (!quiet) message("Downloading: ", basename(out_paths[i]))
    utils::download.file(
      url = download_urls[i],
      destfile = out_paths[i],
      mode = "wb",
      quiet = quiet
    )
  }
  
  invisible(out_paths)
}
