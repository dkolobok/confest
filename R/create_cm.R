#' A custommodel-class constructor
#'
#' Takes all the information needed and returns a cm-class object
#' @param data_ode a data.frame with ode data or a name of xlsx file with data
#' @param data_expl_list a list of data.frames with explicit data or a name of xlsx file with data (each data.frame in a separate sheet)
#' @param explfunlist a list of explicit functions (length equal to length of data_expl_list)
#' @param implfunlist a list of implicit functions (which need an ode solution)
#' @param cname a name of the compiled c file to be loaded
#' @param ode_names a vector of names of y's for ode system
#' @param ode_parnames a vector of names of parameters for ode system
#' @param fit_start a named vector of initial values for fitting
#' @keywords custommodel
#' @export
#' @examples
#' cat_function()
#'
#'

create_cm <- function(data_ode = NULL, data_expl_list = NULL,
                      explfunlist = NULL, implfunlist = NULL,
                      cname = NULL, ode_names = NULL,
                      ode_parnames = NULL, fit_start = NULL) {
  cm <- list()
  class(cm) <- 'custommodel'
  if (!is.null(data_ode)) if (class(data_ode) == 'data.frame') {
    cm$data_ode <- data_ode
  } else if (class(try({
    wb <- loadWorkbook(data_ode)
    cm$data_ode <- readWorksheet(wb, sheet = getSheets(wb))
  })) != 'data.frame') stop('wrong ode data')
  if (!is.null(data_expl_list)) if (class(data_expl_list) == 'list') {
    cm$data_expl_list <- data_expl_list
  } else if (class(try({
    wb <- loadWorkbook(data_expl_list)
    temp <- readWorksheet(wb, sheet = getSheets(wb))
    if (class(temp) == 'data.frame') temp <- list(temp)
    cm$data_expl_list <- temp
  })) != 'list') stop('wrong explicit data')

  cm$explfunlist <- explfunlist
  cm$implfunlist <- implfunlist

  if (!is.null(cname)) {
    cm$cname <- cname
    dyn.load(paste('data/', cname, .Platform$dynlib.ext, sep = ""))
  }
  if (!is.null(ode_names)) cm$ode_names <- ode_names
  if (!is.null(ode_parnames)) cm$ode_parnames <- ode_parnames
  if (!is.null(fit_start)) cm$fit_start <- fit_start
  cm$scalevec <- fit_start
  cm$timeseq <- unique(as.numeric(as.character(cm$data_ode$time)))
  cm
}
