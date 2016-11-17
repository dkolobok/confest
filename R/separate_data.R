#separates experimental table (cm$data_ode) into a list of data.frames
#with the same experimental conditions (only time differs)
#' Data separator
#'
#' gets a data.frame, divides it into a list of data.frames so that each data.frame has unique experimental conditions
#' @param dat data.frame to be separated
#' @keywords likelihood
#' @export
#' @examples
#' separate_data(dat)
#'
#'
#

separate_data <- function(dat) {
  #determine column indeces
  state_ind <- grep('iv_', names(dat))
  state <- dat[state_ind]
  names(state) <- gsub('iv_', '', names(state))
  time_ind <- which(names(dat) == 'time')
  #  time <- dat['time']
  yname_ind <- which(names(dat) == 'yname')
  yvalue_ind <- which(names(dat) == 'yvalue')
  noise_ind <- which(names(dat) == 'noise')
  parms <- dat[-c(state_ind, time_ind, yname_ind, yvalue_ind, noise_ind)]
  #  parms <- dat[-c(state_ind, time_ind, yname_ind, yvalue_ind)]
  #  browser()
  #  split(dat, dat[,-c(time_ind, yname_ind, yvalue_ind)])
  split(dat, dat[,-c(time_ind, yname_ind, yvalue_ind, noise_ind)])
}

assign_all <- function(dat, cm) {
  dat1 <- as.character(unlist(dat[1,]))
  names(dat1) <- names(dat)
  #determine column indeces
  state_ind <- grep('init_', names(dat1))
  state <- rep(0, length(cm$ode_state_names))
  names(state) <- cm$ode_state_names
  state1 <- dat1[state_ind]
  names(state1) <- gsub('init_', '', names(state1))
  state[names(state1)] <- state1
  time_ind <- which(names(dat1) == 'time')
  time <- dat1['time']
  yname_ind <- which(names(dat1) == 'yname')
  yvalue_ind <- which(names(dat1) == 'yvalue')
  noise_ind <- which(names(dat) == 'noise')
  parms <- dat1[-c(state_ind, time_ind, yname_ind, yvalue_ind, noise_ind)]
  list(state = state, parms = parms, time = time)
}

