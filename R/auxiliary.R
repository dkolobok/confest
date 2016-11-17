#character or factor to numeric preseving the names
ch2n <- function(x) {
  y <- as.numeric(as.character(x))
  names(y) <- names(x)
  y
}
