convertPairToString <- function(x, col.1, col.2) {
	gsub(" ", "", paste(x[col.1], x[col.2], sep = "_"))
}
