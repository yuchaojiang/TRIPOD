convertTrioToString <- function(x, col.1, col.2, col.3){
	gsub(" ", "", paste(x[col.1], x[col.2], x[col.3], sep = "_"))
}
