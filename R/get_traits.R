#' Plant traits from online databases
#'
#' Get plant height data from online databases.
#' This function requires Internet connection.
#'
#' @param X A vector of species names
#' @param database The database to search for traits. Available databases include:
#' \itemize{
#'   \item \dQuote{eflora_US} = Flora of North America.
#'   \item \dQuote{FNA} = Beta version of Flora of North America.
#'   \item \dQuote{plantzafrica} = Plants of southern Africa.
#'   \item \dQuote{wikipedia} = Wikipedia free online encyclopedia.
#' }
#'
#' @rdname get_trait
#' @importFrom xml2 read_html xml_text
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return A dataframe containing the traits data for the species of interest.
#' @examples
#' sp <- c("Graptopetalum bartramii", "Salix raupii",
#'        "Berberis harrisoniana", "Quercus arkansana", "Aloe ferox")
#' get_trait(sp, database = "FNA")
#' @export
get_trait <- function(X, database = "FNA") {

  fna <- "http://beta.floranorthamerica.org/"
  efl_US <- "http://www.efloras.org/florataxon.aspx?flora_id=1&taxon_id="
  wiki <- "https://en.wikipedia.org/wiki/"
  pza <- "http://pza.sanbi.org/"

  if (length(X) > 50L && interactive()) {
    spo <- progress(X, function(x) {
      tryCatch({
        if (database=="wikipedia") {
          web <- paste0(wiki, gsub(" ", "_", as.character(x)))
          page = xml2::read_html(web)
          nat <- xml2::xml_text(page); nat <-  gsub("[ ]+", " ", nat)
          id <- "NA"
          nat <- sub(".*(Description)", "\\1", nat)
          nat <- trimws(nat, "left", "\\D")

          nat <- gsub("(high).*|(tall).*|(height).*|(long).*|(length).*","", nat)[1]
          nat <- sub("\\ and |\\ to ", "-", nat)
          nat <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")

          nat <- gsub('[^ -~]', '-', nat)
          max <- sub('.*\\-|.*(to)', '', nat)
          unit <- sub(".*? ", "", max); unit <- gsub("[^[:alnum:]]", "", unit)
          max <- gsub('([0-9]+) .*', '\\1', max)
          max <- as.numeric(gsub("[^[:digit:].]", "", max))
          min <- sub('\\-.*|(to).*', '', nat)
          min <- trimws(min, "both", "\\D")
          min <- as.numeric(gsub("[^[:digit:].]", "", min))
          res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
        } else if (database=="plantzafrica"){
          x <- tolower(x)
          web <- paste0(pza, gsub(" ", "-", as.character(x)))
          page = xml2::read_html(web)
          nat <- xml2::xml_text(page); nat <-  gsub("[ ]+", " ", nat)
          id <- "NA"
          nat <- sub(".*(Description)", "\\1", nat)
          nat <- trimws(nat, "left", "\\D")

          nat <- gsub("(high).*|(tall).*|(height).*|(long).*|(length).*","", nat)[1]
          nat <- sub("\\ and |\\ to ", "-", nat)
          nat <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")

          nat <- gsub('[^ -~]', '-', nat)
          max <- sub('.*\\-|.*(to)', '', nat)
          unit <- sub(".*? ", "", max); unit <- gsub("[^[:alnum:]]", "", unit)
          max <- gsub('([0-9]+) .*', '\\1', max)
          max <- as.numeric(gsub("[^[:digit:].]", "", max))
          min <- sub('\\-.*|(to).*', '', nat)
          min <- trimws(min, "both", "\\D")
          min <- as.numeric(gsub("[^[:digit:].]", "", min))
          res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
        } else {
          url <- switch(database,
                        FNA = paste0(fna, as.character(gsub(" ", "_", x))),
                        eflora_US = paste0(efl_US, as.character(x)))
          page <- xml2::read_html(url)
          nat <- xml2::xml_text(page); nat <-  gsub("[ ]+", " ", nat)
          id <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")
          one <- sub("Leaves.*", "\\1", nat); one <- sub(".*\\\t|.*\\\n", "\\1", one)
          nat <- trimws(one, "left", "\\D")

          nat <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")
          nat <- gsub('[^ -~]', '-', nat)
          max <- sub('.*\\-|.*(to)', '', nat)
          unit <- sub(".*? ", "", max); unit <- gsub("[^[:alnum:]]", "", unit)
          max <- gsub('([0-9]+) .*', '\\1', max)
          max <- as.numeric(gsub("[^[:digit:].]", "", max))
          min <- sub('\\-.*|(to).*', '', nat)
          min <- trimws(min, "both", "\\D")
          min <- as.numeric(gsub("[^[:digit:].]", "", min))
          res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
        }

        return(res)
      }, error = function(e) {NA})
    })
  } else {
    spo <- lapply(X, function(x) {
      tryCatch({
        if (database=="wikipedia") {
          web <- paste0(wiki, gsub(" ", "_", as.character(x)))
          page = xml2::read_html(web)
          nat <- xml2::xml_text(page); nat <-  gsub("[ ]+", " ", nat)
          id <- "NA"
          nat <- sub(".*(Description)", "\\1", nat)
          nat <- trimws(nat, "left", "\\D")

          nat <- gsub("(high).*|(tall).*|(height).*|(long).*|(length).*","", nat)[1]
          nat <- sub("\\ and |\\ to ", "-", nat)
          nat <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")

          nat <- gsub('[^ -~]', '-', nat)
          max <- sub('.*\\-|.*(to)', '', nat)
          unit <- sub(".*? ", "", max); unit <- gsub("[^[:alnum:]]", "", unit)
          max <- gsub('([0-9]+) .*', '\\1', max)
          max <- as.numeric(gsub("[^[:digit:].]", "", max))
          min <- sub('\\-.*|(to).*', '', nat)
          min <- trimws(min, "both", "\\D")
          min <- as.numeric(gsub("[^[:digit:].]", "", min))
          res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
        } else if (database=="plantzafrica"){
          x <- tolower(x)
          web <- paste0(pza, gsub(" ", "-", as.character(x)))
          page = xml2::read_html(web)
          nat <- xml2::xml_text(page); nat <-  gsub("[ ]+", " ", nat)
          id <- "NA"
          nat <- sub(".*(Description)", "\\1", nat)
          nat <- trimws(nat, "left", "\\D")

          nat <- gsub("(high).*|(tall).*|(height).*|(long).*|(length).*","", nat)[1]
          nat <- sub("\\ and |\\ to ", "-", nat)
          nat <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")

          nat <- gsub('[^ -~]', '-', nat)
          max <- sub('.*\\-|.*(to)', '', nat)
          unit <- sub(".*? ", "", max); unit <- gsub("[^[:alnum:]]", "", unit)
          max <- gsub('([0-9]+) .*', '\\1', max)
          max <- as.numeric(gsub("[^[:digit:].]", "", max))
          min <- sub('\\-.*|(to).*', '', nat)
          min <- trimws(min, "both", "\\D")
          min <- as.numeric(gsub("[^[:digit:].]", "", min))
          res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
        } else {
          url <- switch(database,
                        FNA = paste0(fna, as.character(gsub(" ", "_", x))),
                        eflora_US = paste0(efl_US, as.character(x)))
          page <- xml2::read_html(url)
          nat <- xml2::xml_text(page); nat <-  gsub("[ ]+", " ", nat)
          id <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")
          one <- sub("Leaves.*", "\\1", nat); one <- sub(".*\\\t|.*\\\n", "\\1", one)
          nat <- trimws(one, "left", "\\D")

          nat <- paste(unlist(strsplit(nat, split = "\\s+"))[1:2], collapse = " ")
          nat <- gsub('[^ -~]', '-', nat)
          max <- sub('.*\\-|.*(to)', '', nat)
          unit <- sub(".*? ", "", max); unit <- gsub("[^[:alnum:]]", "", unit)
          max <- gsub('([0-9]+) .*', '\\1', max)
          max <- as.numeric(gsub("[^[:digit:].]", "", max))
          min <- sub('\\-.*|(to).*', '', nat)
          min <- trimws(min, "both", "\\D")
          min <- as.numeric(gsub("[^[:digit:].]", "", min))
          res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
        }

        return(res)
      }, error = function(e) {NA})
    })
  }
  rr <- spo[!is.na(spo)]
  r <- data.frame(do.call(rbind, rr))
  r
}




