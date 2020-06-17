JS <- function(X) {
  X <- gsub(" ","_",X)
  js <- paste("// Braz_fl.js

  // Create a webpage object
  var page = require('webpage').create();

  // Include the File System module for writing to files
  var fs = require('fs');

  // Specify source and path to output file
  var url  = 'http://servicos.jbrj.gov.br/flora/search/",X,"'
  var path = 'Braz_fl.html'

  page.open(url, function (status) {
    var content = page.content;
    fs.write(path,content,'w')
    phantom.exit();
  });",sep="")
  writeLines(js, "Braz_fl.js")
}
#' Plant traits from online databases
#'
#' Get plant height data from online databases.
#' This function requires Internet connection.
#' If the desired database uses java script to store the data
#' (e.g. Flora of Brazil),
#' the function requires phantomjs (.exe for windows) stored in
#' the workind directory.
#'
#' @param X A vector of species names
#' @param database The database to search for traits. Available databases include:
#' \itemize{
#'   \item \dQuote{eflora_US} = Flora of North America.
#'   \item \dQuote{FNA} = Beta version of Flora of North America.
#'   \item \dQuote{plantzafrica} = Plants of southern Africa.
#'   \item \dQuote{wikipedia} = Wikipedia free online encyclopedia.
#'   \item \dQuote{reflora} = Flora of Brazil
#' }
#'
#' @rdname get_trait
#' @importFrom xml2 read_html xml_text
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
        } else if (database=="reflora") {
          if (!(file.exists("phantomjs") | file.exists("phantomjs.exe"))) {
            warning("PhantomJS does not exist; please download it from https://phantomjs.org/download.html and save it in your working directory")
          } else {
            JS(x)
            system("./phantomjs Braz_fl.js")
            txt <- xml2::read_html("Braz_fl.html")
            txt <- xml2::xml_text(txt)
            id <- "NA"
            txt <-  gsub("[ ]+", " ", txt)
            txt <- sub(".*(Description)", "\\1", txt)

            one <- sub(".*Stem:", "\\1", txt)
            one <- sub(";.* ", "", one)
            one <- gsub("administrator.*","",one)
            one <- strsplit(one, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)

            mxun <- one[[1]][length(one[[1]])]
            max <- gsub('([0-9]+) .*', '\\1', mxun)
            max <- as.numeric(gsub("[^[:digit:].]", "", max))

            unit <- sub("\\).*", "", sub(".*\\(", "", mxun))
            min <- one[[1]][length(one[[1]])-1]
            min <- trimws(min, "both", "\\D")
            min <- as.numeric(gsub("[^[:digit:].]", "", min))
            res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
          }
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
        } else if (database=="reflora") {
          if(!(file.exists("phantomjs") | file.exists("phantomjs.exe"))) {
            warning("PhantomJS does not exist; please download it from https://phantomjs.org/download.html and save it in your working directory")
          } else {
            JS(x)
            system("./phantomjs Braz_fl.js")
            txt <- xml2::read_html("Braz_fl.html")
            txt <- xml2::xml_text(txt)
            id <- "NA"
            txt <-  gsub("[ ]+", " ", txt)
            txt <- sub(".*(Description)", "\\1", txt)

            one <- sub(".*Stem:", "\\1", txt)
            one <- sub(";.* ", "", one)
            one <- gsub("administrator.*","",one)
            one <- strsplit(one, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)

            mxun <- one[[1]][length(one[[1]])]
            max <- gsub('([0-9]+) .*', '\\1', mxun)
            max <- as.numeric(gsub("[^[:digit:].]", "", max))

            unit <- sub("\\).*", "", sub(".*\\(", "", mxun))
            min <- one[[1]][length(one[[1]])-1]
            min <- trimws(min, "both", "\\D")
            min <- as.numeric(gsub("[^[:digit:].]", "", min))
            res <- as.data.frame(cbind(species=x, taxa=id, min, max, unit))
          }
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




