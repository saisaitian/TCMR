
#' translateDrugID
#' @param query the identifer of compounds in other database
#' @param from the first
#' @param to  the second
#' @return a data.frame
#' @export
#' @details Interface to CTS (http://cts.fiehnlab.ucdavis.edu/) for metabolite identifier translation between
#'  > 200 of the most common biological databases including: Chemical Name, InChIKey, PubChem CID,
#'  ChemSpider, BioCyc, ChEBI, CAS, HMDB, KEGG and LipidMAPS.
#' @examples
#' query<-c("C15973","C00026")
#' from<-"KEGG"
#' to<-"PubChem CID"
#' translateDrugID(query,from,to)



translateDrugID = function(query, from, to) {

  if(is.null(query)){
    stop("query should not be NULL.\n")
  }

  if(is.list(query)){
    stop("query should be a vector.\n")
  }

  message(paste0("\n", ">>> Translating ",from,"  to  ", to))

  out <- vector('list',length(query))

  for (i in seq_along(query)) {

    display.progress(index = i, totalN = length(query))

    name <- query[i]

    out[[i]] <- single_translate(name, from, to)

    #out[[i]] <- ifelse(is.null(out[[i]])==TRUE ,NA,out[[i]])

  }
  names(out) <- query

  res <- plyr::ldply(out, data.frame)

  names(res) <- c(from,to)

  return(res)

}

single_translate = function(query, from, to){
  query_url <-  paste0('https://cts.fiehnlab.ucdavis.edu/rest/convert/', URLencode(from), '/', URLencode(to), '/', URLencode(query))
  res <- unlist(jsonlite::fromJSON(RCurl::getURL(query_url))$results)
  res <- ifelse(is.null(res)==TRUE ,NA,res)
  res <- ifelse(length(res)>1,res[1],res)
  return(res)
}


display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  if(index/totalN==1){
    cat("\n")
  }

}
