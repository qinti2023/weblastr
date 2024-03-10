library(httr)
library(stringr)
library(progress)
library(tidyverse)
library(data.table)
library(XML)
library(RCurl)
library(rentrez)


#' Countdown seconds
#'
#' @param seconds int
#'
#' @return a letter on console
#' @export
#'
#' @examples countdown(30)
countdown <- function(seconds) {
  for (i in seconds:0) {
    countdown_str <- sprintf("\rRetry after: %02d seconds", i)
    cat(countdown_str)
    flush.console()
    Sys.sleep(1)
  }
  cat("\n")
}

#' get blast result from web
#'
#' @param seq strings
#' @param program strings
#' @param database strings
#' @param org strings
#'
#' @return "data frame"
#' @export
#'
#' @examples get_blast_result("MLRACQLSGVTAAAQSCLCGKFVLRPLRPCRRYSTSGSSGLTTGKIAGAGLLFVGGGIGGTILYAKWDSHFRESVEKTIPYSDKLFEMVLGPAAY","tblastn","nr","Homo sapeins")
get_blast_result= function(seq,program,database,org){
  cat("###### post the blast task ######\n")
  params <- list(
    CMD = "Put",
    PROGRAM = program,
    DATABASE = database,
    QUERY = seq,
    Organism = org
  )

response <- POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi", body = params, encode = "form")
response_content <- content(response, "text",encoding = "UTF-8")
rid <- str_extract(response_content, "(?<=RID = )\\w+")
url <- paste0("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=", rid)
found_result = F
while (!found_result) {
  response <- GET(url)
  content <- content(response, "text",encoding = "UTF-8")
  if(grepl("Sequences producing significant alignments:",content)){
    found_result <- TRUE
    cat("Success\n")
  } else if(grepl("error",content)){
    found_result <- TRUE
    cat("Not Found Results\n")
  }else {
    cat("Did not find 'result'. Waiting for 1 minutes before retrying...\n")
    countdown(60)
    }
  }
pattern <- "(?s)Sequences producing significant alignments:(.*?)ALIGNMENTS"
matches <- regmatches(content, regexec(pattern, content, perl=TRUE))
hit_string = matches[[1]][1]
if(!is.na(hit_string)){
hit_df = data.frame(raw_string = unlist(strsplit(hit_string,"\n")), stringsAsFactors = FALSE)
hit_df =hit_df[3:(nrow(hit_df)-2),,drop=F]

cat("###### get ncbi summary ######\n")
get_summary = function(string){
    spilt_string = unlist(strsplit(string,"\\s"))
    string_len = length(spilt_string)
    Ref_ID = spilt_string[1]
    Ident = spilt_string[string_len-5]
    E_value = spilt_string[string_len-7]
    if(Ident == ""){
      Ident = spilt_string[string_len-6]
      E_value = spilt_string[string_len-8]
    }
    result_list = list(Ref_ID = Ref_ID,
                       Ident = Ident,
                       E_value = E_value)
    return(result_list)
}
get_summary = Vectorize(get_summary)

hit_final = get_summary(hit_df$raw_string)|>
    t()|>
    as.data.table()|>
    mutate(NCBI_url = paste0("https://www.ncbi.nlm.nih.gov/nuccore/",Ref_ID))
}else{
  hit_final = NULL
}}


get_des <- function(url, url_all, len){
    pos <- match(url, url_all)
    cat(sprintf("%d/%d running get_des at: %s \n", pos, len, url))
    try_once = function(url) {
      tryCatch({
        getURL(url) %>%
          htmlParse(asText = TRUE) %>%
          getNodeSet('//*[@id="maincontent"]/div/div[5]/div[1]/h1') %>%
          sapply(xmlValue) %>%
          keep(~ .x != "") %>%
          str_replace_all("(\\n )+","")
      }, error=function(e) {
        cat("Error in URL:", url, "\n", e$message, "\n")
        return(NULL)
      })
    }

    result = try_once(url)

    if (is.null(result)) {
      cat("Retrying after 2 seconds...\n")
      Sys.sleep(2)
      result = try_once(url)
    }

    if (is.null(result)) {
      return(NA)
    } else {
      return(result)

    }
  }

#' Mutate a column about description base on Genebank ID
#'
#' @param df "data frome"
#'
#' @return "data frome"
#' @export
#'
#' @examples get_cds_df(data_frame)
get_des_df <- function(df) {
    df %>%
      mutate(Description= map(df$NCBI_url, ~possibly(get_des, NA)(.x, df$NCBI_url, n())))
  }

get_cds <- function(id, id_all, len){
    pos <- match(id, id_all)
    cat(sprintf("%d/%d running get_cds at: %s \n", pos, len, id))
    ncbi_gb <- tryCatch({
      entrez_fetch(db = "nuccore", id = id, rettype = "gb", retmode = "text")
    }, error=function(e) {
      cat("Error in ID:", id, "\n", e$message, "\n")
      return(NA)
    })
    if (is.na(ncbi_gb)) return(NA)
    pattern_cds <- "(?s)/translation=(.*?)ORIGIN"
    matches_cds <- regmatches(ncbi_gb, regexec(pattern_cds, ncbi_gb, perl=TRUE))
    if (length(matches_cds) < 1 || length(matches_cds[[1]]) < 2) {
      cat("NOT Found CDS :", id, "\n")
      return(NA)
    }

    cds_string = matches_cds[[1]][2]
    cds = cds_string %>%
      unlist() %>%
      str_replace_all("[^A-Z]", "")
    return(cds)
  }
#' Mutate a column about aa sequence of CDS base on Genebank ID
#'
#' @param df "data frame"
#'
#' @return "data frame"
#' @export
#'
#' @examples get_cds_df(data_frame)
get_cds_df <- function(df) {
    df %>%
      mutate(CDS_aa = map(df$Ref_ID, ~possibly(get_cds, NA)(.x, df$Ref_ID, n())))
  }


#' get all results from blast.ncbi.nlm.nih.gov
#'
#' @param seq strings
#' @param program strings
#' @param database strings
#' @param org strings
#'
#' @return "data frame"
#' @export
#'
#' @examples test= rblast_web("MLRACQLSGVTAAAQSCLCGKFVLRPLRPCRRYSTSGSSGLTTGKIAGAGLLFVGGGIGGTILYAKWDSHFRESVEKTIPYSDKLFEMVLGPAAYNVPLPKKMMMMMMMMMMMMMMMMMMMMMMMMMMM","tblastn","nr","Homo sapeins")
rblast_web = function(seq,program,database,org){
  hit_final = get_blast_result(seq,program,database,org)
  if(!is.null(hit_final)){
    cat("###### get description and cds_aa ######\n")
    hit_final = hit_final |>
    get_des_df() |>
    get_cds_df()
    cat("done")
  }
  return(hit_final)
}







