# Dependencies
library(dplyr)
library(rentrez)
library(stringr)
library(ggplot2)


## Clean up the n-end characters from a string 
clean_end <- function(str, x){
  
  ### str, string vector
  ### x, number of characters
  
  str = as.character(str)
  str = substr(str, start = 1, stop = nchar(str)-x)
  str
}


## Extract the end n-chars from a string
get_end <- function(str, x) {
  
  ### str, string vector
  ### x, number of characters
  
  text = as.character(str)
  text = substr(text, start = nchar(text)-x+1, stop = nchar(text))
  text
  
}

## Keep 1 hit per match. Otherwise we would be counting twice, 
# a single hit with two matches in the seq.
one_match <- function(df) {
  
  ## Create new df with number of sequence counts  
  ### df, dataset
  
  acc = df[[1]]
  
  # remove last 2 characters from query.acc.
  # new column with the tidy acc. 
  df$query = clean_end(acc,2)
  
  # index for query.acc with ".2"
  index = which(get_end(df[[1]],2) == ".2")
  
  # only if there are ".2" sequences
  if(length(index) > 0) {
    # remove rows with ".2"
    df = df[-index,]
  }
  
  ## Clean up everything beyond '.' in a string 
  df$query = sub("\\..*", "", df$query)
  
  ## Select columns
  df <- df %>% 
    select(query, reference.acc, score, identity, query.start:reference.end, 
           query.strand:query.length) 
  
  df
}  

## Filter by score
by_score <- function(df, score) {
  df = df[df$score >= score, ]
}

## Filter by identity
by_identity <- function(df, identity) {
  df = df[df$identity >= identity, ]
}


## Create df with number of sequence counts 
getCounts <- function(df) {
  
  ### df, dataset
  
  # rename col names 
  df = df %>% 
    select(query, reference.acc) %>% 
    rename("Query" = "query", "Reference" = "reference.acc")
  
  
  # select columns
  df = df %>% 
    group_by(Query, Reference) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(Count)) %>%
    ungroup()
  
  
  df
  
  
  }


## 'getSize' takes a vector containing NCBI SRA accessions and returns and vector 
## with the run size (Mb) 

getSize <- function(ids) {
  
  sizes = vector(mode = "numeric")
  
  for(id in ids) {
    term =  paste0(id, "[ACCN]")
    run = entrez_search(db = "sra", term = term )
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    
    # biological reps. are shown between ";"
    bioreps = str_split(x, ";")
    # loop over bio reps to find the right run 
    for(i in seq_along(1:length(bioreps[[1]]))) {
      if(str_detect(bioreps[[1]][i], id)) {
        # extract range defined by "total bases" and "load_done"
        size = str_sub(x, start =str_locate(x, "total_bases=")[2]+2, 
                       end =str_locate(x, "load_done")[1]-3)
      }
    }
    
    size_mega = as.numeric(size)/1e6
    sizes <- c(sizes, size_mega)
    
  }
  sizes
  
  }

## Usage 
# ids <- c("DRR071071", "ERR1912953")
# getSize(ids)


## `gatherSize` takes a dataset with SRA accessions as first columns and returns 
## a vector containing the the run size (Mb) for each accession

gatherSize <-  function(df) {
  
  ## df, contains the SRA Query in the first column
  
  # get unique SRA queries from the Query column
  accessions <- unique(df[[1]])   
  
  # getSizes (Mb) per each SRA query
  runsize <- getSize(accessions)         
  
  # create a vector match for each query
  mindex = match(df[[1]], accessions)
  
  # create a vector with size for each SRA query (= row)
  squery = c()
  for(m in mindex) {
    squery = c(squery, runsize[m])
  }
  
  # return the vector
  squery
}

