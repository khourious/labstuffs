# https://lpmor22.github.io/ | 2025-06-04

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", dependencies = TRUE)
library("pacman")

if (!requireNamespace("read.gb", quietly = TRUE))
  devtools::install_github("lpmor22/read.gb")

p_load(countrycode, parallel, read.gb, seqinr, this.path, writexl)

setwd(dirname(this.path()))

input <- "Zika_GenBank_after_2019-11-30.gb"

n_threads <- detectCores() - 1
cluster <- makeCluster(n_threads)
clusterEvalQ(cluster, {library(read.gb)})

gb <- readLines(input)
gb_start <- grep("^LOCUS", gb)
gb_end <- grep("^//$", gb)
gb_temp <- mapply(function(x,y) gb[x:y], gb_start, gb_end)

gb_rec <- gb_temp[sapply(gb_temp, function(rec) {
  locus_line <- rec[1]
  !grepl(" rc", locus_line) || as.numeric(gsub(".* (\\d+) rc.*", "\\1", locus_line)) > 2
})]

start_time <- Sys.time()
read_gb <- unlist(parLapply(cluster, gb_rec, function(part) {
  tmpfile <- tempfile()
  writeLines(part, tmpfile)
  gb <- read.gb(tmpfile, DNA = TRUE, Type = "full", Source = "File")
  file.remove(tmpfile)
  return(gb)}), recursive = FALSE)
end_time <- Sys.time()
end_time - start_time
stopCluster(cluster)

extract_gb <- function(x) {
  locus <- strsplit(x[[1]], " ")[[1]]
  data.frame(accession = locus[1], length = locus[2])}
df <- do.call(rbind, lapply(extract.gb(read_gb, "LOCUS"), extract_gb))

extract_gb <- function(x) {
  location <- x$source$Location
  qualifier <- x$source$Qualifier
  var <- c("organism",
           "isolate",
           "isolation_source",
           "host",
           "country",
           "collection_date",
           "note",
           "serotype")
  data <- sapply(var, function(x) {
    if (x %in% location) {
      if (x == "host") {
        host_raw <- qualifier[which(location == x)]
        host_edited <- sub(";.*", "", host_raw)
        return(host_edited)}
      else if (x == "collection_date") {
        date_raw <- qualifier[which(location == x)]
        date_edited <- gsub(",.*", "", date_raw)
        date_edited <- gsub("\\s+", "", date_edited)
        return(date_edited)}
      else {
        return(qualifier[which(location == x)])}}
    else {
      return(NA)}})
  return(data.frame(t(data)))}
temp <- do.call(rbind, lapply(extract.gb(read_gb, "source"), extract_gb))
df <- cbind(df, temp)

extract_gb <- function(x) {
  return(data.frame(fasta = paste(x$ORIGIN, collapse = '')))}
temp <- do.call(rbind, lapply(extract.gb(read_gb, "ORIGIN"), extract_gb))
df <- cbind(df, temp)
rownames(df) <- NULL

names(df)[which(names(df) == "country")] <- "country_original"
df$country <- sub(":.*", "", df$country_original)
df$country_code <- countrycode(df$country, origin = "country.name", destination = "iso3c", warn = FALSE)

names(df)[which(names(df) == "collection_date")] <- "collection_date_original"
date_edited <- function(date) {
  if (grepl("^\\d{2}-\\w{3}-\\d{4}$", date)) {
    return(format(as.Date(date, "%d-%b-%Y"), "%Y-%m-%d"))}
  else if (grepl("^\\w{3}-\\d{4}$", date)) {
    return(paste(format(as.Date(paste("01", date, sep = "-"), "%d-%b-%Y"), "%Y-%m"), "XX", sep = "-"))}
  else if (grepl("^\\d{4}$", date)) {
    return(paste(date, "XX", "XX", sep = "-"))}
  else {
    return(date)}}
df$collection_date <- sapply(df$collection_date_original, date_edited)

names(df)[names(df) == "serotype"] <- "genotype"
get_genotype <- function(note) {
  if (grepl("genotype:", note)) {
    genotype <- strsplit(note, ": ")[[1]][2]
    genotype <- sub(",.*", "", genotype)
    return(genotype)}
  else {
    return(NA)}}
df$genotype <- ifelse(is.na(df$genotype), sapply(df$note, get_genotype), df$genotype)
df <- df[, !names(df) %in% "note"]
df <- df[, c("accession",
             "length",
             "isolate",
             "isolation_source",
             "host",
             "organism",
             "genotype",
             "country_original",
             "country",
             "country_code",
             "collection_date_original",
             "collection_date",
             "fasta")]
df$isolate <- sapply(df$isolate, paste, collapse = "")
df$isolation_source <- sapply(df$isolation_source, paste, collapse = "")
df$host <- sapply(df$host, paste, collapse = "")
df$organism <- sapply(df$organism, paste, collapse = "")
df$genotype <- sapply(df$genotype, paste, collapse = "")
df$country_original <- sapply(df$country_original, paste, collapse = "")
df$collection_date_original <- sapply(df$collection_date_original, paste, collapse = "")

write_xlsx(df, sub("\\.gb$", ".xlsx", input))
fasta_sequence <- as.list(as.character(df$fasta))
fasta_header <- as.character(df$accession)

write.fasta(sequences = fasta_sequence, names = fasta_header, file.out = sub("\\.gb$", ".fasta", input))