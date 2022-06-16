source("r/header.R")

# new way to get off web:
uri <- "https://api.zotero.org/users/5814932/collections/MF4CQ9M4/items"
.parms <- list(
  key = "ED0yHWd4Qtd6ereWuiLgYN2p", 
  collection = "MF4CQ9M4",
  limit = 1000,
  format = "bibtex",
  uri = "https://api.zotero.org/users/5814932/collections/MF4CQ9M4/items"
)
res <- httr::GET(uri, query = .parms)
res <- httr::content(res, as = "text", encoding = "UTF-8")
write(res, file = "ms/stomata-independence.bib", append = FALSE)
