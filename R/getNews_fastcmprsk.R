getNews_fastcmprsk <- function(...) {
  newsfile <- file.path(system.file(package = "fastcmprsk"), "NEWS")
  file.show(newsfile)
}
