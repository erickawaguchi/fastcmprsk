.onAttach <- function(libname, pkgname) {
  pkg.version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                            fields = "Version")
  packageStartupMessage(paste("\n",
                              paste(pkgname,
                              pkg.version,
                              "is installed!"),
                              "\n",
                              "\n",
                              "Type getNews_fastcmprsk() to see changes and updates to the package.",
                              "\n",
                              "Go to https://github.com/erickawaguchi/fastcmprsk/ for latest updates!",
                              "\n",
                              "\n"))
}
