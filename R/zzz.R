.onAttach <- function(libname, pkgname) {

  mylib <- dirname(system.file(package = "AFR"))
  ver <- packageDescription("AFR", lib.loc = mylib)["Version"]
  txt <- c("\n",
           paste(sQuote("AFR"), "version:", ver),
           "\n",
           paste(sQuote("AFR"),
                 "is a package for banking sector analysis",
                 "and easier interpretation of statistical functions."),
           "\n",
           paste("See",
                 sQuote("library(help=\"AFR\")"),
                 "for details."),
           "\n",
            paste("Please fill in the form", sQuote("https://forms.gle/anRQ8WZbTQ5hQetg9"),
                  "It is necessary for the analysis of the users of our package.", "\n"))

  if(interactive() || getOption("verbose"))
  packageStartupMessage(paste(strwrap(txt,
                                        indent = 4,
                                        exdent = 4),
                                collapse = "\n"))


}

