.First.lib <-
    function(libname, pkgname, where)
    library.dynam("beadarray", pkgname, libname)

.Last.lib <-
    function(libpath)
    dyn.unload(file.path(libpath,
                         "libs",
                         paste("beadarray",
                               .Platform$"dynlib.ext",
                               sep = "")))

