.onLoad <- function(lib, pkg)
{
  library.dynam("richest", pkg, lib)
  .Call('init_pari',PACKAGE="richest")
}
.Last.lib <- function(lib, pkg)
{
  library.dynam.unload("richest", pkg, lib)
  .Call('close_pari',PACKAGE="richest")
}
