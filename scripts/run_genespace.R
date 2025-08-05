library(GENESPACE)

gpar = init_genespace(
         wd = paste0(getwd(),"/gs"),
         nCores = 16,
         path2orthofinder = "/usr/local/bin/orthofinder",
         path2diamond =  "/usr/local/bin/diamond",
         path2mcscanx = "/usr/local/bin/")
out = run_genespace(gsParam = gpar)

