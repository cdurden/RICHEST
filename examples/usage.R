library(richest)
cps = readCPS(file="data/sargassoSeaOTU_CPS_unique.txt")
est = richest(cps,seq(690,10000,25),"Lijou")
plot(slot(est,"sampleSizes"),slot(est,"estimates"))
est = richest(cps,seq(690,10000,25),"PNPMLE")
points(slot(est,"sampleSizes"),slot(est,"estimates"))
est = richest(cps,seq(690,1400,10),"GT")
points(slot(est,"sampleSizes"),slot(est,"estimates"))

spc = as(cps,"richestSPC")
est = richest(spc,seq(690,10000,25),"Lijou")
plot(slot(est,"sampleSizes"),slot(est,"estimates"))
est = richest(spc,seq(690,10000,25),"PNPMLE")
points(slot(est,"sampleSizes"),slot(est,"estimates"))
est = richest(spc,seq(690,1400,10),"GT")
points(slot(est,"sampleSizes"),slot(est,"estimates"))

