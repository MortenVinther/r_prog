runs<-NULL
my.ps<-c(1L,2L)

pn=1L
source(file.path(rsms.root.prog,'batchInclude.R'))


pn=2L
source(file.path(rsms.root.prog,'batchInclude.R'))
pn=3L
source(file.path(rsms.root.prog,'batchInclude.R'))

aa<-lapply(runs,function(x) x)
cleanup()
plotCompareRunSummary(Type=c("compSummaryConf","compSummary")[2],showSpecies=1:12,
                      inpRdata=aa,
                      labels=aa,
                      outFormat=c('screen','pdf','png')[1],
                      longSpNames=FALSE, fileLabel='SS_')
