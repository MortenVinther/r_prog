# make default rsms.dat (control) file
batch_default_configuration <-function(outfile='rsms.dat',dir=data.path,writeConrol=TRUE) {
  
  a<-RSMS.control(first.year=1974,last.year=2022,no.species=27,last.season=4,max.age.all=10,
                  no.other.predators=15,no.VPA.predators=5,
                  species.names=c("FUL", "GLT", "HEG", "KTW", "GBG", "GNT", "PUF", "RAZ", "RAJ", "GUR", "WHM", "NHM", "GSE", "HBP", "HKE", "COD", "WHG", "HAD", "POK", "MAC", "HER", "NSA", "SSA", "NOP", "SPR", "PLE", "SOL"),
                  species.names.long=c("Fulmar","Guillemot","Her.Gull","Kittiwake","GBB.Gull","Gannet","Puffin","Razorbill","A.radiata","G.gurnards","W.horse.mac","N.horse.mac","Grey.seal","H.porpoise","Hake","Cod","Whiting","Haddock","Saithe","Mackerel","Herring","N.sandeel","S.sandeel","Nor.pout","Sprat","Plaice","Sole")
     )
  
  nonPrey<-c('POK','MAC','PLE','SOL')
  a@species.info[nonPrey,'prey']<-0L
  
  firstAgeF1<-c('COD','MAC','SSA','SPR','PLE','SOL')
  a@species.info[firstAgeF1,'first-age F>0']<-1L
  firstAgeF3<-('POK')
  a@species.info[firstAgeF3,'first-age F>0']<-3L
  a@species.info[,'last-age'] <-as.integer(
    # FUL GLT HEG KTW GBG GNT PUF RAZ RAJ GUR WHM NHM GSE HBP HKE COD WHG HAD POK MAC HER NSA SSA,NOP SPR PLE SOL
    c(1,    1,  1,  1,  1,  1,  1,  1,  3,  4,  3,  6,  1,  1,  9, 10,  6, 10, 10, 10,  8,  4,  4,  3,  3, 10, 10))
  
  no.other.predators<-sum(a@species.info[,'predator']==2)
  
  
  bySpAge<-list(
    c(1,2,4), #COD   
    c(0,1,3),   #WHG   
    c(0,1,2), #HAD   
    c(3,4), #POK   
    c(1,2), #MAC   
    c(0,1), #HER   
    c(0,1), #NSA
    c(0,1), #SSA   
    c(0,1), #NOP
    c(0,1), #SPR
    c(0,1), #PLE
    c(0,1)  #SOL 
  )
  a@catch.s2.group<-bySpAge
  
   
                 #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@SSB.R<-as.integer(c(0,    0,    0,    1,    0,    1,    1,    1,    2,    2,    0,    0) )
  
  
  bySpAge<-list(
    c(0,1), #COD   
    c(0,1),   #WHG   
    c(0,1), #HAD   
    c(0,3), #POK   
    c(0,1), #MAC   
    c(0,1), #HER   
    c(0,1), #NSA
    c(0,1), #SSA   
    c(0,1), #NOP
    c(0,1), #SPR
    c(0,1), #PLE
    c(0,1)  #SOL 
  )
  a@keyVarLogN<-bySpAge
  
  bySpAge<-list(
    c(1,2,3), #COD   
    c(0,1,2),     #WHG   
    c(0,1,2),   #HAD   
    c(3,4,5),   #POK   
    c(1),   #MAC   
    c(0,1,2,4), #HER   
    c(0,1),   #NSA
    c(1,2),   #SSA   
    c(0,1,2), #NOP
    c(1),     #SPR
    c(1,2,3), #PLE
    c(1,2,3)  #SOL 
  );
  a@keyLogFsta<-bySpAge
  
  
  #                      COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@use.rho<-as.integer(c(1,    1,    0,    0,    0,    1,    0,    0,    0,    0,    1,    1) )
  
  avgF<-matrix(as.integer(c(
    ## first and last age in calculation of average F by species (option avg.F.ages)
  2, 4,  # COD 
  2, 5,  # WHG 
  2, 4,  # HAD 
  4, 7,  # POK 
  4, 8,  # MAC 
  2, 6,  # HER 
  1, 2,  # NSA 
  1, 2,  # SSA 
  1, 2,  # NOP 
  1, 2,  # SPR 
  2, 6,  # PLE 
  2, 6  # SOL 
  )),ncol=2, byrow = TRUE)
  dimnames(avgF)<-dimnames(a@avg.F.ages)
  a@avg.F.ages<-avgF
  
  if (writeConrol) write.RSMS.control(a,file=file.path(data.path,outfile))
  invisible(a)
}
