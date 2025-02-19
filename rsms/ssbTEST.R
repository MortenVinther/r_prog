
SSB_R<-function(s,y,a=1) {
  
  switch( as.character(stockRecruitmentModelCode[s]), 
    
    "1"=  rec<-rec_loga[s]+log(ssb[s,y-recAge])-exp(rec_logb[s])*ssb[s,y-recAge],
    "2"=  rec<-rec_loga[s]+log(ssb[s,y-recAge])-log(1+exp(rec_logb[s])*ssb[s,y-recAge]),
    "3"=  rec<-rec_loga[s],
    stop(paste0("SR model code ",stockRecruitmentModelCode[s]," not recognized"))
  )
}
 

#   case 61: // Hockey stick
#   // Type log_level = rec_pars(0);
#   // Type log_blim = rec_pars(1);
#   // Type a = thisSSB - exp(log_blim);
#   // Type b = 0.0;
#   // Type cut = 0.5 * (a+b+CppAD::abs(a-b)); // max(a,b)
#   predN = rec_pars(0) - rec_pars(1) +
     log(thisSSB - (0.5 * ((thisSSB - exp(rec_pars(1)))+Type(0.0)+CppAD::abs((thisSSB - exp(rec_pars(1)))-Type(0.0)))));
#   break;
# 

rec_loga<-log(1000)
rec_logb<-log(1000000) # log_blim
SSB<-exp(rec_logb)
SSB<-seq(0.2,2,0.1)*SSB
s=1

rec<-  rec_loga[s] - rec_logb[s] + log(SSB - (0.5 * (SSB - exp(rec_logb[s])+abs(SSB - exp(rec_logb[s])))))

plot(x=SSB,y=exp(rec))

for (i %in% seq(0.2,2,0.1))
log(rec)
stockRecruitmentModelCode<-c(1,2,3,9)


rec_loga            4.29          0.249   -0.143       73.2    520     1 COD                            1      -9      -9
521 rec_loga            4.19          0.203    0.237       65.8    521     2 HER                            2      -9      -9
522 rec_loga            6.55          0.117    0.0230      -9      522     3 NSA                            3      -9      -9
523 rec_loga            6.43          0.0549   0.0456      -9      523     4 SSA                            4      -9      -9
524 rec_logb          -11.6           0.637    0.123        0      524     1 COD                            1      -9      -9
525 rec_logb          -13.3           0.306   -0.134        0      525     2 HER                            2      -9      -9
526 rec_logb          -19.2           5.21     0.0368       0      526     3 NSA                            3      -9      -9
527 rec_logb          -18.3           3.94   



    "d"= cat("Subtraction =", val1 - val2),  
    "r"= cat("Division = ", val1 / val2),  
    "s"= cat("Multiplication =", val1 * val2),
    "m"= cat("Modulus =", val1 %% val2),
    "p"= cat("Power =", val1 ^ val2)
  )  
  
  
  if (stockRecruitmentModelCode[s]==0 | !recruitYears[s,y]){    ## straight RW
    rec = logN[[s]][a, y-1]
  } else {
    if (stockRecruitmentModelCode[s]==1){ ## Ricker
      rec<-rec_loga[s]+log(ssb[s,y-recAge])-exp(rec_logb[s])*ssb[s,y-recAge]
    } else {
      if(stockRecruitmentModelCode[s]==2){  ## B&H
        rec<-rec_loga[s]+log(ssb[s,y-recAge])-log(1+exp(rec_logb[s])*ssb[s,y-recAge])
      } else {
        if(stockRecruitmentModelCode[s]==3){  ## GM
          rec<-rec_loga[s]
        } else {
          stop(paste0("SR model code ",stockRecruitmentModelCode[s]," not recognized"))
        }
      }
    }
  }
  return(rec)
}

result = switch(  
  val3,  
  "a"= cat("Addition =", val1 + val2),  
  "d"= cat("Subtraction =", val1 - val2),  
  "r"= cat("Division = ", val1 / val2),  
  "s"= cat("Multiplication =", val1 * val2),
  "m"= cat("Modulus =", val1 %% val2),
  "p"= cat("Power =", val1 ^ val2)
)  

