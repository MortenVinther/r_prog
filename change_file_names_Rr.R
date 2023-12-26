

ff<-dir(prog.path,pattern="*.R",recursive=T,full.names=TRUE)

for (f in ff) {
 #f<-"/home/morten/SMS/r_prog/makeAllGraphs_NorthSea.R"
  a<-readLines(con=f)
  pat1<-'\\.r")'
  pat2<-"\\.r')"
  pat3<-"SMS\\.dat"
  ll1<-grep(pat1 ,a)
  ll2<-grep(pat2,a)
  ll3<-grep(pat3,a)
  ll<-sort(unique(c(ll1,ll2,ll3)))
  if (length(ll)> 0) {
    cat(f,ll,'\n')
    print(a[ll])
  }

}
