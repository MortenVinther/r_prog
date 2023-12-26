library(knitr)

setwd(file.path('C:','_C_drev','SMS-git','Documentation'))
pandoc('SMS_description_2023.xml',format='Rmd')


pnd<-file.path("C:","Program Files","RStudio","resources","app","bin","quarto","bin","tools","pandoc.exe")
file.exists(pnd)

system(paste0(pnd, "-h")) 

system(paste0(pnd, "-f docx -t markdown SMS_description_2023.docx"))

       
