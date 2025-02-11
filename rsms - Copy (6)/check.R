## Process noise: Include process noise (option incl.process.noise)
#    1= Both recruits and older ages have process noise estimated
#    2= Recruits only (N at age older than recruits in the first year become parameters)
#    3= Do not include process noise at all (N in the first year and recruits become parameters)


inp_all$data$doProcessNoise
inp_all$data$doProcessN_any
inp_all$data$doProcessN_old
inp_all$data$doProcessN_none

inp_all$data$doProcessNoise
inp_all$data$logNfirstYparamfromTo
inp_all$parameters$logNfirstYparam;length(inp_all$parameters$logNfirstYparam)

inp_all$data$doProcessNoise
inp_all$data$logNrecruitParamfromTo
inp_all$parameters$ogNrecruitParam; length(inp_all$parameters$logNrecruitParam)


#####
inp$data$doProcessNoise
inp$data$doProcessN_any
inp$data$doProcessN_old
inp$data$doProcessN_none

inp$data$doProcessNoise
inp$data$logNfirstYparamfromTo
inp$parameters$logNfirstYparam;length(inp$parameters$logNfirstYparam)

inp$data$doProcessNoise
inp$data$logNrecruitParamfromTo
inp$parameters$ogNrecruitParam; length(inp$parameters$logNrecruitParam)

