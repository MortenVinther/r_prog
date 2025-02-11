#    1= include for for both recruits and older ages
#    2= Include for recruits only (N at age older than recruits in first the year become parameters)
#    3= Include for older ages only (recruits become parameters)
#    4= Do not include process noise at all (N in first year and  recruits become parameters)
inp_all$data$inclProcessNoiseinp_all$data$logNfirstYparamfromTo
inp_all$parameters$logNfirstYparam;length(inp_all$parameters$logNfirstYparam)

inp_all$data$inclProcessNoise
inp_all$data$logNrecruitParamfromTo
inp_all$parameters$ogNrecruitParam; length(inp_all$parameters$ogNrecruitParam)

#####

inp$data$inclProcessNoise
inp$data$logNfirstYparamfromTo
inp$parameters$logNfirstYparam;length(inp$parameters$logNfirstYparam)

inp$data$inclProcessNoise
inp$data$logNrecruitParamfromTo
inp$parameters$ogNrecruitParam;  length(inp$parameters$ogNrecruitParam)
