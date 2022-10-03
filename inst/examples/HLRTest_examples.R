

y84 <- as.matrix(hamilton84GNP$GNP_logdiff)


Y <- y84


hlrt_control  <- list(ix = 1, 
                      iv = 0, 
                      iq = 0,
                      gn = 20, 
                      pn = 17,
                      pn_step_size = 0.05,
                      nwband = 4, 
                      reps = 1000)

st <- Sys.time()
hlrt <- HLRTest(Y, p = 4, control = hlrt_control)
end <- Sys.time() - st

end
hlrt$`Test-Stat`
hlrt$params
hlrt$`p-value`