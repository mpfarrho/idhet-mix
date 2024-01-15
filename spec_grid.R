grid.hor <- 60
grid.mix <- c(TRUE)
grid.pub <- c(TRUE)
grid.J <- c(30)
grid.Q <- c(2)
grid.sv <- c(TRUE)
grid.svm <- c(FALSE)
grid.svf <- c(TRUE)
grid.lev <- c(TRUE)
grid.exo <- c(TRUE)

grid.full <- expand.grid("mix"=grid.mix,
                         "pub"=grid.pub,
                         "J"=grid.J,
                         "Q"=grid.Q,
                         "h"=0:grid.hor,
                         "sv"=grid.sv,
                         "svm"=grid.svm,
                         "svf"=grid.svf,
                         "lev"=grid.lev,
                         "exo"=grid.exo)
grid.full <- grid.full[!(grid.full$sv==FALSE & grid.full$svm!=FALSE),]
grid.full <- grid.full[!(grid.full$sv==FALSE & grid.full$lev!=FALSE),]
grid.full <- grid.full[!(grid.full$mix==FALSE & grid.full$J!=30),]

# further subset
grid.full <- grid.full[
            ((grid.full$mix == TRUE  & grid.full$pub == TRUE  & grid.full$J == 30 & grid.full$Q == 2 &  grid.full$sv == TRUE  & grid.full$svm == FALSE & grid.full$svf == TRUE  & grid.full$lev == TRUE  & grid.full$exo == TRUE)  +     # baseline 
             (grid.full$mix == TRUE  & grid.full$pub == TRUE  & grid.full$J == 30 & grid.full$Q == 2 &  grid.full$sv == FALSE & grid.full$svm == FALSE & grid.full$svf == TRUE  & grid.full$lev == FALSE  & grid.full$exo == TRUE)  +    # R1
             (grid.full$mix == TRUE  & grid.full$pub == TRUE  & grid.full$J == 30 & grid.full$Q == 2 &  grid.full$sv == FALSE & grid.full$svm == FALSE & grid.full$svf == FALSE & grid.full$lev == FALSE  & grid.full$exo == TRUE)  +    # R2
             (grid.full$mix == TRUE  & grid.full$pub == TRUE  & grid.full$J == 30 & grid.full$Q == 3 &  grid.full$sv == TRUE  & grid.full$svm == FALSE & grid.full$svf == TRUE  & grid.full$lev == TRUE  & grid.full$exo == TRUE)  +     # R3
             (grid.full$mix == TRUE  & grid.full$pub == TRUE  & grid.full$J == 30 & grid.full$Q == 2 &  grid.full$sv == TRUE  & grid.full$svm == FALSE & grid.full$svf == TRUE  & grid.full$lev == TRUE  & grid.full$exo == FALSE) +     # R4
             (grid.full$mix == TRUE  & grid.full$pub == TRUE  & grid.full$J == 10 & grid.full$Q == 2 &  grid.full$sv == TRUE  & grid.full$svm == FALSE & grid.full$svf == TRUE  & grid.full$lev == TRUE  & grid.full$exo == TRUE)  +     # varying J
             (grid.full$mix == FALSE & grid.full$pub == TRUE  & grid.full$J == 30 & grid.full$Q == 2 &  grid.full$sv == TRUE  & grid.full$svm == FALSE & grid.full$svf == TRUE  & grid.full$lev == TRUE  & grid.full$exo == TRUE)) > 0,] # no mix

rownames(grid.full) <- 1:nrow(grid.full)

# algorithm setup
p <- 1 # lags
nsave <- 9000
nburn <- 3000
thinfac <- 3

# credible sets
conf <- c(0.05,0.16,0.50,0.84,0.95)

# excluded variables
sl.excl <- c("gold_usd_d","oil_brent_d",
             "bofaml_us_hyld_oas_d","bofaml_us_aa_oas_d","bofaml_us_bbb_oas_d",
             "gs1_d","gs10_d","usbkeven5y_d","vixcls_d","nasdaqcom_d",
             "eurinflswap1y_d","eurinflswap2y_d","eurinflswap3y_d","eurinflswap4y_d","eurinflswap5y_d","eurinflswap10y_d")


