######################################################
######################################################
# 
######################################################
######################################################

# Edited from https://github.com/jugne/sRanges-material/blob/main/validation/estimate_all_params_ext_dna/validationPlots.R 
# Inefficient code to plot results for four intervals
# If applying to a greater number of intervals, it would be worth writing
# some functions to reduce repeated code.
library(ggplot2)
library(coda)
options(digits=7)
library(dplyr)

# clear workspace
rm(list = ls())

test_vals <- function(prob, true=0, estimated=c()){
    lower <- HPDinterval(as.mcmc(estimated), prob=prob)[1]
    upper <- HPDinterval(as.mcmc(estimated), prob=prob)[2]
    test <- c(as.numeric(true >= lower & true <= upper))
  return(test)
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Set the directory to the directory of the file
wd<-"/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/graphing"
setwd(wd)
dir.create("figures")

# read in the true rates
true_rates_file <- "/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/true_rates.csv"
true.rates <- read.table(true_rates_file, header=TRUE, sep=",")

n<-200
remove_rows<-c()

ints = 4

diversificationRate_1 = data.frame(true=numeric(n), estimated=numeric(n),
                                 upper=numeric(n), lower=numeric(n),
                                 rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                 rel.err.mode=numeric(n), cv=numeric(n),
                                 hpd.lower=numeric(n), hpd.upper=numeric(n),
                                 hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                 test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                 test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                 test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                 test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                 test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                 test_90=numeric(n), test_95=numeric(n),
                                 test_100=numeric(n))

diversificationRate_2 = data.frame(true=numeric(n), estimated=numeric(n),
                                 upper=numeric(n), lower=numeric(n),
                                 rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                 rel.err.mode=numeric(n), cv=numeric(n),
                                 hpd.lower=numeric(n), hpd.upper=numeric(n),
                                 hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                 test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                 test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                 test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                 test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                 test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                 test_90=numeric(n), test_95=numeric(n),
                                 test_100=numeric(n))

 diversificationRate_3 = data.frame(true=numeric(n), estimated=numeric(n),
                                 upper=numeric(n), lower=numeric(n),
                                 rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                 rel.err.mode=numeric(n), cv=numeric(n),
                                 hpd.lower=numeric(n), hpd.upper=numeric(n),
                                 hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                 test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                 test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                 test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                 test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                 test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                 test_90=numeric(n), test_95=numeric(n),
                                 test_100=numeric(n))

diversificationRate_4 = data.frame(true=numeric(n), estimated=numeric(n),
                                 upper=numeric(n), lower=numeric(n),
                                 rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                 rel.err.mode=numeric(n), cv=numeric(n),
                                 hpd.lower=numeric(n), hpd.upper=numeric(n),
                                 hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                 test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                 test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                 test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                 test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                 test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                 test_90=numeric(n), test_95=numeric(n),
                                 test_100=numeric(n))                                

turnover_1 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                      lower=numeric(n), rel.err.median=numeric(n),
                      rel.err.mean=numeric(n), rel.err.mode=numeric(n),
                      cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n),
                      hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                      test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                      test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                      test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                      test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                      test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                      test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

turnover_2 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                      lower=numeric(n), rel.err.median=numeric(n),
                      rel.err.mean=numeric(n), rel.err.mode=numeric(n),
                      cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n),
                      hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                      test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                      test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                      test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                      test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                      test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                      test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

turnover_3 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                      lower=numeric(n), rel.err.median=numeric(n),
                      rel.err.mean=numeric(n), rel.err.mode=numeric(n),
                      cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n),
                      hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                      test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                      test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                      test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                      test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                      test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                      test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

turnover_4 = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                      lower=numeric(n), rel.err.median=numeric(n),
                      rel.err.mean=numeric(n), rel.err.mode=numeric(n),
                      cv=numeric(n), hpd.lower=numeric(n), hpd.upper=numeric(n),
                      hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                      test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                      test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                      test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                      test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                      test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                      test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))



samplingProportion_1 = data.frame(true=numeric(n), estimated=numeric(n),
                                upper=numeric(n), lower=numeric(n),
                                rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                rel.err.mode=numeric(n), cv=numeric(n),
                                hpd.lower=numeric(n), hpd.upper=numeric(n),
                                hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

samplingProportion_2 = data.frame(true=numeric(n), estimated=numeric(n),
                                upper=numeric(n), lower=numeric(n),
                                rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                rel.err.mode=numeric(n), cv=numeric(n),
                                hpd.lower=numeric(n), hpd.upper=numeric(n),
                                hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

samplingProportion_3 = data.frame(true=numeric(n), estimated=numeric(n),
                                upper=numeric(n), lower=numeric(n),
                                rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                rel.err.mode=numeric(n), cv=numeric(n),
                                hpd.lower=numeric(n), hpd.upper=numeric(n),
                                hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

samplingProportion_4 = data.frame(true=numeric(n), estimated=numeric(n),
                                upper=numeric(n), lower=numeric(n),
                                rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                rel.err.mode=numeric(n), cv=numeric(n),
                                hpd.lower=numeric(n), hpd.upper=numeric(n),
                                hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

samplingAtPresentProb = data.frame(true=numeric(n), estimated=numeric(n),
                                   upper=numeric(n), lower=numeric(n),
                                   rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                   rel.err.mode=numeric(n), cv=numeric(n),
                                   hpd.lower=numeric(n), hpd.upper=numeric(n),
                                   hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                   test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                   test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                   test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                   test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                   test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                   test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

origin = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                    lower=numeric(n), rel.err.median=numeric(n), rel.err.mean=numeric(n),
                    rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n),
                    hpd.upper=numeric(n), hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                    test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                    test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                    test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                    test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                    test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                    test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

mrca = data.frame(true=numeric(n), estimated=numeric(n), upper=numeric(n),
                  lower=numeric(n), rel.err.median=numeric(n), rel.err.mean=numeric(n),
                  rel.err.mode=numeric(n), cv=numeric(n), hpd.lower=numeric(n),
                  hpd.upper=numeric(n), hpd.rel.width=numeric(n), test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                  test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                  test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                  test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                  test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                  test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                  test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))

samplingAtPresentProb = data.frame(true=numeric(n), estimated=numeric(n),
                                    upper=numeric(n), lower=numeric(n),
                                    rel.err.median=numeric(n), rel.err.mean=numeric(n),
                                    rel.err.mode=numeric(n), cv=numeric(n),
                                    hpd.lower=numeric(n), hpd.upper=numeric(n),
                                    hpd.rel.width=numeric(n), test_0=numeric(n),
                                    test_0=numeric(n), test_5=numeric(n), test_10=numeric(n),
                                    test_15=numeric(n), test_20=numeric(n), test_25=numeric(n),
                                    test_30=numeric(n), test_35=numeric(n), test_40=numeric(n),
                                    test_45=numeric(n), test_50=numeric(n), test_55=numeric(n),
                                    test_60=numeric(n), test_65=numeric(n), test_70=numeric(n),
                                    test_75=numeric(n), test_80=numeric(n), test_85=numeric(n),
                                    test_90=numeric(n), test_95=numeric(n), test_100=numeric(n))




post_ess<- numeric(n)
for (i in 1:200){
   cat("i =",i, " ")
  if (!(file.exists(paste0("../logs2/", (i-1), ".log")))){
    remove_rows = cbind(remove_rows, i)
    cat("FILE DOES NOT EXIST ", (i-1), "\n")
    next
  }
  t=read.table(paste0("../logs2/", (i-1), ".log"), header=TRUE, sep="\t")
  
  # print(head(t))
  # used.rates = true.rates[which(true.rates$X==i),];
  used.rates <- NULL
  used.rates$birthRate_1 <- true.rates$lambda.1[i]
  used.rates$birthRate_2 <- true.rates$lambda.2[i]
  used.rates$deathRate_1 <- true.rates$mu.1[i]
  used.rates$deathRate_2 <- true.rates$mu.2[i]
  used.rates$samplingRate_1 <- true.rates$psi.1[i]
  used.rates$samplingRate_2 <- true.rates$psi.2[i]

    used.rates$birthRate_3 <- true.rates$lambda.3[i]
  used.rates$birthRate_4 <- true.rates$lambda.4[i]
  used.rates$deathRate_3 <- true.rates$mu.3[i]
  used.rates$deathRate_4 <- true.rates$mu.4[i]
  used.rates$samplingRate_3 <- true.rates$psi.3[i]
  used.rates$samplingRate_4 <- true.rates$psi.4[i]

  used.rates$rho <- true.rates$rho[i]
  used.rates$div_rate_1 <- true.rates$div_rate.1[i]
    used.rates$div_rate_2 <- true.rates$div_rate.2[i]
  used.rates$turnover_1 <- true.rates$turnover.1[i]
  used.rates$turnover_2 <- true.rates$turnover.2[i]
  used.rates$sampling_prop_1 <- true.rates$sampling_prop.1[i]
  used.rates$sampling_prop_2 <- true.rates$sampling_prop.2[i]

    used.rates$div_rate_3 <- true.rates$div_rate.3[i]
    used.rates$div_rate_4 <- true.rates$div_rate.4[i]
  used.rates$turnover_3 <- true.rates$turnover.3[i]
  used.rates$turnover_4 <- true.rates$turnover.4[i]
  used.rates$sampling_prop_3 <- true.rates$sampling_prop.3[i]
  used.rates$sampling_prop_4 <- true.rates$sampling_prop.4[i]

  used.rates$mrca <- true.rates$mrca[i]
  used.rates$origin <- true.rates$origin[i]
  
  # take a 10% burnin
  t <- t[-seq(1,ceiling(length(t$netDiversification.1)/10)), ]
  
  if (nrow(t)==0){
    remove_rows = cbind(remove_rows, i)
    next
  }
  ess<-  effectiveSize(as.mcmc(t))
  post_ess[i] = as.numeric(ess["posterior"])
  
  if(post_ess[i]<100){
    remove_rows = cbind(remove_rows, i)
    next
  }
  
  
  diversificationRate_1$true[i] <- used.rates$div_rate_1
  diversificationRate_1$estimated[i] <- median(t$netDiversification.1)
  diversificationRate_1$upper[i] <- quantile(t$netDiversification.1,0.975)
  diversificationRate_1$lower[i] <- quantile(t$netDiversification.1,0.025)
  
  
  diversificationRate_1$rel.err.median[i] <- abs(median(t$netDiversification.1)-used.rates$div_rate_1)/used.rates$div_rate_1
  diversificationRate_1$rel.err.mean[i] <-abs(mean(t$netDiversification.1)-used.rates$div_rate_1)/used.rates$div_rate_1
  diversificationRate_1$rel.err.mode[i] <- abs(Mode(t$netDiversification.1)-used.rates$div_rate_1)/used.rates$div_rate_1
  diversificationRate_1$cv[i] <- sqrt(exp(sd(log(t$netDiversification.1))**2)-1)
  diversificationRate_1$hpd.lower[i] <- HPDinterval(as.mcmc(t$netDiversification.1))[1]
  diversificationRate_1$hpd.upper[i] <- HPDinterval(as.mcmc(t$netDiversification.1))[2]
  diversificationRate_1$hpd.rel.width[i] <- (diversificationRate_1$hpd.upper[i]-diversificationRate_1$hpd.lower[i])/used.rates$div_rate_1
  diversificationRate_1[i,(ncol(diversificationRate_1)-21+1):ncol(diversificationRate_1)] <- lapply(seq(0.,1,0.05), test_vals, diversificationRate_1$true[i], t$netDiversification.1)
  
  
  diversificationRate_2$true[i] <- used.rates$div_rate_2
  diversificationRate_2$estimated[i] <- median(t$netDiversification.2)
  diversificationRate_2$upper[i] <- quantile(t$netDiversification.2,0.975)
  diversificationRate_2$lower[i] <- quantile(t$netDiversification.2,0.025)
  
  
  diversificationRate_2$rel.err.median[i] <- abs(median(t$netDiversification.2)-used.rates$div_rate_2)/used.rates$div_rate_2
  diversificationRate_2$rel.err.mean[i] <-abs(mean(t$netDiversification.2)-used.rates$div_rate_2)/used.rates$div_rate_2
  diversificationRate_2$rel.err.mode[i] <- abs(Mode(t$netDiversification.2)-used.rates$div_rate_2)/used.rates$div_rate_2
  diversificationRate_2$cv[i] <- sqrt(exp(sd(log(t$netDiversification.2))**2)-1)
  diversificationRate_2$hpd.lower[i] <- HPDinterval(as.mcmc(t$netDiversification.2))[1]
  diversificationRate_2$hpd.upper[i] <- HPDinterval(as.mcmc(t$netDiversification.2))[2]
  diversificationRate_2$hpd.rel.width[i] <- (diversificationRate_2$hpd.upper[i]-diversificationRate_2$hpd.lower[i])/used.rates$div_rate_2
  diversificationRate_2[i,(ncol(diversificationRate_2)-21+1):ncol(diversificationRate_2)] <- lapply(seq(0.,1,0.05), test_vals, diversificationRate_2$true[i], t$netDiversification.2)
  
    diversificationRate_3$true[i] <- used.rates$div_rate_3
  diversificationRate_3$estimated[i] <- median(t$netDiversification.3)
  diversificationRate_3$upper[i] <- quantile(t$netDiversification.3,0.975)
  diversificationRate_3$lower[i] <- quantile(t$netDiversification.3,0.025)
  
  
  diversificationRate_3$rel.err.median[i] <- abs(median(t$netDiversification.3)-used.rates$div_rate_3)/used.rates$div_rate_3
  diversificationRate_3$rel.err.mean[i] <-abs(mean(t$netDiversification.3)-used.rates$div_rate_3)/used.rates$div_rate_3
  diversificationRate_3$rel.err.mode[i] <- abs(Mode(t$netDiversification.3)-used.rates$div_rate_3)/used.rates$div_rate_3
  diversificationRate_3$cv[i] <- sqrt(exp(sd(log(t$netDiversification.3))**2)-1)
  diversificationRate_3$hpd.lower[i] <- HPDinterval(as.mcmc(t$netDiversification.3))[1]
  diversificationRate_3$hpd.upper[i] <- HPDinterval(as.mcmc(t$netDiversification.3))[2]
  diversificationRate_3$hpd.rel.width[i] <- (diversificationRate_3$hpd.upper[i]-diversificationRate_3$hpd.lower[i])/used.rates$div_rate_3
  diversificationRate_3[i,(ncol(diversificationRate_3)-21+1):ncol(diversificationRate_3)] <- lapply(seq(0.,1,0.05), test_vals, diversificationRate_3$true[i], t$netDiversification.3)
  
    diversificationRate_4$true[i] <- used.rates$div_rate_4
  diversificationRate_4$estimated[i] <- median(t$netDiversification.4)
  diversificationRate_4$upper[i] <- quantile(t$netDiversification.4,0.975)
  diversificationRate_4$lower[i] <- quantile(t$netDiversification.4,0.025)
  
  
  diversificationRate_4$rel.err.median[i] <- abs(median(t$netDiversification.4)-used.rates$div_rate_4)/used.rates$div_rate_4
  diversificationRate_4$rel.err.mean[i] <-abs(mean(t$netDiversification.4)-used.rates$div_rate_4)/used.rates$div_rate_4
  diversificationRate_4$rel.err.mode[i] <- abs(Mode(t$netDiversification.4)-used.rates$div_rate_4)/used.rates$div_rate_4
  diversificationRate_4$cv[i] <- sqrt(exp(sd(log(t$netDiversification.4))**2)-1)
  diversificationRate_4$hpd.lower[i] <- HPDinterval(as.mcmc(t$netDiversification.4))[1]
  diversificationRate_4$hpd.upper[i] <- HPDinterval(as.mcmc(t$netDiversification.4))[2]
  diversificationRate_4$hpd.rel.width[i] <- (diversificationRate_4$hpd.upper[i]-diversificationRate_4$hpd.lower[i])/used.rates$div_rate_4
  diversificationRate_4[i,(ncol(diversificationRate_4)-21+1):ncol(diversificationRate_4)] <- lapply(seq(0.,1,0.05), test_vals, diversificationRate_4$true[i], t$netDiversification.4)
  

  turnover_1$true[i] <- used.rates$turnover_1
  turnover_1$estimated[i] <- median(t$turnOver.1)
  turnover_1$upper[i] <- quantile(t$turnOver.1,0.975)
  turnover_1$lower[i] <- quantile(t$turnOver.1,0.025)
  
  
  turnover_1$rel.err.median[i] <- abs(median(t$turnOver.1)-used.rates$turnover_1)/used.rates$turnover_1
  turnover_1$rel.err.mean[i] <-abs(mean(t$turnOver.1)-used.rates$turnover_1)/used.rates$turnover_1
  turnover_1$rel.err.mode[i] <- abs(Mode(t$turnOver.1)-used.rates$turnover_1)/used.rates$turnover_1
  turnover_1$cv[i] <- sqrt(exp(sd(log(t$turnOver.1))**2)-1)
  turnover_1$hpd.lower[i] <- HPDinterval(as.mcmc(t$turnOver.1))[1]
  turnover_1$hpd.upper[i] <- HPDinterval(as.mcmc(t$turnOver.1))[2]
  turnover_1$hpd.rel.width[i] <- (turnover_1$hpd.upper[i]-turnover_1$hpd.lower[i])/used.rates$turnover_1
  turnover_1[i,(ncol(turnover_1)-21+1):ncol(turnover_1)] <- lapply(seq(0.,1,0.05), test_vals, turnover_1$true[i], t$turnOver.1)
  
  turnover_2$true[i] <- used.rates$turnover_2
  turnover_2$estimated[i] <- median(t$turnOver.2)
  turnover_2$upper[i] <- quantile(t$turnOver.2,0.975)
  turnover_2$lower[i] <- quantile(t$turnOver.2,0.025)
  
  
  turnover_2$rel.err.median[i] <- abs(median(t$turnOver.2)-used.rates$turnover_2)/used.rates$turnover_2
  turnover_2$rel.err.mean[i] <-abs(mean(t$turnOver.2)-used.rates$turnover_2)/used.rates$turnover_2
  turnover_2$rel.err.mode[i] <- abs(Mode(t$turnOver.2)-used.rates$turnover_2)/used.rates$turnover_2
  turnover_2$cv[i] <- sqrt(exp(sd(log(t$turnOver.2))**2)-1)
  turnover_2$hpd.lower[i] <- HPDinterval(as.mcmc(t$turnOver.2))[1]
  turnover_2$hpd.upper[i] <- HPDinterval(as.mcmc(t$turnOver.2))[2]
  turnover_2$hpd.rel.width[i] <- (turnover_2$hpd.upper[i]-turnover_2$hpd.lower[i])/used.rates$turnover_2
  turnover_2[i,(ncol(turnover_2)-21+1):ncol(turnover_2)] <- lapply(seq(0.,1,0.05), test_vals, turnover_2$true[i], t$turnOver.2)
  

    turnover_3$true[i] <- used.rates$turnover_3
  turnover_3$estimated[i] <- median(t$turnOver.3)
  turnover_3$upper[i] <- quantile(t$turnOver.3,0.975)
  turnover_3$lower[i] <- quantile(t$turnOver.3,0.025)
  
  
  turnover_3$rel.err.median[i] <- abs(median(t$turnOver.3)-used.rates$turnover_3)/used.rates$turnover_3
  turnover_3$rel.err.mean[i] <-abs(mean(t$turnOver.3)-used.rates$turnover_3)/used.rates$turnover_3
  turnover_3$rel.err.mode[i] <- abs(Mode(t$turnOver.3)-used.rates$turnover_3)/used.rates$turnover_3
  turnover_3$cv[i] <- sqrt(exp(sd(log(t$turnOver.3))**2)-1)
  turnover_3$hpd.lower[i] <- HPDinterval(as.mcmc(t$turnOver.3))[1]
  turnover_3$hpd.upper[i] <- HPDinterval(as.mcmc(t$turnOver.3))[2]
  turnover_3$hpd.rel.width[i] <- (turnover_3$hpd.upper[i]-turnover_3$hpd.lower[i])/used.rates$turnover_3
  turnover_3[i,(ncol(turnover_3)-21+1):ncol(turnover_3)] <- lapply(seq(0.,1,0.05), test_vals, turnover_3$true[i], t$turnOver.3)
  
    turnover_4$true[i] <- used.rates$turnover_4
  turnover_4$estimated[i] <- median(t$turnOver.4)
  turnover_4$upper[i] <- quantile(t$turnOver.4,0.975)
  turnover_4$lower[i] <- quantile(t$turnOver.4,0.025)
  
  
  turnover_4$rel.err.median[i] <- abs(median(t$turnOver.4)-used.rates$turnover_4)/used.rates$turnover_4
  turnover_4$rel.err.mean[i] <-abs(mean(t$turnOver.4)-used.rates$turnover_4)/used.rates$turnover_4
  turnover_4$rel.err.mode[i] <- abs(Mode(t$turnOver.4)-used.rates$turnover_4)/used.rates$turnover_4
  turnover_4$cv[i] <- sqrt(exp(sd(log(t$turnOver.4))**2)-1)
  turnover_4$hpd.lower[i] <- HPDinterval(as.mcmc(t$turnOver.4))[1]
  turnover_4$hpd.upper[i] <- HPDinterval(as.mcmc(t$turnOver.4))[2]
  turnover_4$hpd.rel.width[i] <- (turnover_4$hpd.upper[i]-turnover_4$hpd.lower[i])/used.rates$turnover_4
  turnover_4[i,(ncol(turnover_4)-21+1):ncol(turnover_4)] <- lapply(seq(0.,1,0.05), test_vals, turnover_4$true[i], t$turnOver.4)

  samplingProportion_1$true[i] <- used.rates$sampling_prop_1
  samplingProportion_1$estimated[i] <- median(t$samplingProportion.1)
  samplingProportion_1$upper[i] <- quantile(t$samplingProportion.1,0.975)
  samplingProportion_1$lower[i] <- quantile(t$samplingProportion.1,0.025)
  
  
  samplingProportion_1$rel.err.median[i] <- abs(median(t$samplingProportion.1)-used.rates$sampling_prop_1)/used.rates$sampling_prop_1
  samplingProportion_1$rel.err.mean[i] <-abs(mean(t$samplingProportion.1)-used.rates$sampling_prop_1)/used.rates$sampling_prop_1
  samplingProportion_1$rel.err.mode[i] <- abs(Mode(t$samplingProportion.1)-used.rates$sampling_prop_1)/used.rates$sampling_prop_1
  samplingProportion_1$cv[i] <- sqrt(exp(sd(log(t$samplingProportion.1))**2)-1)
  samplingProportion_1$hpd.lower[i] <- HPDinterval(as.mcmc(t$samplingProportion.1))[1]
  samplingProportion_1$hpd.upper[i] <- HPDinterval(as.mcmc(t$samplingProportion.1))[2]
  samplingProportion_1$hpd.rel.width[i] <- (samplingProportion_1$hpd.upper[i]-samplingProportion_1$hpd.lower[i])/used.rates$sampling_prop_1
  samplingProportion_1[i,(ncol(samplingProportion_1)-21+1):ncol(samplingProportion_1)] <- lapply(seq(0.,1,0.05), test_vals, samplingProportion_1$true[i], t$samplingProportion.1)
  
  samplingProportion_2$true[i] <- used.rates$sampling_prop_2
  samplingProportion_2$estimated[i] <- median(t$samplingProportion.2)
  samplingProportion_2$upper[i] <- quantile(t$samplingProportion.2,0.975)
  samplingProportion_2$lower[i] <- quantile(t$samplingProportion.2,0.025)
  
  
  samplingProportion_2$rel.err.median[i] <- abs(median(t$samplingProportion.2)-used.rates$sampling_prop_2)/used.rates$sampling_prop_2
  samplingProportion_2$rel.err.mean[i] <-abs(mean(t$samplingProportion.2)-used.rates$sampling_prop_2)/used.rates$sampling_prop_2
  samplingProportion_2$rel.err.mode[i] <- abs(Mode(t$samplingProportion.2)-used.rates$sampling_prop_2)/used.rates$sampling_prop_2
  samplingProportion_2$cv[i] <- sqrt(exp(sd(log(t$samplingProportion.2))**2)-1)
  samplingProportion_2$hpd.lower[i] <- HPDinterval(as.mcmc(t$samplingProportion.2))[1]
  samplingProportion_2$hpd.upper[i] <- HPDinterval(as.mcmc(t$samplingProportion.2))[2]
  samplingProportion_2$hpd.rel.width[i] <- (samplingProportion_2$hpd.upper[i]-samplingProportion_2$hpd.lower[i])/used.rates$sampling_prop_2
  samplingProportion_2[i,(ncol(samplingProportion_2)-21+1):ncol(samplingProportion_2)] <- lapply(seq(0.,1,0.05), test_vals, samplingProportion_2$true[i], t$samplingProportion.2)
  
  samplingProportion_3$true[i] <- used.rates$sampling_prop_3
  samplingProportion_3$estimated[i] <- median(t$samplingProportion.3)
  samplingProportion_3$upper[i] <- quantile(t$samplingProportion.3,0.975)
  samplingProportion_3$lower[i] <- quantile(t$samplingProportion.3,0.025)
  
  
  samplingProportion_3$rel.err.median[i] <- abs(median(t$samplingProportion.3)-used.rates$sampling_prop_3)/used.rates$sampling_prop_3
  samplingProportion_3$rel.err.mean[i] <-abs(mean(t$samplingProportion.3)-used.rates$sampling_prop_3)/used.rates$sampling_prop_3
  samplingProportion_3$rel.err.mode[i] <- abs(Mode(t$samplingProportion.3)-used.rates$sampling_prop_3)/used.rates$sampling_prop_3
  samplingProportion_3$cv[i] <- sqrt(exp(sd(log(t$samplingProportion.3))**2)-1)
  samplingProportion_3$hpd.lower[i] <- HPDinterval(as.mcmc(t$samplingProportion.3))[1]
  samplingProportion_3$hpd.upper[i] <- HPDinterval(as.mcmc(t$samplingProportion.3))[2]
  samplingProportion_3$hpd.rel.width[i] <- (samplingProportion_3$hpd.upper[i]-samplingProportion_3$hpd.lower[i])/used.rates$sampling_prop_3
  samplingProportion_3[i,(ncol(samplingProportion_3)-21+1):ncol(samplingProportion_3)] <- lapply(seq(0.,1,0.05), test_vals, samplingProportion_3$true[i], t$samplingProportion.3)
  
    samplingProportion_4$true[i] <- used.rates$sampling_prop_4
  samplingProportion_4$estimated[i] <- median(t$samplingProportion.4)
  samplingProportion_4$upper[i] <- quantile(t$samplingProportion.4,0.975)
  samplingProportion_4$lower[i] <- quantile(t$samplingProportion.4,0.025)
  
  
  samplingProportion_4$rel.err.median[i] <- abs(median(t$samplingProportion.4)-used.rates$sampling_prop_4)/used.rates$sampling_prop_4
  samplingProportion_4$rel.err.mean[i] <-abs(mean(t$samplingProportion.4)-used.rates$sampling_prop_4)/used.rates$sampling_prop_4
  samplingProportion_4$rel.err.mode[i] <- abs(Mode(t$samplingProportion.4)-used.rates$sampling_prop_4)/used.rates$sampling_prop_4
  samplingProportion_4$cv[i] <- sqrt(exp(sd(log(t$samplingProportion.4))**2)-1)
  samplingProportion_4$hpd.lower[i] <- HPDinterval(as.mcmc(t$samplingProportion.4))[1]
  samplingProportion_4$hpd.upper[i] <- HPDinterval(as.mcmc(t$samplingProportion.4))[2]
  samplingProportion_4$hpd.rel.width[i] <- (samplingProportion_4$hpd.upper[i]-samplingProportion_4$hpd.lower[i])/used.rates$sampling_prop_4
  samplingProportion_4[i,(ncol(samplingProportion_4)-21+1):ncol(samplingProportion_4)] <- lapply(seq(0.,1,0.05), test_vals, samplingProportion_4$true[i], t$samplingProportion.4)
  

  samplingAtPresentProb$true[i] <- used.rates$rho
  samplingAtPresentProb$estimated[i] <- median(t$rho)
  samplingAtPresentProb$upper[i] <- quantile(t$rho,0.975)
  samplingAtPresentProb$lower[i] <- quantile(t$rho,0.025)
  
  
  samplingAtPresentProb$rel.err.median[i] <- abs(median(t$rho)-used.rates$rho)/used.rates$rho
  samplingAtPresentProb$rel.err.mean[i] <-abs(mean(t$rho)-used.rates$rho)/used.rates$rho
  samplingAtPresentProb$rel.err.mode[i] <- abs(Mode(t$rho)-used.rates$rho)/used.rates$rho
  samplingAtPresentProb$cv[i] <- sqrt(exp(sd(log(t$rho))**2)-1)
  samplingAtPresentProb$hpd.lower[i] <- HPDinterval(as.mcmc(t$rho))[1]
  samplingAtPresentProb$hpd.upper[i] <- HPDinterval(as.mcmc(t$rho))[2]
  samplingAtPresentProb$hpd.rel.width[i] <- (samplingAtPresentProb$hpd.upper[i]-samplingAtPresentProb$hpd.lower[i])/used.rates$rho
  samplingAtPresentProb[i,(ncol(samplingAtPresentProb)-21+1):ncol(samplingAtPresentProb)] <- lapply(seq(0.,1,0.05), test_vals, samplingAtPresentProb$true[i], t$rho)
  
  
  
  origin$true[i] <- used.rates$origin
  origin$estimated[i] <- median(t$origin)
  origin$upper[i] <- quantile(t$origin,0.975)
  origin$lower[i] <- quantile(t$origin,0.025)
  
  
  origin$rel.err.median[i] <- abs(median(t$origin)-used.rates$origin)/used.rates$origin
  origin$rel.err.mean[i] <-abs(mean(t$origin)-used.rates$origin)/used.rates$origin
  origin$rel.err.mode[i] <- abs(Mode(t$origin)-used.rates$origin)/used.rates$origin
  origin$cv[i] <- sqrt(exp(sd(log(t$origin))**2)-1)
  origin$hpd.lower[i] <- HPDinterval(as.mcmc(t$origin))[1]
  origin$hpd.upper[i] <- HPDinterval(as.mcmc(t$origin))[2]
  origin$hpd.rel.width[i] <- (origin$hpd.upper[i]-origin$hpd.lower[i])/used.rates$origin
  origin[i,(ncol(origin)-21+1):ncol(origin)] <- lapply(seq(0.,1,0.05), test_vals, origin$true[i], t$origin)
  
  
  mrca$true[i] <- used.rates$mrca
  mrca$estimated[i] <- median(t$TreeHeight)
  mrca$upper[i] <- quantile(t$TreeHeight,0.975)
  mrca$lower[i] <- quantile(t$TreeHeight,0.025)
  
  
  mrca$rel.err.median[i] <- abs(median(t$TreeHeight)-used.rates$mrca)/used.rates$mrca
  mrca$rel.err.mean[i] <-abs(mean(t$TreeHeight)-used.rates$mrca)/used.rates$mrca
  mrca$rel.err.mode[i] <- abs(Mode(t$TreeHeight)-used.rates$mrca)/used.rates$mrca
  mrca$cv[i] <- sqrt(exp(sd(log(t$TreeHeight))**2)-1)
  mrca$hpd.lower[i] <- HPDinterval(as.mcmc(t$TreeHeight))[1]
  mrca$hpd.upper[i] <- HPDinterval(as.mcmc(t$TreeHeight))[2]
  mrca$hpd.rel.width[i] <- (mrca$hpd.upper[i]-mrca$hpd.lower[i])/used.rates$mrca
  mrca[i,(ncol(mrca)-21+1):ncol(mrca)] <- lapply(seq(0.,1,0.05), test_vals, mrca$true[i], t$TreeHeight)
  
}

save.image(file="figures/validationPlotEnvironment.RData")
# load(file="figures/validationPlotEnvironment.RData")

if (!is.null(remove_rows)){
  diversificationRate_1 <- diversificationRate_1[-remove_rows, ]
    diversificationRate_2 <- diversificationRate_2[-remove_rows, ]
  turnover_1 <- turnover_1[-remove_rows, ]
  turnover_2 <- turnover_2[-remove_rows, ]
  mrca <- mrca[-remove_rows, ]
  origin <- origin[-remove_rows, ]
  samplingAtPresentProb <- samplingAtPresentProb[-remove_rows, ]
  samplingProportion_1 <- samplingProportion_1[-remove_rows, ]
    samplingProportion_2 <- samplingProportion_2[-remove_rows, ]

      diversificationRate_3 <- diversificationRate_3[-remove_rows, ]
    diversificationRate_4 <- diversificationRate_4[-remove_rows, ]
  turnover_3 <- turnover_3[-remove_rows, ]
  turnover_4 <- turnover_4[-remove_rows, ]
    samplingProportion_3 <- samplingProportion_1[-remove_rows, ]
    samplingProportion_4 <- samplingProportion_2[-remove_rows, ]

}


d1 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate_1[,(ncol(diversificationRate_1)-21+1):ncol(diversificationRate_1)], sum)/length(diversificationRate_1$test_0))

l <- length(diversificationRate_1$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd1 <- data.frame(p=cc, vals=unlist(diversificationRate_1[,(ncol(diversificationRate_1)-21+1):ncol(diversificationRate_1)]))

p.qq_div_rate_1 <- ggplot(d1)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_div_rate_1 <- p.qq_div_rate_1+
  stat_summary(data=dd1, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 1")
ggsave(plot=p.qq_div_rate_1,paste("figures/qq_diversificationRate_1", ".pdf", sep=""),width=6, height=5)

x.min_diversificationRate = min(diversificationRate_1$true)
x.max_diversificationRate = max(diversificationRate_1$true)

y.min_diversificationRate = min(diversificationRate_1$lower)
y.max_diversificationRate = max(diversificationRate_1$upper)

lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate_1 <- ggplot(diversificationRate_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.diversificationRate_1,paste("figures/diversificationRate_1", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_1_log <- ggplot(diversificationRate_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.diversificationRate_1_log,paste("figures/diversificationRate_1_logScale", ".pdf", sep=""),width=6, height=5)



## div 2
d2 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate_2[,(ncol(diversificationRate_2)-21+1):ncol(diversificationRate_2)], sum)/length(diversificationRate_2$test_0))

l <- length(diversificationRate_2$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd2 <- data.frame(p=cc, vals=unlist(diversificationRate_2[,(ncol(diversificationRate_2)-21+1):ncol(diversificationRate_2)]))

p.qq_div_rate_2 <- ggplot(d2)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_div_rate_2 <- p.qq_div_rate_2+
  stat_summary(data=dd2, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 2")


ggsave(plot=p.qq_div_rate_2,paste("figures/qq_diversificationRate_2", ".pdf", sep=""),width=6, height=5)


x.min_diversificationRate = min(diversificationRate_2$true)
x.max_diversificationRate = max(diversificationRate_2$true)

y.min_diversificationRate = min(diversificationRate_2$lower)
y.max_diversificationRate = max(diversificationRate_2$upper)

lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate_2 <- ggplot(diversificationRate_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()


ggsave(plot=p.diversificationRate_2,paste("figures/diversificationRate_2", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_2_log <- ggplot(diversificationRate_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.diversificationRate_2_log,paste("figures/diversificationRate_2_logScale", ".pdf", sep=""),width=6, height=5)

# div 3
d3 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate_3[,(ncol(diversificationRate_3)-21+1):ncol(diversificationRate_3)], sum)/length(diversificationRate_3$test_0))

l <- length(diversificationRate_3$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd3 <- data.frame(p=cc, vals=unlist(diversificationRate_3[,(ncol(diversificationRate_3)-21+1):ncol(diversificationRate_3)]))

p.qq_div_rate_3 <- ggplot(d3)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_div_rate_3 <- p.qq_div_rate_3+
  stat_summary(data=dd3, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 3")
ggsave(plot=p.qq_div_rate_3,paste("figures/qq_diversificationRate_3", ".pdf", sep=""),width=6, height=5)

x.min_diversificationRate = min(diversificationRate_3$true)
x.max_diversificationRate = max(diversificationRate_3$true)

y.min_diversificationRate = min(diversificationRate_3$lower)
y.max_diversificationRate = max(diversificationRate_3$upper)

lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate_3 <- ggplot(diversificationRate_3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.diversificationRate_3,paste("figures/diversificationRate_3", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_log <- ggplot(diversificationRate_3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.diversificationRate_log,paste("figures/diversificationRate_logScale", ".pdf", sep=""),width=6, height=5)


# div 4
d4 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate_4[,(ncol(diversificationRate_4)-21+1):ncol(diversificationRate_4)], sum)/length(diversificationRate_4$test_0))

l <- length(diversificationRate_4$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd4 <- data.frame(p=cc, vals=unlist(diversificationRate_4[,(ncol(diversificationRate_4)-21+1):ncol(diversificationRate_4)]))

x.min_diversificationRate = min(diversificationRate_4$true)
x.max_diversificationRate = max(diversificationRate_4$true)

y.min_diversificationRate = min(diversificationRate_4$lower)
y.max_diversificationRate = max(diversificationRate_4$upper)

lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate_4 <- ggplot(diversificationRate_4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.diversificationRate_4,paste("figures/diversificationRate_4", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_4_log <- ggplot(diversificationRate_4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.diversificationRate_4_log,paste("figures/diversificationRate_logScale", ".pdf", sep=""),width=6, height=5)

p.qq_div_rate_4 <- ggplot(d4)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()


p.qq_div_rate_4 <- p.qq_div_rate_4+
  stat_summary(data=dd4, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 4")
ggsave(plot=p.qq_div_rate_4,paste("figures/qq_diversificationRate_4", ".pdf", sep=""),width=6, height=5)

# combined div plot

d_combined <- bind_rows(d1, d2, d3, d4) # 
dd_combined <- bind_rows(dd1,dd2, dd3, dd4) #  


p.qq_div_rate_combined <- ggplot(d_combined)+
  geom_abline(intercept = 0, color="#BBBBBB", linewidth=1)+
   xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()

legend_colors <- c("4" = "#4477AA", "1" = "#228833", "2" = "#EE6677", "3" = "#AA3377")

p.qq_div_rate_combined <- p.qq_div_rate_combined +
                 stat_summary(data=dd4, aes(x=(p+0.008), y=vals, colour="4"),fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
  stat_summary(data=dd1, aes(x=(p-0.008), y=vals, colour="1"),fun.data = "mean_cl_boot",  geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
                 stat_summary(data=dd2, aes(x=(p-0.0025), y=vals, colour="2"), fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1,  linewidth=0.65) +
                                stat_summary(data=dd3, aes(x=(p+0.003), y=vals, colour="3"),fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
                 scale_color_manual(values = legend_colors, name = "Interval") + 
                 theme(legend.position = "inside", legend.position.inside=c(0.85, 0.25), legend.background=element_rect(fill="white", colour=NA)) # +
 # ggtitle("Diversification")
  
ggsave(plot=p.qq_div_rate_combined,paste("figures/qq_combined_div", ".pdf", sep=""),width=6, height=5)


########## turnover ################


d1 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover_1[,(ncol(turnover_1)-21+1):ncol(turnover_1)], sum)/length(turnover_1$test_0))

l <- length(turnover_1$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd1 <- data.frame(p=cc, vals=unlist(turnover_1[,(ncol(turnover_1)-21+1):ncol(turnover_1)]))

p.qq_turnover_1 <- ggplot(d1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_turnover_1 <- p.qq_turnover_1 +
  stat_summary(data=dd1, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Turnover 1")

ggsave(plot=p.qq_turnover_1 ,paste("figures/qq_turnover_1", ".pdf", sep=""),width=6, height=5)

x.min_turnover = min(turnover_1$true)
x.max_turnover = max(turnover_1$true)

y.min_turnover = min(turnover_1$lower)
y.max_turnover = max(turnover_1$upper)


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover_1 <- ggplot(turnover_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.turnover_1,paste("figures/turnover_1", ".pdf", sep=""),width=6, height=5)

p.turnover_1_log <- ggplot(turnover_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.turnover_1_log,paste("figures/turnover_1_logScale", ".pdf", sep=""),width=6, height=5)

# turnover 2

d2 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover_2[,(ncol(turnover_2)-21+1):ncol(turnover_2)], sum)/length(turnover_2$test_0))

l <- length(turnover_2$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd2 <- data.frame(p=cc, vals=unlist(turnover_2[,(ncol(turnover_2)-21+1):ncol(turnover_2)]))

p.qq_turnover_2 <- ggplot(d2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_turnover_2 <- p.qq_turnover_2 +
  stat_summary(data=dd2, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Turnover 2")

ggsave(plot=p.qq_turnover_2 ,paste("figures/qq_turnover_2", ".pdf", sep=""),width=6, height=5)

x.min_turnover = min(turnover_2$true)
x.max_turnover = max(turnover_2$true)

y.min_turnover = min(turnover_2$lower)
y.max_turnover = max(turnover_2$upper)


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover_2 <- ggplot(turnover_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.turnover_2,paste("figures/turnover_2", ".pdf", sep=""),width=6, height=5)

p.turnover_2_log <- ggplot(turnover_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.turnover_2_log,paste("figures/turnover_2_logScale", ".pdf", sep=""),width=6, height=5)


# turnover 3
d3 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover_3[,(ncol(turnover_3)-21+1):ncol(turnover_3)], sum)/length(turnover_3$test_0))

l <- length(turnover_3$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd3 <- data.frame(p=cc, vals=unlist(turnover_3[,(ncol(turnover_3)-21+1):ncol(turnover_3)]))

p.qq_turnover_3 <- ggplot(d3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_turnover_3 <- p.qq_turnover_3 +
  stat_summary(data=dd3, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Turnover 3")

ggsave(plot=p.qq_turnover_3 ,paste("figures/qq_turnover_3", ".pdf", sep=""),width=6, height=5)

x.min_turnover = min(turnover_3$true)
x.max_turnover = max(turnover_3$true)

y.min_turnover = min(turnover_3$lower)
y.max_turnover = max(turnover_3$upper)


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover_3 <- ggplot(turnover_3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.turnover_3,paste("figures/turnover_3", ".pdf", sep=""),width=6, height=5)

p.turnover_3_log <- ggplot(turnover_3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.turnover_3_log,paste("figures/turnover_3_logScale", ".pdf", sep=""),width=6, height=5)

# turnover 4
d4 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover_4[,(ncol(turnover_4)-21+1):ncol(turnover_4)], sum)/length(turnover_4$test_0))

l <- length(turnover_4$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd4 <- data.frame(p=cc, vals=unlist(turnover_4[,(ncol(turnover_4)-21+1):ncol(turnover_4)]))

p.qq_turnover_4 <- ggplot(d4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_turnover_4 <- p.qq_turnover_4 +
  stat_summary(data=dd4, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Turnover")

ggsave(plot=p.qq_turnover_4 ,paste("figures/qq_turnover_4", ".pdf", sep=""),width=6, height=5)

x.min_turnover = min(turnover_4$true)
x.max_turnover = max(turnover_4$true)

y.min_turnover = min(turnover_4$lower)
y.max_turnover = max(turnover_4$upper)


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover_4 <- ggplot(turnover_4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.turnover_4,paste("figures/turnover_4", ".pdf", sep=""),width=6, height=5)

p.turnover_4_log <- ggplot(turnover_4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.turnover_4_log,paste("figures/turnover_4_logScale", ".pdf", sep=""),width=6, height=5)



# combined turnover plot

d_combined <- bind_rows(d1, d2, d3, d4) # 
dd_combined <- bind_rows(dd1,dd2, dd3, dd4) #  


p.qq_turnover_combined <- ggplot(d_combined)+
  geom_abline(intercept = 0, color="#BBBBBB", linewidth=1)+
   xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()

legend_colors <- c("4" = "#4477AA", "1" = "#228833", "2" = "#EE6677", "3" = "#AA3377")

p.qq_turnover_combined <- p.qq_turnover_combined +
                 stat_summary(data=dd4, aes(x=(p+0.008), y=vals, colour="4"),fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
  stat_summary(data=dd1, aes(x=(p-0.008), y=vals, colour="1"),fun.data = "mean_cl_boot",  geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
                 stat_summary(data=dd2, aes(x=(p-0.0025), y=vals, colour="2"), fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1,  linewidth=0.65) +
                                stat_summary(data=dd3, aes(x=(p+0.003), y=vals, colour="3"),fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
                 scale_color_manual(values = legend_colors, name = "Interval") + 
                 theme(legend.position = "inside", legend.position.inside=c(0.85, 0.25), legend.background=element_rect(fill="white", colour=NA)) # +
 # ggtitle("Turnover")
  
ggsave(plot=p.qq_turnover_combined,paste("figures/qq_combined_turnover", ".pdf", sep=""),width=6, height=5)

########## mrca 2 ################


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(mrca[,(ncol(mrca)-21+1):ncol(mrca)], sum)/length(mrca$test_0))

l <- length(mrca$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(mrca[,(ncol(mrca)-21+1):ncol(mrca)]))

x.min_mrca = min(mrca$true)
x.max_mrca = max(mrca$true)

y.min_mrca = min(mrca$lower)
y.max_mrca = max(mrca$upper)


lim.min = min(x.min_mrca, y.min_mrca)
lim.max = max(x.max_mrca, y.max_mrca)

p.mrca <- ggplot(mrca)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()


ggsave(plot=p.mrca,paste("figures/mrca", ".pdf", sep=""),width=6, height=5)

p.mrca_log <- ggplot(mrca)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.mrca_log,paste("figures/mrca_logScale", ".pdf", sep=""),width=6, height=5)


########## origin ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(origin[,(ncol(origin)-21+1):ncol(origin)], sum)/length(origin$test_0))

l <- length(origin$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(origin[,(ncol(origin)-21+1):ncol(origin)]))


p.qq_origin <- ggplot(d)+
  geom_abline(intercept = 0, color="#BBBBBB", linewidth = 1)+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()
p.qq_origin <- p.qq_origin+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "#2D9689") # +
  # ggtitle("Origin")

ggsave(plot=p.qq_origin,paste("figures/qq_origin", ".pdf", sep=""),width=6, height=5)

x.min_origin = min(origin$true)
x.max_origin = max(origin$true)

y.min_origin = min(origin$lower)
y.max_origin = max(origin$upper)


lim.min = min(x.min_origin, y.min_origin)
lim.max = max(x.max_origin, y.max_origin)

p.origin <- ggplot(origin)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()


ggsave(plot=p.origin,paste("figures/origin", ".pdf", sep=""),width=6, height=5)

p.origin_log <- ggplot(origin)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.origin_log,paste("figures/origin_logScale", ".pdf", sep=""),width=6, height=5)


########## samplingAtPresentProb ################

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingAtPresentProb[,(ncol(samplingAtPresentProb)-21+1):ncol(samplingAtPresentProb)], sum)/length(samplingAtPresentProb$test_0))


l <- length(samplingAtPresentProb$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(samplingAtPresentProb[,(ncol(samplingAtPresentProb)-21+1):ncol(samplingAtPresentProb)]))

p.qq_samplingAtPresentProb <- ggplot(d)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="#BBBBBB", linewidth = 1)+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") +
  theme_minimal()

p.qq_samplingAtPresentProb <- p.qq_samplingAtPresentProb+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "#2D9689") # +
  # ggtitle("Probability of sampling at present")

ggsave(plot=p.qq_samplingAtPresentProb,paste("figures/qq_samplingAtPresentProb", ".pdf", sep=""),width=6, height=5)

x.min_samplingAtPresentProb = min(samplingAtPresentProb$true)
x.max_samplingAtPresentProb = max(samplingAtPresentProb$true)

y.min_samplingAtPresentProb = min(samplingAtPresentProb$lower)
y.max_samplingAtPresentProb = max(samplingAtPresentProb$upper)


lim.min = min(x.min_samplingAtPresentProb, y.min_samplingAtPresentProb)
lim.max = max(x.max_samplingAtPresentProb, y.max_samplingAtPresentProb)

p.samplingAtPresentProb <- ggplot(samplingAtPresentProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()


ggsave(plot=p.samplingAtPresentProb,paste("figures/samplingAtPresentProb", ".pdf", sep=""),width=6, height=5)

p.samplingAtPresentProb_log <- ggplot(samplingAtPresentProb)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingAtPresentProb_log,paste("figures/samplingAtPresentProb_logScale", ".pdf", sep=""),width=6, height=5)


########## samplingProportion ################

d1 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion_1[,(ncol(samplingProportion_1)-21+1):ncol(samplingProportion_1)], sum)/length(samplingProportion_1$test_0))

l <- length(samplingProportion_1$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd1 <- data.frame(p=cc, vals=unlist(samplingProportion_1[,(ncol(samplingProportion_1)-21+1):ncol(samplingProportion_1)]))


p.qq_samplingProportion_1 <- ggplot(d1)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") +
  theme_minimal()

p.qq_samplingProportion_1 <- p.qq_samplingProportion_1+
  stat_summary(data=dd1, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Sampling proportion 1")

ggsave(plot=p.qq_samplingProportion_1,paste("figures/qq_samplingProportion_1", ".pdf", sep=""),width=6, height=5)

x.min_samplingProportion_1 = min(samplingProportion_1$true)
x.max_samplingProportion_1 = max(samplingProportion_1$true)

y.min_samplingProportion_1 = min(samplingProportion_1$lower)
y.max_samplingProportion_1 = max(samplingProportion_1$upper)


lim.min = min(x.min_samplingProportion_1, y.min_samplingProportion_1)
lim.max = max(x.max_samplingProportion_1, y.max_samplingProportion_1)

p.samplingProportion_1 <- ggplot(samplingProportion_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.samplingProportion_1,paste("figures/samplingProportion_1", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_1_log <- ggplot(samplingProportion_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingProportion_1_log,paste("figures/samplingProportion_1_logScale", ".pdf", sep=""),width=6, height=5)

# sampling proportion 2

d2 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion_2[,(ncol(samplingProportion_2)-21+1):ncol(samplingProportion_2)], sum)/length(samplingProportion_2$test_0))

l <- length(samplingProportion_2$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd2 <- data.frame(p=cc, vals=unlist(samplingProportion_2[,(ncol(samplingProportion_2)-21+1):ncol(samplingProportion_2)]))


p.qq_samplingProportion_2 <- ggplot(d2)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") +
  theme_minimal()

p.qq_samplingProportion_2 <- p.qq_samplingProportion_2+
  stat_summary(data=dd2, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Sampling proportion 2")

ggsave(plot=p.qq_samplingProportion_2,paste("figures/qq_samplingProportion_2", ".pdf", sep=""),width=6, height=5)

x.min_samplingProportion_2 = min(samplingProportion_2$true)
x.max_samplingProportion_2 = max(samplingProportion_2$true)

y.min_samplingProportion_2 = min(samplingProportion_2$lower)
y.max_samplingProportion_2 = max(samplingProportion_2$upper)


lim.min = min(x.min_samplingProportion_2, y.min_samplingProportion_2)
lim.max = max(x.max_samplingProportion_2, y.max_samplingProportion_2)

p.samplingProportion_2 <- ggplot(samplingProportion_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.samplingProportion_2,paste("figures/samplingProportion_2", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_2_log <- ggplot(samplingProportion_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingProportion_2_log,paste("figures/samplingProportion_2_logScale", ".pdf", sep=""),width=6, height=5)

# sampling proportion 3
d3 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion_3[,(ncol(samplingProportion_3)-21+1):ncol(samplingProportion_3)], sum)/length(samplingProportion_3$test_0))

l <- length(samplingProportion_3$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd3 <- data.frame(p=cc, vals=unlist(samplingProportion_3[,(ncol(samplingProportion_3)-21+1):ncol(samplingProportion_3)]))


p.qq_samplingProportion_3 <- ggplot(d3)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") +
  theme_minimal()

p.qq_samplingProportion_3 <- p.qq_samplingProportion_3+
  stat_summary(data=dd3, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Sampling proportion 3")

ggsave(plot=p.qq_samplingProportion_3,paste("figures/qq_samplingProportion_3", ".pdf", sep=""),width=6, height=5)

x.min_samplingProportion_3 = min(samplingProportion_3$true)
x.max_samplingProportion_3 = max(samplingProportion_3$true)

y.min_samplingProportion_3 = min(samplingProportion_3$lower)
y.max_samplingProportion_3 = max(samplingProportion_3$upper)


lim.min = min(x.min_samplingProportion_3, y.min_samplingProportion_3)
lim.max = max(x.max_samplingProportion_3, y.max_samplingProportion_3)

p.samplingProportion_3 <- ggplot(samplingProportion_3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.samplingProportion_3,paste("figures/samplingProportion_3", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_3_log <- ggplot(samplingProportion_3)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingProportion_3_log,paste("figures/samplingProportion_3_logScale", ".pdf", sep=""),width=6, height=5)

# sampling proportion 4
d4 <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion_4[,(ncol(samplingProportion_4)-21+1):ncol(samplingProportion_4)], sum)/length(samplingProportion_4$test_0))

l <- length(samplingProportion_4$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd4 <- data.frame(p=cc, vals=unlist(samplingProportion_4[,(ncol(samplingProportion_4)-21+1):ncol(samplingProportion_4)]))


p.qq_samplingProportion_4 <- ggplot(d4)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("True recovery proportion") +
  theme_minimal()

p.qq_samplingProportion_4 <- p.qq_samplingProportion_4+
  stat_summary(data=dd4, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Sampling proportion 4")

ggsave(plot=p.qq_samplingProportion_4,paste("figures/qq_samplingProportion_4", ".pdf", sep=""),width=6, height=5)

x.min_samplingProportion_4 = min(samplingProportion_4$true)
x.max_samplingProportion_4 = max(samplingProportion_4$true)

y.min_samplingProportion_4 = min(samplingProportion_4$lower)
y.max_samplingProportion_4 = max(samplingProportion_4$upper)


lim.min = min(x.min_samplingProportion_4, y.min_samplingProportion_4)
lim.max = max(x.max_samplingProportion_4, y.max_samplingProportion_4)

p.samplingProportion_4 <- ggplot(samplingProportion_4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()

ggsave(plot=p.samplingProportion_4,paste("figures/samplingProportion_4", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_4_log <- ggplot(samplingProportion_4)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))

ggsave(plot=p.samplingProportion_4_log,paste("figures/samplingProportion_4_logScale", ".pdf", sep=""),width=6, height=5)


# combined sampling proportion plot

d_combined <- bind_rows(d1, d2, d3, d4) # 
dd_combined <- bind_rows(dd1,dd2, dd3, dd4) #  


p.qq_sampling_prop_combined <- ggplot(d_combined)+
  geom_abline(intercept = 0, color="#BBBBBB", linewidth=1)+
   xlab("Credible interval width") +
  ylab("True recovery proportion") + 
  theme_minimal()

legend_colors <- c("4" = "#4477AA", "1" = "#228833", "2" = "#EE6677", "3" = "#AA3377")

p.qq_sampling_prop_combined <- p.qq_sampling_prop_combined +
                 stat_summary(data=dd4, aes(x=(p+0.008), y=vals, colour="4"),fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
  stat_summary(data=dd1, aes(x=(p-0.008), y=vals, colour="1"),fun.data = "mean_cl_boot",  geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
                 stat_summary(data=dd2, aes(x=(p-0.0025), y=vals, colour="2"), fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1,  linewidth=0.65) +
                                stat_summary(data=dd3, aes(x=(p+0.003), y=vals, colour="3"),fun.data = "mean_cl_boot", geom="linerange",
               fun.args=list(conf.int = .95, B = 1000), alpha=1, linewidth=0.65) +
                 scale_color_manual(values = legend_colors, name = "Interval") + 
                 theme(legend.position = "inside", legend.position.inside=c(0.85, 0.25), legend.background=element_rect(fill="white", colour=NA))  # +
  
ggsave(plot=p.qq_sampling_prop_combined,paste("figures/qq_combined_sampling_proportion", ".pdf", sep=""),width=6, height=5)


######### tables #############

hpd.test = data.frame("Parameter"=c('diversificationRate_1', 'diversificationRate_2', 'diversificationRate_3', 'diversificationRate_4',
                                    'turnover_1', 'turnover_2', 'turnover_3', 'turnover_4',
                                    'samplingProportion_1', 'samplingProportion_2', 'samplingProportion_3', 'samplingProportion_4',
                                    'samplingAtPresentProb',
                                    'mrca',
                                    'origin'),
                      "HPD coverage"=c(mean(diversificationRate_1$test_95), mean(diversificationRate_2$test_95), 
                                      mean(diversificationRate_3$test_95), mean(diversificationRate_4$test_95),
                                       mean(turnover_1$test_95), mean(turnover_2$test_95),
                                       mean(turnover_3$test_95), mean(turnover_4$test_95),
                                       mean(samplingProportion_1$test_95),   mean(samplingProportion_2$test_95),
                                        mean(samplingProportion_3$test_95),   mean(samplingProportion_4$test_95),
                                       mean(samplingAtPresentProb$test_95),
                                       mean(mrca$test_95),
                                       mean(origin$test_95)))
write.csv(hpd.test, file = paste("figures/HPD_test",".csv", sep="" ))


hpd.width = data.frame("Parameter"=c('diversificationRate_1', 'diversificationRate_2', 'diversificationRate_3', 'diversificationRate_4',
                                    'turnover_1', 'turnover_2', 'turnover_3', 'turnover_4',
                                    'samplingProportion_1', 'samplingProportion_2', 'samplingProportion_3', 'samplingProportion_4',
                                    'samplingAtPresentProb',
                                    'mrca',
                                    'origin'),
                       "Average relative HPD width"=c(mean(diversificationRate_1$hpd.rel.width), mean(diversificationRate_2$hpd.rel.width),
                       mean(diversificationRate_3$hpd.rel.width), mean(diversificationRate_4$hpd.rel.width),
                                                      mean(turnover_1$hpd.rel.width), mean(turnover_2$hpd.rel.width),
                                                      mean(turnover_3$hpd.rel.width), mean(turnover_4$hpd.rel.width),
                                                      mean(samplingProportion_1$hpd.rel.width),  mean(samplingProportion_2$hpd.rel.width),
                                                      mean(samplingProportion_3$hpd.rel.width),  mean(samplingProportion_4$hpd.rel.width),
                                                      mean(samplingAtPresentProb$hpd.rel.width),
                                                      mean(mrca$hpd.rel.width),
                                                      mean(origin$hpd.rel.width)))
write.csv(hpd.width, file = paste("figures/HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('diversificationRate_1', 'diversificationRate_2', 'diversificationRate_3', 'diversificationRate_4',
                                    'turnover_1', 'turnover_2', 'turnover_3', 'turnover_4',
                                    'samplingProportion_1', 'samplingProportion_2', 'samplingProportion_3', 'samplingProportion_4',
                                    'samplingAtPresentProb',
                                    'mrca',
                                    'origin'),
                "Average Coefficient of Variation"=c(mean(diversificationRate_1$cv), mean(diversificationRate_2$cv), mean(diversificationRate_3$cv), mean(diversificationRate_4$cv),
                                                     mean(turnover_1$cv),  mean(turnover_2$cv), mean(turnover_3$cv),  mean(turnover_4$cv),
                                                     mean(samplingProportion_1$cv), mean(samplingProportion_2$cv), mean(samplingProportion_3$cv), mean(samplingProportion_4$cv),
                                                     mean(samplingAtPresentProb$cv),
                                                     mean(mrca$cv),
                                                     mean(origin$cv)))
write.csv(cv, file = paste("figures/CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('diversificationRate_1', 'diversificationRate_2', 'diversificationRate_3', 'diversificationRate_4',
                                    'turnover_1', 'turnover_2', 'turnover_3', 'turnover_4',
                                    'samplingProportion_1', 'samplingProportion_2', 'samplingProportion_3', 'samplingProportion_4',
                                    'samplingAtPresentProb',
                                    'mrca',
                                    'origin'),
                     "Relative Error at Medians"=c(mean(diversificationRate_1$rel.err.median), mean(diversificationRate_2$rel.err.median),
                     mean(diversificationRate_3$rel.err.median), mean(diversificationRate_4$rel.err.median),
                                 mean(turnover_1$rel.err.median), mean(turnover_2$rel.err.median), mean(turnover_3$rel.err.median), mean(turnover_4$rel.err.median),
                                 mean(samplingProportion_1$rel.err.median), mean(samplingProportion_2$rel.err.median),
                                  mean(samplingProportion_3$rel.err.median), mean(samplingProportion_4$rel.err.median),
                                 mean(samplingAtPresentProb$rel.err.median),
                                 mean(mrca$rel.err.median),
                                 mean(origin$rel.err.median)))
write.csv(medians, file = paste("figures/rel_error_median",".csv", sep="" ))

means = data.frame("Parameter"=c('diversificationRate_1', 'diversificationRate_2', 'diversificationRate_3', 'diversificationRate_4',
                                    'turnover_1', 'turnover_2', 'turnover_3', 'turnover_4',
                                    'samplingProportion_1', 'samplingProportion_2', 'samplingProportion_3', 'samplingProportion_4',
                                    'samplingAtPresentProb',
                                    'mrca',
                                    'origin'),
                   "Relative Error at Means"=c(mean(diversificationRate_1$rel.err.mean), mean(diversificationRate_2$rel.err.mean),
                   mean(diversificationRate_3$rel.err.mean), mean(diversificationRate_4$rel.err.mean),
                                                 mean(turnover_1$rel.err.mean), mean(turnover_2$rel.err.mean),
                                                  mean(turnover_3$rel.err.mean), mean(turnover_4$rel.err.mean),
                                                 mean(samplingProportion_1$rel.err.mean), mean(samplingProportion_2$rel.err.mean),
                                                  mean(samplingProportion_3$rel.err.mean), mean(samplingProportion_4$rel.err.mean),
                                                 mean(samplingAtPresentProb$rel.err.mean),
                                                 mean(mrca$rel.err.mean),
                                                 mean(origin$rel.err.mean)))
write.csv(means, file = paste("figures/rel_error_mean",".csv", sep=""))

