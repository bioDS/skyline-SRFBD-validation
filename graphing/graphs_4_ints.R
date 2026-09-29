######################################################
######################################################
# 
######################################################
######################################################
library(ggplot2)
library(coda)
library(data.table)
options(digits=7)
library("Hmisc")


# clear workspace
rm(list = ls())

valid = 0
removed = 0

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
#wd<-"/Users/jugne/Documents/Source/beast2.7/sRanges-material/validation/estimate_all_params_priors/"
setwd("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/graphing")
dir.create("figures")

# system("bash renamelogs.sh")
# system("bash movelogs.sh")

# read in the true rates
true_rates_file <- "../adjusted_params"
true.rates <- fread(true_rates_file, sep=" ", fill = TRUE)
######################################################
######################################################
# 
######################################################
######################################################
library(ggplot2)
library(coda)
library(data.table)
options(digits=7)

n<-200
remove_rows<-c()
ints = 4

diversificationRate = data.frame(true=numeric(ints*n), estimated=numeric(ints*n),
                                 upper=numeric(ints*n), lower=numeric(ints*n),
                                 rel.err.median=numeric(ints*n), rel.err.mean=numeric(ints*n),
                                 rel.err.mode=numeric(ints*n), cv=numeric(ints*n),
                                 hpd.lower=numeric(ints*n), hpd.upper=numeric(ints*n),
                                 hpd.rel.width=numeric(ints*n), test_0=numeric(ints*n), test_5=numeric(ints*n), test_10=numeric(ints*n),
                                 test_15=numeric(ints*n), test_20=numeric(ints*n), test_25=numeric(ints*n),
                                 test_30=numeric(ints*n), test_35=numeric(ints*n), test_l0=numeric(ints*n),
                                 test_l5=numeric(ints*n), test_50=numeric(ints*n), test_55=numeric(ints*n),
                                 test_60=numeric(ints*n), test_65=numeric(ints*n), test_70=numeric(ints*n),
                                 test_75=numeric(ints*n), test_80=numeric(ints*n), test_85=numeric(ints*n),
                                 test_90=numeric(ints*n), test_95=numeric(ints*n),
                                 test_100=numeric(ints*n))

turnover = data.frame(true=numeric(ints*n), estimated=numeric(ints*n), upper=numeric(ints*n),
                      lower=numeric(ints*n), rel.err.median=numeric(ints*n),
                      rel.err.mean=numeric(ints*n), rel.err.mode=numeric(ints*n),
                      cv=numeric(ints*n), hpd.lower=numeric(ints*n), hpd.upper=numeric(ints*n),
                      hpd.rel.width=numeric(ints*n), test_0=numeric(ints*n), test_5=numeric(ints*n), test_10=numeric(ints*n),
                      test_15=numeric(ints*n), test_20=numeric(ints*n), test_25=numeric(ints*n),
                      test_30=numeric(ints*n), test_35=numeric(ints*n), test_40=numeric(ints*n),
                      test_45=numeric(ints*n), test_50=numeric(ints*n), test_55=numeric(ints*n),
                      test_60=numeric(ints*n), test_65=numeric(ints*n), test_70=numeric(ints*n),
                      test_75=numeric(ints*n), test_80=numeric(ints*n), test_85=numeric(ints*n),
                      test_90=numeric(ints*n), test_95=numeric(ints*n), test_100=numeric(ints*n))

samplingProportion = data.frame(true=numeric(ints*n), estimated=numeric(ints*n),
                                upper=numeric(ints*n), lower=numeric(ints*n),
                                rel.err.median=numeric(ints*n), rel.err.mean=numeric(ints*n),
                                rel.err.mode=numeric(ints*n), cv=numeric(ints*n),
                                hpd.lower=numeric(ints*n), hpd.upper=numeric(ints*n),
                                hpd.rel.width=numeric(ints*n), test_0=numeric(ints*n), test_5=numeric(ints*n), test_10=numeric(ints*n),
                                test_15=numeric(ints*n), test_20=numeric(ints*n), test_25=numeric(ints*n),
                                test_30=numeric(ints*n), test_35=numeric(ints*n), test_40=numeric(ints*n),
                                test_45=numeric(ints*n), test_50=numeric(ints*n), test_55=numeric(ints*n),
                                test_60=numeric(ints*n), test_65=numeric(ints*n), test_70=numeric(ints*n),
                                test_75=numeric(ints*n), test_80=numeric(ints*n), test_85=numeric(ints*n),
                                test_90=numeric(ints*n), test_95=numeric(ints*n), test_100=numeric(ints*n))


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

rsamplingAtPresentProb = data.frame(true=numeric(n), estimated=numeric(n),
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
  #  cat("i =",i, " ")
  if (!(file.exists(paste0("../logs/", (i-1), ".log")))){
    remove_rows = cbind(remove_rows, i)
    cat("FILE DOES NOT EXIST ", (i-1), "\n")
    next
  }
  t=read.table(paste0("../logs/", (i-1), ".log"), header=TRUE, sep="\t")
  
  
  used.rates = NULL
  used.rates$birthRate = true.rates$birth[(ints*i-(ints-1)):(ints*i)]
  used.rates$deathRate = true.rates$death[(ints*i-(ints-1)):(ints*i)]
  used.rates$samplingRate = true.rates$sampling[(ints*i-(ints-1)):(ints*i)]
  used.rates$rho = true.rates$rho[(ints*i)]
  used.rates$div_rate = true.rates$d[(ints*i-(ints-1)):(ints*i)]
  used.rates$turnover = true.rates$v[(ints*i-(ints-1)):(ints*i)]
  used.rates$sampling_prop = true.rates$s[(ints*i-(ints-1)):(ints*i)]
  used.rates$mrca = true.rates$mrca[(ints*i)]
  used.rates$origin = true.rates$origin[(ints*i)]
  
  
  # take a 10% burnin
  t <- t[-seq(1,ceiling(length(t$deathRate)/10)), ]
  
  if (nrow(t)==0){
    remove_rows = cbind(remove_rows, i)
    next
  }
  ess<-  effectiveSize(as.mcmc(t))
  post_ess = as.numeric(ess["posterior"])
  
  if(post_ess<100){
    remove_rows = cbind(remove_rows, i)
    cat("post ess for ", (i-1), " is too small (", post_ess, ")\n")
    next
  }
  valid = valid + 1
  origin$true[i] <- used.rates$origin
  mrca$true[i] <- used.rates$mrca
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


  for (j in 1:ints){
    diversificationRate$true[ints*i+j-ints] <- used.rates$div_rate[j]
    turnover$true[ints*i+j-ints] <- used.rates$turnover[j]
    samplingProportion$true[ints*i+j-ints] <- used.rates$sampling_prop[j]
    
    
    divcol = 4 + j
    turncol = 4 + j + ints 
    sampcol = 4 + j + 2*ints
    
    diversificationRate$estimated[ints*i+j-ints] <- median(t[[divcol]])
    diversificationRate$upper[ints*i+j-ints] <- quantile(t[[divcol]],0.975)
    diversificationRate$lower[ints*i+j-ints] <- quantile(t[[divcol]],0.025)
    
    
    diversificationRate$rel.err.median[ints*i+j-ints] <- abs(median(t[[divcol]])-used.rates$div_rate[j])/used.rates$div_rate[j]
    diversificationRate$rel.err.mean[ints*i+j-ints] <-abs(mean(t[[divcol]])-used.rates$div_rate[j])/used.rates$div_rate[j]
    diversificationRate$rel.err.mode[ints*i+j-ints] <- abs(Mode(t[[divcol]])-used.rates$div_rate[j])/used.rates$div_rate[j]
    diversificationRate$cv[ints*i+j-ints] <- sqrt(exp(sd(log(t[[divcol]]))**2)-1)
    diversificationRate$hpd.lower[ints*i+j-ints] <- HPDinterval(as.mcmc(t[[divcol]]))[1]
    diversificationRate$hpd.upper[ints*i+j-ints] <- HPDinterval(as.mcmc(t[[divcol]]))[2]
    diversificationRate$hpd.rel.width[ints*i+j-ints] <- (diversificationRate$hpd.upper[ints*i+j-ints]-diversificationRate$hpd.lower[ints*i+j-ints])/used.rates$div_rate[j]
    diversificationRate[(ints*i+j-ints),(ncol(diversificationRate)-21+1):ncol(diversificationRate)] <- lapply(seq(0.,1,0.05), test_vals, diversificationRate$true[ints*i+j-ints], t[[divcol]])
    
    
    turnover$estimated[ints*i+j-ints] <- median(t[[turncol]])
    turnover$upper[ints*i+j-ints] <- quantile(t[[turncol]],0.975)
    turnover$lower[ints*i+j-ints] <- quantile(t[[turncol]],0.025)
    
    
    turnover$rel.err.median[ints*i+j-ints] <- abs(median(t[[turncol]])-used.rates$turnover[j])/used.rates$turnover[j]
    turnover$rel.err.mean[ints*i+j-ints] <-abs(mean(t[[turncol]])-used.rates$turnover[j])/used.rates$turnover[j]
    turnover$rel.err.mode[ints*i+j-ints] <- abs(Mode(t[[turncol]])-used.rates$turnover[j])/used.rates$turnover[j]
    turnover$cv[ints*i+j-ints] <- sqrt(exp(sd(log(t[[turncol]]))**2)-1)
    turnover$hpd.lower[ints*i+j-ints] <- HPDinterval(as.mcmc(t[[turncol]]))[1]
    turnover$hpd.upper[ints*i+j-ints] <- HPDinterval(as.mcmc(t[[turncol]]))[2]
    turnover$hpd.rel.width[ints*i+j-ints] <- (turnover$hpd.upper[ints*i+j-ints]-turnover$hpd.lower[ints*i+j-ints])/used.rates$turnover[j]
    turnover[(ints*i+j-ints),(ncol(turnover)-21+1):ncol(turnover)] <- lapply(seq(0.,1,0.05), test_vals, turnover$true[ints*i+j-ints], t[[turncol]])
    
    samplingProportion$estimated[ints*i+j-ints] <- median(t[[sampcol]])
    samplingProportion$upper[ints*i+j-ints] <- quantile(t[[sampcol]],0.975)
    samplingProportion$lower[ints*i+j-ints] <- quantile(t[[sampcol]],0.025)
    
    
    samplingProportion$rel.err.median[ints*i+j-ints] <- abs(median(t[[sampcol]])-used.rates$sampling_prop[j])/used.rates$sampling_prop[j]
    samplingProportion$rel.err.mean[ints*i+j-ints] <-abs(mean(t[[sampcol]])-used.rates$sampling_prop[j])/used.rates$sampling_prop[j]
    samplingProportion$rel.err.mode[ints*i+j-ints] <- abs(Mode(t[[sampcol]])-used.rates$sampling_prop[j])/used.rates$sampling_prop[j]
    samplingProportion$cv[ints*i+j-ints] <- sqrt(exp(sd(log(t[[sampcol]]))**2)-1)
    samplingProportion$hpd.lower[ints*i+j-ints] <- HPDinterval(as.mcmc(t[[sampcol]]))[1]
    samplingProportion$hpd.upper[ints*i+j-ints] <- HPDinterval(as.mcmc(t[[sampcol]]))[2]
    samplingProportion$hpd.rel.width[ints*i+j-ints] <- (samplingProportion$hpd.upper[ints*i+j-ints]-samplingProportion$hpd.lower[ints*i+j-ints])/used.rates$sampling_prop[j]
    
    samplingProportion[(ints*i+j-ints),(ncol(samplingProportion)-21+1):ncol(samplingProportion)] <- lapply(seq(0.,1,0.05), test_vals, samplingProportion$true[(ints*i+j-ints)], t[[sampcol]])
    
  }
}

save.image(file="figures/validationPlotEnvironment.RData")

if (!is.null(remove_rows)){
  extended_rows = c(remove_rows*ints - 1, remove_rows*ints)
  diversificationRate <- diversificationRate[-(extended_rows), ]
  turnover <- turnover[-(extended_rows), ]
  samplingProportion <- samplingProportion[-(extended_rows), ]
  mrca <- mrca[-remove_rows, ]
  origin <- origin[-remove_rows, ]
  samplingAtPresentProb <- samplingAtPresentProb[-remove_rows, ]
}

#DIVERSIFICATION 1

l <- length(diversificationRate$test_0)/ints
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}

indices = c(TRUE, FALSE, FALSE, FALSE)

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)], sum)/(length(diversificationRate$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)]))

p.qq_div_rate_1 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_div_rate_1 <- p.qq_div_rate_1 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 1")



x.min_diversificationRate = min(diversificationRate$true[indices])
x.max_diversificationRate = max(diversificationRate$true[indices])

y.min_diversificationRate = min(diversificationRate$lower[indices])
y.max_diversificationRate = max(diversificationRate$upper[indices])


lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate1 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("Diversification rate 1")
ggsave(plot=p.diversificationRate1,paste("figures/diversificationRate1", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_log1 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Diversification rate 1")

p.qq_div_rate_1
p.diversificationRate1
p.diversificationRate_log1
ggsave(plot=p.diversificationRate1,paste("figures/diversificationRate1", ".pdf", sep=""),width=6, height=5)
ggsave(plot=p.qq_div_rate_1,paste("figures/qq_div1", ".pdf", sep=""),width=6, height=5)

ggsave(plot=p.diversificationRate_log1,paste("figures/diversificationRate1_logScale", ".pdf", sep=""),width=6, height=5)


#DIVERSIFICATION 2

indices = c(FALSE, TRUE, FALSE, FALSE)



d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)], sum)/(length(diversificationRate$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)]))

p.qq_div_rate_2 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_div_rate_2 <- p.qq_div_rate_2 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 2")


x.min_diversificationRate = min(diversificationRate$true[indices])
x.max_diversificationRate = max(diversificationRate$true[indices])

y.min_diversificationRate = min(diversificationRate$lower[indices])
y.max_diversificationRate = max(diversificationRate$upper[indices])


lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate2 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  ggtitle("Diversification rate 2")
ggsave(plot=p.diversificationRate2,paste("figures/diversificationRate2", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_log2 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Diversification rate 2")

ggsave(plot=p.qq_div_rate_2,paste("figures/qq_div2", ".pdf", sep=""),width=6, height=5)
p.qq_div_rate_2
p.diversificationRate2
p.diversificationRate_log2
ggsave(plot=p.diversificationRate2,paste("figures/diversificationRate2", ".pdf", sep=""),width=6, height=5)

ggsave(plot=p.diversificationRate_log2,paste("figures/diversificationRate2_logScale", ".pdf", sep=""),width=6, height=5)


#DIVERSIFICATION 3

indices = c(FALSE, FALSE, TRUE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)], sum)/(length(diversificationRate$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)]))

p.qq_div_rate_3 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_div_rate_3 <- p.qq_div_rate_3 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 3")


x.min_diversificationRate = min(diversificationRate$true[indices])
x.max_diversificationRate = max(diversificationRate$true[indices])

y.min_diversificationRate = min(diversificationRate$lower[indices])
y.max_diversificationRate = max(diversificationRate$upper[indices])


lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate3 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  ggtitle("Diversification rate 3")
ggsave(plot=p.diversificationRate3,paste("figures/diversificationRate3", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_log3 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Diversification rate 3")

ggsave(plot=p.qq_div_rate_3,paste("figures/qq_div3", ".pdf", sep=""),width=6, height=5)
p.qq_div_rate_3
p.diversificationRate3
p.diversificationRate_log3
ggsave(plot=p.diversificationRate3,paste("figures/diversificationRate3", ".pdf", sep=""),width=6, height=5)

ggsave(plot=p.diversificationRate_log3,paste("figures/diversificationRate3_logScale", ".pdf", sep=""),width=6, height=5)


#DIVERSIFICATION 4

indices = c(FALSE, FALSE, FALSE, TRUE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)], sum)/(length(diversificationRate$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(diversificationRate[indices,(ncol(diversificationRate)-21+1):ncol(diversificationRate)]))

p.qq_div_rate_4 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_div_rate_4 <- p.qq_div_rate_4 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Diversification rate 4")


x.min_diversificationRate = min(diversificationRate$true[indices])
x.max_diversificationRate = max(diversificationRate$true[indices])

y.min_diversificationRate = min(diversificationRate$lower[indices])
y.max_diversificationRate = max(diversificationRate$upper[indices])


lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
lim.max = max(x.max_diversificationRate, y.max_diversificationRate)

p.diversificationRate4 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  ggtitle("Diversification rate 4")
ggsave(plot=p.diversificationRate4,paste("figures/diversificationRate4", ".pdf", sep=""),width=6, height=5)

p.diversificationRate_log4 <- ggplot(diversificationRate[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Diversification rate 4")

ggsave(plot=p.qq_div_rate_4,paste("figures/qq_div3", ".pdf", sep=""),width=6, height=5)
p.qq_div_rate_4
p.diversificationRate4
p.diversificationRate_log4
ggsave(plot=p.diversificationRate4,paste("figures/diversificationRate4", ".pdf", sep=""),width=6, height=5)

ggsave(plot=p.diversificationRate_log4,paste("figures/diversificationRate4_logScale", ".pdf", sep=""),width=6, height=5)

#unlist(c((samplingProportion[,(ncol(samplingProportion)-21+1):ncol(samplingProportion)])))
# d <- data.frame(p = seq(0.,1,0.05), 
#                 vals=sapply(diversificationRate[,(ncol(diversificationRate)-21+1):ncol(diversificationRate)], sum)/length(diversificationRate$test_0))
# 
# l <- length(diversificationRate$test_0)
# cc <- c()
# for (x in seq(0.00, 1.00, 0.05)){
#   cc <- c(cc, rep(x, l))
# }
# dd <- data.frame(p=cc, vals=unlist(diversificationRate[,(ncol(diversificationRate)-21+1):ncol(diversificationRate)]))
# 
# p.qq_div_rate <- ggplot(d)+
#   geom_abline(intercept = 0, color="red")+
#   geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
#   ylab("Truth recovery proportion") + 
#   theme_minimal()
# p.qq_div_rate <- p.qq_div_rate+
#   stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
#                fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
#   ggtitle("Diversification rate")
# 
# # p.qq_div_rate <- p.qq_div_rate+
# #   stat_summary(data=dd, aes(x=p, y=vals), fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
# #   ggtitle("Diversification rate")
# 
# ggsave(plot=p.qq_div_rate,paste("figures/qq_diversificationRate", ".pdf", sep=""),width=6, height=5)
# 
# 
# x.min_diversificationRate = min(diversificationRate$true)
# x.max_diversificationRate = max(diversificationRate$true)
# 
# y.min_diversificationRate = min(diversificationRate$lower)
# y.max_diversificationRate = max(diversificationRate$upper)
# 
# 
# lim.min = min(x.min_diversificationRate, y.min_diversificationRate)
# lim.max = max(x.max_diversificationRate, y.max_diversificationRate)
# 
# p.diversificationRate <- ggplot(diversificationRate)+
#   geom_abline(intercept = 0, color="red", linetype="dashed")+
#   geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
#   geom_point(aes(x=true, y=estimated), size=2) +
#   theme_minimal()
# # scale_y_log10(limits=c(lim.min, lim.max)) +
# # scale_x_log10(limits=c(lim.min, lim.max))
# 
# 
# ggsave(plot=p.diversificationRate,paste("figures/diversificationRate", ".pdf", sep=""),width=6, height=5)
# 
# p.diversificationRate_log <- ggplot(diversificationRate)+
#   geom_abline(intercept = 0, color="red", linetype="dashed")+
#   geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
#   geom_point(aes(x=true, y=estimated), size=2) +
#   theme_minimal() +
#   scale_y_log10(limits=c(lim.min, lim.max)) +
#   scale_x_log10(limits=c(lim.min, lim.max))
# 
# ggsave(plot=p.diversificationRate_log,paste("figures/diversificationRate_logScale", ".pdf", sep=""),width=6, height=5)

# dt_rel.err.median.diversificationRate = data.frame(param='diversificationRate', value=diversificationRate$rel.err.median)
# p_rel_err_diversificationRate<- ggplot(dt_rel.err.median.diversificationRate, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.median.diversificationRate$value)$stats[c(1, 5)]*1.05)
# ggsave(plot=p_rel_err_diversificationRate,paste("figures/diversificationRate_rel_err_median", ".pdf", sep=""),width=6, height=5)
# 
# dt_rel.err.mean.diversificationRate = data.frame(param='diversificationRate', value=diversificationRate$rel.err.mean)
# p_rel_err_mean_diversificationRate<- ggplot(dt_rel.err.mean.diversificationRate, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_rel.err.mean.diversificationRate$value)$stats[c(1, 5)]*1.05)
# ggsave(plot=p_rel_err_mean_diversificationRate,paste("figures/diversificationRate_rel_err_mean", ".pdf", sep=""),width=6, height=5)
# 
# dt_hpd.rel.width.diversificationRate = data.frame(param='diversificationRate', value=diversificationRate$hpd.rel.width)
# p_hpd.rel.width.diversificationRate <- ggplot(dt_hpd.rel.width.diversificationRate, aes(x=param, y=value)) +   geom_boxplot() + coord_cartesian(ylim = boxplot.stats(dt_hpd.rel.width.diversificationRate$value)$stats[c(1, 5)]*1.05)
# ggsave(plot=p_hpd.rel.width.diversificationRate,paste("figures/diversificationRate_rel_hpd_width", ".pdf", sep=""),width=6, height=5)

########## turnover ################


#Turnover 1

l <- length(turnover$test_0)/ints
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}

indices = c(TRUE, FALSE, FALSE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover[indices,(ncol(turnover)-21+1):ncol(turnover)], sum)/(length(turnover$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(turnover[indices,(ncol(turnover)-21+1):ncol(turnover)]))

p.qq_turnover_1 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_turnover_1 <- p.qq_turnover_1 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("Turnover 1")

ggsave(plot=p.qq_turnover_1,paste("figures/qq_turnover1", ".pdf", sep=""),width=6, height=5)



x.min_turnover = min(turnover$true[indices])
x.max_turnover = max(turnover$true[indices])

y.min_turnover = min(turnover$lower[indices])
y.max_turnover = max(turnover$upper[indices])


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover1 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("Turnover 1")
ggsave(plot=p.turnover1,paste("figures/turnover1", ".pdf", sep=""),width=6, height=5)

p.turnover_log1 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Turnover rate 1")
ggsave(plot=p.turnover_log1,paste("figures/p.turnover_log1", ".pdf", sep=""),width=6, height=5)

p.qq_turnover_1
p.turnover1
p.turnover_log1

# TURNOVER 2

d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover[c(FALSE, TRUE),(ncol(turnover)-21+1):ncol(turnover)], sum)/(length(turnover$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(turnover[c(FALSE, TRUE),(ncol(turnover)-21+1):ncol(turnover)]))

p.qq_turnover_2 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_turnover_2 <- p.qq_turnover_2 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue")  +
  ggtitle("Turnover 2")

ggsave(plot=p.qq_turnover_2,paste("figures/qq_turnover2", ".pdf", sep=""),width=6, height=5)


indices = c(FALSE, TRUE, FALSE)

x.min_turnover = min(turnover$true[indices])
x.max_turnover = max(turnover$true[indices])

y.min_turnover = min(turnover$lower[indices])
y.max_turnover = max(turnover$upper[indices])


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover2 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("Turnover 2")
ggsave(plot=p.turnover2,paste("figures/turnover2", ".pdf", sep=""),width=6, height=5)

p.turnover_log2 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Turnover 2")
ggsave(plot=p.turnover_log2,paste("figures/p.turnover_log2", ".pdf", sep=""),width=6, height=5)

p.qq_turnover_2
p.turnover2
p.turnover_log2

# TURNOVER 3

indices = c(FALSE, FALSE, TRUE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover[indices,(ncol(turnover)-21+1):ncol(turnover)], sum)/(length(turnover$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(turnover[indices,(ncol(turnover)-21+1):ncol(turnover)]))

p.qq_turnover_3 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_turnover_3 <- p.qq_turnover_3 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue")  +
  ggtitle("Turnover 3")

ggsave(plot=p.qq_turnover_3,paste("figures/qq_turnover3", ".pdf", sep=""),width=6, height=5)



x.min_turnover = min(turnover$true[indices])
x.max_turnover = max(turnover$true[indices])

y.min_turnover = min(turnover$lower[indices])
y.max_turnover = max(turnover$upper[indices])


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover3 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("Turnover 3")
ggsave(plot=p.turnover3,paste("figures/turnover3", ".pdf", sep=""),width=6, height=5)

p.turnover_log3 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Turnover 3")
ggsave(plot=p.turnover_log3,paste("figures/p.turnover_log2", ".pdf", sep=""),width=6, height=5)

p.qq_turnover_3
p.turnover3
p.turnover_log3

# TURNOVER 4

indices = c(FALSE, FALSE, FALSE, TRUE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(turnover[indices,(ncol(turnover)-21+1):ncol(turnover)], sum)/(length(turnover$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(turnover[indices,(ncol(turnover)-21+1):ncol(turnover)]))

p.qq_turnover_4 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_turnover_4 <- p.qq_turnover_4 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue")  +
  ggtitle("Turnover 4")

ggsave(plot=p.qq_turnover_4,paste("figures/qq_turnover4", ".pdf", sep=""),width=6, height=5)



x.min_turnover = min(turnover$true[indices])
x.max_turnover = max(turnover$true[indices])

y.min_turnover = min(turnover$lower[indices])
y.max_turnover = max(turnover$upper[indices])


lim.min = min(x.min_turnover, y.min_turnover)
lim.max = max(x.max_turnover, y.max_turnover)

p.turnover4 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("Turnover 4")
ggsave(plot=p.turnover4,paste("figures/turnover3", ".pdf", sep=""),width=6, height=5)

p.turnover_log4 <- ggplot(turnover[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("Turnover 4")
ggsave(plot=p.turnover_log4,paste("figures/p.turnover_log2", ".pdf", sep=""),width=6, height=5)

p.qq_turnover_4
p.turnover4
p.turnover_log4



########## mrca ################
# 
d <- data.frame(p = seq(0.,1,0.05),
                vals=sapply(mrca[,(ncol(mrca)-21+1):ncol(mrca)], sum)/length(mrca$test_0))


l <- length(mrca$test_0)
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}
dd <- data.frame(p=cc, vals=unlist(mrca[,(ncol(mrca)-21+1):ncol(mrca)]))


p.qq_mrca <- ggplot(d)+
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") +
  theme_minimal()

p.qq_mrca <- p.qq_mrca+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Origin")

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
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


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
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_origin <- p.qq_origin+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("Origin")

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
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


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
  geom_abline(aes(xmin=0, ymin=0), intercept = 0, color="red", linetype="dashed")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") +
  theme_minimal()

p.qq_samplingAtPresentProb <- p.qq_samplingAtPresentProb+
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot", colour = "blue") +
  ggtitle("samplingAtPresentProb ")

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
# scale_y_log10(limits=c(lim.min, lim.max)) +
# scale_x_log10(limits=c(lim.min, lim.max))


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

#sampling proportion 1

l <- length(samplingProportion$test_0)/ints
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l))
}

indices = c(TRUE, FALSE, FALSE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)], sum)/(length(samplingProportion$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)]))

p.qq_samplingProportion_1 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_samplingProportion_1 <- p.qq_samplingProportion_1 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("sampling proportion 1")

ggsave(plot=p.qq_samplingProportion_1,paste("figures/qq_samplingProportion1", ".pdf", sep=""),width=6, height=5)



x.min_samplingProportion = min(samplingProportion$true[indices])
x.max_samplingProportion = max(samplingProportion$true[indices])

y.min_samplingProportion = min(samplingProportion$lower[indices])
y.max_samplingProportion = max(samplingProportion$upper[indices])


lim.min = min(x.min_samplingProportion, y.min_samplingProportion)
lim.max = max(x.max_samplingProportion, y.max_samplingProportion)

p.samplingProportion1 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("sampling proportion 1")
ggsave(plot=p.samplingProportion1,paste("figures/samplingProportion1", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_log1 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("sampling proportion rate 1")
ggsave(plot=p.samplingProportion_log1,paste("figures/p.samplingProportion_log1", ".pdf", sep=""),width=6, height=5)

p.qq_samplingProportion_1
p.samplingProportion1
p.samplingProportion_log1

# sampling proportion 2

l2 <- length(samplingProportion$test_0)/ints
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l2))
}

indices = c(FALSE, TRUE, FALSE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)], sum)/(length(samplingProportion$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)]))

p.qq_samplingProportion_2 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_samplingProportion_2 <- p.qq_samplingProportion_2 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("sampling proportion 2")

ggsave(plot=p.qq_samplingProportion_2,paste("figures/qq_samplingProportion2", ".pdf", sep=""),width=6, height=5)



x.min_samplingProportion = min(samplingProportion$true[indices])
x.max_samplingProportion = max(samplingProportion$true[indices])

y.min_samplingProportion = min(samplingProportion$lower[indices])
y.max_samplingProportion = max(samplingProportion$upper[indices])


lim.min = min(x.min_samplingProportion, y.min_samplingProportion)
lim.max = max(x.max_samplingProportion, y.max_samplingProportion)

p.samplingProportion2 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("sampling proportion 2")
ggsave(plot=p.samplingProportion2,paste("figures/samplingProportion2", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_log2 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("sampling proportion rate 2")
ggsave(plot=p.samplingProportion_log2,paste("figures/p.samplingProportion_log2", ".pdf", sep=""),width=6, height=5)

p.qq_samplingProportion_2
p.samplingProportion2
p.samplingProportion_log2

# sampling proportion 3

l2 <- length(samplingProportion$test_0)/ints
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l2))
}

indices = c(FALSE, FALSE, TRUE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)], sum)/(length(samplingProportion$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)]))

p.qq_samplingProportion_3 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_samplingProportion_3 <- p.qq_samplingProportion_2 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("sampling proportion 2")

ggsave(plot=p.qq_samplingProportion_3,paste("figures/qq_samplingProportion3", ".pdf", sep=""),width=6, height=5)



x.min_samplingProportion = min(samplingProportion$true[indices])
x.max_samplingProportion = max(samplingProportion$true[indices])

y.min_samplingProportion = min(samplingProportion$lower[indices])
y.max_samplingProportion = max(samplingProportion$upper[indices])


lim.min = min(x.min_samplingProportion, y.min_samplingProportion)
lim.max = max(x.max_samplingProportion, y.max_samplingProportion)

p.samplingProportion3 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("sampling proportion 3")
ggsave(plot=p.samplingProportion3,paste("figures/samplingProportion3", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_log3 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("sampling proportion rate 3")
ggsave(plot=p.samplingProportion_log3,paste("figures/p.samplingProportion_log3", ".pdf", sep=""),width=6, height=5)

p.qq_samplingProportion_3
p.samplingProportion3
p.samplingProportion_log3
# sampling proportion 3

l2 <- length(samplingProportion$test_0)/ints
cc <- c()
for (x in seq(0.00, 1.00, 0.05)){
  cc <- c(cc, rep(x, l2))
}

indices = c(FALSE, FALSE, TRUE, FALSE)


d <- data.frame(p = seq(0.,1,0.05), 
                vals=sapply(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)], sum)/(length(samplingProportion$test_0)/ints))
dd <- data.frame(p=cc, vals=unlist(samplingProportion[indices,(ncol(samplingProportion)-21+1):ncol(samplingProportion)]))

p.qq_samplingProportion_4 <- ggplot(d)+
  geom_abline(intercept = 0, color="red")+
  geom_point(aes(x=p, y=vals), size=2) + xlab("Credible interval width") +
  ylab("Truth recovery proportion") + 
  theme_minimal()
p.qq_samplingProportion_4 <- p.qq_samplingProportion_4 +
  stat_summary(data=dd, aes(x=p, y=vals),fun.data = "mean_cl_boot",
               fun.args=list(conf.int = .95, B = 1000), colour = "blue") +
  ggtitle("sampling proportion 4")

ggsave(plot=p.qq_samplingProportion_4,paste("figures/qq_samplingProportion4", ".pdf", sep=""),width=6, height=5)



x.min_samplingProportion = min(samplingProportion$true[indices])
x.max_samplingProportion = max(samplingProportion$true[indices])

y.min_samplingProportion = min(samplingProportion$lower[indices])
y.max_samplingProportion = max(samplingProportion$upper[indices])


lim.min = min(x.min_samplingProportion, y.min_samplingProportion)
lim.max = max(x.max_samplingProportion, y.max_samplingProportion)

p.samplingProportion4 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal()+
  ggtitle("sampling proportion 4")
ggsave(plot=p.samplingProportion4,paste("figures/samplingProportion4", ".pdf", sep=""),width=6, height=5)

p.samplingProportion_log4 <- ggplot(samplingProportion[indices,])+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) +
  theme_minimal() +
  scale_y_log10(limits=c(lim.min, lim.max)) +
  scale_x_log10(limits=c(lim.min, lim.max))+
  ggtitle("sampling proportion rate 4")
ggsave(plot=p.samplingProportion_log4,paste("figures/p.samplingProportion_log3", ".pdf", sep=""),width=6, height=5)

p.qq_samplingProportion_4
p.samplingProportion4
p.samplingProportion_log4


######### tables #############

hpd.test = data.frame("Parameter"=c('diversificationRate1', 'diversificationRate2', 'diversificationRate3', 'diversificationRate4', 
                                    'turnover1', 'turnover2', 'turnover3', 'turnover4', 
                                    'samplingProportion1', 'samplingProportion2', 'samplingProportion3', 'samplingProportion4', 
                                    'samplingAtPresentProb', 
                                    'mrca',
                                    'origin'),
                      "HPD coverage"=c(mean(diversificationRate$test_95[c(TRUE, FALSE,FALSE, FALSE)]), mean(diversificationRate$test_95[c(FALSE, TRUE, FALSE, FALSE)]), 
                      mean(diversificationRate$test_95[c(FALSE, FALSE, TRUE, FALSE)]), mean(diversificationRate$test_95[c(FALSE, FALSE, FALSE, TRUE)]), 
                                       mean(turnover$test_95[c(TRUE, FALSE, FALSE, FALSE)]), mean(turnover$test_95[c(FALSE, TRUE, FALSE, FALSE)]), 
                                       mean(turnover$test_95[c(FALSE, FALSE, TRUE, FALSE)]), mean(turnover$test_95[c(FALSE, FALSE, FALSE, TRUE)]), 
                                       mean(samplingProportion$test_95[c(TRUE, FALSE, FALSE, FALSE)]), mean(samplingProportion$test_95[c(FALSE, TRUE, FALSE, FALSE)]), 
                                       mean(samplingProportion$test_95[c(FALSE, FALSE, TRUE, FALSE)]), mean(samplingProportion$test_95[c(FALSE, FALSE, FALSE, TRUE)]),
                                       mean(samplingAtPresentProb$test_95),
                                       mean(mrca$test_95),
                                       mean(origin$test_95)))
write.csv(hpd.test, file = paste("figures/HPD_test",".csv", sep="" ))


hpd.width = data.frame("Parameter"=c('diversificationRate1', 'diversificationRate2', 'diversificationRate3', 'diversificationRate4', 
                                     'turnover1', 'turnover2', 'turnover3', 'turnover4', 
                                     'samplingProportion1', 'samplingProportion2', 'samplingProportion3', 'samplingProportion4', 
                                     'samplingAtPresentProb',
                                     'mrca',
                                     'origin'),
                       "Average relative HPD width"=c(mean(diversificationRate$hpd.rel.width[c(TRUE, FALSE, FALSE, FALSE)]), mean(diversificationRate$hpd.rel.width[c(FALSE, TRUE, FALSE, FALSE)]), 
                       mean(diversificationRate$hpd.rel.width[c(FALSE, FALSE, TRUE, FALSE)]), mean(diversificationRate$hpd.rel.width[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                      mean(turnover$hpd.rel.width[c(TRUE, FALSE, FALSE, FALSE)]), mean(turnover$hpd.rel.width[c(FALSE, TRUE, FALSE, FALSE)]), 
                                                      mean(turnover$hpd.rel.width[c(FALSE, FALSE, TRUE, FALSE)]), mean(turnover$hpd.rel.width[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                      mean(samplingProportion$hpd.rel.width[c(TRUE, FALSE, FALSE, FALSE)]), mean(samplingProportion$hpd.rel.width[c(FALSE, TRUE, FALSE, FALSE)]), 
                                                      mean(samplingProportion$hpd.rel.width[c(FALSE, FALSE, TRUE, FALSE)]), mean(samplingProportion$hpd.rel.width[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                      mean(samplingAtPresentProb$hpd.rel.width),
                                                      mean(mrca$hpd.rel.width),
                                                      mean(origin$hpd.rel.width)))
write.csv(hpd.width, file = paste("figures/HPD_rel_width_score_",".csv", sep="" ))

cv = data.frame("Parameter"=c('diversificationRate1', 'diversificationRate2', 'diversificationRate3', 'diversificationRate4', 
                              'turnover1', 'turnover2', 'turnover3', 'turnover4', 
                              'samplingProportion1', 'samplingProportion2', 'samplingProportion3', 'samplingProportion4', 
                              'samplingAtPresentProb',
                              'mrca',
                              'origin'),
                "Average Coefficient of Variation"=c(mean(diversificationRate$cv[c(TRUE, FALSE, FALSE, FALSE)]), mean(diversificationRate$cv[c(FALSE, TRUE, FALSE, FALSE)]), 
                mean(diversificationRate$cv[c(FALSE, FALSE, TRUE, FALSE)]), mean(diversificationRate$cv[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                     mean(turnover$cv[c(TRUE, FALSE, FALSE, FALSE)]), mean(turnover$cv[c(FALSE, TRUE, FALSE, FALSE)]), 
                                                     mean(turnover$cv[c(FALSE, FALSE, TRUE, FALSE)]), mean(turnover$cv[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                     mean(samplingProportion$cv[c(TRUE, FALSE, FALSE)]), mean(samplingProportion$cv[c(FALSE, TRUE, FALSE)]), 
                                                     mean(samplingProportion$cv[c(FALSE, FALSE, TRUE, FALSE)]), mean(samplingProportion$cv[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                     mean(samplingAtPresentProb$cv),
                                                     mean(mrca$cv),
                                                     mean(origin$cv)))
write.csv(cv, file = paste("figures/CV_score",".csv", sep="" ))

medians = data.frame("Parameter"=c('diversificationRate1', 'diversificationRate2', 'diversificationRate3', 'diversificationRate4', 
                                   'turnover1', 'turnover2', 'turnover3', 'turnover4', 
                                   'samplingProportion1', 'samplingProportion2', 'samplingProportion3', 'samplingProportion4',
                                   'samplingAtPresentProb',
                                   'mrca',
                                   'origin'),
                     "Relative Error at Medians"=c(mean(diversificationRate$rel.err.median[c(TRUE, FALSE, FALSE, FALSE)]), mean(diversificationRate$rel.err.median[c(FALSE, TRUE, FALSE, FALSE)]), 
                                                  mean(diversificationRate$rel.err.median[c(FALSE, FALSE, TRUE, FALSE)]), mean(diversificationRate$rel.err.median[c(FALSE, FALSE, FALSE, TRUE)]), 
                                                   mean(turnover$rel.err.median[c(TRUE, FALSE, FALSE, FALSE)]), mean(turnover$rel.err.median[c(FALSE, TRUE, FALSE, FALSE)]), 
                                                    mean(turnover$rel.err.median[c(FALSE, FALSE, TRUE, FALSE)]), mean(turnover$rel.err.median[c(FALSE, FALSE, TRUE, FALSE)]), 
                                                   mean(samplingProportion$rel.err.median[c(TRUE, FALSE, FALSE, FALSE)]), mean(samplingProportion$rel.err.median[c(FALSE, TRUE, FALSE, FALSE)]), 
                                                   mean(samplingProportion$rel.err.median[c(FALSE, FALSE, TRUE, FALSE)]), mean(samplingProportion$rel.err.median[c(FALSE, FALSE, TRUE, FALSE)]), 
                                                   mean(samplingAtPresentProb$rel.err.median),
                                                   mean(mrca$rel.err.median),
                                                   mean(origin$rel.err.median)))
write.csv(medians, file = paste("figures/rel_error_median",".csv", sep="" ))

means = data.frame("Parameter"=c('diversificationRate1', 'diversificationRate2', 'diversificationRate3', 'diversificationRate4', 
                                 'turnover1', 'turnover2', 'turnover3', 'turnover4', 
                                 'samplingProportion1', 'samplingProportion2', 'samplingProportion3', 'samplingProportion4',
                                 'samplingAtPresentProb',
                                 'mrca',
                                 'origin'),
                   "Relative Error at Means"=c(mean(diversificationRate$rel.err.mean[c(TRUE, FALSE, FALSE, FALSE)]), mean(diversificationRate$rel.err.mean[c(FALSE, TRUE, FALSE, FALSE)]), 
                   mean(diversificationRate$rel.err.mean[c(FALSE, FALSE, TRUE, FALSE)]),  mean(diversificationRate$rel.err.mean[c(FALSE, FALSE, FALSE, TRUE)]), 
                                               mean(turnover$rel.err.mean[c(TRUE, FALSE, FALSE, FALSE)]), mean(turnover$rel.err.mean[c(FALSE, TRUE, FALSE, FALSE)]), 
                                               mean(turnover$rel.err.mean[c(FALSE, FALSE, TRUE, FALSE)]),  mean(turnover$rel.err.mean[c(FALSE, FALSE, FALSE, TRUE)]), 
                                               mean(samplingProportion$rel.err.mean[c(TRUE, FALSE, FALSE, FALSE)]), mean(samplingProportion$rel.err.mean[c(FALSE, TRUE, FALSE, FALSE)]), 
                                               mean(samplingProportion$rel.err.mean[c(FALSE, FALSE, TRUE, FALSE)]), mean(samplingProportion$rel.err.mean[c(FALSE, FALSE, FALSE, TRUE)]),
                                               mean(samplingAtPresentProb$rel.err.mean),
                                               mean(mrca$rel.err.mean),
                                               mean(origin$rel.err.mean)))
write.csv(means, file = paste("figures/rel_error_mean",".csv", sep=""))

removed = length(remove_rows)