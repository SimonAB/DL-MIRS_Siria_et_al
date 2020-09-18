# Simulations to ask how many mosquitoes would need to be sampled
# to allow a given difference in age structure between two populations
# (e.g. pre- and post-intervention) to be detected, when age-class is 
# inferred by the the MIRS-CNN method. The effect of enriching the training 
# set of lab-reared mosquitoes using increasing numbers of mosquitoes reared with
# environmental variation (EV) is assessed.
# (Assumptions adapted from on Fig 4 of: http://dx.doi.org/10.12688/wellcomeopenres.15201.2)

# load packages
library(parallel)
library(scales)
library(RColorBrewer)

# clear memory
rm(list = ls())

# get date
date.today <- Sys.Date()

# load "EV variation" (TRUE) or "sampling variation" (FALSE) matrices 
# selecting "EV <- TRUE" will run the simulations to generate Fig 4b,c,d
# and Table S4
# selecting "EV <- FALSE" will run the simulations to generate Fig S2
EV <- TRUE

# assumptions:

# Simulate a population with survival rate of 0.91 (gambiae) or 0.82 (arabiensis)
s <- c(gambiae = 0.91, arabiensis = 0.82)[1]
p <- 1 - s # daily death probability before intervention

# increase in death rate due to two interventions 
# LLIN: 4-fold increase in death rate, starting on day 4 of life
# toxic sugar bait: 3-fold increase in death rate, works from day 1
intervention.tab <- data.frame(effect = c(1, 4, 2), day.active = c(2, 4, 2))
#intervention.tab <- data.frame(effect = c(1, 2, 2), day.active = c(2, 4, 2))
rownames(intervention.tab) <- 
  c("Control",   # natural mortality 
    "LLIN",      # long-lasting insecticide nets (LLIN)
    "ATSB")      # attractive toxic sugar baits (ATSB)
intervention.tab

# maximum lifespan of mosquitoes
n.day <- 20
day <- 1:n.day

# age-classes: 1-4, 5-10, 11+
age.cut <- c(min(day) - 0.5, 4.5, 10.5, max(day) + 0.5)
age.bin <- 
  lapply(2:length(age.cut), function(i) {
    day[sapply(day, function(x) all((x < age.cut[c(i-1, i)]) == c(FALSE, TRUE)))]
  })
names(age.bin) <- apply(sapply(age.bin, range), 2, paste, collapse = "-")
age.bin.num <- 1:length(age.bin)

# age structure in the 3 groups
ageprob.list <-
  lapply(1:nrow(intervention.tab), function(i) {
    death.prob <- 
      c(0, 
        rep(p, intervention.tab$day.active[i] - 2), 
        rep(p *  intervention.tab$effect[i], n.day - intervention.tab$day.active[i] + 1))
    ageprob <- cumprod(1 - death.prob)
    ageprob <- ageprob/sum(ageprob)
    names(ageprob) <- day
    ageprob
  })
names(ageprob.list) <- rownames(intervention.tab)


# convert per-day structure to binned structure
ageprob.bin.list <-
  lapply(1:nrow(intervention.tab), function(i) {
    sapply(age.bin, function(x) sum(ageprob.list[[i]][as.character(x)]))
  })
names(ageprob.bin.list) <- rownames(intervention.tab)

# make plot comparing age structures
cols <- c("black", "blue", "red")
old.par <- par(mar = c(5.1, 4.1, 4.1, 4.1))
ylim <- c(0, ceiling(10 * max(unlist(ageprob.list)))/10)
plot(day, ageprob.list[[1]], type = "n", ylim = ylim,
     xlab = "Mosquito age (days)", ylab = "Proportion in population")
lapply(1:length(ageprob.list), function(i) {
  points(day, ageprob.list[[i]], type = "b", pch = 16, col = cols[i])
})

# add age class proportions
bin.scale <- max(unlist(ageprob.bin.list))/max(unlist(ageprob.list))
lapply(names(age.bin), function(x) {
  lapply(1:nrow(intervention.tab), function(i) {
    lines(cbind(age.bin[[x]], ageprob.bin.list[[i]][x]/bin.scale), col = alpha(cols[i], 0.4), lwd = 4)
  })
})

axis(4, at = pretty(0:1)/bin.scale, labels = pretty(0:1))
mtext("Proportion in population (binned)", side = 4, line = 2.5)


# add legend and title
legend("topright", legend = rownames(intervention.tab), pch = 16, 
       col = cols, lty = 1, bty = "n")
legend("topright", lwd = 4, 
       legend = rownames(intervention.tab),
       col = alpha(cols, 0.4), lty = 1, bty = "n")

par(old.par)


# make bar chart
nice.cols <- c("#e0f3db", "#a8ddb5", "#43a2ca")
ageprob.bin.tab <- do.call("cbind", ageprob.bin.list)
pdf(paste0("agestructure.barchart.", date.today, ".pdf"), height = 5/2.54, width = 6/2.54, pointsize = 10)
old.par <- par(mar = c(2.1, 3.1, 0.6, 0.1))
rownames(ageprob.bin.tab) <- paste0(rownames(ageprob.bin.tab), "d")
barplot(ageprob.bin.tab, beside = TRUE, legend.text = TRUE, 
        args.legend = list(x = "topleft", bty = "n", x.intersp = 0.2, inset = -0.03),
        ylab = "", xlab = "", axes = FALSE, ylim = c(0, 1.1 * max(ageprob.bin.tab)),
        col = nice.cols, padj = -1)
mtext("Proportion", side = 2, line = 2)
rownames(ageprob.bin.tab) <- names(ageprob.bin.list[[1]])
axis(2, at = pretty(c(ageprob.bin.tab)), padj = 0.7)
par(old.par)
dev.off()

# read in confusion matrices
if(EV) {
  mat.tab <- read.csv("Confusion_Matrices/confusion_matrices_all_2020-01-22.csv", header = FALSE)
} else {
  mat.tab <- read.csv("Confusion_Matrices/confusion_matrices_0_05_2020-01-22.csv", header = FALSE)
  #for(j in 1:ncol(mat.tab)) mat.tab[, j] <- mean(mat.tab[, j])
}


mat.names <- paste0("r", rep(1:3, each = 3), "c", rep(1:3, 3))
names(mat.tab) <- mat.names
dim(mat.tab)
if(EV) {
  mat.tab$n.tcv <- c(0, 162, 324, 486, 654, 815, 973, 1131, 1294, 1452)
} else {
  mat.tab$n.tcv <- 1:nrow(mat.tab)
}

# confusion matrices explained:
# each row of 
mat.tab 
# contains 9 values from 3x3 confusion matrix
# which defines the accuracy of the MIRS-CNN method in inferring 
# the age of a mosquito. 
# the final column of mat.tab gives the number of mosquitoes from
# the environmental variation (EV) data set that were added to the training data
# to improve the training of the the CNN (convolutional neural network).
# for example, row 6 of mat.tab
mat.tab[6, ]
# gives the confusion matrix where 815 EV mosquitoes were used, 
# and the value 815. turn these back into a confusion matrix:
matrix(unlist(mat.tab[6, mat.names]), ncol = 3, byrow = TRUE)
# the rows represent true age classes. the columns give the probability
# that a mosquito of that age class will be assigned to each of the three 
# age classes. e.g. the probability of a mosquito in the first age class (1-4 days)
# being correctly assigned to that age class by the MIRS-CNN method is:
matrix(unlist(mat.tab[6, mat.names]), ncol = 3, byrow = TRUE)[1, 1]


# make table of assumptions choices (scenarios to simulate)
# try all combinations of 
#  enrichment (degree of enrichment of the training data with EV)
#  sample size (n wild mosquitoes per intervention group)
assumptions <- 
  expand.grid(
    mat.row = 1:nrow(mat.tab),              # which confusion matrix to use
    n = c(20, 50, 100, 150, 200, 250, 300), # sample size from each population
    nsim = 10000,                           # n data sets to simulate per scenario
    stringsAsFactors = FALSE) 

# set random seeds
RNGkind("L'Ecuyer-CMRG")
global.rand.seed <- 782120569 
# https://www.random.org/integers/?num=1&min=0&max=1000000000&col=1&base=10&format=html&rnd=new
# Random Integer Generator
# Here are your random numbers:
#  782120569
# Timestamp: 2020-02-25 13:07:09 UTC
set.seed(global.rand.seed)
assumptions$global.rand.seed <- global.rand.seed
assumptions$rand.seed <- sample(1e9, nrow(assumptions))  

# simulate populations
start.time <- Sys.time()
simres.tab <-
  sapply(1:nrow(assumptions), function(j) {  # loop over scenarios
    set.seed(assumptions$rand.seed[j])
    mc.reset.stream()
    simres.list <-
      mclapply(1:assumptions$nsim[j],       # analyse nsim simulated data sets
               function(i) {
                 # simulate data with true age in days
                 n <- assumptions$n[j]
                 dat <- 
                   do.call("rbind", 
                           lapply(1:nrow(intervention.tab), function(k) {
                             data.frame(intervention = rownames(intervention.tab)[k], 
                                        age = c(day %*% rmultinom(n, 1, ageprob.list[[k]])))
                           }))
                 
                 # bin true age in age classes
                 dat$age.cat <- as.numeric(cut(dat$age, age.cut, labels = names(age.bin)))
                 
                 # apply confusion matrix to give estimated age class
                 mat <- 
                   matrix(unlist(mat.tab[assumptions$mat.row[j], mat.names]), 
                          ncol = 3, byrow = TRUE)
                 dimnames(mat) <- list(names(age.bin), names(age.bin))
                 # check rows sum to 1
                 rowSums(mat)
                 # "estimate" age class by drawing from a multinomial distribution
                 dat$age.cat.est <- 
                   sapply(dat$age.cat, function(a) age.bin.num %*% rmultinom(1, 1, mat[a, ]))
                 
                 # test both interventions against the control population
                 # using wilcoxon-mann-whitney test and chi-squared test
                 # (only chi-squared test was used, ultimately, as this test had greater power)
                 out.list <-
                   lapply(2:nrow(intervention.tab), function(h) {
                     dat.test <- droplevels(dat[dat$intervention %in% rownames(intervention.tab)[c(1, h)], ])
                     # do wilcoxon-mann-whitney test to compare age distributions
                     table(dat.test$age.cat, dat.test$intervention)
                     wil.pow <- wilcox.test(age.cat ~ intervention, data = dat.test)$p.value < 0.05
                     wil.pow.est <- wilcox.test(age.cat.est ~ intervention, data = dat.test)$p.value < 0.05
                     if(is.na(wil.pow)) wil.pow <- 0 
                     if(is.na(wil.pow.est)) wil.pow.est <- 0
                     # do chi-squared test to compare age distributions
                     xtab <- table(factor(dat.test$age.cat.est, 1:3), dat.test$intervention)
                     chi.pow <- chisq.test(table(dat.test$age.cat, dat.test$intervention))$p.value < 0.05
                     chi.pow.est <- chisq.test(xtab[rowSums(xtab) > 0, ])$p.value < 0.05
                     
                     # export test results
                     out <-
                       c(wil.pow = wil.pow, wil.pow.est = wil.pow.est, 
                         chi.pow = chi.pow, chi.pow.est = chi.pow.est, 
                         prop.control = prop.table(xtab, 2)[, rownames(intervention.tab)[1]],
                         prop.intervention = prop.table(xtab, 2)[, rownames(intervention.tab)[h]])
                     names(out) <- paste(names(out), rownames(intervention.tab)[h], sep = ".")
                     out
                   })
                 unlist(out.list)
               }, mc.cores = detectCores() - 1)
    print(paste0(round(100*j/nrow(assumptions)), "% complete"))
    
    # bind results together as a table
    simres <- do.call("rbind.data.frame", simres.list)
    dim(simres)
    names(simres) <- names(simres.list[[1]])
    
    
    # take mean across all nsim simulations, giving power estimates for each scenario
    apply(simres, 2, mean)
    
  })

# bind assumptions table to results
out <- cbind(assumptions, mat.tab[assumptions$mat.row, ], t(simres.tab))
out[, grep("prop\\.", names(out))] <- round(out[, grep("prop\\.", names(out))], 3)
out$mat.row <- NULL

# compare wilcox and chi-squared results
plot(chi.pow.est.LLIN ~ wil.pow.est.LLIN, data = out, xlab = "Wilcoxon", ylab = "Chisq")
points(chi.pow.est.ATSB ~ wil.pow.est.ATSB, data = out, col = "red")
abline(0, 1)
legend("topleft", legend = c("LLIN", "ATSB"), col = 1:2, pch = 1)

# plot power against sample size broken down by enrichment level
lapply(2:nrow(intervention.tab), function(i) {
  gp <- rownames(intervention.tab)[i]
  form <- formula(paste0("wil.pow.est.", gp, " ~ n"))
  form2 <- formula(paste0("wil.pow.", gp, " ~ n"))
  ntcv.lev <- unique(out$n.tcv)
  ntcv.col <- brewer.pal(length(ntcv.lev), "RdYlBu") 
  if(!EV) ntcv.col <- rep(ntcv.col[2], length(ntcv.lev))
  names(ntcv.col) <- ntcv.lev
  powercurve.file <- 
    paste0("agestructure.powercurve.", 
           ifelse(EV, "", "var."),
           names(s), ".", gp, ".", date.today, ".pdf")
  pdf(powercurve.file, height = 7/2.54, width = 8/2.54, pointsize = 10)
  old.par <- par(mar = c(2.6, 2.6, 0.6, 0.2))
  plot(form, data = out, ylim = 0:1, xlim = c(min(out$n), max(out$n) * 1.20^(!EV - 1)), 
       type = "n", ylab = "", xlab = "", axes = FALSE)
  mtext("N per population", 1, line = 1.5)
  mtext("Power", 2, line = 1.5)
  tcl <- -0.3
  axis(2, padj = 0.9, tcl = tcl)
  axis(1, at = unique(out$n), padj = -0.9, tcl = tcl, gap.axis = 0.25)
  box()
  lapply(ntcv.lev, function(ntcv) {
    points(form, data = out[out$n.tcv == ntcv, ], 
           type = "b", pch = 21, bg = ntcv.col[as.character(ntcv)])
  })
  max.power <- tapply(out[, paste0("wil.pow.", gp)], out$n, mean)
  if(EV) {
    lines(as.numeric(names(max.power)), max.power, lty = 3)
    temp <- legend("bottomright", bty = "n",
                   legend = rep(" ", length(ntcv.lev)), pch = 21,
                   text.width = max(strwidth(ntcv.lev)), xjust = 1, yjust = 1,
                   pt.bg = rev(ntcv.col), x.intersp = 0.4, inset = -0.01)
    text(temp$rect$left + temp$rect$w, temp$text$y,
         rev(ntcv.lev), pos = 2)
  }
  #title(rownames(intervention.tab)[i])
  par(old.par)
  dev.off()
})


if(!EV) {
  head(out)
  sapply(out[out$n == 20, grep("wil.pow.est", names(out), value = TRUE)], sd) / 
    sapply(out[out$n == 20, grep("wil.pow.est", names(out), value = TRUE)], function(x) {
      pwr <- mean(x)
      sqrt(pwr * (1-pwr) / unique(out$nsim))
    })
  
}

# write results to csv
out.file <- paste0("agestructure.power.",
                   ifelse(EV, "", "var."),
                   names(s), ".", date.today, ".csv")
write.csv(out, out.file, row.names = FALSE)

print(Sys.time() - start.time)

# post-formatting of results for Table S4

# formatting numbers
# this function is better than round because it doesn't strip off trailing zeroes
library(gdata)
my.format<-
  function(x,ndp=0,na.string="") {
    out<-
      format(round(x,ndp),ns=ndp,scientific=FALSE,just='none')
    out[grep("NA",out)]<-na.string
    trim(out)
  }

if(EV) {
  TableS4 <- 
    read.csv(out.file)[, c("n", 
                           "n.tcv", 
                           "wil.pow.est.LLIN", 
                           "wil.pow.est.ATSB")]
  TableS4$wil.pow.est.LLIN <- paste0(my.format(TableS4$wil.pow.est.LLIN * 100, 1), "%")
  TableS4$wil.pow.est.ATSB <- paste0(my.format(TableS4$wil.pow.est.ATSB * 100, 1), "%")
  write.csv(TableS4, "TableS4.csv", row.names = FALSE)
}

