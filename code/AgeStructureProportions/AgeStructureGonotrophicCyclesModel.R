#' DL-MIRS project
#' 
#' R Script to estimate age structure proportions with 95% CIs based
#' on mosquitoes that have been manually dissected and assigned to three age classes 
#' (nulliparous, one gonotrophic cycle, 2-4 gonotrophic cycles). Do this separately 
#' for An. arabiensis from IHI and An. coluzzi from IRSS.
#'    There is a fourth class, gravid, which is assumed to have the same age 
#' class proportions as non-gravid mosquitoes, so can be ignored.
#' 
#' Paul Johnson
#' 2021-10-29

#' Clear objects from memory
rm(list = ls())

#' Load real data or use demo data?
demo.data <- !file.exists("Total_Summary_dissection2021_corrected.csv")

#' Load or create data
if(demo.data) {
  #' Create demo data
  dat <-
    data.frame(P0 =   c( 50,  70), 
               P1 =   c(100, 200), 
               P234 = c( 20,  40), 
               Location.species = c("IRSS - AC demo", "IHI - AA demo"))
  
} else {
  #' Load data
  dat <- read.csv("Total_Summary_dissection2021_corrected.csv")
  head(dat)
  
  #' Create columns with counts of each age class
  dat$P0 <- dat$DissectedTrainingNulliparous
  dat$P1 <- dat$DissectedTraining1cycle
  dat$P234 <- 
    dat$DissectedTraining2cycles + 
    dat$DissectedTraining3cycles + 
    dat$DissectedTraining4cycles
  
  #' Drop columns that aren't required
  dat$DissectedTrainingNulliparous <- dat$DissectedTraining1cycle <- 
    dat$DissectedTraining2cycles <- dat$DissectedTraining3cycles <- 
    dat$DissectedTraining4cycles <- dat$Total.physiological.age.determined <-
    dat$NotDissectedTestsetAgeUnknown <- dat$DissectedNotTrainingGravid <-
    dat$Date.of.collection <- dat$Date.of.dissection <- NULL
}



#' Much much data from each location?
dat$Location.species <- factor(dat$Location.species)
table(dat$Location.species, useNA = "always")

#' Break down by location-species
results.list <-
  lapply(levels(dat$Location.species), function(loc) {
    
    #' Model the three counts (N_0, N_1, N_234) as multinomial-Dirichlet
    #' The prior distribution of the proportions is Dirichlet:
    #' (P_0, P_1, P_234) ~ Dirichlet(alpha_0, alpha_1, alpha_234)
    #' The posterior is also Dirichlet:
    #' (P_0, P_1, P_234) ~ Dirichlet(N_0 + alpha_0, N_1 + alpha_1, N_234 + alpha_234)
    #' The marginal distribution of each P is Beta, so we can use the Beta quantile
    #' function below to get the 95% CIs.

    #' Observed counts in each age class
    N <- colSums(dat[dat$Location.species == loc, c("P0", "P1", "P234")])
    
    #' Concentration parameter of Dirichlet prior on age class proportions:
    alpha <- rep(1, length(N))
    names(alpha) <- names(N)

    #' Posterior means and 95% CIs from Beta quantile finction
    out <-
      sapply(names(N), function(age) {
        a <- N[age] + alpha[age]
        b <- sum(N) - N[age] + sum(alpha) - alpha[age]
        c(as.vector(N[age]), as.vector(a/(a+b)), qbeta(c(0.025, 0.975), a, b))
      })
    
    rownames(out) <- c("data", "posterior.mean", "ci95lower", "ci95upper")
    out
  })

names(results.list) <- levels(dat$Location.species)

#' Bar plot of estimates and CIs
par(mfrow = c(1, 2))
lapply(names(results.list), function(loc) {
  res <- results.list[[loc]]
  x <- barplot(res["posterior.mean", ], ylim = 0:1)
  arrows(x, res["ci95lower", ], x, res["ci95upper", ], code = 3, angle = 90, length = 0.1,
         ylab = "Proportion")
  title(loc)
})

#' Print results to console
print(results.list)

if(demo.data) warning("This analysis was run using demo data, not real data.")
