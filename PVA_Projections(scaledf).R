#------------------------------------------------------------------------------#
# Population Evaluation Tool (PET) for the North Atlantic Right Whale (NARW)
# Version 1.0

# Structure of input files (actual names may vary):
#   PVA_population_t0.Rdata: nBoot or more samples of starting stage and wound
#      distributions.
#   NARW_posteriors_REPRO.csv: nBoot (or more) samples of parameters from 
#     integrated survival/reproduction model. 
#   NARW_posteriors_MORT.csv: nBoot (or more) samples of fractions of mortality and 
#     wounding probabilities.
#   NARW_food_covariates_1986-2019.Rdata: Historical data on prey availability, not yet 
#     sampled, averaged, or weighted.
#------------------------------------------------------------------------------#
library(zoo);library(abind)
out_drive <- "C:\\temp\\"
library(ggplot2)
plot_by_scenario <- ggplot(summary_data_ntot, aes(x = Year, y = Ntot.mean)) +
  geom_ribbon(aes(ymin = Ntot.lower, ymax = Ntot.upper, fill = Scenario),
              alpha = 0.2, color = NA) +
  geom_line(aes(color = Scenario)) +
  geom_hline(yintercept = starting_abundance,
             linetype = "dashed", color = "black") +
  facet_wrap(~ Scenario, ncol = 2) +
  labs(title = "Population Trajectories by Scenario",
       x = "Year", y = "Total Population Size") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 5, "mm"))

print(plot_by_scenario)
#------------------------------------------------------------------------------#
# Basic settings----
#------------------------------------------------------------------------------#
version <- "1.00"

nBoot <- 5#, number of bootstrap runs
nRep  <- 1   #, number of replications (Monte Carlo loop). 
nT <- 100# number of years

## Reproductive stages----
stages <- c('F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 
            'F9', 'F10+', 'FC', 'FR', 'FW', 
            'M1', 'M2', 'M3', 'M4', 'M5+')
stagesLong <- c('Female Prebreeder 0.5 Years Old', 'Female Prebreeder 1.5 Years Old', 
                'Female Prebreeder 2.5 Years Old', 'Female Prebreeder 3.5 Years Old', 
                'Female Prebreeder 4.5 Years Old', 'Female Prebreeder 5.5 Years Old', 
                'Female Prebreeder 6.5 Years Old', 'Female Prebreeder 7.5 Years Old', 
                'Female Prebreeder 8.5 Years Old', 'Female Prebreeder 9.5+ Years Old', 
                'Female With Calf', 'Female Resting', 'Female Breeder (No Calf)', 
                'Male 0.5 Years Old', 'Male 1.5 Years Old', 'Male 2.5 Years Old', 
                'Male 3.5 Years Old', 'Male 4.5+ Years Old')
nStages <- length(stages)
femStages <- stages[grep("F", stages)]
nFemStages <- length(femStages)
adFemStages <- femStages[-(1:4)]

### Reproductive stage effects struct32ure on B and S----
reproStages <- c(5:10, "W")
nReproStages <- length(reproStages)
print(nReproStages)
dim(survAges <- 1:5)
survStages <- 1:5
nSurvAges <- length(survAges)
nSurvStages <- length(survStages)

## Wound/death states----
woundStates <- c("fine", "ent", "vessel") #, "other")
woundStatesLong <- c("Not Wounded (Or Not Severely)", "Severe Entanglement Wound",
                     "Severe Vessel Wound") #, "Severe Other Wound")
nWoundStates <- length(woundStates)

deadStates <- c("ent", "vessel", "other")
deadStatesLong <- c("Recently Dead Entanglement", "Recently Dead Vessel",
                    "Recently Dead Other")
nDeadStates <- length(deadStates)

liveDeadStates <- c("alive", deadStates)
liveDeadStatesLong <- c("Alive", deadStatesLong)
nLiveDeadStates <- nDeadStates + 1

entangleStages <- c("Calf", "Juvenile", "Adult", "FC", "FR") 
nEntangleStages <- length(entangleStages)
stagesToES <- c(rep(entangleStages[-c(4:5)], c(1, 3, 9)), #Females
                rep(entangleStages[-c(4:5)], c(1, 3, 1))) #Males
# females currently/recently w/ calf are entanglement specific states
stagesToES[stages %in% c("FC","FR")] <- c("FC","FR")

## Extinction thresholds----
thresholds <- c(1, 10, 50, 100)
nThresholds <- length(thresholds)

ceiling_N <- 10000 # Just in case
verbose <- 2 #integer 0 to 5
seed <- 2023 #opportunity to replicate random sampling
set.seed(seed)

#------------------------------------------------------------------------------#
# Read in the functions----
#------------------------------------------------------------------------------#
source("PVA_functions.R")
source("PVA_scenarios.R")

#------------------------------------------------------------------------------#
# Read in parameters, set up starting values, etc.----
#------------------------------------------------------------------------------#

## posterior distributions from reproduction and mortality models----
post_repro <- read.csv("inputs/NARW_posteriors_REPRO.csv")
post_mort  <- read.csv("inputs/NARW_posteriors_MORT.csv")

print(names(post_mort))

## match Greek letter convention for mortality model
mort.coeff <- grep(paste(c("beta.m","beta.i"),collapse="|"),names(post_mort),value=T)
names(post_mort)[match(mort.coeff,names(post_mort))] <- gsub("beta","a",mort.coeff)

## reduce size to number of bootstrap runs per projection----
post_repro <- post_repro[nBootKeep(post_repro),]
post_mort  <-  post_mort[nBootKeep(post_mort),]

## natural mortality (set to 0)----
# based on Dureuil & Froese (2021) and t_max = 268 (bowhead whales)
post_mort$Mu.mO <- Mu.mO_value #rgamma(nrow(post_mort),1,99) # set to specific hazard mortality rate 


## starting stage/state distributions----
load("inputs/PVA_population_t0.Rdata")
N0 <- N0[nBootKeep(N0),]
dimnames(N0) <- list(NULL,stages)
wound0 <- wound0[nBootKeep(wound0),,]
dimnames(wound0) <- list(NULL, stages, woundStates)

## prey index data----
load(file="inputs/NARW_food_covariates_1986-2019.Rdata")

## reproductive parameters----
post_repro$beta.ageW <- (post_repro$p.beta.ageW)
post_repro[, grep("beta.age.first",colnames(post_repro),value=TRUE)] <-
  post_repro[, grep("beta.age.first",colnames(post_repro),value=TRUE)] + post_repro[,"beta.ageW"]
betas <- post_repro[, c(grep("beta.age.first",colnames(post_repro),value=TRUE),
                        "beta.regime2", "beta.regime2W", 
                        "beta.p1", "beta.p2", 
                        "beta.inj"
)]
colnames(betas) <- c(5:10, "W", "b.regime2", "b.regime2W",
                     "b.prey1", "b.prey2", "b.inj") 

## calf loss parameter----
kappa <- post_repro[, "kappa"]*0.5

## survival/mortality parameters----
alphas <- post_mort[,c(
  "Mu.mO","Mu.mE", "Mu.mV",
  "a.mE.age", "a.mE.calf", "a.mE.rest",
  "a.mV.age", "a.mV.calf", "a.mV.rest"
)]

## injury rates/coefficients----
iTheta <- post_mort[,c(
  "Mu.iE","Mu.iV","psiE",
  "a.iE.age", "a.iE.calf", "a.iE.regime2", "a.iE.rest", 
  "a.iV.age", "a.iV.calf", "a.iV.regime2", "a.iV.rest" 
)]

## random effects----
## annual deviations in injury/mortality and reproduction
eps.i <- eps.m <- array(0, c(2, nBoot, nT),
                        list(c("E","V"), NULL, NULL))
eps.r <- eps.i[1,,]
eps <- list(random = list(eps.i,eps.m,eps.r),
            # [NOT USED] for using estimated epsilon values
            period = list(eps.i,eps.m,eps.r))
names(eps[[1]]) <- names(eps[[2]]) <- c("inj","mort","repro")
rm(list = c("eps.i","eps.m","eps.r"))
eps.period <- which(1991:2019 > 2013)

for (i in 1:nBoot) {
  # sample of year-specific deviations
  r.samp <- sample(1:length(eps.period),nT,replace=T)
  
  # injury
  eps[["random"]][["inj"]][1, i, ] <- rnorm(nT, 0, post_mort[i, "sigma.iE.t"])
  eps[["random"]][["inj"]][2, i, ] <- rnorm(nT, 0, post_mort[i, "sigma.iV.t"])
  
  # mortality
  eps[["random"]][["mort"]][1, i, ] <- 0 #rnorm(nT, 0, post_mort[i, "sigma.mE.t"])
  eps[["random"]][["mort"]][2, i, ] <- 0 #rnorm(nT, 0, post_mort[i, "sigma.mV.t"])
  
  # reproduction
  eps[["random"]][["repro"]][i, ] <- rnorm(nT, 0, post_repro[i, "sigma.B.t"])
}
eps.type <- "random"

#------------------------------------------------------------------------------#
# Scenario building----
#------------------------------------------------------------------------------#
load.prey <- FALSE
preyChange_v <- c(0.7,0.8,0.9,1.1,1.2,1.3)

if (!load.prey){
  preyArrays <- list()
  
  # post-2010 prey availability
  preyArrays[["prey post2010"]] <- build.scenario.prey(
    plankton_w_food,
    ref_yrs = list(2010:2019),
    proj_yrs = list(1:nT)
  )
  # historic fluctuations
  preyArrays[["prey historic"]] <- build.scenario.prey(
    plankton_w_food,
    # good vs. bad decade
    #ref_yrs = rep(list(2001:2009,1990:2000),nT/20),
    ref_yrs = list(1990:2009)
    # decadal cycling
    # proj_yrs = split(1:nT,ceiling((1:nT)/10))
  )
  
  # add prey availability due to noise
  for (a in 1:6){
    preyArrays[[a+2]] <- build.scenario.prey(
      plankton_w_food,
      ref_yrs = list(2010:2019),
      proj_yrs = list(1:nT), preyChange = preyChange_v[a]
      
    )
  }
  names(preyArrays)[(1:length(preyChange_v))+2] <- paste0(
    "prey post2010 ", preyChange_v*100,"%")
  save(preyArrays,preyChange_v,file="./inputs/preyArrays_ALLscenarios.Rdata")
}
load(file="./inputs/preyArrays_ALLscenarios.Rdata")

# Multiplying preyArrays by 0 to eliminate prey influence
preyArrays <- lapply(preyArrays, function(arr) {
  arr * 0
})

woundArrays <- list()

#ent.reduce <- c(0,25)
ent.reduce <- c(seq(0,100,by=10),25)
# total entanglement reductions
start <- Sys.time()
for (e in ent.reduce){
  woundArrays[[paste0("reduce entanglement ",e,
                      "%; vessel strike up 0%; speed restrict 0%")]] <- build.scenario.wound(
                        iTheta,
                        eps.i = eps[[eps.type]][["inj"]],
                        iE.change = list(rep((100-e)/100,5))
                      )
}
(end <- Sys.time()-start)

# weak rope coverage (entanglement reductions)
weak.e <- c(50)
for (e in weak.e){
  woundArrays[[paste0("weak rope coverage ",e,
                      "%; vessel strike up 0%; speed restrict 0%")]] <- build.scenario.wound(
                        iTheta,
                        eps.i = eps[[eps.type]][["inj"]],
                        iE.change = list(c(rep((100-(e*.9))/100,2),rep((100-(e*.6))/100,3)))
                      )
}



# all vessel scenarios
ent.reduce <- c(0,25)
vess_t_change <- c(0.7,0,-0.3)
vess_speed <- c(1,0.75,0)

#entanglement reduction
for (e in ent.reduce){        
  #speed restrict
  for (v2 in vess_speed[list(1:3,1:2)[[match(e,ent.reduce)]]]){
    #traffic change
    for (v1 in vess_t_change[list(1:3,1:3,2)[[match(v2,vess_speed)]]]){  
      woundArrays[[paste0(
        "reduce entanglement ",e,"%; vessel strike up ",v1,
        "%; speed restrict ",(1-v2)*100,"%")]] <- build.scenario.wound(
          iTheta,
          eps.i = eps[[eps.type]][["inj"]],
          iE.change = list(rep(1-(e/100),5)), #iE = 0.50 of status quo rate
          iV.change = list(rep(v2,5)),
          iV.change.t = list(rep(1+(v1/100),5))
        )
    }
  }
}

#-----------------------------------#
# count scenarios and name them
#-----------------------------------#

scenario.dat <- expand.grid(
  mortality = c("Low Mortality", "High Mortality"),
  calving = c("Low Fecundity", "High Fecundity"),
  stringsAsFactors = FALSE
)

scenario.dat$ceiling_N <- ceiling_N #?

scenario.dat$names <- apply(scenario.dat[, c("mortality", "calving")], 1, function(x) paste(x, collapse = "; "))
scenarios <- scenario.dat$names
nS <- nrow(scenario.dat)

print(scenario.dat)



#load dans data and assemble 
# Read in your mortality‐posterior matrix
mort_df <- read.csv("/Users/graciesemmens/Downloads/NARW_mortality_rates_1990_2022.csv")
library(dplyr)
library(tidyr)
library(dplyr)

# building the year‐column names
cols_90_00 <- paste0("X", 1990:2000)
cols_15_19 <- paste0("X", 2015:2019)

#computing period means 
mort_periods <- mort_df %>%
  mutate(bootstrap = row_number()) %>%
  mutate(
    mort_early = rowMeans(across(all_of(cols_90_00)), na.rm = TRUE),
    mort_late = rowMeans(across(all_of(cols_15_19)), na.rm = TRUE)
  ) %>%
  select(mort_early, mort_late)



#------------------------------------------------------------------------------#
# Output data structures----
#------------------------------------------------------------------------------#
N <- array(NA, c(nBoot, nRep, nT, nS, nStages),
           dimnames = list(NULL, NULL, NULL, scenarios, stages)
)
Ntot <- array(NA, c(nBoot, nRep, nT, nS))
PQE <- array(0, c(nBoot, nT, nS, nThresholds),
             dimnames = list(NULL, NULL, scenarios, thresholds)
)
Ndead <- array(0, c(nBoot, nRep, nT, nS, nStages, nDeadStates),
               list(NULL, NULL, 1:nT, scenarios, stages, deadStates)
)
Nborn <- array(0, c(nBoot, nRep, nT, nS),
               list(NULL, NULL, 1:nT, scenarios)
)
propEntangled <- array(NA, c(nBoot, nRep, nT, nS))
propStruck <- array(NA, c(nBoot, nRep, nT, nS))
propOther <- array(NA, c(nBoot, nRep, nT, nS))

#------------------------------------------------------------------------------#
# Run simulation----
#------------------------------------------------------------------------------#
library(parallel)
library(foreach)
library(doParallel)


#Enrico's data
load("/Users/graciesemmens/Downloads/New_fecundity_means_GSemmens.RData") # Adjust path
print(calving)
print(str(calving))


#new cluster processing
cl<-makeCluster(n.cores)
clusterExport(cl,varlist=ls())
registerDoParallel(cl)

print(1- mean_surv_early)
print(1-mean_surv_late)
testlate <- .0697
testearly <- .0396
print(mean_surv_late)
print(mean_surv_early)
  
#scenarios using posterior data
for (scenario in start.scenario:nS) {
  if (verbose > 0)
    cat("Running", scenarios[scenario], "scenario\n")
  print(Sys.time())
  
  # Determine natural mortality (Mu.mO) based on the scenario (SINGLE VALUE NOW)
  if (scenario.dat$mortality[scenario] == "Low Mortality") {
    Mu.mO_value <- sample(mort_periods$mort_early, 1)
  } else {
    Mu.mO_value <-sample(mort_periods$mort_late, 1)
  }
  
  alphas$Mu.mO <- Mu.mO_value # Set Mu.mO to the calculated single value
  alphas$Mu.mE <- 0
  alphas$Mu.mV <- 0 
  # fecund scenario (SAMPLING FOR EACH BOOTSTRAP)
  base_calf <- numeric(nBoot) # Initialize base_calf as a numeric vector
  if (scenario.dat$calving[scenario] == "Low Fecundity") {
    for (i in 1:nBoot) {
      base_calf[i] <- sample(calving$calving_late)*1.2
    }
  } else {
    for (i in 1:nBoot) {
      base_calf[i] <- (calving$calving_early)*0.9533
    }
  }
  
  # fill the four non‐zero stages FOR EACH BOOTSTRAP:
  for (i in 1:nBoot) {
    betas[i, "5"] <- 0
    betas[i, "6"] <- 0
    betas[i, "7"] <- 0
    betas[i, "8"] <- 0
    betas[i, "9"] <- base_calf[i] * 0.256
    betas[i, "10"] <- base_calf[i] * 0.336
    betas[i, "W"] <- base_calf[i]
  }
  
  ##simulation
  print(system.time({
    temp <- foreach(i = 1:nBoot) %dopar% {
      print(kappa[i])
      runPVA.par(
        params = list(
          N0 = N0[i, ], betas = as.matrix(betas)[i, ],
          eps.i = eps[[eps.type]][["inj"]][, i, ],
          eps.m = eps[[eps.type]][["mort"]][, i, ],
          alphas = as.matrix(alphas)[i, ],
          woundProb = woundArrays[[1]][i, , , , ],
          wound0 = wound0[i, , ],
          kappa = kappa[i]
        ),
        ceiling_N = scenario.dat$ceiling_N[scenario],
        nT = nT,
        nRep = nRep
      )
    }
  }))
  
  
  PQE[, , scenario, ] <- abind(lapply(1:nBoot, function(x) temp[[x]]$PQE), along = 0)
  N[, , , scenario, ] <- abind(lapply(1:nBoot, function(x) temp[[x]]$N), along = 0)
  Ntot[, , , scenario] <- abind(lapply(1:nBoot, function(x) temp[[x]]$Ntot), along = 0)
  Ndead[, , , scenario, , ] <- abind(lapply(1:nBoot, function(x) temp[[x]]$Ndead), along = 0)
  Nborn[, , , scenario] <- abind(lapply(1:nBoot, function(x) temp[[x]]$Nborn), along = 0)
  propEntangled[, , , scenario] <- abind(lapply(1:nBoot, function(x) temp[[x]]$propEntangled), along = 0)
  propStruck[, , , scenario] <- abind(lapply(1:nBoot, function(x) temp[[x]]$propStruck), along = 0)
  propOther[, , , scenario] <- abind(lapply(1:nBoot, function(x) temp[[x]]$propOther), along = 0)
}
stopCluster(cl)
gc()

#------------------------------------------------------------------------------#
#graph
#------------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)



# Create a data frame for ggplot
plot_data_ntot <- expand.grid(
  Year = 1:nT,
  Scenario = scenarios,
  Boot = 1:nBoot,
  Rep = 1:nRep
)

plot_data_ntot$Ntot <- NA

# Populate the Ntot column
for (b in 1:nBoot) {
  for (r in 1:nRep) {
    for (t in 1:nT) {
      for (s in 1:length(scenarios)) {
        plot_data_ntot$Ntot[plot_data_ntot$Year == t &
                              plot_data_ntot$Scenario == scenarios[s] &
                              plot_data_ntot$Boot == b &
                              plot_data_ntot$Rep == r] <- Ntot[b, r, t, s]
      }
    }
  }
}

# Calculate mean and CI
summary_data_ntot <- plot_data_ntot  %>%
  group_by(Year, Scenario) %>%
  summarize(
    Ntot.mean  = mean(Ntot, na.rm = TRUE),
    Ntot.lower = quantile(Ntot, 0.025, na.rm = TRUE),
    Ntot.upper = quantile(Ntot, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Extract starting abundance (at Year 1)
starting_abundance <- summary_data_ntot$Ntot.mean[summary_data_ntot$Year == 1][1] #Take the first value, as all scenarios start at the same value.

# Create the plot with horizontal line
plot_original <- ggplot(summary_data_ntot, aes(x = Year, y = Ntot.mean, color = Scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = Ntot.lower, ymax = Ntot.upper, fill = Scenario),
              alpha = 0.2, color = NA) +
  labs(title = "Population Trajectories",
       x = "Year",
       y = "Total Population Size") +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 5, "mm")) +
  geom_hline(yintercept = starting_abundance, linetype = "dashed", color = "black")

# Create the zoomed-in plot
plot_zoomed <- plot_original +
  coord_cartesian(ylim = c(0, 2500))

# Display the plots
print(plot_original)
print(plot_zoomed)

##other plots 

final_Ntot <- Ntot[, , nT, ] # Extract total population size at the final year for all boots and scenarios

# Calculate mean and 95% confidence intervals for final population size for each scenario
mean_final_Ntot <- apply(final_Ntot, 2, mean, na.rm = TRUE)
lower_final_Ntot <- apply(final_Ntot, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_final_Ntot <- apply(final_Ntot, 2, quantile, probs = 0.975, na.rm = TRUE)

# Assuming your threshold for PQE is, for example, the 50th threshold (thresholds[3] if thresholds <- c(1, 10, 50, 100))
threshold_index <- which(thresholds == 50) # Find the index for the threshold of 50

final_PQE <- PQE[, nT, , threshold_index] # Extract PQE at the final year for the threshold of 50 for all boots and scenarios

# Calculate mean and 95% confidence intervals for PQE for each scenario
mean_final_PQE <- apply(final_PQE, 2, mean, na.rm = TRUE)
lower_final_PQE <- apply(final_PQE, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_final_PQE <- apply(final_PQE, 2, quantile, probs = 0.975, na.rm = TRUE)

scenario_names <- scenarios # Your scenario names

results_df <- data.frame(
  Scenario = scenario_names,
  Mean_Final_Ntot = mean_final_Ntot,
  Lower_Final_Ntot = lower_final_Ntot,
  Upper_Final_Ntot = upper_final_Ntot,
  Mean_Final_PQE = mean_final_PQE,
  Lower_Final_PQE = lower_final_PQE,
  Upper_Final_PQE = upper_final_PQE
)

print(results_df)
# Assuming Ntot has dimensions (nBoot, nRep, nT, nS) and 'scenarios' is a vector
# with your scenario names (e.g., c("Low Mort; Low Fecund", ...))

##lambda visualizations: 
# Function to calculate approximate lambda for a given scenario index (s)
calculate_approx_lambda <- function(n_tot_array, nT, scenario_index) {
  mean_Ntot_start <- mean(n_tot_array[, , 1, scenario_index])
  mean_Ntot_end <- mean(n_tot_array[, , nT, scenario_index])
  approx_lambda <- (mean_Ntot_end / mean_Ntot_start)^(1 / (nT - 1))
  return(approx_lambda)
}

# Calculate lambda for each of your scenarios
lambda_scenario <- numeric(length(scenarios))
for (i in 1:length(scenarios)) {
  lambda_scenario[i] <- calculate_approx_lambda(Ntot, nT, i)
  cat(paste("Lambda (", scenarios[i], "):", round(lambda_scenario[i], 3), "\n"))
}

# Analyze the impact of mortality and fecundity (assuming your scenarios are ordered)
# 1: Low Mort, Low Fecund
# 2: High Mort, Low Fecund
# 3: Low Mort, High Fecund
# 4: High Mort, High Fecund

if (length(scenarios) == 4) {
  # Impact of Mortality
  percent_change_lambda_mortality_low_fecund <- ((lambda_scenario[2] - lambda_scenario[1]) / lambda_scenario[1]) * 100
  cat(paste("\nChange in lambda due to High Mortality (vs. Low) at Low Fecundity:", round(percent_change_lambda_mortality_low_fecund, 1), "%\n"))
  
  percent_change_lambda_mortality_high_fecund <- ((lambda_scenario[4] - lambda_scenario[3]) / lambda_scenario[3]) * 100
  cat(paste("Change in lambda due to High Mortality (vs. Low) at High Fecundity:", round(percent_change_lambda_mortality_high_fecund, 1), "%\n"))
  
  # Impact of Fecundity
  percent_change_lambda_fecundity_low_mort <- ((lambda_scenario[3] - lambda_scenario[1]) / lambda_scenario[1]) * 100
  cat(paste("Change in lambda due to High Fecundity (vs. Low) at Low Mortality:", round(percent_change_lambda_fecundity_low_mort, 1), "%\n"))
  
  percent_change_lambda_fecundity_high_mort <- ((lambda_scenario[4] - lambda_scenario[2]) / lambda_scenario[2]) * 100
  cat(paste("Change in lambda due to High Fecundity (vs. Low) at High Mortality:", round(percent_change_lambda_fecundity_high_mort, 1), "%\n"))
  
  # Overall Sensitivity Statement
  avg_mortality_effect <- mean(abs(c(percent_change_lambda_mortality_low_fecund, percent_change_lambda_mortality_high_fecund)))
  avg_fecundity_effect <- mean(abs(c(percent_change_lambda_fecundity_low_mort, percent_change_lambda_fecundity_high_mort)))
  
  if (avg_mortality_effect > avg_fecundity_effect) {
    relative_sensitivity <- avg_mortality_effect / avg_fecundity_effect
    cat(paste("\nOverall, the approximate population growth rate (lambda) appears to be more sensitive to mortality than fecundity by a factor of approximately", round(relative_sensitivity, 1), ".\n"))
  } else if (avg_fecundity_effect > avg_mortality_effect) {
    relative_sensitivity <- avg_fecundity_effect / avg_mortality_effect
    cat(paste("\nOverall, the approximate population growth rate (lambda) appears to be more sensitive to fecundity than mortality by a factor of approximately", round(relative_sensitivity, 1), ".\n"))
  } else {
    cat("\nOverall, the approximate population growth rate (lambda) appears to be similarly sensitive to mortality and fecundity.\n")
  }
} else {
  cat("\nNote: The sensitivity analysis comparing mortality and fecundity is set up assuming you have four scenarios in the order: Low Mort/Low Fecund, High Mort/Low Fecund, Low Mort/High Fecund, High Mort/High Fecund. Please adjust the indexing if your scenarios are ordered differently.\n")
}

library(ggplot2)
library(ggplot2)

ggplot(results_df, aes(x = Scenario, y = Mean_Final_Ntot, fill = Scenario)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower_Final_Ntot, ymax = Upper_Final_Ntot), width = 0.2) +
  labs(title = "Mean Projected Population Size at Year 100",
       x = "Scenario (Mortality; Fecundity)",
       y = "Mean Total Population Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") # Removed legend as Scenario is on x-axis



## some testing to make it look like the real popylation growth 
# Function to calculate annual average growth/decline rate
calculate_annual_rate <- function(df, start_year, end_year) {
  start_pop <- df %>%
    filter(Year == start_year) %>%
    pull(mean)
  end_pop <- df %>%
    filter(Year == end_year) %>%
    pull(mean)
  n_years <- end_year - start_year
  if (length(start_pop) > 0 && length(end_pop) > 0 && n_years > 0) {
    proportional_change <- end_pop / start_pop
    annual_rate <- (proportional_change^(1 / n_years)) - 1
    return(annual_rate)
  } else {
    return(NA)
  }
}

# Calculate historical annual average growth/decline rates
growth_1990_2010 <- calculate_annual_rate(pop_df, 1990, 2010)
decline_2010_2018 <- calculate_annual_rate(pop_df, 2010, 2018)

print(paste("Historical Annual Average Growth Rate (1990-2010):", round(growth_1990_2010 * 100, 2), "%"))
print(paste("Historical Annual Average Decline Rate (2010-2018):", round(decline_2010_2018 * 100, 2), "%"))

# --- Comparison to PVA Scenarios ---
# Scenario 1 (Low Mortality, Low Fecundity) for 1990-2010
pva_growth_1990_2010_s1 <- lambda_scenario[1] - 1
cat(paste("PVA Approximate Annual Growth Rate (Scenario 1):",
          round(pva_growth_1990_2010_s1 * 100, 2), "%\n"))

# Scenario 2 (High Mortality, Low Fecundity) for 2010-2018
pva_decline_2010_2018_s2 <- lambda_scenario[2] - 1
cat(paste("PVA Approximate Annual Decline Rate (Scenario 2):",
          round(pva_decline_2010_2018_s2 * 100, 2), "%\n"))

# Scenario 3 (Low Mortality, High Fecundity) for 1990-2010
pva_growth_1990_2010_s3 <- lambda_scenario[3] - 1
cat(paste("PVA Approximate Annual Growth Rate (Scenario 3):",
          round(pva_growth_1990_2010_s3 * 100, 2), "%\n"))

# Scenario 4 (High Mortality, High Fecundity) for 2010-2018
pva_decline_2010_2018_s4 <- lambda_scenario[4] - 1
cat(paste("PVA Approximate Annual Decline Rate (Scenario 4):",
          round(pva_decline_2010_2018_s4 * 100, 2), "%\n"))

# --- 1. Compare Historical vs. Scenario 1 (1990-2010) ---
cat("\n1) Comparing 1990-2010 (Historical vs. Scenario 1):\n")
if (!is.na(growth_1990_2010)) {
  cat(paste("   Historical:", round(growth_1990_2010 * 100, 2), "% annual growth\n"))
  cat(paste("   PVA (Scenario 1):", round(pva_growth_1990_2010_s1 * 100, 2), "% annual growth\n"))
  difference_growth1 <- pva_growth_1990_2010_s1 - growth_1990_2010
  if (abs(difference_growth1) < 0.01) {
    cat("   PVA growth is roughly hitting the historical target.\n")
  } else if (difference_growth1 > 0) {
    cat(paste("   PVA growth is", round(difference_growth1 * 100, 2), "% higher than historical.\n"))
  } else {
    cat(paste("   PVA growth is", round(abs(difference_growth1) * 100, 2), "% lower than historical.\n"))
  }
} else {
  cat("   Historical growth rate for 1990-2010 could not be calculated (missing data).\n")
}

# --- 2. Compare Historical vs. Scenario 2 (2010-2018) ---
cat("\n2) Comparing 2010-2018 (Historical vs. Scenario 2):\n")
if (!is.na(decline_2010_2018)) {
  cat(paste("   Historical:", round(decline_2010_2018 * 100, 2), "% annual decline\n"))
  cat(paste("   PVA (Scenario 2):", round(pva_decline_2010_2018_s2 * 100, 2), "% annual decline\n"))
  difference_decline2 <- pva_decline_2010_2018_s2 - decline_2010_2018
  if (abs(difference_decline2) < 0.01) {
    cat("   PVA decline is roughly hitting the historical target.\n")
  } else if (difference_decline2 > 0) {
    cat(paste("   PVA decline is", round(difference_decline2 * 100, 2), "% less severe than historical.\n"))
  } else {
    cat(paste("   PVA decline is", round(abs(difference_decline2) * 100, 2), "% more severe than historical.\n"))
  }
} else {
  cat("   Historical decline rate for 2010-2018 could not be calculated (missing data).\n")
}

# --- 3. Compare Historical vs. Scenario 3 (1990-2010) ---
cat("\n3) Comparing 1990-2010 (Historical vs. Scenario 3):\n")
if (!is.na(growth_1990_2010)) {
  cat(paste("   Historical:", round(growth_1990_2010 * 100, 2), "% annual growth\n"))
  cat(paste("   PVA (Scenario 3):", round(pva_growth_1990_2010_s3 * 100, 2), "% annual growth\n"))
  difference_growth3 <- pva_growth_1990_2010_s3 - growth_1990_2010
  if (abs(difference_growth3) < 0.01) {
    cat("   PVA growth is roughly hitting the historical target.\n")
  } else if (difference_growth3 > 0) {
    cat(paste("   PVA growth is", round(difference_growth3 * 100, 2), "% higher than historical.\n"))
  } else {
    cat(paste("   PVA growth is", round(abs(difference_growth3) * 100, 2), "% lower than historical.\n"))
  }
} else {
  cat("   Historical growth rate for 1990-2010 could not be calculated (missing data).\n")
}

# --- 4. Compare Historical vs. Scenario 4 (2010-2018) ---
cat("\n4) Comparing 2010-2018 (Historical vs. Scenario 4):\n")
if (!is.na(decline_2010_2018)) {
  cat(paste("   Historical:", round(decline_2010_2018 * 100, 2), "% annual decline\n"))
  cat(paste("   PVA (Scenario 4):", round(pva_decline_2010_2018_s4 * 100, 2), "% annual decline\n"))
  difference_decline4 <- pva_decline_2010_2018_s4 - decline_2010_2018
  if (abs(difference_decline4) < 0.01) {
    cat("   PVA decline is roughly hitting the historical target.\n")
  } else if (difference_decline4 > 0) {
    cat(paste("   PVA decline is", round(difference_decline4 * 100, 2), "% less severe than historical.\n"))
  } else {
    cat(paste("   PVA decline is", round(abs(difference_decline4) * 100, 2), "% more severe than historical.\n"))
  }
} else {
  cat("   Historical decline rate for 2010-2018 could not be calculated (missing data).\n")
}
