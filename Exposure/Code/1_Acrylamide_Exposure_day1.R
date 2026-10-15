# =============================================================================
# Title:  Acrylamide_Exposure
# Autor:  Trine Husøy
# =============================================================================


# =============================================================================
# Activate packages
# =============================================================================

library(ggplot2)
library(ggpubr)
library(writexl)
#library(MASS)
#library(survival)
library(fitdistrplus)
#library(Rcpp)
library(openxlsx)
library(tidyverse)
library(kableExtra)


# =============================================================================
# Get work directory and organise results
# =============================================================================

# Run the script with the working directory set to the Exposure folder of the
# repository (the folder that contains Code/ and Data/).
HOME <- getwd()
if (!dir.exists(file.path(HOME, "Data")) || !dir.exists(file.path(HOME, "Code"))) {
  stop("Set the working directory to the Exposure folder (containing Code/ and Data/).")
}

# Create a results folder named date_commit, so the results can be tied to the
# committed version of the script. "-dirty" is added if the script has
# uncommitted changes.
git_hash <- tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) "nogit", warning = function(w) "nogit")
git_dirty <- tryCatch(
  length(system2("git", c("status", "--porcelain", "--", "Code"), stdout = TRUE, stderr = FALSE)) > 0,
  error = function(e) FALSE, warning = function(w) FALSE)
run_id <- paste0(Sys.Date(), "_", git_hash, if (git_dirty) "-dirty" else "")

newday <- file.path(HOME, "Results", run_id)
dir.create(newday, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# Read in data
# =============================================================================

ConcAcrylamide <- read_csv2("Data/NO_2014_akrylamid_Trine_240522.csv")
FoodIntakeBothDays <- read_csv2("Data/1_inntak_per_måltid_BothDays_300721.csv")
SexWeight <- read_csv2("Data/1-sex_weight.csv")


# =============================================================================
# Functions
# =============================================================================

# Calculate the location (Loc) and scale (shape) parameters of the log normal distribution

Lognorm_Loc <- function(x){
  log(
    x$mean^2/sqrt(x$sd^2+x$mean^2)
  )
}

Lognorm_shape <- function(x){
  sqrt(
    log(
      1 + x$sd^2/x$mean^2
    )
  )
}

# Calculate summary data
Sum <- function(x){
  x %>%
  group_by(FoodCat) %>%
  summarise(
    N = n(),
    mean = mean(Conc, na.rm=TRUE),
    sd=sd(Conc, na.rm=TRUE),
    min=min(Conc, na.rm=TRUE),
    P05=quantile(Conc, .05, na.rm=TRUE),
    P50=quantile(Conc, .50, na.rm=TRUE),
    P95=quantile(Conc, .95, na.rm=TRUE),
    max=max(Conc, na.rm=TRUE)
  )
}

Sum_1 <- function(x){
  x %>%
  group_by(IDkode) %>%
  summarise(
    N = n(),
    mean = mean(value, na.rm=TRUE),
    sd=sd(value, na.rm=TRUE),
    min=min(value, na.rm=TRUE),
    P05=quantile(value, .05, na.rm=TRUE),
    P50=quantile(value, .50, na.rm=TRUE),
    P95=quantile(value, .95, na.rm=TRUE),
    max=max(value, na.rm=TRUE)
  )
}

# The MC_Sim function performs the Monte Carlo (MC) simulation.
# It creates an empty array (x) and, for each iteration, food category and
# person, draws an independent concentration from the lognormal distribution
# defined by Loc and Shape (z) and multiplies it with the person's amount (y).
# The MC is repeated 1000 times.
#
# To have the same value in different trials, addition of a "seed" using
# "set.seed(123)" is needed, where "123" could be any number, but always the same.

# y: intake table with IDkode, BW and one column per food category (g/day)
# z: concentration summary with FoodCat, Loc and Shape
# Food categories are matched by name (not by column position), and the number
# of categories/persons is taken from the data.
MC_Sim <- function(y, z, MC = 1000, seed = 123){
  set.seed(seed)
  cats <- z$FoodCat
  stopifnot(all(cats %in% names(y)))
  amounts <- as.matrix(y[, cats])
  x <- array(dim = c(nrow(y), length(cats), MC),
             dimnames = list(y$IDkode, cats, NULL))
  for (u in 1:MC){
    for (i in seq_along(cats)){
      # independent concentration draw for each person (n = number of persons)
      x[, i, u] <- amounts[, i] * rlnorm(n = nrow(y), z$Loc[i], z$Shape[i])
    }
  }
  return(x)
}


# =============================================================================
# Cleaning of data
# =============================================================================

# --- Sex and weight -----------------------------------------------------------

SexWeight[] <- SexWeight %>%
  lapply(function(x) as.numeric(x))


# --- Acrylamid concentrations in food ----------------------------------------

# The measurements from the Food Authorities are from 2014, 2020 and 2022 in
# the units microg/kg, which is the same as ng/g.

ConcAcrylamide <- ConcAcrylamide %>% # ug/kg = ng/g - select the variables you want to use
  dplyr::select(FoodCat,Conc)

stopifnot(!anyNA(ConcAcrylamide$Conc)) # no NAs in the current file (the old replace_na() call was never assigned, so it did nothing)

# Make a wide table with food categories as columns
ConcAcrylamide_wide <- ConcAcrylamide %>% mutate(observations = 1:n()) # define each entry as individual observations. Otherwise spread will not work
ConcAcrylamide_wide <- spread(ConcAcrylamide_wide, FoodCat, Conc) # A lot of NAs, but will be removed when making a vector

Biscuits <- as.vector(ConcAcrylamide_wide$Biscuits)
Biscuits <- Biscuits[!is.na(Biscuits)] # remove the NA's

Crisps <- as.vector(ConcAcrylamide_wide$Crisps)
Crisps <- Crisps[!is.na(Crisps)] # remove the NA's

Crispbread <- as.vector(ConcAcrylamide_wide$Crispbread)
Crispbread <- Crispbread[!is.na(Crispbread)] # remove the NA's


# --- Check the distribution of the acrylamide concentrations (several distributions) ---

# Biscuits
# lognormal
fit.lnorm <- fitdist(Biscuits, "lnorm")
summary(fit.lnorm)$est
plot(fit.lnorm)

# gamma
fit.gamma <- fitdist(Biscuits, "gamma")
summary(fit.gamma)$est
#plot(fit.gamma)

# weibull
fit.weibull <- fitdist(Biscuits, "weibull")
summary(fit.weibull)$est
#plot(fit.weibull)

# Crisps
# lognormal
fit.lnorm <- fitdist(Crisps, "lnorm")
summary(fit.lnorm)$est
plot(fit.lnorm)

# gamma
fit.gamma <- fitdist(Crisps, "gamma")
summary(fit.gamma)$est
#plot(fit.gamma)

# weibull
fit.weibull <- fitdist(Crisps, "weibull")
summary(fit.weibull)$est
#plot(fit.weibull)

# vanlig normalfordeling
fit.norm <- fitdist(Crisps, "norm")
summary(fit.norm)$est
#plot(fit.norm)

# Crispbread
# lognormal
fit.lnorm <- fitdist(Crispbread, "lnorm")
summary(fit.lnorm)$est
plot(fit.lnorm)

# gamma
fit.gamma <- fitdist(Crispbread, "gamma")
summary(fit.gamma)$est
#plot(fit.gamma)

# vanlig normalfordeling
fit.norm <- fitdist(Crispbread, "norm")
summary(fit.norm)$est
#plot(fit.norm)


# --- Make summary data of the acrylamide concentrations ----------------------

SumConcAcrylamide <- Sum(ConcAcrylamide)

knitr::kable(SumConcAcrylamide, digits = 1) %>% kable_styling(font_size = 10, bootstrap_options = c("striped"))

# Delete the FoodCat that will not be used in the exposure assessment later
SumConcAcrylamide <- SumConcAcrylamide[!SumConcAcrylamide$FoodCat == "Mixed", ]
SumConcAcrylamide <- SumConcAcrylamide[!SumConcAcrylamide$FoodCat == "CoffeeBean", ]
SumConcAcrylamide <- SumConcAcrylamide[!SumConcAcrylamide$FoodCat == "Crisps_veg", ]

# Calculate the location and shape parameter from mean and SD defining the lognormal distribution
SumConcAcrylamide$Loc <- SumConcAcrylamide %>% Lognorm_Loc()
SumConcAcrylamide$Shape <- SumConcAcrylamide %>% Lognorm_shape()

# Save the file
write.xlsx(SumConcAcrylamide,
           file = file.path(newday,"SumConcAcrylamide.xlsx"),
           colNames = TRUE, borders = "rows"
)


# --- Food intake -------------------------------------------------------------

# The food intake from the diaries are in g/day
names(FoodIntakeBothDays)[names(FoodIntakeBothDays) == 'IDKode'] <- 'IDkode' # renames IDKode to IDkode. my input file has "IDKode", which =!"IDkode", like the rest of the files

FoodIntakeBothDays <- FoodIntakeBothDays %>% dplyr::select(-Meal, -Time, -Name)

FoodIntakeDay1 <- FoodIntakeBothDays[FoodIntakeBothDays$Day==1,]

FoodIntakeDay1 <- FoodIntakeDay1 %>% dplyr::select(-Day)
FoodIntakeDay1$Amount <- FoodIntakeDay1$Amount %>% as.numeric()
FoodIntakeDay1[is.na(FoodIntakeDay1)] <- 0

# Sum all food intake for each food category for all meals for each individual
SumFoodIntakeDay1 <- FoodIntakeDay1 %>%
  group_by(IDkode, FoodCat) %>%
  summarize(Amount = sum(Amount))

# Make a wide table from the long table
SumFoodIntakeDay1 <- SumFoodIntakeDay1 %>% pivot_wider(names_from = FoodCat, values_from = Amount)
SumFoodIntakeDay1[is.na(SumFoodIntakeDay1)] <- 0

#SumFoodIntakeDay1 <- dplyr::select(SumFoodIntakeDay1, -c(Mussel))

# Merge with the file SexWeight
SumFoodIntakeDay1 <- merge(SumFoodIntakeDay1, SexWeight, by="IDkode", all = TRUE) # TH: included all = TRUE. Then the non-matching IDs will end up as NAs. (This is then "18, 82, 83 and 120 as expected)
SumFoodIntakeDay1 <- na.omit(SumFoodIntakeDay1) # Remove rows with NAs

# Select the variables for which we have acrylamide concentrations
SumFoodIntakeDay1 <- dplyr::select(SumFoodIntakeDay1, c(IDkode, BW = con_weight, Biscuits, Bread, Cereal, Coffee, Crispbread, Crisps, FrenchFries, PotatoFried, PotatoTortilla))

knitr::kable(head(SumFoodIntakeDay1), digits = 2) %>% kable_styling(font_size = 10, bootstrap_options = c("striped"))


# =============================================================================
# Probabilistic exposure estimates from the diaries
# =============================================================================

Acrylamide_Day1 <- MC_Sim(SumFoodIntakeDay1, SumConcAcrylamide, MC = 1000) / 1000 # g x ng/g = ng; divide by 1000 to get ug/day

print(Acrylamide_Day1[1,,7])

# Add acrylamide exposure (in ug) from all food categories to a total exposure
# per person from 1000 iterations
Acrylamide_Day1_total <- as_tibble(apply(Acrylamide_Day1, 3, rowSums), .name_repair = "minimal")
names(Acrylamide_Day1_total) <- as.character(1:1000)

Acrylamide_Day1_total$IDkode <- SumFoodIntakeDay1$IDkode # Get IDs

# Make a long table of the total exposure to make summary data for each ID
Acrylamide_sum <- pivot_longer(Acrylamide_Day1_total, cols = -IDkode, names_to = "MC", values_to = "value")

# Make summary data of total acrylamide intake for each ID based on the diaries
SumAcrylamideDay1 <- Sum_1(Acrylamide_sum) # ug/day
knitr::kable(head(SumAcrylamideDay1), digits = 2) %>% kable_styling(font_size = 10, bootstrap_options = c("striped"))

# Save the file
write.xlsx(Acrylamide_Day1_total,
           file = file.path(newday, "Acrylamide_Day1_total.xlsx"),
           colNames = TRUE, borders = "rows"
)

write.xlsx(SumAcrylamideDay1,
           file = file.path(newday, "SumAcrylamideDay1.xlsx"),
           colNames = TRUE, borders = "rows"
)


# =============================================================================
# COMMENTED OUT SECTIONS (not used in current analysis)
# =============================================================================

# --- Calculate the contribution from the different food groups from the diaries ---

# Take the average over all MC for each food category for each individual.

# Acrylamide_Diary_FoodCat <- as_tibble(matrix(NA, nrow = nrow(Acrylamide_Diary), ncol = 9)) #ev: mine does not have 124 rows
# for (i in 1:9){
#   Acrylamide_Diary_FoodCat[,i] <- rowSums(Acrylamide_Diary[,i,], na.rm = TRUE)/1000
# }
#
# ColName_FoodCat <- names(SumFoodIntakeFFQ[,3:11])  # The columns named with food categories
# ColName_FoodCat
# colnames(Acrylamide_Diary_FoodCat) <- ColName_FoodCat
#
# Acrylamide_Diary_FoodCat$IDkode <- SumFoodIntakeBothDays[,1]
#
# Acrylamide_Diary_FoodCat <- gather(Acrylamide_Diary_FoodCat, FoodCat, Conc, Biscuits:PotatoTortilla, factor_key = TRUE)


# --- Make summary data for each food category ---

# SumAcrylamide_Diary_FoodCat <- Sum(Acrylamide_Diary_FoodCat)


# Save the file
# write.xlsx(SumAcrylamide_Diary_FoodCat,
#            file = file.path(newday,"SumAcrylamide_Diary_FoodCat.xlsx"),
#            colNames = TRUE, borders = "rows"
# )


# --- Box-plot of all IDs per food category ---

# Plot_acrylamide_Diary_FoodCat <- ggplot(data = Acrylamide_Diary_FoodCat, aes(y=Conc, x=FoodCat)) +
#   geom_boxplot(fill="grey")+
#   theme_minimal()+
#   scale_y_log10(limits=c(0.5,60), labels = scales::comma)+
#   theme(axis.text.x = element_text(angle = 90, size = 20), axis.title = element_text(size = 24),
#         axis.text.y = element_text(size = 20),legend.position= "none")+
#   xlab("Food category")+
#   ylab("Acrylamide exposure (ug/day)")+
#   annotate("text", x=9, y=40, label="Diary", size=9)
#
# Plot_acrylamide_Diary_FoodCat
#
# ggsave(filename=file.path(newday, "Plot_acrylamide_Diary_FoodCat.jpeg"),
#        device = NULL,
#        width=NA,
#        height=NA,
#        units="mm")


# --- Make cumulative density plots ---

# Estimated exposure from adducts compared to the mean exposure from the diaries and FFQ
# #ev: i do not have the same length for B, C & D as stated here (128) so i've changed it to count the length of the df
#
# A = c(Exposure_all$ADDUCT, Exposure_all$mean_diary, Exposure_all$mean_FFQ)
# B = replicate(nrow(Exposure_all),"Adduct")
# C = replicate(nrow(Exposure_all),"mean_diary")
# D = replicate(nrow(Exposure_all),"mean_FFQ")
# E = c(B,C,D)
#
# Acrylamide_mean <- as.data.frame(matrix(NA,nrow = nrow(Exposure_all)*3, ncol = 2))
# Acrylamide_mean[1] <- A
# Acrylamide_mean[2] <- E
#
# names(Acrylamide_mean)[1] <- "Acrylamide"
# names(Acrylamide_mean)[2] <- "Exposure"
#
# Acrylamide_mean <- as.data.frame(Acrylamide_mean)


# Plot_Acrylamide_mean <- ggplot(data = Acrylamide_mean, aes(Acrylamide)) +
#   geom_line(stat = "ecdf", size = 0.8, aes(linetype=Exposure))+
#   theme_minimal()+
#   scale_x_continuous(limits=c(0.001, 1.5))+
#   theme(axis.text.x = element_text(angle = 90, size = 20), axis.title = element_text(size = 20),
#         axis.text.y = element_text(size = 20), legend.key.size = unit(1, "cm"), legend.text = element_text(size =20), legend.title = element_blank())+
#   xlab("Acrylamide (ug/kg bw/day)")+
#   ylab("Cumulative probability")
#
# Plot_Acrylamide_mean
#
# ggsave(filename=file.path(newday,"Plot_Acrylamide_mean.jpeg"),
#        device = NULL,
#        width=NA,
#        height=NA,
#        units="mm")


# Estimated exposure from adducts compared to the median (P50) exposure from the diaries and FFQ
# #ev: again, i have different length...
# G = c(Exposure_all$ADDUCT, Exposure_all$P50_diary, Exposure_all$P50_FFQ)
# H = replicate(nrow(Exposure_all),"Adduct")
# I = replicate(nrow(Exposure_all),"P50_diary")
# J = replicate(nrow(Exposure_all),"P50_FFQ")
# K = c(H,I,J)
#
# Acrylamide_P50 <- as.data.frame(matrix(NA,nrow = nrow(Exposure_all)*3, ncol = 2))
# Acrylamide_P50[1] <- G
# Acrylamide_P50[2] <- K
#
# names(Acrylamide_P50)[1] <- "Acrylamide"
# names(Acrylamide_P50)[2] <- "Exposure"
#
# Acrylamide_P50 <- as.data.frame(Acrylamide_P50)


# Plot_Acrylamide_P50 <- ggplot(data = Acrylamide_P50, aes(Acrylamide)) +
#   geom_line(stat = "ecdf", size = 0.8, aes(linetype=Exposure))+
#   theme_minimal()+
#   scale_x_continuous(limits=c(0.001, 1.3))+
#   theme(axis.text.x = element_text(angle = 90, size = 20), axis.title = element_text(size = 20),
#         axis.text.y = element_text(size = 20), legend.key.size = unit(1, "cm"), legend.text = element_text(size =20), legend.title = element_blank())+
#   xlab("Acrylamide (ug/kg bw/day)")+
#   ylab("Cumulative probability")
#
# Plot_Acrylamide_P50
#
# ggsave(filename=file.path(newday,"Plot_Acrylamide_P50.jpeg"),
#        device = NULL,
#        width=NA,
#        height=NA,
#        units="mm")


# --- Bland Altman Plots ---

# Adducts versus P50 from the diaries
# Exposure_all <- Exposure_all %>% mutate(avg_P50_diary = rowMeans(Exposure_all %>% dplyr::select("ADDUCT", "P50_diary"))) # calculate means of the two measurements for each ID
# Exposure_all$diff_P50_diary <- Exposure_all$ADDUCT - Exposure_all$P50_diary # calculate difference in measurements for each ID
#
# Exposure_all <- merge(Exposure_all, SexWeight, by="IDkode", all=TRUE) %>%
#   na.omit() %>%
#   as_tibble()
#
# Exposure_all$smoking <- ifelse(Exposure_all$cat_Smoking==3, "Smoking", "Not smoking") # Make a variable gender to use in the Bland-Altman plot
#
# mean_diff_diary <- mean(Exposure_all$diff_P50_diary) # calculate average difference
# mean_diff_diary
# mean_sd_diary <- sd(Exposure_all$diff_P50_diary)
# mean_sd_diary
# lower_diff_diary <- mean_diff_diary - 1.96*sd(Exposure_all$diff_P50_diary) # calculate lower 95% conf int
# lower_diff_diary
# upper_diff_diary <- mean_diff_diary + 1.96*sd(Exposure_all$diff_P50_diary) # calculate upper 95% conf int
# upper_diff_diary
#
# Plot_BA_adduct_DiaryP50 <- ggplot(Exposure_all, aes(x = avg_P50_diary, y = diff_P50_diary)) +
#   geom_point(shape = 21, size=3, aes(fill = factor(smoking))) +
#   scale_fill_grey(start = 0.5, end = .9)+
#   annotate("rect", xmin = 0, xmax=0.75,  ymin = mean_diff_diary-mean_sd_diary, ymax = mean_diff_diary+mean_sd_diary, alpha = .1, fill = "grey" )+
#   geom_hline(yintercept = mean_diff_diary) +
#   geom_hline(yintercept = lower_diff_diary, color = "black", linetype="dashed") +
#   geom_hline(yintercept = upper_diff_diary, color = "black", linetype="dashed") +
#   ggtitle("Acrylamide adduct versus P50 from diary") +
#   ylab("Difference Between Measurements") +
#   xlab("Average Measurement")+
#   theme( panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "light grey", size = 0.2), plot.title = element_text(size = 22, face = "bold"), axis.text=element_text(size=22), axis.title=element_text(size=22),
#         legend.title = element_blank())
#
# Plot_BA_adduct_DiaryP50
#
# ggsave(filename=file.path(newday, "Plot_BA_adduct_DiaryP50.jpeg"),
#        device = NULL,
#        width=NA,
#        height=NA,
#        units="mm")
