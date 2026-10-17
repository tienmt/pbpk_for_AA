# =============================================================================
# Title:  Acrylamide_Exposure
# Autor:  Trine Husøy
# =============================================================================


# =============================================================================
# Activate packages
# =============================================================================

library(fitdistrplus)
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
# Figures
# =============================================================================

# --- Box-plot of all IDs per food category -----------------------------------

# Mean over the MC iterations for each person and food category (ug/day)
Acrylamide_Day1_FoodCat <- as_tibble(apply(Acrylamide_Day1, c(1, 2), mean), .name_repair = "minimal") %>%
  mutate(IDkode = SumFoodIntakeDay1$IDkode) %>%
  pivot_longer(cols = -IDkode, names_to = "FoodCat", values_to = "Conc") %>%
  filter(Conc > 0) # persons who did not eat the food category (cannot be shown on a log scale)

Plot_acrylamide_Day1_FoodCat <- ggplot(Acrylamide_Day1_FoodCat, aes(x = FoodCat, y = Conc, fill = FoodCat)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  # fixed y-axis range with room above the highest box for the label
  scale_y_log10(limits = c(0.05, 500), breaks = c(0.1, 1, 10, 100), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        legend.position = "none") +
  xlab("Food category") +
  ylab("Acrylamide exposure (ug/day)") +
  annotate("text", x = -Inf, y = Inf, label = "Diary, day 1", hjust = -0.1, vjust = 1.5, size = 9)

Plot_acrylamide_Day1_FoodCat

ggsave(filename = file.path(newday, "Plot_acrylamide_Day1_FoodCat.jpeg"),
       plot = Plot_acrylamide_Day1_FoodCat, width = 250, height = 200, units = "mm")


# --- Cumulative distribution of total exposure per kg body weight ------------

# Total exposure (ug/kg bw/day) for each person and MC iteration
Acrylamide_Day1_total_bw <- sweep(as.matrix(Acrylamide_Day1_total[, as.character(1:1000)]),
                                  1, SumFoodIntakeDay1$BW, "/")

# Person-level P05, P50 and P95 over the MC iterations
Acrylamide_Day1_bw <- tibble(
  IDkode = SumFoodIntakeDay1$IDkode,
  P05 = apply(Acrylamide_Day1_total_bw, 1, quantile, 0.05),
  P50 = apply(Acrylamide_Day1_total_bw, 1, quantile, 0.50),
  P95 = apply(Acrylamide_Day1_total_bw, 1, quantile, 0.95)
) %>%
  filter(P95 > 0) %>% # persons with no intake of any of the food categories on day 1
  pivot_longer(cols = -IDkode, names_to = "Percentile", values_to = "Acrylamide")

Plot_Acrylamide_Day1_ecdf <- ggplot(Acrylamide_Day1_bw, aes(Acrylamide, colour = Percentile)) +
  geom_line(stat = "ecdf", linewidth = 1.2) +
  scale_colour_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 20), axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20), legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 20), legend.title = element_blank()) +
  xlab("Acrylamide (ug/kg bw/day)") +
  ylab("Cumulative probability")

Plot_Acrylamide_Day1_ecdf

ggsave(filename = file.path(newday, "Plot_Acrylamide_Day1_ecdf.jpeg"),
       plot = Plot_Acrylamide_Day1_ecdf, width = 250, height = 200, units = "mm")


# =============================================================================
# Persons with highest and lowest exposure
# =============================================================================

# Ranked by the mean over the MC iterations of total exposure (ug/kg bw/day).
# The lowest are selected among persons with exposure > 0 (persons with no
# intake of any of the food categories on day 1 have exposure = 0).
Exposure_per_person <- tibble(
  IDkode = SumFoodIntakeDay1$IDkode,
  BW = SumFoodIntakeDay1$BW,
  mean = rowMeans(Acrylamide_Day1_total_bw),
  P05 = apply(Acrylamide_Day1_total_bw, 1, quantile, 0.05),
  P50 = apply(Acrylamide_Day1_total_bw, 1, quantile, 0.50),
  P95 = apply(Acrylamide_Day1_total_bw, 1, quantile, 0.95)
)

Exposure_extremes <- bind_rows(
  Exposure_per_person %>% slice_max(mean, n = 3) %>% mutate(Group = "Highest"),
  Exposure_per_person %>% filter(mean > 0) %>% slice_min(mean, n = 3) %>% mutate(Group = "Lowest")
) %>%
  relocate(Group)

knitr::kable(Exposure_extremes, digits = 3) %>% kable_styling(font_size = 10, bootstrap_options = c("striped"))

# Contains person-level data: saved in Results/ only (not tracked by git)
write.xlsx(Exposure_extremes,
           file = file.path(newday, "Exposure_extremes_ug_kg_bw.xlsx"),
           colNames = TRUE, borders = "rows"
)
