# =============================================================================
# Title:  Acrylamide_Exposure_day2
# Autor:  Trine Husøy
#
# Acrylamide exposure per meal on day 2 for the individuals with spot urine
# measurements. The meal-level exposures are input to the PBPK model.
#
# Note: only the food categories with acrylamide concentration data are included
# (see "Food categories" below). Intake of acrylamide from other foods is not
# estimated, so the exposure per meal is an UNDERESTIMATE.
# =============================================================================


# =============================================================================
# Activate packages
# =============================================================================

library(openxlsx)
library(tidyverse)


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

# The urine file is only used to identify the individuals with spot urine
# measurements. The urine values are not used here. The first line of the file
# is an extra header, and the IDs have inconsistent leading zeros (009, 40).
Urine <- read_csv2("Data/ADS-20-0109b_NIPH_urine_20201201.csv", skip = 1,
                   col_types = cols(.default = col_character()),
                   name_repair = "minimal")
Urine_IDs <- sort(unique(as.numeric(na.omit(Urine[[1]]))))


# =============================================================================
# Cleaning of data
# =============================================================================

SexWeight[] <- SexWeight %>%
  lapply(function(x) as.numeric(x))

names(FoodIntakeBothDays)[names(FoodIntakeBothDays) == 'IDKode'] <- 'IDkode'

stopifnot(!anyNA(ConcAcrylamide$Conc))


# =============================================================================
# Acrylamide concentrations in food (same as in the day 1 script)
# =============================================================================

# The measurements are in microg/kg, which is the same as ng/g.

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

SumConcAcrylamide <- ConcAcrylamide %>%
  dplyr::select(FoodCat, Conc) %>%
  group_by(FoodCat) %>%
  summarise(mean = mean(Conc, na.rm = TRUE), sd = sd(Conc, na.rm = TRUE)) %>%
  filter(!FoodCat %in% c("Mixed", "CoffeeBean", "Crisps_veg")) # not used in the exposure assessment

SumConcAcrylamide$Loc <- SumConcAcrylamide %>% Lognorm_Loc()
SumConcAcrylamide$Shape <- SumConcAcrylamide %>% Lognorm_shape()

# Food categories with acrylamide concentration data
FoodCat_conc <- SumConcAcrylamide$FoodCat


# =============================================================================
# Food intake day 2 for the individuals with urine measurements
# =============================================================================

# All items eaten on day 2 by the selected individuals (g)
FoodIntakeDay2 <- FoodIntakeBothDays %>%
  filter(Day == 2, IDkode %in% Urine_IDs) %>%
  mutate(Amount = replace_na(as.numeric(Amount), 0))

# Check that all selected individuals are in the diaries and have a body weight
stopifnot(all(Urine_IDs %in% FoodIntakeDay2$IDkode))
stopifnot(all(Urine_IDs %in% SexWeight$IDkode))

# One row per meal (individual x meal x time). Meals without any food category
# with concentration data are kept, with zero exposure, so that the time series
# for the PBPK model is complete. The meals are sorted by time within each
# individual (the meal numbers in the diaries are not always in time order);
# Meal is the original meal number, Time is hhmm, Clock hh:mm and Hours the
# time of day in hours.
Meals <- FoodIntakeDay2 %>%
  distinct(IDkode, Meal, Time) %>%
  arrange(IDkode, Time) %>%
  left_join(dplyr::select(SexWeight, IDkode, BW = con_weight), by = "IDkode") %>%
  mutate(Clock = sprintf("%02d:%02d", Time %/% 100, Time %% 100),
         Hours = Time %/% 100 + (Time %% 100) / 60)
stopifnot(!anyDuplicated(Meals[, c("IDkode", "Meal")])) # one time per meal

# Items in the food categories with concentration data
Items <- FoodIntakeDay2 %>%
  filter(FoodCat %in% FoodCat_conc) %>%
  left_join(dplyr::select(SumConcAcrylamide, FoodCat, Loc, Shape), by = "FoodCat")

cat("Individuals:", length(Urine_IDs), "\n",
    "Meals:", nrow(Meals), "\n",
    "Items on day 2:", nrow(FoodIntakeDay2), "of which with concentration data:", nrow(Items), "\n")


# =============================================================================
# Monte Carlo simulation per meal
# =============================================================================

# For each item and iteration an independent concentration is drawn from the
# lognormal distribution of the food category (as for day 1). Amount (g) x
# concentration (ng/g) = ng; divided by 1000 to get ug. The items are summed
# to the exposure per meal.

MC <- 1000
set.seed(123)

Exposure_items <- Items$Amount / 1000 *
  matrix(rlnorm(nrow(Items) * MC, rep(Items$Loc, MC), rep(Items$Shape, MC)),
         nrow = nrow(Items), ncol = MC) # ug, items x MC

# Sum the items to meals (rows in the same order as Meals, zero for meals without items)
meal_key <- paste(Meals$IDkode, Meals$Meal)
Exposure_meal <- matrix(0, nrow = nrow(Meals), ncol = MC, dimnames = list(meal_key, NULL))
item_sum <- rowsum(Exposure_items, paste(Items$IDkode, Items$Meal))
Exposure_meal[rownames(item_sum), ] <- item_sum

# Total exposure per person on day 2 (all meals) for each iteration
Exposure_day <- rowsum(Exposure_meal, Meals$IDkode)

# Summary over the MC iterations
Summarise_MC <- function(m) {
  tibble(
    mean = rowMeans(m),
    sd = apply(m, 1, sd),
    P05 = apply(m, 1, quantile, 0.05),
    P50 = apply(m, 1, quantile, 0.50),
    P95 = apply(m, 1, quantile, 0.95)
  )
}

SumAcrylamideDay2_meal <- bind_cols(Meals, Summarise_MC(Exposure_meal)) # ug per meal

SumAcrylamideDay2_total <- bind_cols(
  tibble(IDkode = as.numeric(rownames(Exposure_day))),
  Summarise_MC(Exposure_day)) %>% # ug/day
  left_join(dplyr::select(SexWeight, IDkode, BW = con_weight), by = "IDkode") %>%
  mutate(mean_ug_kg = mean / BW) # ug/kg bw/day

print(head(SumAcrylamideDay2_meal, 10))
print(SumAcrylamideDay2_total)


# =============================================================================
# Additional calculation: meals merged within a time window
# =============================================================================

# Items registered close in time are merged into one larger meal. The window
# starts at the first meal of a group: a meal belongs to the group if it is at
# most Window_min minutes after the first meal of the group, otherwise a new
# group is started. The total per person is the same as above; only the time
# resolution changes. The dose time of a merged meal is the time of the first
# meal (Hours); Hours_wmean is the mean time weighted by the mean exposure of
# the meals (equal to the first time if all meals have zero exposure).

Window_min <- 30

Assign_group <- function(minutes, window) {
  group <- integer(length(minutes))
  g <- 0
  start <- -Inf
  for (i in seq_along(minutes)) { # minutes must be sorted
    if (minutes[i] - start > window) {
      g <- g + 1
      start <- minutes[i]
    }
    group[i] <- g
  }
  group
}

Meals_win <- Meals %>%
  mutate(Mean_ug = SumAcrylamideDay2_meal$mean) %>%
  group_by(IDkode) %>%
  mutate(Group = Assign_group(Hours * 60, Window_min)) %>%
  ungroup()

Meals_merged <- Meals_win %>%
  group_by(IDkode, Group) %>%
  summarise(
    N_meals = n(),
    Meals = paste(Meal, collapse = ","),
    BW = first(BW),
    Time = first(Time),
    Clock = first(Clock),
    # Hours_wmean and Hours_last are calculated before Hours is overwritten by the first time
    Hours_wmean = if (sum(Mean_ug) > 0) weighted.mean(Hours, Mean_ug) else first(Hours),
    Hours_last = last(Hours),
    Hours = first(Hours),
    .groups = "drop"
  )

# Sum the exposure of the meals in each group (rows in the same order as Meals_merged)
Exposure_merged <- rowsum(Exposure_meal, paste(Meals_win$IDkode, Meals_win$Group))
Exposure_merged <- Exposure_merged[paste(Meals_merged$IDkode, Meals_merged$Group), , drop = FALSE]

SumAcrylamideDay2_merged <- bind_cols(Meals_merged, Summarise_MC(Exposure_merged)) # ug per merged meal

# The total per person must be unchanged by merging
stopifnot(isTRUE(all.equal(
  as.numeric(rowsum(Exposure_merged, Meals_merged$IDkode)),
  as.numeric(Exposure_day))))

cat("Meals:", nrow(Meals), "-> merged meals (window", Window_min, "min):", nrow(Meals_merged), "\n")


# =============================================================================
# Save the results (individual-level data: saved in Results/ only, not tracked by git)
# =============================================================================

write.xlsx(SumAcrylamideDay2_meal,
           file = file.path(newday, "SumAcrylamideDay2_meal.xlsx"),
           colNames = TRUE, borders = "rows"
)

write.xlsx(SumAcrylamideDay2_total,
           file = file.path(newday, "SumAcrylamideDay2_total.xlsx"),
           colNames = TRUE, borders = "rows"
)

write.xlsx(SumAcrylamideDay2_merged,
           file = file.path(newday, paste0("SumAcrylamideDay2_merged_", Window_min, "min.xlsx")),
           colNames = TRUE, borders = "rows"
)

# All MC iterations per meal (ug), with the meal table, as input to the PBPK model
saveRDS(list(meals = Meals, exposure_ug = Exposure_meal),
        file = file.path(newday, "Acrylamide_Day2_meal_MC.rds"))

# The same for the meals merged within the time window
saveRDS(list(meals = Meals_merged, window_min = Window_min, exposure_ug = Exposure_merged),
        file = file.path(newday, paste0("Acrylamide_Day2_merged_", Window_min, "min_MC.rds")))
