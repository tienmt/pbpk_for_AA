# Acrylamide exposure (input to the PBPK model)

Probabilistic estimate of dietary acrylamide exposure (ug/day) for the EuroMix
human study, based on 24-hour food diaries and measured acrylamide
concentrations in Norwegian foods.

## Script

`Code/1_Acrylamide_Exposure_day1.R` - exposure for day 1 (matches the 24-hour
urine collection of day 1).

Method: for each person, food category and Monte Carlo iteration (1000), an
independent concentration is drawn from a lognormal distribution (fitted from
mean and SD of the measured concentrations) and multiplied with the person's
daily intake of that food category. Only the food categories with acrylamide
concentration data are included (Biscuits, Bread, Cereal, Coffee, Crispbread,
Crisps, FrenchFries, PotatoFried, PotatoTortilla), so the exposure is an
**underestimate** of total dietary exposure.

Run the script with the working directory set to this `Exposure/` folder.
Results are written to `Results/<date>_<commit>/` (`-dirty` is added if the
script has uncommitted changes). `Results/` is not version controlled.

## Data

| File | In repository |
|------|---------------|
| `Data/NO_2014_akrylamid_Trine_240522.csv` - acrylamide concentrations in food (ug/kg), 2014-2022 | yes |
| `Data/1_inntak_per_måltid_BothDays_300721.csv` - food intake per meal, both days | **no** (human study data) |
| `Data/1-sex_weight.csv` - sex, body weight and smoking per participant | **no** (human study data) |

The human study data are kept locally only and are excluded via
`.git/info/exclude`. Contact the author for access.
