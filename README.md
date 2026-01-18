## antipsychotics_dementia
Code for data analysis using time-varying antipsychotic exposure - Cox PH model


## Study Overview
The goal of this project is to estimate the time‑varying hazard of death following antipsychotic initiation in a dementia cohort. The analysis uses:

- Multiple imputation (MICE) for missing baseline covariates

- Discrete time intervals (3‑month windows) for time‑varying covariates

- Episode splitting at the exact date of antipsychotic initiation

- Delayed entry methods

- Pooling across imputed datasets using Rubin’s rules

The mock datasets allow users to run the full pipeline without accessing protected health information.

**Data Description**

**1. baseline_mock**
A synthetic baseline dataset containing: demographic variables (age, sex, education), cognitive and psychiatric measures (MMSE, BDI‑II), 
comorbidities and medication use, key dates (distorted data, converted to synthetic shifted dates)

**2. long_mock** 
A synthetic longitudinal dataset containing: repeated measures of time‑varying covariates, medication exposures, follow‑up intervals, synthetic dates aligned with baseline_mock. 

### Privacy and Ethical Considerations

No real patient‑level data are included in this repository.
All mock data are synthetic, non‑identifiable, and cannot be linked to any real individual.
The mock datasets are intended solely for code reproducibility, not for scientific inference.

### Contact
For questions about the analysis or data structure, please contact:
Augusto Ferraris, MD, MPH - UW, Seattle, WA | Email: aferra@uw.edu
