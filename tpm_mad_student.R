# =============================================================================
# TPM Estimation: Minimum Absolute Deviations (MAD)
# Lee, Judge & Zellner (1970)
# =============================================================================

# Load functions (run once; requires internet connection)
source(r"(C:\Users\KRUPESH SK\OneDrive\Documents\tpm_mad_functions.R)")

# -----------------------------------------------------------------------------
# STEP 1: Load your data
#   - Rows = years, Columns = crops, Values = crop areas
# -----------------------------------------------------------------------------
your_data <- read.csv("C:/Users/KRUPESH SK/OneDrive/Documents/markov-chain-data.csv", row.names = 1)

# -----------------------------------------------------------------------------
# STEP 2: Estimate TPM
# -----------------------------------------------------------------------------
result <- estimate_tpm_mad(your_data)

# -----------------------------------------------------------------------------
# STEP 3: View results
# -----------------------------------------------------------------------------
print_tpm_results(result)
