import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

# --------------------------------------------------
# Experimental constants
# --------------------------------------------------

TOTAL_VOLUME_ML = 5.00
SAMPLE_VOLUME_ML = 0.50 
LOADED_MASS_MG = 0.107296137 

CALIBRATION_SLOPE = 662.29 
CALIBRATION_INTERCEPT = 0.025 

# --------------------------------------------------
# Load experimental data
# --------------------------------------------------

data = pd.read_csv("release_data.csv")

required_columns = ["time_min", "absorbance"]

if not all(column in data.columns for column in required_columns):
    raise ValueError("CSV file is missing required columns.")

# --------------------------------------------------
# Data processing
# --------------------------------------------------

data["time_hours"] = data["time_min"] / 60.0

# Calculate concentration in mg/mL
data["concentration"] = (
    data["absorbance"]-CALIBRATION_INTERCEPT
    )/CALIBRATION_SLOPE

# Clean up the concentration column by eliminating any unrealistic negative values caused by baseline noise in the UV-Vis machine.
# Downstream functions like np.log() or np.sqrt() will break or output NaN errors if given negative concentration inputs.
data["concentration"] = data["concentration"].apply(lambda x: max(0, x)) # .apply() loops through every single number in the concentration column and transforms it using the function provided inside the parentheses.
# lambda x: creates a small, anonymous inline function where x represents each individual concentration value in that column, one at a time.

cumulative_mass =[] # Total mg of drug released up to each time point.
prev_sum = 0.0

for C in data["concentration"]:
    # Calculate the corrected total mass (M_n) at the current time point using two parts:
    # the mass of drug currently remaining inside the release vessel and
    # the total mass of drug that was already removed in all previous sampling steps combined
    M_n = (C * TOTAL_VOLUME_ML) + (SAMPLE_VOLUME_ML * prev_sum)
    cumulative_mass.append(M_n)
    prev_sum += C # Updates the running total of concentrations for the next iteration.

data["cumulative_mass_mg"] = cumulative_mass
data["percent_released"] = (data["cumulative_mass_mg"] / LOADED_MASS_MG) * 100

# Filter non-zero points for logarithmic/sqrt model fitting.
valid_data = data[data["time_hours"] > 0].copy()
t = valid_data["time_hours"].values # .values extracts the raw numerical data out of a pandas Series (or DataFrame column) and converts it into a NumPy array.
Mt_Minf = (valid_data["percent_released"] / 100.0).values # Kinetic equations require fractional release values rather than whole percentage integers.

# Filter for Korsmeyer-Peppas (Mt/Minf < 0.6)
kp_mask = (Mt_Minf > 0) & (Mt_Minf <= 0.6)
t_kp = t[kp_mask] # Contains only the time points where kp_mask is True
Mt_Minf_kp = Mt_Minf[kp_mask]

# 3. Kinetic Model Fits
# Zero Order: Mt/Minf = q0 + k0 * t
res_zero = linregress(t, Mt_Minf) # The slope gives k0 and the intercept gives q0

# First Order: ln(1 - Mt/Minf) = -k1 * t
res_first = linregress(t, np.log(np.maximum(1 - Mt_Minf, 1e-5))) # Use np.maximum to avoid log(0) issues, the slope gives -k1

# Higuchi: Mt/Minf = kH * sqrt(t)
res_higuchi = linregress(np.sqrt(t_kp), Mt_Minf_kp) # The slope gives kH

# Korsmeyer-Peppas: log(Mt/Minf) = log(k) + n * log(t)
res_kp = linregress(np.log(t_kp), np.log(Mt_Minf_kp)) # The slope gives n and the intercept gives log(k)

# 4. Display Results
print("--- Kinetic Model Fitting Results ---")
print(f"Zero-Order R²:      {res_zero.rvalue**2:.4f}")
print(f"First-Order R²:     {res_first.rvalue**2:.4f}")
print(f"Higuchi R²:         {res_higuchi.rvalue**2:.4f}")
print(f"Korsmeyer-Peppas R²: {res_kp.rvalue**2:.4f} (Release exponent n = {res_kp.slope:.4f})")

# 5. Plotting Cumulative Release Curve
plt.figure(figsize=(8, 5))
plt.plot(data["time_hours"], data["percent_released"], label="Experimental Data")
plt.xlabel("Time (hours)")
plt.ylabel("Cumulative Drug Release (%)")
plt.title("Pevonedistat Release from Cellulose Particles")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("pevonedistat_release_curve.png")
plt.show()