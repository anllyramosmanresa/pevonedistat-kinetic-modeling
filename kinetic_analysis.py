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

# The 0.5 mL sample was diluted 10x before UV-Vis measurement.
# This converts the measured concentration back to the
# concentration in the original release vessel.
DILUTION_FACTOR = 10.0

# --------------------------------------------------
# Loading experimental data
# --------------------------------------------------

data = pd.read_csv("release_data.csv")

required_columns = ["time_min", "absorbance"]

if not all(column in data.columns for column in required_columns):
    raise ValueError("CSV file is missing required columns.")

# --------------------------------------------------
# Data processing
# --------------------------------------------------

data["time_hours"] = data["time_min"] / 60.0

# Calculates concentration in mg/mL
data["concentration_raw_mg_ml"] = (
    data["absorbance"]-CALIBRATION_INTERCEPT
    )/CALIBRATION_SLOPE * DILUTION_FACTOR

# Cleans up the concentration column by eliminating any unrealistic negative values caused by baseline noise in the UV-Vis machine.
# Downstream functions like np.log() or np.sqrt() will break or output NaN errors if given negative concentration inputs.
# Keep the raw values above, but use non-negative values for cumulative-release calculations.
data["concentration"] = (
    data["concentration_raw_mg_ml"].clip(lower=0)
)

# Calculates cumulative released mass
cumulative_mass =[] 
previous_concentration_sum = 0.0

for concentration in data["concentration"]:
    # Calculate the corrected total mass at the current time point using two parts:
    # 1. The mass of drug currently remaining inside the release vessel.
    current_mass = concentration * TOTAL_VOLUME_ML
    # 2. The total mass of drug that was already removed in all previous sampling steps combined.
    previously_removed_mass = SAMPLE_VOLUME_ML * previous_concentration_sum

    total_mass = current_mass + previously_removed_mass
    cumulative_mass.append(total_mass)
    previous_concentration_sum += concentration # Updates the running total of concentrations for the next iteration.

data["cumulative_mass_mg"] = cumulative_mass

# This describes release relative to the original loaded mass.
data["percent_released"] = (
    data["cumulative_mass_mg"] / LOADED_MASS_MG
    ) * 100

# --------------------------------------------------
# Data Preparation for kinetic models
# --------------------------------------------------

# M_infinity is taken as the MAXIMUM experimentally observed release amount.
M_INF_MG = data["cumulative_mass_mg"].max()

M_INF_TIME_MIN = data.loc[
    data["cumulative_mass_mg"].idxmax(), "time_min"]

M_INF_PERCENT = (M_INF_MG / LOADED_MASS_MG) * 100

# Display M_infinity and related information
print("\n--- Experimental M_infinity ---")
print("-" * 40)

print(f"M_infinity = {M_INF_MG:.8f} mg")

print(f"Maximum observed release time = "
    f"{M_INF_TIME_MIN:.0f} min")

print(f"M_infinity = {M_INF_PERCENT:.2f}% "
    f"of initially loaded mass")

# Data filtering for kinetic model fitting
valid_data = data[data["time_hours"] > 0].copy()

t = valid_data["time_hours"].to_numpy() #.to_numpy() converts a Pandas Series into a clean, pure NumPy array.

fraction_released = (valid_data["cumulative_mass_mg"].to_numpy() / M_INF_MG)

valid_data["Mt_Minf"] = fraction_released # Kinetic equations require fractional release values rather than whole percentage integers.

# --------------------------------------------------
# Kinetic model fitting
# --------------------------------------------------

# Zero Order
# Mt/Minf = q0 + k0 * t
res_zero = linregress(t, fraction_released) # The slope gives k0 and the intercept gives q0

# First Order
# ln(1 - Mt/Minf) = -k1 * t
first_order_mask = (fraction_released < 1) & (fraction_released > 0) 

t_first_order = t[first_order_mask]
release_first = fraction_released[first_order_mask]
res_first = linregress(
    t_first_order,
    np.log(1 - release_first)) # The slope gives -k1

# Filter for Korsmeyer-Peppas and Higuchi models to only include data points where 0 < Mt/Minf <= 0.6
early_release_mask = (
    (fraction_released > 0)
    & (fraction_released <= 0.60))

t_early = t[early_release_mask]
release_early = fraction_released[early_release_mask]

# Higuchi
# Mt/Minf = kH * sqrt(t)
res_higuchi = linregress(np.sqrt(t_early), release_early) # The slope gives kH

# Korsmeyer-Peppas
# log(Mt/Minf) = log(k) + n * log(t)
res_kp = linregress(np.log(t_early), np.log(release_early)) # The slope gives n and the intercept gives log(k)

# --------------------------------------------------
# Results Display
# --------------------------------------------------

# Rate Constants Extraction
k0 = res_zero.slope
k1 = -res_first.slope
kH = res_higuchi.slope
k_kp = np.exp(res_kp.intercept)  # Converted back from log(k)

print("\n--- Complete Kinetic Parameters ---")
print("-" * 35)
print(f"Zero-Order:           R² = {res_zero.rvalue**2:.4f},  k0 = {k0:.6f} h⁻¹")
print(f"First-Order:          R² = {res_first.rvalue**2:.4f},  k1 = {k1:.6f} h⁻¹")
print(f"Higuchi:              R² = {res_higuchi.rvalue**2:.4f},  kH = {kH:.6f} h⁻⁰·⁵")
print(f"Korsmeyer-Peppas:     R² = {res_kp.rvalue**2:.4f},  k_KP = {k_kp:.6f},  n = {res_kp.slope:.4f}")

# --------------------------------------------------
# Primary Plot
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.plot(
    data["time_hours"], 
    data["percent_released"], 
    "o-",
    label="Experimental Data"
    )
plt.xlabel("Time (hours)")
plt.ylabel("Cumulative Drug Release (%)")
plt.title("Pevonedistat Release from Cellulose Particles")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()

plt.savefig("pevonedistat_release_curve.png", dpi=300)
plt.show()

# --------------------------------------------------
# 4-Panel Subplot
# --------------------------------------------------

t_smooth = np.linspace(t.min(), t.max(), 200)
t_early_smooth = np.linspace(t_early.min(), t_early.max(), 200)

fig, axs = plt.subplots(2, 2, figsize=(11, 10))

# Zero-Order
axs[0, 0].plot(t, fraction_released, "bo", label="Data")
axs[0, 0].plot(t_smooth, res_zero.intercept + res_zero.slope * t_smooth, "r--", label="Fit")
axs[0, 0].set_title(f"Zero-Order (R² = {res_zero.rvalue**2:.4f})")
axs[0, 0].set_xlabel("Time (h)")
axs[0, 0].set_ylabel("Mt / Minf")
axs[0, 0].grid(True, linestyle="--", alpha=0.5)
axs[0, 0].legend()

# First-Order
axs[0, 1].plot(t_first_order, np.log(1 - release_first), "bo", label="Data")
axs[0, 1].plot(t_smooth, res_first.intercept + res_first.slope * t_smooth, "r--", label="Fit")
axs[0, 1].set_title(f"First-Order (R² = {res_first.rvalue**2:.4f})")
axs[0, 1].set_xlabel("Time (h)")
axs[0, 1].set_ylabel("ln(1 - Mt / Minf)")
axs[0, 1].grid(True, linestyle="--", alpha=0.5)
axs[0, 1].legend()

# Higuchi
axs[1, 0].plot(np.sqrt(t_early), release_early, "bo", label="Data")
axs[1, 0].plot(np.sqrt(t_early_smooth), res_higuchi.intercept + res_higuchi.slope * np.sqrt(t_early_smooth), "r--", label="Fit")
axs[1, 0].set_title(f"Higuchi (R² = {res_higuchi.rvalue**2:.4f})")
axs[1, 0].set_xlabel("√Time (h^0.5)")
axs[1, 0].set_ylabel("Mt / Minf")
axs[1, 0].grid(True, linestyle="--", alpha=0.5)
axs[1, 0].legend()

# Korsmeyer-Peppas
axs[1, 1].plot(np.log(t_early), np.log(release_early), "bo", label="Data")
axs[1, 1].plot(np.log(t_early_smooth), res_kp.intercept + res_kp.slope * np.log(t_early_smooth), "r--", label="Fit")
axs[1, 1].set_title(f"Korsmeyer-Peppas (R² = {res_kp.rvalue**2:.4f}, n = {res_kp.slope:.4f})")
axs[1, 1].set_xlabel("ln(Time (h))")
axs[1, 1].set_ylabel("ln(Mt / Minf)")
axs[1, 1].grid(True, linestyle="--", alpha=0.5)
axs[1, 1].legend()

plt.tight_layout(pad=3.0)
plt.savefig("kinetic_models_fit.png", dpi=300)
plt.show()

# --------------------------------------------------
# Saves Kinetic Parameter Results to CSV
# --------------------------------------------------

results_df = pd.DataFrame({
    "Model": ["Zero-Order", "First-Order", "Higuchi", "Korsmeyer-Peppas"],
    "R2_Value": [res_zero.rvalue**2, res_first.rvalue**2, res_higuchi.rvalue**2, res_kp.rvalue**2],
    "Rate_Constant_k": [res_zero.slope, -res_first.slope, res_higuchi.slope, np.exp(res_kp.intercept)],
    "Exponent_n": [np.nan, np.nan, np.nan, res_kp.slope]
})

# Export to CSV
results_df.to_csv("kinetic_parameters_summary.csv", index=False) # index=False prevents pandas from writing its automatic row numbers
print("\nExported summary table to 'kinetic_parameters_summary.csv'")