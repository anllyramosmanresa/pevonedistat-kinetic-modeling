"""
Controlled Drug-Release Kinetic Analysis
=========================================
 
Workflow:
    1. Load release data (raw time/absorbance measurements)
    2. Validate the data
    3. Convert absorbance to concentration using a calibration curve
    4. Calculate cumulative drug release, correcting for sampling losses
    5. Fit standard release kinetic models (Zero-order, First-order,
       Higuchi, Korsmeyer-Peppas)
    6. Compare model fits using R²
    7. Plot the cumulative release profile
 
Data status
-----------
The dataset embedded in `load_data()` below is REAL experimental data
(Pevonedistat release from cellulose particles), not simulated data.
Simulated/test data should be kept in a separate file/folder and clearly
labelled as such -- see the "Testing" section at the bottom of this file.
"""
