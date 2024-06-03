# Constants
molar_mass_NaCl = 58.44  # g/mol

# Given data
concentration_HCl = 0.52  # M
volume_HCl = 80  # mL
concentration_NaOH = 0.86  # M
volume_NaOH = None  # To be calculated

# a) Balanced reaction
# HCl + NaOH -> NaCl + H2O

# b) Calculate the volume of NaOH needed
# Moles of HCl = concentration * volume
moles_HCl = concentration_HCl * (volume_HCl / 1000)  # Convert mL to L

# Since the reaction is 1:1, moles of NaOH needed = moles of HCl
moles_NaOH_needed = moles_HCl

# Volume of NaOH needed = moles / concentration
volume_NaOH_needed = moles_NaOH_needed / concentration_NaOH  # in L
volume_NaOH_needed_mL = volume_NaOH_needed * 1000  # Convert to mL

# c) Calculate the mass of NaCl formed
# Moles of NaCl formed = moles of HCl reacted (1:1 ratio)
moles_NaCl_formed = moles_HCl

# Mass of NaCl = moles * molar mass
mass_NaCl = moles_NaCl_formed * molar_mass_NaCl

# d) Reaction with 40 mL of each solution
volume_HCl_experiment = 40  # mL
volume_NaOH_experiment = 40  # mL

moles_HCl_experiment = concentration_HCl * (volume_HCl_experiment / 1000)
moles_NaOH_experiment = concentration_NaOH * (volume_NaOH_experiment / 1000)

# Determine the limiting reagent
if moles_HCl_experiment > moles_NaOH_experiment:
    limiting_reagent = 'NaOH'
    moles_remaining = moles_HCl_experiment - moles_NaOH_experiment
else:
    limiting_reagent = 'HCl'
    moles_remaining = moles_NaOH_experiment - moles_HCl_experiment

# Output results with rounding to 3 decimals
results = {
    "volume_NaOH_needed_mL": round(volume_NaOH_needed_mL, 3),
    "mass_NaCl_formed_g": round(mass_NaCl, 3),
    "limiting_reagent": limiting_reagent,
    "moles_remaining": round(moles_remaining, 3)
}

print(results)
