# Constantes
molar_mass_NaCl = 58.44  # g/mol

# Inputs
concentration_HCl = float(input("Ingrese la concentración de HCl (M): "))
volume_HCl = float(input("Ingrese el volumen de HCl (mL): "))
concentration_NaOH = float(input("Ingrese la concentración de NaOH (M): "))

# Reacción balanceada
# HCl + NaOH -> NaCl + H2O

# Calcular el volumen de NaOH necesario
# Moles de HCl = concentración * volumen
moles_HCl = concentration_HCl * (volume_HCl / 1000)  # Convertir mL a L

# Dado que la reacción es 1:1, moles de NaOH necesarias = moles de HCl
moles_NaOH_needed = moles_HCl

# Volumen de NaOH necesario = moles / concentración
volume_NaOH_needed = moles_NaOH_needed / concentration_NaOH  # en L
volume_NaOH_needed_mL = volume_NaOH_needed * 1000  # Convertir a mL

# Calcular la masa de NaCl formada
# Moles de NaCl formadas = moles de HCl reaccionadas (relación 1:1)
moles_NaCl_formed = moles_HCl

# Masa de NaCl = moles * masa molar
mass_NaCl = moles_NaCl_formed * molar_mass_NaCl

# Experimento con 40 mL de cada solución
volume_experiment = float(input("Ingrese experimento de NaOH (mL): "))
volume_experiment2 = float(input("Ingrese experimento de NaOH (mL): "))

volume_HCl_experiment = volume_experiment  # mL
volume_NaOH_experiment = volume_experiment2  # mL

moles_HCl_experiment = concentration_HCl * (volume_HCl_experiment / 1000)
moles_NaOH_experiment = concentration_NaOH * (volume_NaOH_experiment / 1000)

# Determinar el reactivo limitante
if moles_HCl_experiment > moles_NaOH_experiment:
    limiting_reagent = 'NaOH'
    moles_remaining = moles_HCl_experiment - moles_NaOH_experiment
else:
    limiting_reagent = 'HCl'
    moles_remaining = moles_NaOH_experiment - moles_HCl_experiment

# Resultados con redondeo a 3 decimales
results = {
    "volume_NaOH_needed_mL": round(volume_NaOH_needed_mL, 3),
    "mass_NaCl_formed_g": round(mass_NaCl, 3),
    "limiting_reagent": limiting_reagent,
    "moles_remaining": round(moles_remaining, 3)
}

print(results)
