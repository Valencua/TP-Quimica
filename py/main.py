import re
from sympy import Matrix, lcm 
import sympy as sp
import tkinter as tk
from tkinter import messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Variables para calculos de ecuaciones
x = sp.symbols('x')
R = 0.0821 
P, V, T = sp.symbols('P V T')

# Diccionario de masas molares de los elementos

MASAS_MOLARES = {
    "H": 1.01,
    "He": 4.00,
    "Li": 6.94,
    "Be": 9.01,
    "B": 10.81,
    "C": 12.01,
    "N": 14.01,
    "O": 16.00,
    "F": 19.00,
    "Ne": 20.18,
    "Na": 22.99,
    "Mg": 24.31,
    "Al": 26.98,
    "Si": 28.09,
    "P": 30.97,
    "S": 32.07,
    "Cl": 35.45,
    "K": 39.10,
    "Ar": 39.95,
    "Ca": 40.08,
    "Sc": 44.96,
    "Ti": 47.87,
    "V": 50.94,
    "Cr": 52.00,
    "Mn": 54.94,
    "Fe": 55.85,
    "Co": 58.93,
    "Ni": 58.69,
    "Cu": 63.55,
    "Zn": 65.38,
    "Ga": 69.72,
    "Ge": 72.63,
    "As": 74.92,
    "Se": 78.96,
    "Br": 79.90,
    "Kr": 83.80,
    "Rb": 85.47,
    "Sr": 87.62,
    "Y": 88.91,
    "Zr": 91.22,
    "Nb": 92.91,
    "Mo": 95.95,
    "Tc": 98.00,
    "Ru": 101.07,
    "Rh": 102.91,
    "Pd": 106.42,
    "Ag": 107.87,
    "Cd": 112.41,
    "In": 114.82,
    "Sn": 118.71,
    "Sb": 121.76,
    "Te": 127.60,
    "I": 126.90,
    "Xe": 131.29,
    "Cs": 132.91,
    "Ba": 137.33,
    "La": 138.91,
    "Ce": 140.12,
    "Pr": 140.91,
    "Nd": 144.24,
    "Pm": 145.00,
    "Sm": 150.36,
    "Eu": 151.96,
    "Gd": 157.25,
    "Tb": 158.93,
    "Dy": 162.50,
    "Ho": 164.93,
    "Er": 167.26,
    "Tm": 168.93,
    "Yb": 173.04,
    "Lu": 174.97,
    "Hf": 178.49,
    "Ta": 180.95,
    "W": 183.84,
    "Re": 186.21,
    "Os": 190.23,
    "Ir": 192.22,
    "Pt": 195.08,
    "Au": 196.97,
    "Hg": 200.59,
    "Tl": 204.38,
    "Pb": 207.2,
    "Bi": 208.98,
    "Th": 232.04,
    "Pa": 231.04,
    "U": 238.03
}

RESULT_EQ_BALANCEADA = ""

def calcular_masa_molar(compound):

    pattern = r'([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)'
    stack = []
    total_molar_mass = 0

    for match in re.finditer(pattern, compound):
        element, count, open_paren, close_paren, close_count = match.groups()

        if element:
            count = int(count) if count else 1
            stack.append(MASAS_MOLARES[element] * count)

        elif open_paren:
            stack.append('(')

        elif close_paren:
            temp_mass = 0
            while stack[-1] != '(':
                temp_mass += stack.pop()
            stack.pop()  # remove '(' from stack
            count = int(close_count) if close_count else 1
            stack.append(temp_mass * count)

    total_molar_mass = sum(stack)
    return total_molar_mass


def addToMatrix(element, index, count, side, elementMatrix, elementList):
    if(index == len(elementMatrix)):
       elementMatrix.append([])
       for x in elementList:
            elementMatrix[index].append(0)
    if(element not in elementList):
        elementList.append(element)
        for i in range(len(elementMatrix)):
            elementMatrix[i].append(0)
    column=elementList.index(element)
    elementMatrix[index][column]+=count*side
    
def findElements(segment,index, multiplier, side, elementMatrix, elementList):
    elementsAndNumbers=re.split('([A-Z][a-z]?)',segment)
    i=0
    while(i<len(elementsAndNumbers)-1):
          i+=1
          if(len(elementsAndNumbers[i])>0):
            if(elementsAndNumbers[i+1].isdigit()):
                count=int(elementsAndNumbers[i+1])*multiplier
                addToMatrix(elementsAndNumbers[i], index, count, side, elementMatrix, elementList)
                i+=1
            else:
                addToMatrix(elementsAndNumbers[i], index, multiplier, side, elementMatrix, elementList)        
    
def compoundDecipher(compound, index, side, elementMatrix, elementList):
    segments=re.split('(\([A-Za-z0-9]*\)[0-9]*)',compound)    
    for segment in segments:
        if segment.startswith("("):
            segment=re.split('\)([0-9]*)',segment)
            multiplier=int(segment[1])
            segment=segment[0][1:]
        else:
            multiplier=1
        findElements(segment, index, multiplier, side, elementMatrix, elementList)

def balance_equation():
    elementMatrix = []  
    elementList = []    
    reactants = reactants_entry.get()
    products = products_entry.get()
    
    reactants=reactants.replace(' ', '').split("+")
    products=products.replace(' ', '').split("+")
    
    for i in range(len(reactants)):
        compoundDecipher(reactants[i],i,1, elementMatrix, elementList)
    for i in range(len(products)):
        compoundDecipher(products[i],i+len(reactants),-1, elementMatrix, elementList)
        
    elementMatrix = Matrix(elementMatrix)
    elementMatrix = elementMatrix.transpose()
    solution=elementMatrix.nullspace()[0]
    multiple = lcm([val.q for val in solution])
    solution = multiple*solution
    coEffi=solution.tolist()
    output=""
    for i in range(len(reactants)):
        output+=str(coEffi[i][0])+reactants[i]
        if i<len(reactants)-1:
            output+=" + "
    output+=" -> "
    for i in range(len(products)):
        output+=str(coEffi[i+len(reactants)][0])+products[i]
        if i<len(products)-1:
            output+=" + "
            
    result_window = tk.Toplevel(root)
    result_window.title("Equivalencia")
    result_window.geometry("500x250")
    result_window.configure(background="#F2EFE8")

    result_label = tk.Label(result_window, text="Ecuación Balanceada:", bg="#F2EFE8", fg="#8F788B", font=("Luckiest Guy", 14, "bold"))
    result_label.pack(padx=10, pady=(20,5), anchor="center")
    
    
    result_equation = tk.Label(result_window, text=output, bg="#F2EFE8", fg="#8F788B", font=("Arial", 12, "bold"))
    result_equation.pack(padx=10, pady=(5,20), anchor="center")
    # es donde guarda la ecuación balanceada
    RESULT_EQ_BALANCEADA = output

    return_button = tk.Button(result_window, text="Volver", command=result_window.destroy, bg="#EED0EA", fg="#8F788B", bd=3, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=5)
    return_button.pack(padx=10, pady=(10,20), anchor="center")
    
    # Display stoichiometric calculations
    masses_window = tk.Toplevel(root)
    masses_window.title("Cálculos Estequiométricos")
    masses_window.geometry("500x250")
    masses_window.configure(background="#F2EFE8")

    masses_label = tk.Label(masses_window, text="Masas Molares y Cálculos:", bg="#F2EFE8", fg="#8F788B", font=("Luckiest Guy", 14, "bold"))
    masses_label.pack(padx=10, pady=(20,5), anchor="center")

    stoichiometry_text = "Masas Molares:\n"
    masses = {}
    for i in range(len(reactants)):
        #molar_mass = calculate_molar_mass(reactants[i])
        molar_mass = calcular_masa_molar(reactants[i])
        masses[reactants[i]] = molar_mass
        stoichiometry_text += f"{reactants[i]}: {molar_mass} g/mol\n"
    for i in range(len(products)):
        #molar_mass = calculate_molar_mass(products[i])
        molar_mass = calcular_masa_molar(products[i])
        stoichiometry_text += f"{products[i]}: {molar_mass} g/mol\n"

    stoichiometry_label = tk.Label(masses_window, text=stoichiometry_text, bg="#F2EFE8", fg="#8F788B", font=("Arial", 12, "bold"))
    stoichiometry_label.pack(padx=10, pady=(5,20), anchor="center")

    return_button = tk.Button(masses_window, text="Volver", command=masses_window.destroy, bg="#EED0EA", fg="#8F788B", bd=3, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=5)
    return_button.pack(padx=10, pady=(10,20), anchor="center")

    reactivolimitante_window = tk.Toplevel(root)
    reactivolimitante_window.title("Cálculo de reactivo limitante")
    reactivolimitante_window.geometry("900x250")
    reactivolimitante_window.configure(background="#F2EFE8")
    result_label_eq = tk.Label(reactivolimitante_window, text=RESULT_EQ_BALANCEADA, bg="#F2EFE8", fg="#8F788B", font=("Luckiest Guy", 14, "bold"))
    result_label_eq.pack(padx=10, pady=(20,5), anchor="center")
    equation = RESULT_EQ_BALANCEADA
    
    limiting_reagent, max_moles_product, excess_reagents = find_limiting_reagent(equation, masses)
    label_result_react_limitante = f"El límite de reactivo es {limiting_reagent} y este puede producir un maximo de {max_moles_product} moles de producto. \n Y al reactivo en exceso {list(excess_reagents.keys())[0]} le sobran {list(excess_reagents.values())[0]} moles"

    result_label_calculo_rl = tk.Label(reactivolimitante_window, text=label_result_react_limitante, bg="#F2EFE8", fg="#8F788B", font=("Luckiest Guy", 14, "bold"))
    result_label_calculo_rl.pack(padx=10, pady=(20,5), anchor="center")

def calculate_moles(mass, molar_mass):
    return mass / molar_mass

def parse_equation(equation):
    reactants_part, products_part = equation.split('->')
    reactants = reactants_part.split('+')
    products = products_part.split('+')
    
    reactants = [r.strip() for r in reactants]
    products = [p.strip() for p in products]
    
    return reactants, products

def extract_compound_data(compound):
    pattern = r"(\d*)([A-Za-z0-9\(\)]+)"
    matches = re.match(pattern, compound.strip())
    
    coefficient = int(matches.group(1)) if matches.group(1) else 1
    formula = matches.group(2)
    
    return coefficient, formula

def find_limiting_reagent(equation, masses):
    reactants, products = parse_equation(equation)
    
    reactant_data = [extract_compound_data(reactant) for reactant in reactants]
    product_data = [extract_compound_data(product) for product in products]
    
    # Calculate molar masses for each reactant and product
    molar_masses = {formula: calcular_masa_molar(formula) for _, formula in reactant_data + product_data}
    
    # Check if all reactants are in the masses dictionary
    for _, formula in reactant_data:
        if formula not in masses:
            raise KeyError(f"Mass for reactant {formula} not provided.")
    
    # Calculate moles for each reactant using the provided masses
    moles_reagents = {formula: calculate_moles(masses[formula], molar_masses[formula]) for _, formula in reactant_data}
    
    # Determine the stoichiometric ratios
    reagent_ratios = {formula: coef for coef, formula in reactant_data}
    
    # Find the limiting reagent
    limiting_reagent = None
    max_moles_product = float('inf')
    
    for formula, moles in moles_reagents.items():
        max_moles_for_reagent = moles / reagent_ratios[formula]
        if max_moles_for_reagent < max_moles_product:
            max_moles_product = max_moles_for_reagent
            limiting_reagent = formula
    
    # Calculate the amount of excess reagent left
    excess_reagents = {}
    for formula, moles in moles_reagents.items():
        if formula != limiting_reagent:
            moles_used = reagent_ratios[formula] * max_moles_product
            excess_reagents[formula] = moles - moles_used
    
    return limiting_reagent, max_moles_product, excess_reagents

# Funciones para el calculo de las ecuaciones de gases ideales
def leer_valor(mensaje):
    while True:
        try:
            valor = input(mensaje)
            if valor == "":
                return None
            else:
                return float(valor)
        except ValueError:
            print("Por favor, ingrese un valor numérico o deje en blanco para omitir.")
def gases_compleja(nombre, simbolo_quimico,presion,volumen,moles,temperatura):
    reactivos = []
    R = 0.0821  
    P, V, n, T = sp.symbols('P V n T')
    ecuacion = sp.Eq(P * V, n * R * T)
   
    #print(f"Ingrese los valores del reactivo")
    #nombre = input("Ingrese el nombre del reactivo: ")
    #simbolo_quimico = input("Ingrese el símbolo químico del elemento: ")
    #presion = leer_valor("Ingrese el valor de la presión (en atm): ")
    #volumen = leer_valor("Ingrese el valor del volumen (en L): ")
    #moles = leer_valor("Ingrese la cantidad de moles: ")
    #temperatura = leer_valor("Ingrese el valor de la temperatura (en K): ")
    reactivos.append((nombre, simbolo_quimico, presion, volumen, moles, temperatura))


    if presion is None:
        solucion = sp.solve(ecuacion.subs({V: volumen, n: moles, T: temperatura}), P)
    elif volumen is None:
        solucion = sp.solve(ecuacion.subs({P: presion, n: moles, T: temperatura}), V)
    elif moles is None:
        solucion = sp.solve(ecuacion.subs({P: presion, V: volumen, T: temperatura}), n)
    elif temperatura is None:
        solucion = sp.solve(ecuacion.subs({P: presion, V: volumen, n: moles}), T)
    else:
        #print("Todos los valores son conocidos, no hay nada que resolver.")
        return "Todos los valores son conocidos, no hay nada que resolver."
    #print(solucion)
    return solucion
def gases_simple(presion1, volumen1, moles1, temperatura1, presion2, volumen2, moles2, temperatura2):
    reactivos = {}

    # print("Ingrese los valores del primer estado del reactivo")
    # presion1 = leer_valor("Ingrese el valor de la presión (en atm): ")
    # volumen1 = leer_valor("Ingrese el valor del volumen (en L): ")
    # temperatura1 = leer_valor("Ingrese el valor de la temperatura (en K): ")

    # print("Ingrese los valores del segundo estado del reactivo")
    # presion2 = leer_valor("Ingrese el valor de la presión (en atm): ")
    # volumen2 = leer_valor("Ingrese el valor del volumen (en L): ")
    # temperatura2 = leer_valor("Ingrese el valor de la temperatura (en K): ")

    reactivos['estado1'] = {'presion1': presion1, 'volumen1': volumen1, 'moles1': moles1, 'temperatura1': temperatura1}
    reactivos['estado2'] = {'presion2': presion2, 'volumen2': volumen2, 'moles2': moles2, 'temperatura2': temperatura2}

    if presion1 is None:
        ecuacion = sp.Eq(x * volumen1 / moles1 * temperatura1, presion2 * volumen2 / moles2 * temperatura2)
    elif volumen1 is None:
        ecuacion = sp.Eq(presion1 * x / moles1 * temperatura1, presion2 * volumen2 / moles2 * temperatura2)
    elif moles1 is None:
        ecuacion = sp.Eq(presion1 * volumen1 / x * temperatura1, presion2 * volumen2 / moles2 * temperatura2)
    elif temperatura1 is None:
        ecuacion = sp.Eq(presion1 * volumen1 / moles1 * x, presion2 * volumen2 / moles2 * temperatura2)
    elif presion2 is None:
        ecuacion = sp.Eq(presion1 * volumen1 / moles1 * temperatura1, x * volumen2 / moles2 * temperatura2)
    elif volumen2 is None:
        ecuacion = sp.Eq(presion1 * volumen1 / moles1 * temperatura1, presion2 * x / moles2 * temperatura2)
    elif moles2 is None:
        ecuacion = sp.Eq(presion1 * volumen1 / moles1 * temperatura1, presion2 * volumen2 / x * temperatura2)
    elif temperatura2 is None:
        ecuacion = sp.Eq(presion1 * volumen1 / moles1 * temperatura1, presion2 * volumen2 / moles2 * x)
    else:
        #print("Todos los valores son conocidos, no hay nada que resolver.")
        ecuacion = None
        return "Todos los valores son conocidos, no hay nada que resolver."
        
    solucion = sp.solve(ecuacion, x)

    return solucion

def call_gases_compleja():
    nombre = gases_input_reactivo.get()
    print(nombre)
    simbolo_quimico = gases_input_reactivo_s.get()
    print(simbolo_quimico)
    presion = gases_input_presion.get()
    if presion == "": 
        presion = None
    volumen = gases_input_volumen.get()
    if volumen == "":
        volumen = None
    moles = gases_input_moles.get()
    if moles == "":
        moles = None
    temperatura = gases_input_temperatura.get()
    if temperatura == "":
        temperatura = None

    resultado_gases_compleja = gases_compleja(nombre,simbolo_quimico,presion,volumen, moles, temperatura)
    gases_complejos_label_resultado.config(text=resultado_gases_compleja)

def call_gases_simple():
    temperatura1 = gases_simples_input_temperatura_inicial.get()
    if temperatura1 == "":
        temperatura1 = None
    else : 
        temperatura1 = float(temperatura1)

    temperatura2 = gases_simples_input_temperatura_final.get()
    if temperatura2 == "":
        temperatura2 = None
    else :
        temperatura2 = float(temperatura2)

    presion1 = gases_simples_input_presion_inicial.get()
    if presion1 == "":
        presion1 = None
    else:
        presion1 = float(presion1)
    
    presion2 = gases_simples_input_presion_final.get()
    if presion2 == "":
        presion2 = None
    else:
        presion2 = float(presion2)

    volumen1 = gases_simples_input_volumen_inicial.get()
    if volumen1 == "":
        volumen1 = None
    else:
        volumen1 = float(volumen1)

    volumen2 = gases_simples_input_volumen_final.get()
    if volumen2 == "":
        volumen2 = None
    else:
        volumen2 = float(volumen2)
    moles1 = gases_simples_input_moles_inicial.get()
    if moles1 == "":
        moles1 = None
    else : 
        moles1 = float(moles1)
    moles2 = gases_simples_input_moles_final.get()
    if moles2 == "":
        moles2 = None
    else : 
        moles2 = float(moles2)

    resultado_gases_simple = gases_simple(presion1, volumen1, moles1, temperatura1, presion2, volumen2, moles2, temperatura2)
    gases_simples_label_resultado.config(text=resultado_gases_simple)

root = tk.Tk()
root.title("Balanceador de ecuaciones químicas")
root.configure(background="#F2EFE8")
root.geometry("500x600")
reactants_label = tk.Label(root, text="Reactivos:", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
reactants_entry = tk.Entry(root, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
products_label = tk.Label(root, text="Productos:", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
products_entry = tk.Entry(root, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
balance_button = tk.Button(root, text="Balancear", command=balance_equation, bg="#EED0EA", fg="#8F788B", bd=5, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=10)
balance_button.config(highlightbackground="#F2EFE8")  
reactants_label.place(relx=0.5, rely=0.3, anchor="center", y=-30)
reactants_entry.place(relx=0.5, rely=0.4, anchor="center", y=0) 
products_label.place(relx=0.5, rely=0.5, anchor="center", y=30)
products_entry.place(relx=0.5, rely=0.6, anchor="center", y=60)
balance_button.place(relx=0.5, rely=0.8, anchor="center", y=70)

# Ventana secundaria para el calculador de gases ideales complejos
ventana_gases_complejos = tk.Tk()
ventana_gases_complejos.title("Calulador de gases ideales complejos")
ventana_gases_complejos.configure(background="#F2EFE8")
ventana_gases_complejos.geometry("700x600")

# Configurar la rejilla para que los widgets se expandan con la ventana
ventana_gases_complejos.columnconfigure(0, weight=1)
ventana_gases_complejos.columnconfigure(1, weight=3)

# Crear y colocar los labels y entries en la rejilla
gases_label_reactivo = tk.Label(ventana_gases_complejos, text="Ingrese el nombre de la sustancia", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_label_reactivo.grid(row=0, column=0, sticky="e", padx=5, pady=5)
gases_input_reactivo = tk.Entry(ventana_gases_complejos, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_input_reactivo.grid(row=0, column=1, padx=5, pady=5)

gases_label_reactivo_s = tk.Label(ventana_gases_complejos, text="Ingrese el símbolo químico de la sustancia", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_label_reactivo_s.grid(row=1, column=0, sticky="e", padx=5, pady=5)
gases_input_reactivo_s = tk.Entry(ventana_gases_complejos, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_input_reactivo_s.grid(row=1, column=1, padx=5, pady=5)

gases_label_presion = tk.Label(ventana_gases_complejos, text="Ingrese el valor de la presión (en atm)", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_label_presion.grid(row=2, column=0, sticky="e", padx=5, pady=5)
gases_input_presion = tk.Entry(ventana_gases_complejos, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_input_presion.grid(row=2, column=1, padx=5, pady=5)

gases_label_volumen = tk.Label(ventana_gases_complejos, text="Ingrese el valor de volumen (en L)", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_label_volumen.grid(row=3, column=0, sticky="e", padx=5, pady=5)
gases_input_volumen = tk.Entry(ventana_gases_complejos, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_input_volumen.grid(row=3, column=1, padx=5, pady=5)

gases_label_moles = tk.Label(ventana_gases_complejos, text="Ingrese la cantidad de moles", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_label_moles.grid(row=4, column=0, sticky="e", padx=5, pady=5)
gases_input_moles = tk.Entry(ventana_gases_complejos, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_input_moles.grid(row=4, column=1, padx=5, pady=5)

gases_label_temperatura = tk.Label(ventana_gases_complejos, text="Ingrese el valor de la temperatura", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_label_temperatura.grid(row=5, column=0, sticky="e", padx=5, pady=5)
gases_input_temperatura = tk.Entry(ventana_gases_complejos, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_input_temperatura.grid(row=5, column=1, padx=5, pady=5)

gases_complejos_button = tk.Button(ventana_gases_complejos, text="Calcular", command=call_gases_compleja, bg="#EED0EA", fg="#8F788B", bd=5, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=10)
gases_complejos_button.grid(row=6, column=0, columnspan=3, pady=20)

gases_complejos_label_resultado = tk.Label(ventana_gases_complejos, text="RESULTADO", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_complejos_label_resultado.grid(row=7, column=0, columnspan=3, pady=20)

ventana_gases_simples = tk.Tk()
ventana_gases_simples.title("Calculador de gases ideales simples")
ventana_gases_simples.configure(background="#F2EFE8")
ventana_gases_simples.geometry("700x500")

# Configurar la rejilla para que los widgets se expandan con la ventana
ventana_gases_simples.columnconfigure(0, weight=1)
ventana_gases_simples.columnconfigure(1, weight=1)
ventana_gases_simples.columnconfigure(2, weight=1)

# Crear y colocar los labels y entries en la rejilla
gases_simples_label_reactivo_nombre = tk.Label(ventana_gases_simples, text="Si el valor es desconocido y constante, insertar 1 en inicial y final. Si es desconocido dejarlo vacío.", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_reactivo_nombre.grid(row=0, column=0, columnspan=3, padx=5, pady=5)

gases_simples_label_reactivo_simbolo = tk.Label(ventana_gases_simples, text="Ingrese el símbolo químico de la sustancia", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_reactivo_simbolo.grid(row=1, column=0, sticky="e", padx=5, pady=5)
gases_simples_input_reactivo_simbolo = tk.Entry(ventana_gases_simples, width=40, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_reactivo_simbolo.grid(row=1, column=1, columnspan=2, padx=5, pady=5)

# Agregar los títulos "Inicial" y "Final"
titulo_inicial = tk.Label(ventana_gases_simples, text="Inicial", font=("Luckiest Guy", 12, "bold"), bg="#F2EFE8")
titulo_inicial.grid(row=2, column=1, padx=5, pady=5)

titulo_final = tk.Label(ventana_gases_simples, text="Final", font=("Luckiest Guy", 12, "bold"), bg="#F2EFE8")
titulo_final.grid(row=2, column=2, padx=5, pady=5)

# Crear y colocar los labels y entries para presión, volumen y temperatura en la rejilla
gases_simples_label_presion = tk.Label(ventana_gases_simples, text="Presión (atm)", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_presion.grid(row=3, column=0, sticky="e", padx=5, pady=5)
gases_simples_input_presion_inicial = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_presion_inicial.grid(row=3, column=1, padx=5, pady=5)
gases_simples_input_presion_final = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_presion_final.grid(row=3, column=2, padx=5, pady=5)

gases_simples_label_volumen = tk.Label(ventana_gases_simples, text="Volumen (L)", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_volumen.grid(row=4, column=0, sticky="e", padx=5, pady=5)
gases_simples_input_volumen_inicial = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_volumen_inicial.grid(row=4, column=1, padx=5, pady=5)
gases_simples_input_volumen_final = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_volumen_final.grid(row=4, column=2, padx=5, pady=5)

gases_simples_label_temperatura = tk.Label(ventana_gases_simples, text="Temperatura (K)", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_temperatura.grid(row=5, column=0, sticky="e", padx=5, pady=5)
gases_simples_input_temperatura_inicial = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_temperatura_inicial.grid(row=5, column=1, padx=5, pady=5)
gases_simples_input_temperatura_final = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_temperatura_final.grid(row=5, column=2, padx=5, pady=5)

gases_simples_label_moles = tk.Label(ventana_gases_simples, text="Moles (mol)", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_moles.grid(row=6, column=0, sticky="e", padx=5, pady=5)
gases_simples_input_moles_inicial = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_moles_inicial.grid(row=6, column=1, padx=5, pady=5)
gases_simples_input_moles_final = tk.Entry(ventana_gases_simples, width=20, bd=2, relief="groove", font=("Luckiest Guy", 12))
gases_simples_input_moles_final.grid(row=6, column=2, padx=5, pady=5)

# Colocar el botón "Calcular" en una nueva fila
gases_simples_button = tk.Button(ventana_gases_simples, text="Calcular", command=call_gases_simple, bg="#EED0EA", fg="#8F788B", bd=5, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=10)
gases_simples_button.grid(row=7, column=0, columnspan=3, pady=20)

# Colocar el label de resultado en una nueva fila
gases_simples_label_resultado = tk.Label(ventana_gases_simples, text="RESULTADO", bg="#EEE1D0", bd=2, relief="groove", padx=10, pady=10, font=("Luckiest Guy", 12))
gases_simples_label_resultado.grid(row=8, column=0, columnspan=3, pady=20)


# Iniciar el loop principal de la aplicación
root.mainloop()
ventana_gases_complejos.mainloop()
ventana_gases_simples.mainloop()