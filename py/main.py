import re
from sympy import Matrix, lcm 
import tkinter as tk
from tkinter import messagebox
import pandas as pd
import numpy as np

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
    multiplier = 1

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

# def calculate_molar_mass(compound):
#     elements_and_numbers = re.split('([A-Z][a-z]?)', compound)
#     molar_mass = 0
#     i = 0
#     while i < len(elements_and_numbers) - 1:
#         i += 1
#         if len(elements_and_numbers[i]) > 0:
#             element = elements_and_numbers[i]
#             if elements_and_numbers[i + 1].isdigit():
#                 count = int(elements_and_numbers[i + 1])
#                 molar_mass += atomic_masses[element] * count
#                 i += 1
#             else:
#                 molar_mass += atomic_masses[element]
#     return molar_mass

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

# Prueba calcular masa molar



# root = tk.Tk()
# root.title("Calculo del reactivo limitante")

# root.configure(background="#F2EFE8")
# root.geometry("500x600")

root.mainloop()