import re
from sympy import Matrix, lcm 
import tkinter as tk
from tkinter import messagebox
import pandas as pd
import numpy as np

# Diccionario de masas molares de los elementos
atomic_masses = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.95,
    'Tc': 98, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41,
    'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29,
    'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
    'Pm': 145, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.50,
    'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.04, 'Lu': 174.97, 'Hf': 178.49,
    'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08,
    'Au': 196.97, 'Hg': 200.59, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Th': 232.04,
    'Pa': 231.04, 'U': 238.03, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 'Bk': 247,
    'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262, 'Rf': 267,
    'Db': 270, 'Sg': 271, 'Bh': 270, 'Hs': 277, 'Mt': 276, 'Ds': 281, 'Rg': 282,
    'Cn': 285, 'Nh': 286, 'Fl': 289, 'Mc': 290, 'Lv': 293, 'Ts': 294, 'Og': 294
}

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

def calculate_molar_mass(compound):
    elements_and_numbers = re.split('([A-Z][a-z]?)', compound)
    molar_mass = 0
    i = 0
    while i < len(elements_and_numbers) - 1:
        i += 1
        if len(elements_and_numbers[i]) > 0:
            element = elements_and_numbers[i]
            if elements_and_numbers[i + 1].isdigit():
                count = int(elements_and_numbers[i + 1])
                molar_mass += atomic_masses[element] * count
                i += 1
            else:
                molar_mass += atomic_masses[element]
    return molar_mass

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
    for i in range(len(reactants)):
        molar_mass = calculate_molar_mass(reactants[i])
        stoichiometry_text += f"{reactants[i]}: {molar_mass} g/mol\n"
    for i in range(len(products)):
        molar_mass = calculate_molar_mass(products[i])
        stoichiometry_text += f"{products[i]}: {molar_mass} g/mol\n"

    stoichiometry_label = tk.Label(masses_window, text=stoichiometry_text, bg="#F2EFE8", fg="#8F788B", font=("Arial", 12, "bold"))
    stoichiometry_label.pack(padx=10, pady=(5,20), anchor="center")

    return_button = tk.Button(masses_window, text="Volver", command=masses_window.destroy, bg="#EED0EA", fg="#8F788B", bd=3, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=5)
    return_button.pack(padx=10, pady=(10,20), anchor="center")

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


root.mainloop()