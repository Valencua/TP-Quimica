import re
from sympy import Matrix, lcm 
import tkinter as tk
from tkinter import messagebox

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
    result_window.title("Resultado")
    result_window.geometry("500x250")
    result_window.configure(background="#F2EFE8")

    result_label = tk.Label(result_window, text="Ecuación Balanceada:", bg="#F2EFE8", fg="#8F788B", font=("Luckiest Guy", 14, "bold"))
    result_label.pack(padx=10, pady=(20,5), anchor="center")

    result_equation = tk.Label(result_window, text=output, bg="#F2EFE8", fg="#8F788B", font=("Arial", 12, "bold"))
    result_equation.pack(padx=10, pady=(5,20), anchor="center")

    return_button = tk.Button(result_window, text="Volver", command=result_window.destroy, bg="#EED0EA", fg="#8F788B", bd=3, relief="groove", font=("Luckiest Guy", 12, "bold"), padx=10, pady=5)
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