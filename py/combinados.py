import tkinter as tk
from tkinter import messagebox
from chemlib import Compound, Reaction

# Inicialización de diccionarios
reactivos = {}
productos = {}

# Constante de gas ideal
r = 0.082

# Funciones de cálculo (las funciones definidas en el prompt)

def get_by_MM():
    compuestos = [x for x in reactivos.keys() if reactivos[x]["MM"] is not None and reactivos[x]['V'] is not None] + [x for x in productos.keys() if productos[x]["MM"] is not None and productos[x]['V'] is not None]
    if len(compuestos) == 0: return
    
    for compuesto in compuestos:
        if compuesto in reactivos: reactivos[compuesto]["mol"] = reactivos[compuesto]["V"] * reactivos[compuesto]["MM"]
        else: productos[compuesto]["mol"] = productos[compuesto]["V"] * productos[compuesto]["MM"]
    return

def get_by_MMol():
    compuestos = [x for x in reactivos.keys() if reactivos[x]["mol"] is not None and reactivos[x]['MM'] is not None] + [x for x in productos.keys() if productos[x]["mol"] is not None and productos[x]['MM'] is not None]

    if len(compuestos) == 0: return
    
    for compuesto in compuestos:
        if compuesto in reactivos: reactivos[compuesto]["V"] = reactivos[compuesto]["mol"] / reactivos[compuesto]["MM"]
        else: productos[compuesto]["V"] = productos[compuesto]["mol"] / productos[compuesto]["MM"]
    return

def get_PVT():
    compuestos = [*reactivos, *productos]

    for compuesto in compuestos:
        P = reactivos[compuesto]["P"] if compuesto in reactivos.keys() else productos[compuesto]["P"] if compuesto in productos.keys() else None
        V = reactivos[compuesto]["V"] if compuesto in reactivos.keys() else productos[compuesto]["V"] if compuesto in productos.keys() else None
        T = reactivos[compuesto]["T"] if compuesto in reactivos.keys() else productos[compuesto]["T"] if compuesto in productos.keys() else None
        n = reactivos[compuesto]["mol"] if compuesto in reactivos.keys() else productos[compuesto]["mol"] if compuesto in productos.keys() else None

        if [P,V,T].count(None) > 1: return False        

        if P is None: 
            P = n * T * r / V
            if compuesto in reactivos: reactivos[compuesto]["P"] = P
            else: productos[compuesto]["P"] = P
        
        if V is None: 
            V = n * T * r / P
            if compuesto in reactivos: reactivos[compuesto]["V"] = V
            else: productos[compuesto]["V"] = V
    
        if T is None: 
            T = n * r / P * V
            if compuesto in reactivos: reactivos[compuesto]["T"] = T
            else: productos[compuesto]["T"] = T

    return True

def get_masses(reaction):
    compuestos = [x for x in reactivos.keys() if reactivos[x]["gr"] is not None or reactivos[x]['mol'] is not None] + [x for x in productos.keys() if productos[x]["mol"] is not None or productos[x]['mol'] is not None]
    if len(compuestos) == 0:
        return False
    
    compuesto = compuestos[0]

    listaCompuestos = [*reactivos, *productos]

    if compuesto in reactivos: 
        if reactivos[compuesto]["mol"] is not None:
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, moles = reactivos[compuesto]["mol"])
        else: 
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, grams = reactivos[compuesto]["gr"])
    else:
        if productos[compuesto]["mol"] is not None:
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, moles = productos[compuesto]["mol"])
        else: 
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, grams = productos[compuesto]["gr"])

    for x in listaCompuestos:
        if x in reactivos: 
            reactivos[x]["mol"] = listasa[listaCompuestos.index(x)]['moles']
            reactivos[x]["gr"] = listasa[listaCompuestos.index(x)]['grams']
        elif x in productos: 
            productos[x]["mol"] = listasa[listaCompuestos.index(x)]['moles']
            productos[x]["gr"] = listasa[listaCompuestos.index(x)]['grams']

    return get_PVT()

def get_mass_by_PVT(reaction):
    compuestos = [x for x in reactivos.keys() if reactivos[x]["P"] is not None and reactivos[x]["V"] is not None and reactivos[x]["T"] is not None] + [x for x in productos.keys() if productos[x]["P"] is not None and productos[x]["V"] is not None and productos[x]["T"] is not None]
    if len(compuestos) == 0:
        return False
    
    compuesto = compuestos[0]

    P = reactivos[compuesto]["P"] if compuesto in reactivos else productos[compuesto]["P"]
    V = reactivos[compuesto]["V"] if compuesto in reactivos else productos[compuesto]["V"]
    T = reactivos[compuesto]["T"] if compuesto in reactivos else productos[compuesto]["T"]
    
    if P is None or V is None or T is None: 
        return False
    
    mol = P * V / (r * T)

    if compuesto in reactivos: 
        reactivos[compuesto]["mol"] = mol
        get_masses(reaction)
    else: 
        productos[compuesto]["mol"] = mol
        get_masses(reaction)

    return True

def lim_reag(reaction):
    if len(reactivos) > 1:
        a = [x for x in reactivos if reactivos[x]["mol"] is not None]
        b = [reactivos[x]["mol"] for x in a]
        if len(a) == len(reactivos):
            l = reaction.limiting_reagent(*b, mode = 'moles')
            return l
        return "No se puede calcular el reactivo limitante"
    return "No hay reactivo limitante en una reaccion de un único reactivo"

def agregar_reactivo():
    c = entry_reactivo.get()
    mol = entry_mol_reac.get()
    gr = entry_gr_reac.get()
    P = entry_P_reac.get()
    V = entry_V_reac.get()
    T = entry_T_reac.get()
    MM = entry_MM_reac.get()
    
    if c:
        reactivos[Compound(c).formula] = {
            "mol": None if mol == '' else float(mol),
            "gr": None if gr == '' else float(gr),
            "P": None if P == '' else float(P),
            "V": None if V == '' else float(V),
            "T": None if T == '' else float(T),
            "MM": None if MM == '' else float(MM),
        }
        limpiar_entries_reac()
        actualizar_reaccion()

def agregar_producto():
    c = entry_producto.get()
    mol = entry_mol_prod.get()
    gr = entry_gr_prod.get()
    P = entry_P_prod.get()
    V = entry_V_prod.get()
    T = entry_T_prod.get()
    MM = entry_MM_prod.get()
    
    if c:
        productos[Compound(c).formula] = {
            "mol": None if mol == '' else float(mol),
            "gr": None if gr == '' else float(gr),
            "P": None if P == '' else float(P),
            "V": None if V == '' else float(V),
            "T": None if T == '' else float(T),
            "MM": None if MM == '' else float(MM),
        }
        limpiar_entries_prod()
        actualizar_reaccion()

def balancear():
    react = " + ".join([*reactivos])
    prod = " + ".join([*productos])
    
    try:
        reaccion = Reaction.by_formula(react + "-->" + prod)
        reaccion.balance()
    except Exception as e:
        messagebox.showerror("Error", "Reacción química no válida\n" + str(e))
        return
    
    get_mass_by_PVT(reaccion)
    get_by_MM()
    limitante = lim_reag(reaccion)
    get_masses(reaccion)
    get_by_MMol()
    get_PVT()
    
    

    P = productos
    R = reactivos

    messagebox.showinfo("Resultado", f"Reactivo limitante: {limitante}\nReactivos: {R}\nProductos: {P}")

def limpiar_entries_reac():
    entry_reactivo.delete(0, tk.END)
    entry_mol_reac.delete(0, tk.END)
    entry_gr_reac.delete(0, tk.END)
    entry_P_reac.delete(0, tk.END)
    entry_V_reac.delete(0, tk.END)
    entry_T_reac.delete(0, tk.END)
    entry_MM_reac.delete(0, tk.END)

def limpiar_entries_prod():
    entry_producto.delete(0, tk.END)
    entry_mol_prod.delete(0, tk.END)
    entry_gr_prod.delete(0, tk.END)
    entry_P_prod.delete(0, tk.END)
    entry_V_prod.delete(0, tk.END)
    entry_T_prod.delete(0, tk.END)
    entry_MM_prod.delete(0, tk.END)

def actualizar_reaccion():
    react = " + ".join([*reactivos])
    prod = " + ".join([*productos])
    reaccion_var.set(f"{react} ---> {prod}")

# Configuración de la interfaz gráfica
root = tk.Tk()
root.title("Balanceador de Reacciones Químicas")

reaccion_var = tk.StringVar()
reaccion_var.set("Reacción: ")

# Sección de Reactivos
frame_reactivos = tk.LabelFrame(root, text="Reactivos")
frame_reactivos.grid(row=0, column=0, padx=10, pady=10)

tk.Label(frame_reactivos, text="Reactivo").grid(row=0, column=0)
entry_reactivo = tk.Entry(frame_reactivos)
entry_reactivo.grid(row=0, column=1)

tk.Label(frame_reactivos, text="Moles").grid(row=1, column=0)
entry_mol_reac = tk.Entry(frame_reactivos)
entry_mol_reac.grid(row=1, column=1)

tk.Label(frame_reactivos, text="Gramos").grid(row=2, column=0)
entry_gr_reac = tk.Entry(frame_reactivos)
entry_gr_reac.grid(row=2, column=1)

tk.Label(frame_reactivos, text="Presión").grid(row=3, column=0)
entry_P_reac = tk.Entry(frame_reactivos)
entry_P_reac.grid(row=3, column=1)

tk.Label(frame_reactivos, text="Volumen").grid(row=4, column=0)
entry_V_reac = tk.Entry(frame_reactivos)
entry_V_reac.grid(row=4, column=1)

tk.Label(frame_reactivos, text="Temperatura").grid(row=5, column=0)
entry_T_reac = tk.Entry(frame_reactivos)
entry_T_reac.grid(row=5, column=1)

tk.Label(frame_reactivos, text="Concentración Molar").grid(row=6, column=0)
entry_MM_reac = tk.Entry(frame_reactivos)
entry_MM_reac.grid(row=6, column=1)

btn_agregar_reac = tk.Button(frame_reactivos, text="Agregar Reactivo", command=agregar_reactivo)
btn_agregar_reac.grid(row=7, columnspan=2)

# Sección de Productos
frame_productos = tk.LabelFrame(root, text="Productos")
frame_productos.grid(row=0, column=1, padx=10, pady=10)

tk.Label(frame_productos, text="Producto").grid(row=0, column=0)
entry_producto = tk.Entry(frame_productos)
entry_producto.grid(row=0, column=1)

tk.Label(frame_productos, text="Moles").grid(row=1, column=0)
entry_mol_prod = tk.Entry(frame_productos)
entry_mol_prod.grid(row=1, column=1)

tk.Label(frame_productos, text="Gramos").grid(row=2, column=0)
entry_gr_prod = tk.Entry(frame_productos)
entry_gr_prod.grid(row=2, column=1)

tk.Label(frame_productos, text="Presión").grid(row=3, column=0)
entry_P_prod = tk.Entry(frame_productos)
entry_P_prod.grid(row=3, column=1)

tk.Label(frame_productos, text="Volumen").grid(row=4, column=0)
entry_V_prod = tk.Entry(frame_productos)
entry_V_prod.grid(row=4, column=1)

tk.Label(frame_productos, text="Temperatura").grid(row=5, column=0)
entry_T_prod = tk.Entry(frame_productos)
entry_T_prod.grid(row=5, column=1)

tk.Label(frame_productos, text="Concentración Molar").grid(row=6, column=0)
entry_MM_prod = tk.Entry(frame_productos)
entry_MM_prod.grid(row=6, column=1)

btn_agregar_prod = tk.Button(frame_productos, text="Agregar Producto", command=agregar_producto)
btn_agregar_prod.grid(row=7, columnspan=2)

# Botón para balancear
btn_balancear = tk.Button(root, text="Calcular y Balancear", command=balancear)
btn_balancear.grid(row=1, columnspan=2, pady=10)

# Mostrar reacción en tiempo real
lbl_reaccion = tk.Label(root, textvariable=reaccion_var)
lbl_reaccion.grid(row=2, columnspan=2, pady=10)

root.mainloop()
