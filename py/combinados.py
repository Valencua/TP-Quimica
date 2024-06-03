from chemlib import Compound, Reaction
compound = dict["c":Compound, "mol":float, "gr": float, "P": float, "V": float, "T": float]
reactivos: dict[str: compound] = {}
productos: dict[str: compound] = {}

r = 0.082

def get_PVT():
    pass

def get_masses(reaction: Reaction, cantidad:float, type:str, compuesto:int = 0): #en base a una reacion, una cantidad (de moles o gr) y un compuesto puede sacar el resto de gr y moles de los otros compuestos.
    #esta funcion habria que cambiarla para que no acepte mas data que la reaccion. deberia funcionar checkeando si exsiten compuestos con n (mol) o gr que tengan valores definidos
    listaCompuestos= [*reactivos, *productos]

    if type == "moles": 
        listasa = reaction.get_amounts(compuesto, moles = cantidad)
    elif type == "grams":
        listasa = reaction.get_amounts(compuesto, grams = cantidad)
    else: return False



    for x in range(0, len(listaCompuestos)):

        if listaCompuestos[x] in reactivos.keys(): 
            reactivos[listaCompuestos[x]]["mol"] = listasa[x]['moles']
            reactivos[listaCompuestos[x]]["gr"] = listasa[x]['grams']
        else: 
            productos[listaCompuestos[x]]["mol"] = listasa[x]['moles']
            productos[listaCompuestos[x]]["gr"] = listasa[x]['grams']
    

    print("\n\n",reactivos, "\n\n", productos, "\n\n", listasa, "\n\n")
    return




def get_mass_by_PVT(reaction: Reaction): #saca la masa de un compuesto en base a la presion, volumen y temperatura... y r que es una constante
    compuestos = [x for x in reactivos.keys() if reactivos[x]["P"] is not None and reactivos[x]["V"] is not None and reactivos[x]["T"] is not None] + [x for x in productos.keys() if productos[x]["P"] is not None and productos[x]["V"] is not None and productos[x]["T"] is not None]
    if len(compuestos) == 0: 
        return
    
    compuesto = compuestos[0] if len(compuestos) > 0 else None

    P = reactivos[compuesto]["P"] if compuesto in reactivos.keys() else productos[compuesto]["P"] if compuesto in productos.keys() else None
    V = reactivos[compuesto]["V"] if compuesto in reactivos.keys() else productos[compuesto]["V"] if compuesto in productos.keys() else None
    T = reactivos[compuesto]["T"] if compuesto in reactivos.keys() else productos[compuesto]["T"] if compuesto in productos.keys() else None
    
    if P is None or V is None or T is None: 
        return False
    
    r = 0.082
    mol = P * V / (r * T)
    print(P, V, T, mol)

    if compuesto in reactivos.keys(): 
        reactivos[compuesto]["mol"] = mol
        get_masses(reaction, mol, "moles", compuestos.index(compuesto)+1)
    else: 
        productos[compuesto]["mol"] = mol
        get_masses(reaction, mol, "moles", compuestos.index(compuesto)+1)


    return True




def main(): #main

    match input("\n1.Ingresar reactivo \n2.Ingresar producto \n3.Calcular \n"): #Input de datos. podría estar mejor? si. Tengo sueño? tmb
         
        case "1":
            c   = input("Ingrese el reactivo: ")
            mol = input("Ingrese la cantidad de moles: ")
            gr  = input("Ingrese la cantidad de gramos: ")
            P   = input("Ingrese la presion: ")
            V   = input("Ingrese el volumen: ")
            T   = input("Ingrese la temperatura: ")
            
            reactivos[Compound(c).formula] = {
                "c":   Compound(c), 
                "mol": None if mol == '' else float(mol),
                "gr":  None if gr == '' else float(gr),
                "P":   None if P == '' else float(P),
                "V":   None if V == '' else float(V),
                "T":   None if T == '' else float(T),
            }
            
            return main()
            
        case "2":
            c   = input("Ingrese el reactivo: ")
            mol = input("Ingrese la cantidad de moles: ")
            gr  = input("Ingrese la cantidad de gramos: ")
            P   = input("Ingrese la presion: ")
            V   = input("Ingrese el volumen: ")
            T   = input("Ingrese la temperatura: ")
            
            productos[Compound(c).formula] = {
                "c":   Compound(c), 
                "mol": None if mol == '' else float(mol),
                "gr":  None if gr == '' else float(gr),
                "P":   None if P == '' else float(P),
                "V":   None if V == '' else float(V),
                "T":   None if T == '' else float(T),
            }

            return main()
            
        case "3":
            pass
        
        case _:
            print("Opcion no Valida")
            return main()

    # Termina el Input - Crear la reaccion quimica y la balancea.

    react = " + ".join(reactivos.keys())
    prod  = " + ".join(productos.keys())

    try:  #Inenta generar (chemlib viene con una clase y si no puede generar esa clase es poque el compuesto no existe) la reaccion quimica mediante la formula
        r = Reaction.by_formula(react + "-->" + prod)
        r.balance()
    except: 
        print("\nReaccion quimica no valida\n")
        return
    
    if not get_mass_by_PVT(r): #aca va una escaera de intento de conseguir la data. se reiteran funciones hasta que se consiga la data necesaria o toda la que se pueda.
        print("\nNo se puede calcular la masa\n")
    pass

main()
    
