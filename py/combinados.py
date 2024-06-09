from chemlib import Compound, Reaction
compound = dict["c":Compound, "mol":float, "gr": float, "P": float, "V": float, "T": float, "CM": float]
reactivos: dict[str: compound] = {}
productos: dict[str: compound] = {}

r = 0.082

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
            if compuesto in [*reactivos]: reactivos[compuesto]["P"] = P
            elif compuesto in [*productos]: productos[compuesto]["P"] = P
       
        if V is None:
            V = n * T * r / P
            #print (V, compuesto, [*reactivos])
            if compuesto in [*reactivos]: reactivos[compuesto]["V"] = V
            elif compuesto in [*productos]: productos[compuesto]["V"] = V
   
        if T is None:
            T = n * r / P * V
            if compuesto in [*reactivos]: reactivos[compuesto]["T"] = T
            elif compuesto in [*productos]: productos[compuesto]["T"] = T


        if n is not None and V is not None:
            CM = n / V
            if compuesto in [*reactivos]: reactivos[compuesto]["CM"] = CM
            elif compuesto in [*productos]: productos[compuesto]["CM"] = CM


    return True






def get_masses(reaction: Reaction): #en base a una reacion, una cantidad (de moles o gr) y un compuesto puede sacar el resto de gr y moles de los otros compuestos.
    compuestos = [x for x in reactivos.keys() if reactivos[x]["gr"] is not None or reactivos[x]['mol'] is not None] + [x for x in productos.keys() if productos[x]["mol"] is not None or productos[x]['mol'] is not None]
    if len(compuestos) == 0:
        return False
   
    compuesto = compuestos[0]
    #print(compuesto, reactivos[compuesto]["mol"])


    #esta funcion habria que cambiarla para que no acepte mas data que la reaccion. deberia funcionar checkeando si exsiten compuestos con n (mol) o gr que tengan valores definidos
   
    listaCompuestos= [*reactivos, *productos]


    if compuesto in [*reactivos]:
        if reactivos[compuesto]["mol"] is not None:
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, moles = reactivos[compuesto]["mol"])


        else: listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, grams = reactivos[compuesto]["gr"])


    else:
        if productos[compuesto]["mol"] is not None:
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, moles = productos[compuesto]["mol"])


        else: listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, grams = productos[compuesto]["gr"])


    for x in listaCompuestos:
        if x in reactivos.keys():
            reactivos[x]["mol"] = listasa[listaCompuestos.index(x)]['moles']
            reactivos[x]["gr"] = listasa[listaCompuestos.index(x)]['grams']
        elif x in productos.keys():
            productos[x]["mol"] = listasa[listaCompuestos.index(x)]['moles']
            productos[x]["gr"] = listasa[listaCompuestos.index(x)]['grams']


    return get_PVT()






def get_mass_by_PVT(reaction: Reaction): #saca la masa de un compuesto en base a la presion, volumen y temperatura... y r que es una constante
    compuestos = [x for x in reactivos.keys() if reactivos[x]["P"] is not None and reactivos[x]["V"] is not None and reactivos[x]["T"] is not None] + [x for x in productos.keys() if productos[x]["P"] is not None and productos[x]["V"] is not None and productos[x]["T"] is not None]
    if len(compuestos) == 0:
        return False
   
    compuesto = compuestos[0] if len(compuestos) > 0 else None


    P = reactivos[compuesto]["P"] if compuesto in reactivos.keys() else productos[compuesto]["P"] if compuesto in productos.keys() else None
    V = reactivos[compuesto]["V"] if compuesto in reactivos.keys() else productos[compuesto]["V"] if compuesto in productos.keys() else None
    T = reactivos[compuesto]["T"] if compuesto in reactivos.keys() else productos[compuesto]["T"] if compuesto in productos.keys() else None
   
    if P is None or V is None or T is None:
        return False
   
    mol = P * V / (r * T)


    if compuesto in reactivos.keys():
        reactivos[compuesto]["mol"] = mol
        get_masses(reaction)
    elif compuesto in productos.keys():
        productos[compuesto]["mol"] = mol
        get_masses(reaction)


    return True






def calcular_reactivo_limitante():
    # Obtener la ecuación balanceada
    react = " + ".join([*reactivos])
    prod = " + ".join([*productos])


    try:
        # Intentar crear la reacción química
        reaccion = Reaction.by_formula(react + "-->" + prod)
        reaccion.balance()
    except:
        print("\nReaccion quimica no valida\n")
        return None


    # Calcular las cantidades proporcionales de moles de productos necesarios para cada reactivo
    cantidades_proporcionales = {}
    for reactivo, datos in reactivos.items():
        if datos["mol"] is not None:
            coeficiente_reactivo = reaccion.coefficients.get(reactivo, 0)
            cantidades_proporcionales[reactivo] = datos["mol"] / coeficiente_reactivo if coeficiente_reactivo != 0 else None


    # Verificar si las cantidades de moles son equivalentes
    mol_equivalentes = all(abs(mol - list(cantidades_proporcionales.values())[0]) < 1e-10 for mol in cantidades_proporcionales.values())


    if mol_equivalentes:
        print("No hay reactivo limitante, las cantidades de moles son equivalentes.")
        return None


    # Obtener la cantidad mínima de moles proporcionales
    min_moles = min(cantidades_proporcionales.values())


    # Determinar el reactivo limitante
    reactivo_limitante = None
    for reactivo, cantidad_proporcional in cantidades_proporcionales.items():
        if cantidad_proporcional == min_moles:
            reactivo_limitante = reactivo
            break


    if reactivo_limitante:
        print(f"El reactivo limitante es: {reactivo_limitante}")
    else:
        print("No se pudo determinar el reactivo limitante.")


    return reactivo_limitante






def print_compuestos(compuestos):
    for formula, data in compuestos.items():
        print(f"Compuesto: {formula}")
        print(f"  - Moles: {data['mol']} mol")
        print(f"  - Gramos: {data['gr']} g")
        print(f"  - Presión: {data['P']} atm")
        print(f"  - Volumen: {data['V']} L")
        print(f"  - Temperatura: {data['T']}° K")
#        print(f"  - Concentración Molar: {data['CM']}")
        if data["V"] is not None and data["mol"] is not None:
            print(f"  - Concentración Molar: {data['CM']} mol/L")
        print()






def main(): #main
    match input("\n1.Ingresar reactivo \n2.Ingresar producto \n3.Calcular \n"): #Input de datos. podría estar mejor? si. Tengo sueño? tmb
         
        case "1":
            c   = input("Ingrese el reactivo: ")
            mol = input("Ingrese la cantidad de moles (mol): ")
            gr  = input("Ingrese la cantidad de gramos (g): ")
            P   = input("Ingrese la presion (atm): ")
            V   = input("Ingrese el volumen (l): ")
            T   = input("Ingrese la temperatura (°K): ")
            CM   = input("Ingrese la concentración molar (mol/L): ")
           
            compuesto_obj = Compound(c)




            reactivos[compuesto_obj.formula] = {
                "c":   compuesto_obj,
                "mol": None if mol == '' else float(mol),
                "gr":  None if gr == '' else float(gr),
                "P":   None if P == '' else float(P),
                "V":   None if V == '' else float(V),
                "T":   None if T == '' else float(T),
                "CM":  None if CM == '' else float(CM),
            }
           
            return main()
           
        case "2":
            c   = input("Ingrese el producto: ")
            mol = input("Ingrese la cantidad de moles: ")
            gr  = input("Ingrese la cantidad de gramos: ")
            P   = input("Ingrese la presion: ")
            V   = input("Ingrese el volumen: ")
            T   = input("Ingrese la temperatura: ")
            CM   = input("Ingrese la concentración molar: ")
           
            compuesto_obj = Compound(c)
           
            productos[Compound(c).formula] = {
                "c":   compuesto_obj,
                "mol": None if mol == '' else float(mol),
                "gr":  None if gr == '' else float(gr),
                "P":   None if P == '' else float(P),
                "V":   None if V == '' else float(V),
                "T":   None if T == '' else float(T),
                "CM":  None if CM == '' else float(CM),
            }




            return main()
           
        case "3":
            pass
       
        case _:
            print("Opcion no Valida")
            return main()




    # Termina el Input - Crear la reaccion quimica y la balancea.


    react = " + ".join([*reactivos])
    prod  = " + ".join([*productos])


    if react == prod:
        pass


    try:  #Inenta generar (chemlib viene con una clase y si no puede generar esa clase es poque el compuesto no existe) la reaccion quimica mediante la formula
        reaccion = Reaction.by_formula(react + "-->" + prod)
        reaccion.balance()
    except:
        print("\nReaccion quimica no valida\n")
        return
   
    get_mass_by_PVT(reaccion) #aca va una escaera de intento de conseguir la data. se reiteran funciones hasta que se consiga la data necesaria o toda la que se pueda.
    get_masses(reaccion)
    get_PVT()
   
    print(f"\nEcuación balanceada: ")
    print(reaccion)
    print("\nReactivos:")
    print_compuestos(reactivos)
    print("\nProductos:")
    print_compuestos(productos)
    calcular_reactivo_limitante()
    return


main()


from chemlib import Compound, Reaction
compound = dict["c":Compound, "mol":float, "gr": float, "P": float, "V": float, "T": float, "CM": float]
reactivos: dict[str: compound] = {}
productos: dict[str: compound] = {}




r = 0.082




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
            if compuesto in [*reactivos]: reactivos[compuesto]["P"] = P
            elif compuesto in [*productos]: productos[compuesto]["P"] = P
       
        if V is None:
            V = n * T * r / P
            print (V, compuesto, [*reactivos])
            if compuesto in [*reactivos]: reactivos[compuesto]["V"] = V
            elif compuesto in [*productos]: productos[compuesto]["V"] = V
   
        if T is None:
            T = n * r / P * V
            if compuesto in [*reactivos]: reactivos[compuesto]["T"] = T
            elif compuesto in [*productos]: productos[compuesto]["T"] = T




        if n is not None and V is not None:
            CM = n / V
            if compuesto in [*reactivos]: reactivos[compuesto]["CM"] = CM
            elif compuesto in [*productos]: productos[compuesto]["CM"] = CM




    return True




def get_masses(reaction: Reaction): #en base a una reacion, una cantidad (de moles o gr) y un compuesto puede sacar el resto de gr y moles de los otros compuestos.
   
    compuestos = [x for x in reactivos.keys() if reactivos[x]["gr"] is not None or reactivos[x]['mol'] is not None] + [x for x in productos.keys() if productos[x]["mol"] is not None or productos[x]['mol'] is not None]
    if len(compuestos) == 0:
        return False
   
    compuesto = compuestos[0]
    print(compuesto, reactivos[compuesto]["mol"])




    #esta funcion habria que cambiarla para que no acepte mas data que la reaccion. deberia funcionar checkeando si exsiten compuestos con n (mol) o gr que tengan valores definidos
   
    listaCompuestos= [*reactivos, *productos]




    if compuesto in [*reactivos]:
        if reactivos[compuesto]["mol"] is not None:
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, moles = reactivos[compuesto]["mol"])




        else: listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, grams = reactivos[compuesto]["gr"])




    else:
        if productos[compuesto]["mol"] is not None:
            listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, moles = productos[compuesto]["mol"])




        else: listasa = reaction.get_amounts(listaCompuestos.index(compuesto) + 1, grams = productos[compuesto]["gr"])








    for x in listaCompuestos:




        if x in reactivos.keys():
            reactivos[x]["mol"] = listasa[listaCompuestos.index(x)]['moles']
            reactivos[x]["gr"] = listasa[listaCompuestos.index(x)]['grams']
        elif x in productos.keys():
            productos[x]["mol"] = listasa[listaCompuestos.index(x)]['moles']
            productos[x]["gr"] = listasa[listaCompuestos.index(x)]['grams']




    return get_PVT()
















def get_mass_by_PVT(reaction: Reaction): #saca la masa de un compuesto en base a la presion, volumen y temperatura... y r que es una constante
    compuestos = [x for x in reactivos.keys() if reactivos[x]["P"] is not None and reactivos[x]["V"] is not None and reactivos[x]["T"] is not None] + [x for x in productos.keys() if productos[x]["P"] is not None and productos[x]["V"] is not None and productos[x]["T"] is not None]
    if len(compuestos) == 0:
        return False
   
    compuesto = compuestos[0] if len(compuestos) > 0 else None




    P = reactivos[compuesto]["P"] if compuesto in reactivos.keys() else productos[compuesto]["P"] if compuesto in productos.keys() else None
    V = reactivos[compuesto]["V"] if compuesto in reactivos.keys() else productos[compuesto]["V"] if compuesto in productos.keys() else None
    T = reactivos[compuesto]["T"] if compuesto in reactivos.keys() else productos[compuesto]["T"] if compuesto in productos.keys() else None
   
    if P is None or V is None or T is None:
        return False
   
    mol = P * V / (r * T)




    if compuesto in reactivos.keys():
        reactivos[compuesto]["mol"] = mol
        get_masses(reaction)
    elif compuesto in productos.keys():
        productos[compuesto]["mol"] = mol
        get_masses(reaction)








    return True






def calcular_reactivo_limitante():
    # Obtener la ecuación balanceada
    react = " + ".join([*reactivos])
    prod = " + ".join([*productos])


    try:
        # Intentar crear la reacción química
        reaccion = Reaction.by_formula(react + "-->" + prod)
        reaccion.balance()
    except:
        print("\nReaccion quimica no valida\n")
        return None


    # Obtener las cantidades proporcionales de moles de productos para cada reactivo
    moles_producidos_por_reactivo = {}
    for reactivo in reactivos:
        mol_reactivo = reactivos[reactivo]["mol"]
        if mol_reactivo is not None:
            coeficiente_reactivo = reaccion.coefficients.get(reactivo, 0)
            moles_producidos_por_reactivo[reactivo] = mol_reactivo / coeficiente_reactivo if coeficiente_reactivo != 0 else None


    # Verificar si las cantidades de moles son equivalentes
    equivalentes = True
    cantidad_anterior = None
    for mol in moles_producidos_por_reactivo.values():
        if cantidad_anterior is not None and abs(mol - cantidad_anterior) > 1e-10:
            equivalentes = False
            break
        cantidad_anterior = mol


    if equivalentes:
        print("No hay reactivo limitante, las cantidades de moles son equivalentes.")
        return None


    # Determinar el reactivo limitante
    reactivo_limitante = min(moles_producidos_por_reactivo, key=moles_producidos_por_reactivo.get)


    print(f"El reactivo limitante es: {reactivo_limitante}")
    return reactivo_limitante






def print_compuestos(compuestos):
    for formula, data in compuestos.items():
        print(f"Compuesto: {formula}")
        print(f"  - Moles: {data['mol']}")
        print(f"  - Gramos: {data['gr']}")
        print(f"  - Presión: {data['P']}")
        print(f"  - Volumen: {data['V']}")
        print(f"  - Temperatura: {data['T']}")
        print(f"  - Concentración Molar: {data['CM']}")
        print()












def main(): #main




    match input("\n1.Ingresar reactivo \n2.Ingresar producto \n3.Calcular \n"): #Input de datos. podría estar mejor? si. Tengo sueño? tmb
         
        case "1":
            c   = input("Ingrese el reactivo: ")
            mol = input("Ingrese la cantidad de moles: ")
            gr  = input("Ingrese la cantidad de gramos: ")
            P   = input("Ingrese la presion: ")
            V   = input("Ingrese el volumen: ")
            T   = input("Ingrese la temperatura: ")
            CM   = input("Ingrese la concentración molar: ")
           
            compuesto_obj = Compound(c)




            reactivos[compuesto_obj.formula] = {
                "c":   compuesto_obj,
                "mol": None if mol == '' else float(mol),
                "gr":  None if gr == '' else float(gr),
                "P":   None if P == '' else float(P),
                "V":   None if V == '' else float(V),
                "T":   None if T == '' else float(T),
                "CM":  compuesto_obj.molar_mass(),
            }
           
            return main()
           
        case "2":
            c   = input("Ingrese el producto: ")
            mol = input("Ingrese la cantidad de moles: ")
            gr  = input("Ingrese la cantidad de gramos: ")
            P   = input("Ingrese la presion: ")
            V   = input("Ingrese el volumen: ")
            T   = input("Ingrese la temperatura: ")
            CM   = input("Ingrese la concentración molar: ")
           
            compuesto_obj = Compound(c)
           
            productos[Compound(c).formula] = {
                "c":   compuesto_obj,
                "mol": None if mol == '' else float(mol),
                "gr":  None if gr == '' else float(gr),
                "P":   None if P == '' else float(P),
                "V":   None if V == '' else float(V),
                "T":   None if T == '' else float(T),
                "CM":  compuesto_obj.molar_mass(),
            }




            return main()
           
        case "3":
            pass
       
        case _:
            print("Opcion no Valida")
            return main()




    # Termina el Input - Crear la reaccion quimica y la balancea.




    react = " + ".join([*reactivos])
    prod  = " + ".join([*productos])




    if react == prod:
        pass








    try:  #Inenta generar (chemlib viene con una clase y si no puede generar esa clase es poque el compuesto no existe) la reaccion quimica mediante la formula
        reaccion = Reaction.by_formula(react + "-->" + prod)
        reaccion.balance()
    except:
        print("\nReaccion quimica no valida\n")
        return
   
    get_mass_by_PVT(reaccion) #aca va una escaera de intento de conseguir la data. se reiteran funciones hasta que se consiga la data necesaria o toda la que se pueda.
    get_masses(reaccion)
    get_PVT()
   
    print("\nReactivos:")
    print_compuestos(reactivos)
    print("\nProductos:")
    print_compuestos(productos)
    calcular_reactivo_limitante()
    return
main()