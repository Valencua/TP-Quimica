from chemlib import Compound, Reaction
compound = dict["c":Compound, "mol":float, "gr": float, "P": float, "V": float, "T": float, "MM": float]
reactivos: dict[str: compound] = {}
productos: dict[str: compound] = {}

def main():

    r = 0.082

    match input("\n1.Ingresar reactivo \n2.Ingresar producto \n3.Calcular \n"):

        case "1":
                c = input("Ingrese el reactivo: ")
                mol = input("Ingrese la cantidad de moles: ")
                gr = input("Ingrese la cantidad de gramos: ")
                P = input("Ingrese la presion: ")
                V = input("Ingrese el volumen: ")
                T = input("Ingrese la temperatura: ")

                reactivos[c] = {
                    "c": Compound(c), 
                    "mol": None if mol == '' else float(mol),
                    "gr": None if gr == '' else float(gr),
                    "P": None if P == '' else float(P),
                    "V": None if V == '' else float(V),
                    "T": None if mol == '' else float(T),
                }

                return main()
        case "2":
                c = input("Ingrese el producto: ")
                mol = input("Ingrese la cantidad de moles: ")
                gr = input("Ingrese la cantidad de gramos: ")
                P = input("Ingrese la presion: ")
                V = input("Ingrese el volumen: ")
                T = input("Ingrese la temperatura: ")

                productos[c] = {
                    "c": Compound(c), 
                    "mol": None if mol == '' else float(mol),
                    "gr": None if gr == '' else float(gr),
                    "P": None if P == '' else float(P),
                    "V": None if V == '' else float(V),
                    "T": None if mol == '' else float(T),
                }

                return main()

        case "3":
            react = " + ".join(reactivos.keys())
            prod = " + ".join(productos.keys())
            r = Reaction.by_formula(react + "-->" + prod)

        case _:
            print("Opcion no Valida")
            return main()

    r.balance()

main()