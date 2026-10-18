""""
Generador de grafo de sufijos

Código adaptado de https://www.studysmarter.es/resumenes/ciencias-de-la-computacion/estructuras-de-datos/arbol-de-sufijos/
"""

class Nodo:
    def __init__(self, sub=""):
        self.sub = sub
        self.hijos = {}
        # (nombre_cadena, indice_inicio)
        self.indices = []

class ArbolSufijos:
    def __init__(self):
        self.raíz = Nodo()

    def añadir_cadena(self, cadena, nombre_id):
        """
        Inserta todos los sufijos. 
        Se suma +1 al índice 'i' para que la posición 0 sea el índice 1.
        """
        for i in range(len(cadena)):
            sufijo = cadena[i:]
            # i + 1 índice a guardar
            self.añadir_sufijo(sufijo, i + 1, nombre_id)

    def añadir_sufijo(self, sufijo, indice, nombre_id):
        """ Añadir sufijo al grafo y reordena las ramas con un único flujo de salida """
        actual = self.raíz
        i = 0
        finalizado = False
        
        while i < len(sufijo) and not finalizado:
            char = sufijo[i]
            if char not in actual.hijos:
                nuevo_nodo = Nodo(sufijo[i:])
                nuevo_nodo.indices.append((nombre_id, indice))
                actual.hijos[char] = nuevo_nodo
                finalizado = True 
                
            else:
                hijo = actual.hijos[char]
                sufijo_previo = hijo.sub
                j = 0
                while j < len(sufijo_previo) and i + j < len(sufijo) and sufijo_previo[j] == sufijo[i+j]:
                    j += 1
                
                if j < len(sufijo_previo):
                    # División de nodo
                    nodo_intermedio = Nodo(sufijo_previo[:j])
                    hijo.sub = sufijo_previo[j:]
                    nodo_intermedio.hijos[hijo.sub[0]] = hijo
                    actual.hijos[char] = nodo_intermedio
                    
                    residuo_sufijo = sufijo[i+j:]
                    if residuo_sufijo:
                        nuevo_hoja = Nodo(residuo_sufijo)
                        nuevo_hoja.indices.append((nombre_id, indice))
                        nodo_intermedio.hijos[residuo_sufijo[0]] = nuevo_hoja
                    else:
                        nodo_intermedio.indices.append((nombre_id, indice))
                    
                    finalizado = True
                else:
                    actual = hijo
                    i += j
        
        if not finalizado:
            actual.indices.append((nombre_id, indice))

    def buscar(self, subcadena):
        """Buscar una subcadena y devolver las posiciones encontradas"""
        actual = self.raíz
        i = 0
        posible = True
        resultado = []

        while i < len(subcadena) and posible:
            char = subcadena[i]
            if char not in actual.hijos:
                posible = False
            else:    
                hijo = actual.hijos[char]
                etiqueta = hijo.sub
                j = 0
                while j < len(etiqueta) and i < len(subcadena) and etiqueta[j] == subcadena[i]:
                    i += 1
                    j += 1
                
                if i < len(subcadena) and j < len(etiqueta):
                    posible = False
                
                actual = hijo

        if posible:
            resultado = self.recorrer_arbol(actual)

        return resultado

    def recorrer_arbol(self, nodo):
        """Recorre el subárbol para recolectar todos los índices almacenados."""
        resultados = list(nodo.indices)
        for hijo in nodo.hijos.values():
            resultados.extend(self.recorrer_arbol(hijo))
        
        return resultados

    def visualizar(self, nodo=None, indent=0):
        """Muestra la estructura del árbol."""
        if nodo is None: 
            nodo = self.raíz
            
        for char in sorted(nodo.hijos.keys()):
            hijo = nodo.hijos[char]
            print("  " * indent + "|--" + hijo.sub + (f" {hijo.indices}" if hijo.indices else ""))
            self.visualizar(hijo, indent + 1)

def ejecutar_busqueda(arbol: ArbolSufijos, patron: str):
    print(f"\nPatrón: '{patron}'")
    resultados = arbol.buscar(patron)

    if resultados:
        for id_cadena, pos in resultados:
            print(f"Encontrado en {id_cadena}, posición: {pos}")
    else:
        print("La cadena no está en ninguna secuencia.")

def main():
    s1 = "AGTGAGT"
    s2 = "ACAGTAGTAT"
    arbol = ArbolSufijos()
    arbol.añadir_cadena(s1 + "$", "S1")
    arbol.añadir_cadena(s2 + "#", "S2")
    print("Arbol")
    arbol.visualizar()
    ejecutar_busqueda(arbol, "AGTA")
    ejecutar_busqueda(arbol, "TAGT")

if __name__ == "__main__":
    main()