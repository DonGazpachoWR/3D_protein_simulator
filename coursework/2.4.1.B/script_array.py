"""
Generador de Array de Sufijos
"""

class ArraySufijos:
    def __init__(self, secuencias):
        """
        Inicializa la estructura del array.
        
        Args:
            secuencias (dict): Diccionario donde la clave es el nombre (S1, S2...) 
                              y el valor es la cadena de ADN.
        """
        self.secuencias = secuencias 
        # (texto_sufijo, id_secuencia, posicion_original)
        self.tabla_auxiliar = []      
        self.as_indices = []        
        self.construir_array()

    def construir_array(self):
        """ Extrae todos los sufijos de todas las secuencias y los ordena """
        lista_sufijos = []
        for nombre, cadena in self.secuencias.items():
            for i in range(len(cadena)):
                # sufijo, identificador y posició
                lista_sufijos.append((cadena[i:], nombre, i + 1))
        
        # Ordenar la tabla 
        self.tabla_auxiliar = sorted(lista_sufijos, key=lambda x: x[0])
        self.as_indices = self.tabla_auxiliar

    def buscar(self, subcadena)-> list:
        """
        Busca una subcadena utilizando búsqueda binaria sobre el Array .
        
        Args:
            subcadena (str): El patrón
            
        Return:
            list: Lista de tuplas (id_secuencia, posicion) con todas las coincidencias.
        """
        resultados = []
        bajo = 0
        alto = len(self.as_indices) - 1
        encontrado = False
        
        # Búsqueda binaria
        while bajo <= alto and not encontrado:
            medio = (bajo + alto) // 2
            sufijo_actual = self.as_indices[medio][0]
            
            if sufijo_actual.startswith(subcadena):
                encontrado = True
                
                # Encontrar todas las coincidencias contiguas hacia atrás
                izq = medio
                while izq >= 0 and self.as_indices[izq][0].startswith(subcadena):
                    resultados.append((self.as_indices[izq][1], self.as_indices[izq][2]))
                    izq -= 1
                
                # Encontrar todas las coincidencias contiguas hacia adelante
                der = medio + 1
                while der < len(self.as_indices) and self.as_indices[der][0].startswith(subcadena):
                    resultados.append((self.as_indices[der][1], self.as_indices[der][2]))
                    der += 1
            
            elif subcadena < sufijo_actual:
                alto = medio - 1
            else:
                bajo = medio + 1
        
        return resultados


def mostrar_array(array: ArraySufijos):
    """Imprime el array"""
    print(f"{'Índice':<8} | {'Sufijo':<12} | {'Origen':<8} | {'Pos'}")
    print("-" * 42)
    for i, item in enumerate(array.as_indices):
        print(f"{i:<8} | {item[0]:<12} | {item[1]:<8} | {item[2]}")


def ejecutar_busqueda(array: ArraySufijos, patrones):
    """Ejecuta la orden de busqueda de patrones en el array"""
    for patron in patrones:
        res = array.buscar(patron)
        print(f"\nPatrón: '{patron}':")
        if res:
            for r in res:
                print(f"Encontrado en {r[0]}, posición: {r[1]}")
        else:
            print("La subcadena no está en la base de datos.")


def main():
    datos = {
        "S1": "AGTGAGT",
        "S2": "ACAGTAG",
        "S3": "TATTTCGA"
    }
    array = ArraySufijos(datos)
    mostrar_array(array)
    patrones = ["AGTA", "ATT"]
    ejecutar_busqueda(array, patrones)
    

if __name__ == "__main__":
    main()