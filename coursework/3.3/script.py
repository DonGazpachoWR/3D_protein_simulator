import numpy as np

def matriz_BWT(secuencia):

    lista1 = []
    lista_nt = []
    lista_i = []
    for iteracion in range(len(secuencia)):
        lista1.append([secuencia.copy(), iteracion])
        secuencia.append(secuencia.pop(0))

        
    lista1.sort()
    # Mostrar matriz BWT
    matriz_0 = np.array([item[0] + [str(item[1])] for item in lista1])
    lista_nt = [lista[0][-1] for lista in lista1]
    lista_i = [lista[-1] for lista in lista1]

    return matriz_0, lista_nt, lista_i



def matriz_BWT_invertida(lista_nt, lista_i):
    # Mostrar matriz BWT invertida
    n = len(lista_nt)

    matriz = np.empty((n, n + 1), dtype = object)
    matriz[:] = ""

    matriz[:, n] = lista_i
    matriz[:, n - 1] = lista_nt
    matriz[:, 0] = sorted(lista_nt)

    for col in range(1,len(lista_nt) - 1):

        col_prev = [''.join(matriz[i,:col]) for i in range(n)]
        
        col_fin = matriz[:, n - 1]
        col_nueva = sorted([fin + prev for prev,fin in zip(col_prev, col_fin)])
        matriz[:,col] = [ item[col] for item in col_nueva ]

    return matriz

# calculo posiciones

def matriz_BWT_minima(lista_nt):

    primero = []
    contador = {}
    for char in sorted(lista_nt):
        contador[char] = contador.get(char, 0) + 1
        primero.append([char, contador[char]])

    ultimo = []
    contador = {}

    for char in lista_nt:
        contador[char] = contador.get(char, 0) + 1
        ultimo.append([char, contador[char]])

    pack = [ prim + ult for prim, ult in zip(primero, ultimo)]

    return pack


# algoritmo busqueda
def BTW_busqueda(query, pack, lista_i):

    char_actual = query[-1]
    top = -1
    bottom = -1

    # Buscar la primera y última aparición del char
    for i, fila in enumerate(pack):
        if fila[0] == char_actual:
            if top == -1:
                top = i  
            bottom = i

    # buscar hacia atrás
    if top != -1:
        i = len(query) - 2
        
        while i >= 0 and top != -1:
            char_buscado = query[i]
            nuevo_top = -1
            nuevo_bottom = -1
            
            # Mirar la columna L en el rango
            for fila_idx in range(top, bottom + 1):
                L_char = pack[fila_idx][2]
                L_count = pack[fila_idx][3]
                
                if L_char == char_buscado:
                    f_idx = 0
                    encontrado = False
                    
                    while f_idx < len(pack) and not encontrado:
                        f_fila = pack[f_idx]
                        if f_fila[0] == L_char and f_fila[1] == L_count:
                            if nuevo_top == -1:
                                nuevo_top = f_idx
                            nuevo_bottom = f_idx
                            encontrado = True  
                        f_idx += 1
                        
            # Actualizamos 
            top = nuevo_top
            bottom = nuevo_bottom
            
            # Bajar el índice
            i -= 1

    # posiciones encontrandas
    solucion = []
    if top != -1 and bottom != -1:
        for fila_idx in range(top, bottom + 1):
            solucion.append(lista_i[fila_idx] + 1)

    return solucion

def main():
    secuencia = "ACTGAACTGACAATCA"
    secuencia += "#"
    secuencia = list(secuencia)
    query = "GACA"

    matriz_0, lista_nt, lista_i = matriz_BWT(secuencia)
    matriz = matriz_BWT_invertida(lista_nt, lista_i)
    pack = matriz_BWT_minima(lista_nt)
    solucion = BTW_busqueda(query, pack, lista_i)

    print("Matriz BWT\n")
    for fila in matriz_0:
        print(fila)
    print("\nMatriz BWT invertida\n")
    print(matriz)
    print("\nSecuencia BWT: ", ''.join(lista_nt))
    print("Posiciones de la query ", query, "en la cadena: ", ','.join(str(pos) for pos in solucion))


if __name__ == "__main__":
    main()





        
    
