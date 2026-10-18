"""
Clase Grafo
Autor: Adrián Berenguer Agustí
Fecha: 12/03/2026
Versión: 1.0

"""
from typing import TextIO
import re
import statistics
import random
import math
from matplotlib import pyplot as plt

class Grafo:
    """
    Grafo compuesto por nodos y arcos (relaciones entre nodos)

    Atributos
    ---------
    nodos : lista
        Lista de objetos nodo
    arcos: lista
        Lista de objetos arco
    
    Métodos 
    --------
    insertar_nodo():
        Añade un nodo a la lista de nodos
    insertar_arco():
        Actualiza A True la relación entre dos nodos
    nodos_adyacentes()
        Devuelve una lista de los nodos adyacentes al grado diana
    listar_nodos():
        Se obtiene lista de los nodos existentes
    listar_arcos():
        Se obtiene lista booleanos sobre las relaciones entre nodos
    abrir_fichero():
        Abre fichero
    extraer_nodos():
        A partir de archivo se cargan los nodos
        del objeto grafo
    extraer_relaciones():
        A partir de un archivo carga arcos entre nodos
    grafo_apartir_fichero():
        A partir de una ruta se obtiene un Grafo rellenado con nodos y arcos
    estadisticas():
        nº de nodos e  interacciones directas
    algoritmo_floyd():
        calcula distancia mínima entre todos los nodos
    mostrar_distancias():
        imprime camino mínimo entre todos los nodos
    d_max():
        muestra la distancia entre los 2 nodos del grafo más alejados entre ellos 
    d_promedio():
        Muestra la distancia promedio de la red
    ordenar_nodos_adyacentes():
        Muestra los nodos ordenados por nº de nodos adyacentes
    ordenar_d_promedio():
        Muestra los nodos ordenados por distancia promedio al resto de nodos
    com_nod_n():
        Calcula com nod de n nodos
    hill_climbing():
        agrupa nodos en n grupos por máximos locales de com nod
    sim_annealing():
        agrupa nodos en n grupos por máximo global de com nod
    mostrar_nodos():
        imprime el nombre de los nodos de una lista de nodos
    graficar():
        genera plot de valor com nod para n iteraciones para y ejecuciones del 
        algoritmo
    """
    class Nodo:
        """
        Objeto nodo que está conectado a otros nodos en el objeto grafo
        
        Atributos
        ----------
        info: str
            información almacenada en el nodo, puede ser un objeto planeta
        visitado: bool
            True si ya ha saido visitado, False si no
        
        Métodos
        --------
        visitar():
            marca un nodo como visitado
        limpiar_visita():
            marca un nodo como no visitado
        fue_visitado():
            devuelve si el nodo se encuentra visitado o no
        """
        def __init__(self, info: str = None) -> None: 
            self.info = info
            self.visitado = False
        
        def visitar(self)-> None:
            """ Marca el nodo como visitado"""
            self.visitado = True
        
        def limpiar_visita(self)-> None:
            """Marca el nodo como no visitado"""
            self.visitado = False
        
        def fue_visitado(self)-> bool:
            """Devuelve si el nodo ha sido visitado o no"""

            return self.visitado
        
    class Arco:
        """
        Relación entre dos nodos

        Atributos
        ----------
        existe: bool
            Si True indica que hay conexión entre los dos nodos
        info: str
            información que pudiera ser descriptiva en el arco
        """
        def __init__(self, min: list = None, max: list = None, 
                     distancia: float = float('inf')) -> None: 
            """min almacena camino más corto entre nodos y max el más largo"""
            self.existe = bool()
            if min:
                self.min = min
            else:  
                self.min = []
            if max:
                self.max = max
            else:
                self.max = []
            self.distancia = distancia
    
    def __init__(self, n_nod: int = 0) -> None:
        self.nodos_ID = [Grafo.Nodo() for x in range(n_nod)]
        self.arcos = [[Grafo.Arco() for y in range(n_nod)] 
                      for x in range(n_nod)]
        self.nodos_nombre = [ x.info for x in self.nodos_ID]
    
    def insertar_nodo(self, nuevo_nodo: Nodo = None) -> None:
        """Añade un nodo a la lista de nodos del objeto grafo"""
        # print("Nuevo nodo nombre: ", nuevo_nodo)
        self.nodos_ID.append(Grafo.Nodo(nuevo_nodo))
        nombre = self.nodos_ID[-1].info
        self.nodos_nombre.append(nombre)
        for fila_existente in self.arcos:
            # Añade una nueva columna
            fila_existente.append(Grafo.Arco()) 
        # Añade nueva fila
        nueva_fila = [Grafo.Arco() for x in range(len(self.nodos_ID))]
        self.arcos.append(nueva_fila)
        # Valor 0 cuando es el nodo consigo mismo
        idx_nuevo = len(self.nodos_ID) - 1
        self.arcos[idx_nuevo][idx_nuevo].distancia = 0
    
    def insertar_arco(self, origen: int, destino: int, info: str = None)-> None: 
        """
        Actualiza como True la relacion entre un nodo origen y otro destino
        Se actualiza en ambos sentidos pues se trata de grafo dirigido
        Args:
            origen: posicion del nodo origen en la lista self.nodos_ID
            destino: posicion del nodo destino en la lista self.nodos_ID
            info: información relevante a añadir en el nodo
        """
        self.arcos[origen][destino].existe = True
        self.arcos[origen][destino].min.append(self.nodos_nombre[destino])
        # print(origen, destino, self.arcos[origen][destino].min)
        # input()
        # distancia 1 entre los 2 nodos
        self.arcos[origen][destino].distancia = 1

    def nodos_adyacentes (self, ori:int)-> list["Arco"]:
        """
        Devuelve lista de posiciones de self.nodos_ID que son adyacentes a la
        posicion del nodo introducido como argumento
        Args:
            ori: posicion del nodo introducido como argumento
        Return:
            lista de nodos adyacentes
        """
        return [destino for destino in range(0,len(self.arcos[ori]))
                if self.arcos[ori][destino].existe]
    
    def listar_arcos(self) -> None:
        """Imprime la relación entre nodos de la lista de arcos"""
        for x in range (0,len(self.arcos)):
            for y in range (0,len(self.arcos)):
                print(self.arcos[x][y].existe)

    def listar_nodos(self)-> None:
        """Imprime por pantalla todos los nodos"""
        for nodo in self.nodos_ID:
            print(nodo.info)

    def abrir_fichero(self, ruta:str)-> TextIO | None:
        """Abre fichero con try
        Return: objetivo tipo archivo o None si no lo pudo abrir
        """
        try:
            archivo = open(ruta, "r", encoding = "utf-8")
        except:
            print("\nError al abrir el archivo")
            archivo = None
        else: print("exito al leer: ", ruta)
        
        return archivo
    
    def extraer_nodos(self, archivo: TextIO)-> None:
        """A partir de archivo con relaciones entre nodos se generan los nodos 
        y arcos sin relaciones"""
        for linea in archivo:
            linea = re.split(r'\s+', linea.strip())
            if linea:
                for nodo in linea: 
                    if nodo not in self.nodos_nombre:
                        self.insertar_nodo(nodo)
                        
    def extraer_relaciones(self, archivo: TextIO)-> None:
        """A partir de archivo con relaciones entre nodos se generan los arcos 
        que determinan relaciones entre nodos"""
        for linea in archivo:
            linea = re.split(r'\s+', linea.strip())
            if linea:
                for nodo1 in linea: 
                    x = self.nodos_nombre.index(nodo1)
                    for nodo2 in linea:
                        # input()
                        # print(nodo1,nodo2)
                        y = self.nodos_nombre.index(nodo2)
                        if nodo1 != nodo2:
                            # print("relacion")
                            # print(self.arcos[x][y].existe)
                            if self.arcos[x][y].existe== False:
                                # print(x,y)
                                self.insertar_arco(x, y) 


    def grafo_apartir_fichero(self, ruta:str) -> bool:
        """
        Lee archivo que contiene relaciones entre nodos
        Genera contenido del grafo

        Args:
            ruta1: ruta del archivo que contiene las relaciones entre nodos
        Return: booleano que indica si se pudo completar el proceso
        """
        archivo = self.abrir_fichero(ruta) 
        ok = False
        if archivo:
            self.extraer_nodos(archivo)
            archivo.seek(0)
            self.extraer_relaciones(archivo)
            archivo.close()
            ok = True

        return ok
    
    def estadisticas(self)-> None:
        """Muestra nº de nodos e  interacciones directas"""
        # Nº de nodos
        print("Nº total de nodos: ", len(self.nodos_ID))
        # Nº de interacciones directas
        suma = sum(arco.existe for nodo in self.arcos for arco in nodo)
        print("Nº de interacciones directas: ", suma)

    def algoritmo_floyd(self):
        """Obtiene la distancia mínima entre todos los nodos"""
        # N es el número de nodos (el tamaño de la matriz)
        n = len(self.nodos_ID)
        # k es el nodo intermedio (el "puente")
        for k in range(0,n):
            # i es el nodo de origen
            for i in range(0,n):
                # j es el nodo de destino
                for j in range(0,n):
                    if i != j:
                        arco = self.arcos[i][j]
                        # si no es nodo adyacente
                        if not arco.existe:
                            arco_puente1 = self.arcos[i][k]
                            arco_puente2 = self.arcos[k][j]
                            # Calculamos si es mejor ir directo (i -> j) 
                            # o pasar por el puente k (i -> k -> j)
                            d_con_puente = arco_puente1.distancia + \
                            arco_puente2.distancia
                            d_existente = arco.distancia
                            if d_existente > d_con_puente:
                                arco.distancia = d_con_puente
                                arco.min = arco_puente1.min + arco_puente2.min

    def d_max(self)-> int:
        """
        Devuelve la distancia entre los 2 nodos del grafo más alejados 
        entre ellos. Si no se ejecuta
        algoritmo_floyd antes da valor 1 por nodos adyacentes
        """
        d_max = 0
        for x in range(len(self.arcos)):
            for y in range(len(self.arcos)):
                if x < y:
                    d_nueva = self.arcos[x][y].distancia
                    if d_nueva > d_max and d_nueva != float('inf'):
                        d_max = d_nueva
        
        return d_max
    
    def d_promedio(self)-> float:
        """Devuelve la distancia promedio de la red"""
        lista = [self.arcos[x][y].distancia 
                 for x in range(len(self.arcos)) 
                 for y in range(len(self.arcos)) if x < y]
        media = statistics.mean(lista)
        print("Dsitancia promedio de la red: ", f"{media:.2f}")
        
        return media

    def mostrar_distancias(self)-> None:
        """muestra camino mínimo entre nodos"""
        print("\nMatriz de distancias mínimas entre nodos\n")
        cabecera = [item if len(item)== 2 else " " + item 
                    for item in self.nodos_nombre]
        print("  ",' '.join(cabecera))
        for x, id in enumerate(self.nodos_ID):
            lista = [" " + str(self.arcos[x][y].distancia) 
                    if x < y else ( " 0" if x == y else " -") 
                    for y in range (len((self.nodos_ID)))]
            if len(id.info) == 1:
                lista.insert(0, " " + id.info)
            else:
                lista.insert(0, id.info)
            print(' '.join(lista))

    def ordenar_nodos_adyacentes(self)-> None:
        """
        Muestra listado de nodos ordenado por número de interacciones
        directar y por orden alfabético
        """
        diccionario = {nodo.info: len(self.nodos_adyacentes(x)) 
                       for x,nodo in enumerate(self.nodos_ID)}
        # print(diccionario)
        # print(statistics.mean(diccionario.values()))
        # print(sorted(diccionario.items(), key = lambda item:item[1]))
        # print("\nmétodo combinado")
        salida = sorted(diccionario.items(), 
                        key = lambda item: (item[1], item[0]))
        for elemento in salida:
            if len(elemento[0]) == 1:
                print(elemento[0], "  : ", elemento[1])
            else:
                print(elemento[0], " : ", elemento[1] )

    def ordenar_d_promedio(self)-> None:
        """Muestra listado de nodos ordenado por distancia promedio al 
        resto de nodos"""
        diccionario = {nodo.info: statistics.mean([self.arcos[x][y].distancia
                        for y in range(len(self.arcos)) if x != y]) 
                        for x, nodo in enumerate(self.nodos_ID)}

        salida = sorted(diccionario.items(), 
                        key = lambda item: (item[1], item[0]))
        print("Nodos ordenados por distancia promedio al resto de nodos")
        for elemento in salida:
            if len(elemento[0])==1:
                print(elemento[0], "  : ", f"{elemento[1]:.2f}")
            else:
                print(elemento[0], " : ", f"{elemento[1]:.2f}")   
    
    def com_nod_n(self, l_grupos_nodos)-> float:
        """Calcula valor com nod de n listas de nodos
        Args:
            listas: lista de listas de nodos. Imprescindible para que funcione
        """
        listas_transformadas = []
        l, s, f1, ss, n, f2 = 0, 0, 0, 0, 0, 0
        for lista in l_grupos_nodos:
            i_n = [self.nodos_ID.index(x) for x in lista]
            listas_transformadas.append(i_n)
            k = len(i_n)
            l += k * (k - 1) / 2
            s += sum(self.arcos[x][y].distancia
                    for a, x in enumerate(i_n) for b, y in enumerate(i_n)
                    if a < b)
        ss = sum(self.arcos[x][y].distancia 
                 for a, i_1 in enumerate(listas_transformadas)
                 for b, i_2 in enumerate(listas_transformadas)
                 if a < b
                 for x in i_1 for y in i_2)
        n = sum(len(i_1) * len(i_2)
                 for a, i_1 in enumerate(listas_transformadas)
                 for b, i_2 in enumerate(listas_transformadas)
                 if a < b)
        f1 = s / l
        f2 = ss / n
        
        return f2 / f1
    
    def agrupacion_aleatoria(self, args: list[int]) -> int:
        """Agrupar aleatoriamente todos los nodos del grafo en n grupos de 
        longitud m
        Args:
            *args: numero de nodos para m grupos
        Return:
            lista de listas de objetos nodos agrupados aleatoriamente
            número de grupos generados
        """
        copia = self.nodos_ID.copy()
        l_grupos_nodos = []
        # Primera agrupación de nodos aleatoria
        for i in args:
            lista1 = []
            for j in range(i):
                lista1.append(random.choice(copia))
                copia.remove(lista1[-1])
            # Añadir lista de nodos a lista general
            l_grupos_nodos.append(lista1)
        num_grupos = len(l_grupos_nodos)

        return l_grupos_nodos, num_grupos
    
    def escoger_grupos(self, num_grupos: int)-> tuple[int, int]:
        """Escoge dos de los n elementos de una lista"""
        x = 0
        y = 0
        while x == y:
            x = random.choice(range(num_grupos))
            y = random.choice(range(num_grupos))
        
        return x, y

    def modificar_agrupacion(self, l_grupos_nodos, grupo: int, posicion: int, nodo: Nodo)-> list:
        """Modifica el nodo contenido en una lista de listas de nodos"""
        l_grupos_nodos[grupo][posicion] = nodo

        return l_grupos_nodos

    def intercambiar_nodos(self, l_grupos_nodos, num_grupos: int)-> tuple[int, 
                            int, int, int, Nodo, Nodo]:
        """
        Intercambia dos nodos de una lista de nodos agrupados.
        Args:
            lista de listas de nodos agrupados
        Return:
            x: primer grupo intercambiado
            x1: posicion del nodo del primer grupo
            y: segundo grupo intercambiado
            y1: posicion del nodo del segundo grupo
            id_x1: Nodo del primero grupo
            id_y1: Nodo del segundo grupo
        """
        # Grupos de nodos escogidos sobre los que interambiar nodos
        g1, g2 = self.escoger_grupos(num_grupos)
        # Escoger nodo para cada grupo de nodos
        p1 = random.choice(range(len(l_grupos_nodos[g1])))
        p2 = random.choice(range(len(l_grupos_nodos[g2])))
        id_x1 = l_grupos_nodos[g1][p1]
        id_y1 = l_grupos_nodos[g2][p2]
        # Intercambiar nodos entre los grupos
        l_grupos_nodos = self.modificar_agrupacion(l_grupos_nodos, g1, p1, id_y1)
        l_grupos_nodos = self.modificar_agrupacion(l_grupos_nodos, g2, p2, id_x1)

        return g1, p1, g2, p2, id_x1, id_y1, l_grupos_nodos

    def hill_climbing(self, *args: int)-> tuple[list[float], float, list[list[Nodo]]]:
        """
        Agrupa por hill climbing el grafo en n agrupaciones. 
        Obtiene maximo local de la función com nod.
        Args:
            *args: numero de nodos para n grupos
        Return:
            lista con el valor de com nod para cada iteración
        """
        l_grupos_nodos, n_grupos = self.agrupacion_aleatoria(args)
        iteracion, comnod_max, paciencia = 0, 0, 0
        resultado = []
        # Maximo de iteraciones sin mejoras
        paciencia_max = len(self.nodos_ID) * 5
        while paciencia_max > paciencia:
            paciencia += 1
            iteracion += 1
            g1, p1, g2, p2, id_x1, id_y1, l_grupos_nodos = self.intercambiar_nodos(l_grupos_nodos, n_grupos)
            # Calcular com nod
            comnod = round(self.com_nod_n(l_grupos_nodos), 5)
            resultado.append(comnod)
            # Evaluar si el nuevo cambio es mejor o no a lo que había
            if comnod > comnod_max:
                comnod_max = comnod
                iteracion_comnod_max = iteracion
                paciencia = 0
            else:
                # Deshacer cambio
                l_grupos_nodos = self.modificar_agrupacion(l_grupos_nodos, g1, p1, id_x1)
                l_grupos_nodos = self.modificar_agrupacion(l_grupos_nodos, g2, p2, id_y1)
        print("iteracion: ", iteracion_comnod_max, " valor: ", comnod_max)

        return resultado, comnod_max, l_grupos_nodos

    def sim_annealing(self, *args: int)-> tuple[list[float], float]:
        """
        Agrupa por simulated annealing el grafo en n agrupaciones. 
        Obtiene maximo global de la función com nod.
        Args:
            *args: numero de nodos para n grupos
        Return:
            lista con el valor de com nod para cada iteración
            valor com nod máximo
            lista de nodos agrupados
        """
        l_grupos_nodos, n_grupos = self.agrupacion_aleatoria(args)
        # Parámetros simulated annealing
        # Valor inicial de temperatura
        t = 1
        t_min = 0.01
        alpha = 0.99999 # enfriamiento
        comnod_max_relativo, iteracion_max_abs, comnod_max_abs = 0, 0, 0
        resultado = []
        while t > t_min:
            g1, p1, g2, p2, id_x1, id_y1, l_grupos_nodos = self.intercambiar_nodos(l_grupos_nodos, n_grupos)
            # Calcular com nod
            comnod = round(self.com_nod_n(l_grupos_nodos), 5)
            resultado.append(comnod)
            # Almacenar valor máximos absoluto de com nod
            if comnod > comnod_max_abs:
                comnod_max_abs = comnod
                iteracion_max_abs = len(resultado)
            # Almacenar valor máximo relativo de com nod
            delta_comnod = comnod - comnod_max_relativo
            if delta_comnod > 0:
                comnod_max_relativo = comnod
            # Acepta peores opciones con p probabilidad
            elif math.e**(delta_comnod / (t * comnod_max_relativo)) > random.random():
                comnod_max_relativo = comnod
            else:
                # Deshacer cambio
                l_grupos_nodos = self.modificar_agrupacion(l_grupos_nodos, g1, p1, id_x1)
                l_grupos_nodos = self.modificar_agrupacion(l_grupos_nodos, g2, p2, id_y1)
            # enfriar función temperatura
            t *= alpha 
        print("iteracion: ", iteracion_max_abs, " valor: ", comnod_max_abs)

        return resultado, comnod_max_abs, l_grupos_nodos

    def mostrar_nodos(self, l_grupos_nodos)-> None:
        """Muestra nombre de los nodos de la lista"""
        for n, grupo in enumerate(l_grupos_nodos):
            print("grupo: ", n)
            lista2 = []
            for nodo in grupo:
                lista2.append(self.nodos_ID[self.nodos_ID.index(nodo)].info)
            print(" ".join(lista2))

    def graficar(self, resultado: list[list[float]], título: str, nombre: str)-> None:
        """Imprime gráfica de valores com nod vs nº de iteraciones"""
        plt.figure()
        for r in resultado:
            plt.plot(range(len(r)), r)
        plt.xlabel("iteración (nº)")
        plt.ylabel("valor com nod (ud.)")
        plt.title(título)
        plt.savefig(nombre + ".png")
        plt.close()

