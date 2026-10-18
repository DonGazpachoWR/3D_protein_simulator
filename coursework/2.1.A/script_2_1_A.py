"""
Generador de grafo y cálculo de estadísticos a partir de un archivo que contiene
las relaciones entre los nodos. El archivo se debe incorporar al ejecutar el 
script. Ejemplo: "python script.py /ruta_relativa"
Autor: Adrián Berenguer
Fecha: 12/03/2026
Autor: Adrián Berenguer
"""
from ClaseGrafo_2_1_A import Grafo
import sys
import os
import statistics

def extraer_directorio(n: int = 1)-> str|None:
    """extrae el directorio con sys.argv con n como parámetro"""
    if len(sys.argv) > n:
        directorio = sys.argv[n]
    else:
        directorio = os.curdir
    
    return directorio


def comparar_algoritmos(alg1:list[float], alg2:list[float] )-> None:
    """Compara com nod maximos de dos tipos de algoritmos"""
    print("com nod medio: ", statistics.mean(alg1), statistics.mean(alg2))
    print("desviación std: ", statistics.stdev(alg1), statistics.stdev(alg2))


def main():
    print(50*"*")
    print("Generador de grafos y estadísticas")
    print(50*"*")
    #input("\nPulse cualquier tecla para iniciar el programa")
    print(50*"*")
    print("Buscando archivo...")
    ruta = extraer_directorio()
    if ruta:
        print("Archivo encontrado con éxito...")
        grafo = Grafo()
        print("Generando grafo...")
        ok = grafo.grafo_apartir_fichero("Algoritmos/github/cubos_26.txt")
        if ok:
            print("Grafo generado con éxito...")
            print("Ejecutando algoritmo de floyd...")
            grafo.algoritmo_floyd()
            print("\nMostrando algunas estadísticas...")
            grafo.estadisticas()
            print("Diámetro del grafo: ", grafo.d_max())
            # input()
            grafo.mostrar_distancias()
            # input()
            grafo.d_promedio()
            # input()
            print("Numero de interacciones directas para cada nodo")
            grafo.ordenar_nodos_adyacentes()
            # input()
            grafo.ordenar_d_promedio()
            # input()
            # Hill climbing
            print("Algoritmo hill climbing para clasificar...")
            algortimo = []
            for n in range(10):
                print(n)
                algortimo.append(grafo.hill_climbing(13,9))
            ordenado = max(algortimo, key = lambda x: x[1])
            grafo.mostrar_nodos(ordenado[2])
            resultado = [ordenado[0]]
            grafo.graficar(resultado, "Algoritmo Hill Climbing", "hill_climbing")
            # Simulated annealing
            print("Algoritmo simulated annealing")
            algortimo = []
            algortimo.append(grafo.sim_annealing(13, 9))
            ordenado = max(algortimo, key = lambda x: x[1])
            grafo.mostrar_nodos(ordenado[2])
            resultado = [ordenado[0]]
            grafo.graficar(resultado, "Algoritmo Simulated Annealing", "SA")

            print("\Programa finalizado con éxito")


if __name__ == "__main__":
    main()