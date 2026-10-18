"""
Generador de grafo y cálculo de estadísticos a partir de un archivo que contiene
las relaciones entre los nodos. El archivo se debe incorporar al ejecutar el 
script. Ejemplo: "python script.py /ruta_relativa"
Autor: Adrián Berenguer
Fecha: 12/03/2026
Autor: Adrián Berenguer
"""
from ClaseGrafo_2_1_B import Grafo
import sys
import os
import statistics
from Bio import Entrez

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


def extraer_función(lista: list[str])-> None:
    """Llama a la API de NCBI para extraer la función de una lista de proteínas
    y las imprime por pantalla
    """
    funciones = []
    Entrez.email = "adrian9420agusti@gmail.com" 
    for gen in lista:
        try:
            # Obtener el ID del gen para realizar la búsqueda
            busqueda_nombre = Entrez.esearch(db="gene", term = gen + "[Gene Name]\
                                      AND Homo Sapiens[Organism]", retmax = "1")
            resultado_nombre = Entrez.read(busqueda_nombre)
            busqueda_nombre.close()
            # Repetir la búsqueda una vez obtenido el ID
            if resultado_nombre["IdList"]:
                gen_id = resultado_nombre["IdList"]
                busqueda_id = Entrez.esummary(db="gene", id = gen_id)
                resultado_id = Entrez.read(busqueda_id)
                busqueda_id.close()
                
            else:
                print("No se encontró en NCBI el gen: ", gen)
        
        except:
            print("Error al buscar el gen ", gen)
        
        else:
            funcion = resultado_id['DocumentSummarySet']['DocumentSummary'][0].get('Summary', 'No se encontró resumen.')
            print("Gen: ", gen)
            print("Función: ", funcion, "\n")
            

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
        ok = grafo.grafo_apartir_fichero("Algoritmos/github/2.1.B/biogrid_26.txt")
        if ok:
            print("Grafo generado con éxito...")
            print("Ejecutando algoritmo de floyd...")
            grafo.algoritmo_floyd()
            print("\nMostrando algunas estadísticas...")
            grafo.estadisticas()
            print("Diámetro del grafo: ", grafo.d_max())
            # input()
            #grafo.mostrar_distancias()
            # input()
            grafo.d_promedio()
            # input()
            print("15 nodos con más interacciones directas")
            lista = grafo.ordenar_nodos_adyacentes()[0:5]
            print("función de las 5 proteínas con más interacciones: \n")
            extraer_función(lista)
            # input()
            print("10 nodos mejor conectados a la red")
            grafo.ordenar_d_promedio()
            # input()
            print("\nPrograma finalizado con éxito")


if __name__ == "__main__":
    main()