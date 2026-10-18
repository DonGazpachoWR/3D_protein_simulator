import json
import os
import itertools
from typing import TextIO
import sys
from pathlib import Path

def abrir_archivo(ruta, modo: str = "r")-> TextIO | None:
    """intenta abrir un archivo en modo lectura o escritura"""
    try:
        archivo = open(ruta, modo)
    
    except:
        print("Error al abrir el archivo ", ruta)
        archivo = None
    
    return archivo

def cerrar_archivo(archivo: TextIO)-> None:
    archivo.close()

def comparador(datos: dict, genoma1: str, genoma2: str, seq: str) -> list[float]:
    """Calcula el cociente de P y la diferencia de ln(P) entre dos genomas."""
    val1 = datos[genoma1][seq]
    val2 = datos[genoma2][seq]
    resultados = [0.0, 0.0]
    
    # p1/p2
    if val2[0] != 0:
        resultados[0] = val1[0] / val2[0]
    else:
        resultados[0] = float('inf')
        
    # ln(p1) - ln(p2)
    resultados[1] = val1[1] - val2[1]
    
    return resultados


def obtener_archivo_y_ruta(entrada: str) -> tuple[str, list[str]]:
    """
    Determina si la entrada es un archivo . o un directorio.
    Devuelve (ruta_directorio, archivo).
    """
    p = Path(entrada) if entrada else Path.cwd()
    salida = []
    # Caso 1: Es un archivo individual
    if p.is_file():
        if p.suffix.lower() == ".json":
            salida = [str(p.parent), [p.name]]
        else:
            print("El archivo ",p.name," no tiene extensión fasta")
            salida =  [str(p.parent), []]
    # Caso 2: No existe
    else:
        salida = [str(p), []]
    
    return salida


def comprobar_y_abrir_json(n: int = 1)-> tuple[TextIO | None, str | None]:
    """Comprueba si el input de la posición n de sys.argv es un json y lo abre"""
    archivo_entrada = None
    ruta = ""
    entrada = extraer_input(n)
    if entrada:
        ruta, nombre_archivo = obtener_archivo_y_ruta(entrada)
        if nombre_archivo and ruta:
            if ruta == ".":
                archivo_entrada = abrir_archivo(nombre_archivo[0])
            else:
                archivo_entrada = abrir_archivo(ruta + nombre_archivo[0])

    return archivo_entrada, ruta


def extraer_input(n: int = 1)-> str:
    """Identifica un input de la terminal"""
    if len(sys.argv) > n:
        salida = sys.argv[n].strip()
    else:
        salida = ""
    
    return salida

def main():
    archivo_json, ruta = comprobar_y_abrir_json(1)
    if archivo_json:
        diccionario = json.load(archivo_json)
        cerrar_archivo(archivo_json)
        secuencias = [[sec for sec in diccionario[d]] for d in diccionario][0]
        if len(secuencias) >= 2:
            l_genomas = list(diccionario.keys())
            for x in range(len(l_genomas)):
                for y in range(len(l_genomas)):
                    if x < y:
                        genoma1 = l_genomas[x]
                        genoma2 = l_genomas[y]
                        print("\nComparación de genoma", genoma1[0:20], " con ", genoma2[0:20], "\n")
                        for secuencia in secuencias:
                            p1 = diccionario[genoma1][secuencia][0]
                            p2 = diccionario[genoma2][secuencia][0]
                            p_log_1 = diccionario[genoma1][secuencia][1]
                            p_log_2 = diccionario[genoma2][secuencia][1]
                            print("Secuencia: ", secuencia)
                            print("Ratio numérico: ", round(p1 / p2, 4), 
                                "ratio log: ", round(p_log_1 - p_log_2, 4))
        else:
            print("Se necesitan al menos 2 genomas en el json para comparar")
    else:
        print("No se encontró el archivo comparar.json")
    


if __name__ == "__main__":
    main()