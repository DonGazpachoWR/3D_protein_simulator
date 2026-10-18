"""
Calculador probabilidad de que una secuencia pertenezca a un genoma
Autor: Adrián Berenguer
Versión: 1.0
Fecha: 19/03/2025
"""
import json
from typing import TextIO
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import os
from pathlib import Path
import numpy as np
import itertools
import math
from numpy.typing import NDArray

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


def extraer_input(n: int = 1)-> str:
    """Identifica un input de la terminal"""
    if len(sys.argv) > n:
        salida = sys.argv[n].strip()
    else:
        salida = ""
    
    return salida

def comprobar_y_abrir_json(n: int = 1)-> tuple[TextIO | None, str | None]:
    """Comprueba si el input de la posición n de sys.argv es un fasta y lo abre"""
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


def multiplicar(lista: list[float])-> float:
    """Multiplica todos los elementos de una lista"""
    resultado = 1
    for f in lista:
        resultado *= f
    
    return resultado

def match_en_genoma(vector: list[float], matriz: NDArray,  seq: str) -> float:
    """Da la probabilidad de encontrar una secuencia en un genoma mediante
    una matriz de transiciones de markov"""
    l_nt = ("A", "C", "G", "T")
    seq = seq.upper()
    # Probabilidad de empezar con el primer nt de la secuecnia
    p = vector[l_nt.index(seq[0])] 
    # Probabilidad de continuar con el resto de nt de la secuencia
    p *= multiplicar([matriz[l_nt.index(seq[x]), l_nt.index(seq[x + 1])] 
                        for x in range(len(seq) - 1) ])
    
    return p


def match_en_genoma_log(vector: list[float], matriz: NDArray, seq: str) -> float:
    """Da el logaritmo de la probabilidad de encontrar una secuencia en un genoma mediante
    una matriz de transiciones de markov"""
    l_nt = ("A", "C", "G", "T")
    seq = seq.upper()
    # Probabilidad de empezar con el primer nt de la secuecnia
    p_log = vector[l_nt.index(seq[0])]
    # Probabilidad de continuar con el resto de nt de la secuencia
    p_log += sum([matriz[l_nt.index(seq[x]), l_nt.index(seq[x + 1])] 
                        for x in range(len(seq) - 1) ])

    return p_log

def abrir_json_comparar(ruta: str) -> dict:
    """Carga el JSON de comparación o crea uno vacío si no existe."""
    if os.path.exists(ruta):
        archivo = abrir_archivo(ruta)
        if archivo:
            try:
                diccionario = json.load(archivo)
            except: 
                diccionario = {}
    cerrar_archivo(archivo)
            
    return diccionario


def main():
    archivo_json, ruta = comprobar_y_abrir_json(2)
    seq = extraer_input(1)
    if archivo_json:
        if seq:
            datos = json.load(archivo_json)
            cerrar_archivo(archivo_json)
            matriz = np.array(datos["transicion"])
            vector = datos["vector"]
            nombre =  datos["nombre"]
            matriz_log = np.array(datos["matriz_log"])
            vector_log = datos["vector_log"]
            p = match_en_genoma(vector, matriz, seq)
            p_log = match_en_genoma_log(vector_log, matriz_log, seq)
            print(70 * "*")
            print("La probabilidad de encontrar la secuencia: ")
            print(seq)  
            print("En el genoma ", nombre) 
            print(f"Es de {p:.2e} y en forma log ", round(p_log, 5))

            # Generar o abrir .json para comparar datos
            ruta_json = os.path.join(ruta, "comparar.json")
            dicc_comparar = abrir_json_comparar(ruta_json)
            # Almacenar dato de probabilidades
            if nombre not in dicc_comparar:
                dicc_comparar[nombre] = {}
            dicc_comparar[nombre][seq] = [p, p_log]
            # Guardar en archivo json
            json_salida = abrir_archivo(ruta_json, "w")
            if json_salida:
                json.dump(dicc_comparar, json_salida, indent = 4)

        else:
            print("No se añadió secuencia")
    else:
        print("No se añadió json ")
        
        
    # el array hay que np.array(json.load(f))


if __name__ == "__main__":
    main()
