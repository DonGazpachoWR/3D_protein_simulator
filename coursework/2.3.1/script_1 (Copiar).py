"""
Obtiene frecuencia de agrupaciones de n nucleótidos a partir de un fasta

"""
from typing import TextIO
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import os
from pathlib import Path
import json
import numpy as np
import itertools
import math


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
        if p.suffix.lower() == ".fa" or p.suffix.lower() == ".fasta":
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
        salida = sys.argv[n]
    else:
        salida = ""
    
    return salida


def comprobar_y_abrir_fasta(n: int = 1)-> tuple[TextIO | None, str | None]:
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


def contar_nt(secuencia: str, n: int = 1, dicc_frecuencias: dict = {})-> None:
    """ Cuenta los nucleótidos agrupados en n grupos de una secuencia
    Args:
        n: Longitud de la agrupación (1=mono, 2=di, 3=tri...)
    Return:
        lista de longitud variable con los contadores para n agrupaciones de 
        nucleótidos
    """
    bases = ['A', 'C', 'G', 'T']
    combinaciones = [''.join(p) for p in itertools.product(bases, repeat = n)]
    # Combinacion de nucleótidos A,C,G,T
    for kmero in combinaciones:
        if kmero not in dicc_frecuencias:
            dicc_frecuencias[kmero] = [0]
    for i in range(len(secuencia) - n + 1):
        kmero = secuencia[i:i+n]
        if kmero in dicc_frecuencias:
            dicc_frecuencias[kmero][0] += 1
        # Por si hay N, Y u otra nomenclatura
        else:
            dicc_frecuencias[kmero] = [1]


def extraer_fasta(archivo: TextIO)-> tuple[dict[list[int, list[int]]], str]:
    """A partir de un archivo extrae un fasta con longitud y contaje de tipo de nt"""
    dicc_fasta = {}
    fasta = SeqIO.read(archivo, "fasta")
    dicc_fasta[fasta.description] = [len(fasta), str(fasta.seq)]

    return dicc_fasta, fasta.description


def frec_relativas(diccionario: dict, l = int)-> None:
    """Calcula frecuencia relativa de nucleótidos a partir de frec. absoluta"""
    for clave in diccionario:
        diccionario[clave].append(round(diccionario[clave][0] / l, 4))


def mostrar_v_estacionario(diccionario: int,)-> tuple[float]:
    """Muestra vector de mono nucleótidos"""
    print("\nVector estacionario del sistema")
    vector = tuple([diccionario[clave][1] for clave in diccionario])
    print(vector)

    return vector


def main():
    archivo_entrada, ruta = comprobar_y_abrir_fasta()
    if archivo_entrada:
        genoma, nombre = extraer_fasta(archivo_entrada)
        cerrar_archivo(archivo_entrada)
        print("Genoma de ", nombre, "\n")
        for n in range(1,4):
            genoma[nombre].append({})
            contar_nt(genoma[nombre][1], n, genoma[nombre][1 + n])
            frec_relativas(genoma[nombre][1 + n], genoma[nombre][0])
            print("\nTabla de frecuencia para agrupación de nucleótidos de ", n, " en ", n, "\n")
            print(f"{'Nt':<8} {'F. Absoluta':>15} {'F. Relativa':>15}")
            for nt_n, valores in genoma[nombre][1 + n].items():
                print(f"{nt_n:<8} {valores[0]:>15,d} {valores[1]:>15.4f}")
            if n == 2:
                contador =  tuple([genoma[nombre][2][clave][0] for clave in genoma[nombre][2]])
                # Matriz de transiciones
                matriz_nt = []
                matriz = []
                k = int(math.sqrt(len(genoma[nombre][1 + n].items())))
                matriz = np.zeros((k,k))
                for x, nt_n in enumerate(genoma[nombre][1 + n]):
                    valor = genoma[nombre][1 + n][nt_n][0]
                    fila = x // k
                    col = x % k
                    valor_matricial = round(valor / contador[fila],4)
                    matriz[fila, col] = valor_matricial
                # Matriz de transiciones logarítmica
                matriz_log = [[math.log(matriz[i, j]) 
                              for j in range(k)] for i in range(k)]
                print("\nMatriz de transiciones para agrupación de nucleótidos de ", n, " en ", n, "\n")
                bases = ['A', 'C', 'G', 'T']
                linea = "      " + "       ".join(bases)
                print(linea)
                for i, fila in enumerate(matriz):
                    valores_fila = "  ".join([f"{val:.4f}" for val in fila])
                    print(f"{bases[i]}  [{valores_fila}]")
                # Imprimir vector estacionario
                vector = mostrar_v_estacionario(genoma[nombre][2])
                vector_log = [math.log(i) for i in vector]
                datos_a_guardar = {
                # El np.array hay que mi_array_numpy.tolist() al guardar en .json
                "transicion": matriz.tolist(),
                "vector": vector,
                "nombre": nombre,
                "matriz_log": matriz_log,
                "vector_log": vector_log  }
        # Guardar datos
        archivo_json = abrir_archivo(os.path.join(ruta, "salida.json"), modo = "w")
        if archivo_json:
            json.dump(datos_a_guardar, archivo_json, indent=4)
            cerrar_archivo(archivo_json)

    else:
        print("No incluyó el fasta al ejecutar el programa")
        




    


if __name__ == "__main__":
    main()