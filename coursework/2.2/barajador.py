"""
Barajador de nucleótidos
Autor: Adrián Berenguer Agustí
Versión: 1.0
Fecha 18/03/2026
"""
from typing import TextIO
import glob
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import os
from pathlib import Path

def abrir_archivo(ruta, modo: str = "r")-> TextIO | None:
    """intenta abrir un archivo en modo lectura o escritura"""
    try:
        archivo = open(ruta, modo)
    
    except:
        print("Error al abrir el archivo")
        archivo = None
    
    return archivo


def cerrar_archivo(archivo: TextIO)-> None:
    archivo.close()


def obtener_archivos_y_ruta(entrada: str) -> tuple[str, list[str]]:
    """
    Determina si la entrada es un archivo .fa o un directorio.
    Devuelve (ruta_directorio, lista_de_archivos).
    """
    p = Path(entrada) if entrada else Path.cwd()
    salida = []
    # Caso 1: Es un archivo individual
    if p.is_file():
        if p.suffix.lower() == ".fa" or p.suffix.lower() == ".fasta":
            return str(p.parent), [p.name]
        else:
            print("El archivo ",p.name," no tiene extensión fasta")
            salida =  [str(p.parent), []]
            
    # Caso 2: Es un directorio
    elif p.is_dir():
        archivos = [f.name for f in p.glob("*.fa" or ".fasta")]
        salida = [str(p), archivos]
    
    # Caso 3: No existe
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


def extraer_multifasta(archivo: TextIO)-> dict[list[int, list[int]]]:
    """A partir de un archivo extrae un fasta o varios con longitud y contaje de tipo de nt"""
    dicc_fastas = {}
    for fasta in SeqIO.parse(archivo, "fasta"):
        contadores = [fasta.seq.count(nt) for nt in ("A", "C", "G", "T")]
        dicc_fastas[fasta.id + fasta.description] = [len(fasta), contadores, fasta.seq]

    return dicc_fastas


def barajar(longitud: int, A: int, C: int, G: int, T: int)-> str:
    """Genera secuencia aleatoria de nucleotidos con diferente porcentaje
    de cada nucleótido"""
    N = longitud - sum((A,C,G,T))
    # Opciones que se pueden escoger
    letras = ("A", "C", "G", "T", "N")
    # Índice que será escogido por random.choices()
    indices = [0, 1, 2, 3, 4]
    # Peso ponderado para cada opción
    p = [A, C, G, T, N]
    secuencia = []
    for n in range(longitud):
        i = random.choices(indices, p)[0]
        secuencia.append(letras[i])
        p[i] -= 1

    return "".join(secuencia)


def barajar_Fisher_Yates(secuencia: str)-> str:
    """Algoritmo Fisher Yates para el barajado"""
    secuencia = [ nt for nt in secuencia]
    longitud = len(secuencia)
    for i in range(0, longitud - 1):
        j = random.choice(range(i, longitud))
        secuencia[i], secuencia[j] = secuencia[j], secuencia[i]
        
    return "".join(secuencia)

def main()-> None:
    print("Iniciando el programa...")
    print("Buscando archivos...")
    entrada = extraer_input(1)
    opcion = extraer_input(2)
    ruta, lista_nombres = obtener_archivos_y_ruta(entrada)
    if lista_nombres:
        print("Archivos encontrados!")
        for nombre in lista_nombres:
            print(50 * "*")
            print("Abriendo archivo ", nombre)
            print(50 * "*")
            ruta_completa = os.path.join(ruta, nombre)
            archivo_entrada = abrir_archivo(ruta_completa)
            if archivo_entrada:
                print("Se logró abrir el archivo ", nombre)
                print("Extrayendo secuencias...")
                dicc_fastas = extraer_multifasta(archivo_entrada)
                cerrar_archivo(archivo_entrada)
                print("Generando archivo de salida...")
                nombre_salida = "Barajado_" + nombre
                ruta_salida = os.path.join(ruta, nombre_salida)
                archivo_salida = abrir_archivo(ruta_salida + nombre, "a")
                if archivo_salida:
                    print("Se logró generar el archivo salida")
                # Generar secuencias aleatorizadas
                    for fasta in dicc_fastas:
                        print(50 * "*")
                        print("Analizando el fasta ", fasta, "del archivo ", nombre)
                        print(50 * "*")
                        longitud = dicc_fastas[fasta][0]
                        A, C, G, T = dicc_fastas[fasta][1]
                        print("tiene una longitud de ", longitud)
                        print("tiene los siguientes nt: A =",A," C=",C, " G=", G, " T=",T)
                        print("Barajando secuencia...")
                        if opcion:
                            if opcion == "1":
                                secuencia = barajar(longitud, A, C, G, T)
                        else:
                            secuencia = barajar_Fisher_Yates(dicc_fastas[fasta][2])
                        contadores = [secuencia.count(nt) for nt in ("A", "C", "G", "T")]
                        A, C, G, T = contadores
                        cabecera = fasta + "_barajado"
                        print("Se barajó ", cabecera, "")
                        print("tiene una longitud de ", longitud)
                        print("tiene los siguientes nt: A =",A," C=",C, " G=", G, " T=",T)
                        # Guardar fichero de secuencias aleatorizadas usando Biopython
                        print("Escribiendo en archivo salida...")
                        secuencia = Seq(secuencia)
                        a_escribir = SeqRecord(secuencia, cabecera)
                        SeqIO.write(a_escribir, archivo_salida, "fasta")
                        print("Archivo escrito correctamente")
                    cerrar_archivo(archivo_salida)
                    print("\nPrograma finalizado con éxito")
    else:
        print("Archivos .fa no encontrados en ", ruta)
            


if __name__ == "__main__":
    main()
