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
import itertools

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


def extraer_multifasta(archivo: TextIO, dicc_fastas: dict = {})-> dict[list[int, list[int]]]:
    """A partir de un archivo extrae un fasta o varios con longitud y contaje de tipo de nt"""
    for fasta in SeqIO.parse(archivo, "fasta"):
        contadores = [fasta.seq.count(nt) for nt in ("A", "C", "G", "T")]
        dicc_fastas[fasta.id] = [len(fasta), contadores, fasta.seq]

    return dicc_fastas


def contar_nt(secuencia: str, n: int = 1)-> list[int]:
    """ Cuenta los nucleótidos agrupados en n grupos de una secuencia
    Args:
        n: Longitud de la agrupación (1=mono, 2=di, 3=tri...)
    Return:
        lista de longitud variable con los contadores para n agrupaciones de 
        nucleótidos
    """
    bases = ['A', 'C', 'G', 'T']
    combinaciones = [''.join(p) for p in itertools.product(bases, repeat = n)]
    contadores = {kmero: 0 for kmero in combinaciones}
    for i in range(len(secuencia) - n + 1):
        kmero = secuencia[i:i+n]
        if kmero in contadores:
            contadores[kmero] += 1
    
    return [contadores[kmero] for kmero in combinaciones]


def complex_n(secuencia: str, n: int = 1)-> float:
    """Calcula valor complex N para una lista de contadores de n nucleótidos"""
    contadores = contar_nt(secuencia, n)
    divisor = len(contadores)
    total = sum(contadores)
    complex_n = sum(abs((contador / total) - (1 / divisor)) 
                   for contador in contadores)
    
    return complex_n


def main()-> None:
    print("Iniciando el programa...")
    entrada = extraer_input(1)
    ruta, lista_nombres = obtener_archivos_y_ruta(entrada)
    if lista_nombres:
        dicc_fastas = {}
        for nombre in lista_nombres:
            print(50 * "*")
            print("Abriendo archivo ", nombre)
            print(50 * "*")
            ruta_completa = os.path.join(ruta, nombre)
            archivo_entrada = abrir_archivo(ruta_completa)
            if archivo_entrada:
                print("Extrayendo secuencias...")
                dicc_fastas = extraer_multifasta(archivo_entrada, dicc_fastas)
                cerrar_archivo(archivo_entrada)
        print("Calculando complex...")
        # Calcular varios complex
        for fasta in dicc_fastas:
            print("Calculando complex de ", fasta)
            secuencia = dicc_fastas[fasta][2]
            for c in range(2,3,1):
                complex = complex_n(secuencia, c)
                dicc_fastas[fasta].append(complex)
        # Mostrar los valores
        #print("Genoma \t longitud (nt) \t complex2 \t complex3 \t complex4 \t complex5 \t " \
        #"\complex 6  \t complex7")
        max_l = max(len(fasta) for fasta in dicc_fastas)
        ancho_col = max_l + 2
        print(f"{'Genoma':<{ancho_col}} {'long':>12} {'c2':>8}")
        #print("Genoma\t\t\tlong\tc2")
        print("-" * (ancho_col + 12 + 8 + 2))
        for fasta in sorted(dicc_fastas.keys()):
            longitud = dicc_fastas[fasta][0]
            linea = f"{fasta:<{ancho_col}} {longitud:>12}"
            for i in range(3, len(dicc_fastas[fasta])):
                valor_c = dicc_fastas[fasta][i]
                linea += f" {valor_c:>8.2f}"
            print(linea)  
            
        print("\nPrograma finalizado con éxito")
    else:
        print("Archivos .fa no encontrados en ", ruta)
            


if __name__ == "__main__":
    main()











