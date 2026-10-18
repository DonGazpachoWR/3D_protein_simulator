"""
python script.py directory/bacteria_genome.fa directory with target fasta
"""
import math
import sys
from Bio import SeqIO
from typing import TextIO
from pathlib import Path
import os
import itertools
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
from sklearn.preprocessing import StandardScaler

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


def extraer_fasta(archivo: TextIO)-> list[str,str]:
    """A partir de un archivo extrae un fasta único"""

    fasta = SeqIO.read(archivo, "fasta")

    return fasta.seq, fasta.description.split('|')[-1].split(',')[0].strip()


def obtener_archivo_y_ruta(entrada: str) -> tuple[str, list[str]]:
    """
    Determina si la entrada es un archivo . o un directorio.
    Devuelve (ruta_directorio, archivo).
    """
    p = Path(entrada) if entrada else Path.cwd()
    salida = []
    # Caso 1: Es un archivo individual
    if p.is_file():
        if p.suffix.lower() in (".fa",".fasta", ".fna") :
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
                archivo_entrada = abrir_archivo(os.path.join(ruta, nombre_archivo[0]))

    return archivo_entrada, ruta


def listar_archivos_ruta(n: int = 2)-> list[str, str]:
    """genera lista con los archivos fasta contenidos en la carpeta"""
    ruta = extraer_input(n)
    if ruta:
        lista = os.listdir(ruta)
    else:
        lista = ""
    
    return lista, ruta


def frecuencia_relativa(seq: str, n: int = 3)-> dict:
    """Genera diccionario con la frecuencia relativa de kmeros de 3 en 3"""

    bases = ("A","C","G","T")
    combinaciones = [''.join(p) for p in itertools.product(bases, repeat = n)]
    # Contar nº de contigs en la secuencia
    d_freq = { tri_nt: seq.count(tri_nt)
        for tri_nt in combinaciones}
    # Nº total de contigs
    contigs = sum(d_freq.values())
    # Frecuencias relativas de cada contig
    d_freq = {tri_nt: value / contigs
        for tri_nt, value in d_freq.items()}
    
    return d_freq


def euclidean_d(dict_virus: dict, dict_bac: dict):
    """Calcula distancia euclidea entre dos diccionarios de frecuencias relativas"""
    # Distancia euclídea
    eu_distance = math.sqrt(sum(abs(x - y)**2 for x,y in zip(dict_virus.values(), dict_bac.values())))

    return eu_distance

def graficar_k_means(x: list[float], y: list[float], nombres: list[str], 
                     titulo: str = "Grafica.png", linea: bool = False):
    """Grafica scatter plot más algoritmo k means"""
    array = np.array(list(zip(x, y)))
    # Escalar los datos porque la diferencia de magnitud entre el eje X e Y es muy grande
    scaler = StandardScaler()
    array = scaler.fit_transform(array)
    # Entrenar k means
    kmeans = KMeans(n_clusters=2, random_state=0, n_init=10)
    agrupacion = kmeans.fit_predict(array)
    # Parámetros gráfica
    plt.figure(figsize=(10, 6)) 
    grafico = plt.scatter(x, y, c=agrupacion)
    # Añadir nombres a los puntos
    for i, nombre in enumerate(nombres):
        plt.annotate(nombre, (x[i] + 0.0005, y[i]), fontsize = 6)

    plt.xlabel("k4_freq")
    plt.ylabel("Phage size (nt)")
    if linea:
        plt.axvline(x = 0.026, color="r", label="thresold", linestyle = ":")
    plt.savefig("k4.png")
    plt.show()
    plt.close()


def main():

    archivo_bac, ruta_bacteria = comprobar_y_abrir_fasta(1)
    lista_rutas, ruta_virus = listar_archivos_ruta()

    if archivo_bac and lista_rutas:

        resultado = {}
        seq_bac, id_bac = extraer_fasta(archivo_bac)
        cerrar_archivo(archivo_bac)

        print("\nBacteria ", id_bac, " encontrada!")
        # Frecuencia relativa para k-meros de 3 y 4
        dic_fq_coli_3 = frecuencia_relativa(seq_bac, 3)
        dic_fq_coli_4 = frecuencia_relativa(seq_bac, 4)

        for ruta in lista_rutas:

            if ruta.endswith(".fa") or ruta.endswith(".fasta"):
                archivo_virus = abrir_archivo(os.path.join(ruta_virus, ruta))

                if archivo_virus:

                    seq_virus, id_virus = extraer_fasta(archivo_virus)
                    print("\nVirus ", id_virus, " encontrado!")
                    cerrar_archivo(archivo_virus)

                    for x,dic_fq_coli in zip((3,4), (dic_fq_coli_3, dic_fq_coli_4)):
                        dic_fq_virus = frecuencia_relativa(seq_virus, x)
                        d_euclidea = euclidean_d(dic_fq_virus, dic_fq_coli)
                        if id_virus not in resultado:
                            resultado[id_virus] = [[x, len(seq_virus), d_euclidea]]
                        else:
                            resultado[id_virus].append([x, len(seq_virus), d_euclidea])

        print("\nDiccionario\n")
        
        for salida in resultado.items():
            print("\n", salida)

        # Extraer datos del diccionario para graficar

        k_3 = [ lista[0][2] for lista in resultado.values()]
        y = [lista[0][1] for lista in resultado.values()]
        k_4 = [ lista[1][2] for lista in resultado.values()]
        nombres = [nombre for nombre in resultado.keys()]

        graficar_k_means(k_3, y, nombres, titulo = "k3.png")
        graficar_k_means(k_4, y, nombres, titulo = "k4.png", linea = True)



if __name__ == "__main__":
    main()
                    



# conda activate redes            
# python Algoritmos/github/3.1.C/script.py Algoritmos/github/3.1.C/e_coli.fna Algoritmos/github/3.1.C


