import sys
import os
import random
longitud = sys.argv[1]
archivo = open(os.curdir, "w")
archivo.write(">Simulado longitud: ", longitud)
ancho_linea = 80
nucleotidos = ("A", "C", "G", "T")
for i in range(longitud // ancho_linea):
    linea = random.choices(nucleotidos, k=ancho_linea)
    archivo.write(linea)
archivo.close()