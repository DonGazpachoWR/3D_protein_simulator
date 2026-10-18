source /home/adrian/aplicaciones/anaconda3/etc/profile.d/conda.sh
conda activate redes
for file in *.fa;  do echo "Procesando $file"; time python barajador.py "$file" 1; done