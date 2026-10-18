for num in 2 3;
do
echo Procesando para nucleótidos de $num en $num
for fasta in *.fa; 
do
echo Procesando genoma "$fasta";
python script_1.py $fasta $num;
# Se genera matriz de transiciones. 
# Muestra Tablas de frecuencia de dinucleótidos y trinucleótidos. 
# Vector de estado inicial. 
grep -v "^>" Secuencias_26.txt | while read -r secuencia; 
do
if [[ -n $secuencia ]]; 
# Script_2py recibe una secuencia y comprueba la probabilidad de pertenencia
then python script_2.py $secuencia "salida.json";
fi;
done;
done;
python script_3.py comparar.json
done;
