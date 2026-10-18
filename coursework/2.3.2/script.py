
"""
Calculador de probabilidades de que unas observaciones sigan una secuencia
determinada de estados mediante modelos ocultos de markov, que intenta
inferir la secuencia de estados.
Autor: Adrián Berenguer
Fecha 31/03/2026
Versión: 1.0
"""
import numpy as np
from numpy.typing import NDArray


def viterbi(l_estados: list, prob_inicial: list, 
            m_trans: NDArray, m_salida: NDArray, observaciones: list,
            p_final: list) -> tuple[list, float]:
    """Utiliza el algoritmo de viterbi para obtener la secuencia de estados
    transitados y la probabilidad de que haya sucedido esa sec. de estados
    Args:
        l_estados: lista de estados ocultos posibles
        prob_inicial: lista de probabilidades de comenzar en cada estado
        m_trans: matriz de probabilidades para mantener el mismo estado o cambiar
        m_salida: probabilidad de que una observación pertenezca a cada estado
        observaciones: lista de observaciones cuya sec. de estados se infiere
        p_final: lista de probabilidades de salir en cada estado
    Return:
        lista de indices de los estados que con mayor probabilidad se ha recorrido
        probabilidad de que haya sucedido esa secuencia de estados
    """
    n_estados = len(l_estados)
    n_obs = len(observaciones)
    # Matriz de probabilidad
    prob = np.zeros((n_obs, n_estados))
    # Matriz que guarda posición del estado anterior más probable
    prev = np.zeros((n_obs, n_estados), dtype = int)
    # Probabilidad de comenzar en cada estado, para este caso 0 menos para el
    # primer estado
    for j in range(n_estados):
        prob[0][j] = prob_inicial[j] * m_salida[j][observaciones[0]]

    # Recorrer a partir de la segunda observacion
    for i in range(1, n_obs):
        # s es el estado actual que alcanzar
        for j in range(n_estados):
            # r es el posible estado anterior. Ver cuál es el mejor paso anterior
            for r in range(n_estados):
                nueva_prob = prob[i - 1][r] * m_trans[r][j] * m_salida[j][observaciones[i]]
                # Cuando es nueva se actualiza el valor y cuál es el estado anterior
                if nueva_prob > prob[i][j]:
                    prob[i][j] = nueva_prob
                    prev[i][j] = r
    # Añadir la probabilidad de salir del sistema
    for s in range(n_estados):
        prob[n_obs - 1][s] = prob[n_obs - 1][s] * p_final[s]
    # Buscamos el índice del estado previo al final con la probabilidad más alta
    max_p_final = -1
    for s in range(n_estados):
        if prob[n_obs - 1][s] > max_p_final:
            max_p_final = prob[n_obs-1][s]
            i_estado_final = s
    ruta = [0] * n_obs       
    ruta[n_obs - 1] = i_estado_final
    # Deshacer camino
    for t in range(n_obs - 2, -1, -1):
        ruta[t] = prev[t + 1][ruta[t + 1]]

    return ruta, max_p_final, prob


def mostrar_m_trans(matriz: list[list], n: int)-> None:
    """Imprime por pantalla la matriz de transiciones"""
    print("\nMatriz de transiciones\n")
    bases = ['Q1', 'Q2', 'Q3', 'Q4']
    nombres_filas = ['C1', 'C2', 'C4', 'C5', 'C6', 'C7']
    linea = " " * n + "          " + "           ".join(bases)
    print(linea)
    for i, fila in enumerate(matriz):
        valores_fila = "  ".join([f"{val:.9f}" for val in fila])
        print(f"{nombres_filas[i]:<{n}}  [{valores_fila}]")


def main():
    # Todos los estados ocultos posibles
    l_estados = ["Q1", "Q2", "Q3", "Q4"] 
    # Variables observadas
    observaciones = [0, 1, 3, 4, 5, 6]
    # Es una secuencia que empieza siempre en el estado Q1
    prob_inicial = np.array([1, 0, 0, 0])
    # Es una secuencia que al salir de Q4 se multplica por 0.6
    p_final = [0.0, 0.0, 0.0, 0.6]
    # probabilidades de mantenerse en el mismo estado o de pasar al siguiente
    matriz_transiciones = np.array([
                        [0.4, 0.6, 0.0, 0.0], 
                        [0.0, 0.8, 0.2, 0.0], 
                        [0.0, 0.0, 0.3, 0.7], 
                        [0.0, 0.0, 0.0, 0.4]])
    # probabilidad de observar cada variables en cada estado
    matriz_salidas = np.array([
                            [0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0],  
                            [0.0, 0.1, 0.2, 0.5, 0.1, 0.1, 0.0],
                            [0.0, 0.0, 0.2, 0.2, 0.4, 0.2, 0.0],
                            [0.0, 0.0, 0.0, 0.1, 0.1, 0.5, 0.3]]) 
    
    posiciones, p_final, prob = viterbi(l_estados, prob_inicial, matriz_transiciones, 
                                  matriz_salidas, observaciones, p_final)
    solución = [l_estados[i] for i in posiciones]
    print("La secuencia de estados más probable es: ", " -> ".join(solución), 
          " con probabilidad ", p_final)
    mostrar_m_trans(prob, 1)


if __name__ == "__main__":
    main()
