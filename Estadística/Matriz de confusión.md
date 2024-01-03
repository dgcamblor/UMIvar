---
field: statistics
---

La matriz de confusión es una tabla de contingencia bidimensional que contiene el número de etiquetas predichas y reales, lo que permite visualizar el rendimiento de una prueba o algoritmo.

![[Pasted image 20240103200602.png]]

## Métricas

Al aplicar una prueba, y respecto de la situación real de la enfermedad, nos podemos encontrar los siguientes casos:

- Verdaderos positivos (TP). El número de casos positivos correctamente clasificados como positivos.
- Verdadero negativo (TN). El número de casos negativos correctamente clasificados como negativos.
- Falso positivo (FP) ([[Error tipo I (FP)]]). El número de casos negativos clasificados incorrectamente como positivos.
- Falso negativo (FN) ([[Error tipo II (FN)]]). El número de casos positivos clasificados incorrectamente como negativos.

Las métricas derivadas de la matriz de confusión son:

- Sensibilidad o recall (sensiblidad): La proporción de casos positivos correctamente clasificados como positivos. Es la probabilidad de que un caso positivo se clasifique correctamente como positivo. Se calcula como $TP/(TP+FN)$.
- Especificidad: Es la proporción de casos negativos correctamente clasificados como negativos. Es la probabilidad de que un caso negativo se clasifique correctamente como negativo. Se calcula como $TN/(TN+FP)$.
- Precisión: Proporción de casos positivos clasificados correctamente como positivos. Es la probabilidad de que un caso positivo se clasifique correctamente como positivo. Se calcula como $TP/(TP+FP)$.
## Otras métricas