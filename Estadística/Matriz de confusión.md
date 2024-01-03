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

A partir de estos casos se pueden calcular las siguientes métricas:

- Sensibilidad. Proporción de individuos enfermos que poseen una prueba positiva. Una prueba sensible es más útil cuando su resultado es negativo, ya que raramente deja de detectar un caso positivo.

$$
Sensibilidad = TP/(TP+FN)
$$

- Especificidad. Proporción de individuos sin la enfermedad que poseen una prueba negativa. Una prueba específica es más útil cuando su resultado es positivo, ya que raramente clasifica un caso negativo como positivo.

$$
Especificidad = TN/(TN+FP)
$$

- Valor Predictivo Positivo (VPP) o precisión. Proporción de individuos con una prueba positiva que tienen la enfermedad.

$$
VPP = TP/(TP+FP)
$$

- Valor Predictivo Negativo (VPN). Proporción de individuos con una prueba negativa que no tienen la enfermedad.

$$
VPN = TN/(TN+FN)
$$

## Reglas

- SENDES: cuando una prueba diagnóstica posee una sensibilidad alta (>= 95%), obtener un resultado negativo descarta el diagnóstico.

- ESPIN: cuando una prueba diagnóstica posee una especificidad alta (>= 95%), obtener un resultado positivo indica el diagnóstico.