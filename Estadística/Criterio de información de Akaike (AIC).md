---
field: statistics
---

El criterio de información de Akaike (AIC) es una medida de la calidad relativa de un modelo estadístico, para un conjunto dado de datos. Como tal, el AIC proporciona un medio para la selección del modelo.

$$
AIC = 2k - 2ln(L)
$$

El AIC se basa en un equilibrio entre la bondad de ajuste del modelo y la complejidad del modelo (el número de parámetros utilizados). Dado un conjunto de modelos, el preferido es aquel con un menor valor de AIC. Recompensa la bondad de ajuste y penaliza el número de parámetros, para desalentar el sobreajuste.

El AIC es útil para la comparativa entre modelos generados de una [[Stepwise regression]].