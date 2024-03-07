---
aliases:
  - Estimador de Kaplan-Meier
---
Estimador de Kaplan-Meier calcula la probabilidad de supervivencia a cada tiempo observado, permitiendo la construcción de curvas de supervivencia.

## Sobre la obtención de un p-valor con K-M

El método de Kaplan-Meier por sí solo no produce un p-valor directamente. Sin embargo, se utiliza en conjunto con otras pruebas estadísticas para comparar las curvas de supervivencia entre grupos. Estas pruebas estadísticas sí proporcionan un p-valor que indica la significancia estadística de las diferencias observadas entre las curvas de supervivencia. Para comparar la distribución de supervivencia entre dos grupos, se puede calcular un test estadístico como:

- [[Test Log-Rank]]
- [[Regresión de Cox]]