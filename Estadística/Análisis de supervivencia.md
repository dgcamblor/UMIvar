---
aliases:
  - Survival analysis
---

El análisis de supervivencia requiere de dos elementos básicos:

- **El tiempo hasta el evento.** Desde una fecha inicial, mide el tiempo (por ejemplo, en días) hasta que ocurre un evento o se deja de seguir al paciente sin que haya ocurrido el evento (censura).
- **El estado del paciente.** Si el evento ocurre, el estado del paciente es 1; si no ocurre, el estado del paciente es 0. Esta variable se emplea como indicadora de censura o no.

El hecho de que los datos puedan estar censurados, lo que significa que para algunas observaciones no se conozca el tiempo exacto en el que ocurre el evento, hace que sea necesario introducir particularidades a la metodología de análisis, y que no sean válidos métodos convencionales para predecir tiempos de supervivencia. #insight 

Algunas de las metodologías más frecuentes son:

- [[Método de Kaplan-Meier]]
- [[Test Log-Rank]]
- [[Regresión de Cox]]
- [[Random Survival Forests]]