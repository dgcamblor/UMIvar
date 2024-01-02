Add mean comparison p-values to a ggplot, such as box blots, dot plots and stripcharts.

```r
stat_compare_means()
```

Use `comparisons = ` to establish the comparisons.

Q-Q plotting the residuals.

```
shapiro.test()
```

El modelo tiene unos residuos. Nos interesa ver si esta desviación tiene los criterios de normalidad. Los residuos son una especie de variable. De eso hay que ver la normalidad, de los residuos y no de la variable dependiente. 

```
fit$residuals
```

Con pocos casos, asumir una distribución normal, aunque te lo diga Shapiro Wilk, es demasiado fuerte. 

>The Shapiro-Wilk tests the Normality of the residuals. The Null Hypothesis is that the residuals are normally distributed. A low P-Value below a given significance level indicates the values are NOT Normally Distributed.

## Kaplan-Meier

0 -> Datos censurados. 

- Censura por la derecha. Antes. Se ha hecho seguimiento y, a partir de ahí, no se sabe qué es lo que pasó.
- Censura por la izquierda. Después.
- Censura de intervalo. El paciente ha muerto entre dos puntos en los que no se sabe con especificidad cuando ha muerto. 

Librerías `survival` y `survminer`.

```r
km1 <- survfit()
```

```
Surv(os.m, os.status) ~ variable.
```

Objeto de supervivencia: primero tiempo y luego estado (EXITUS).

Para recodificar variables 1 a 1:

```
car::recode()
```

```
fit <- survfit(Surv(TimetoOHEmonths, OHE_status) ~ PHESCat
```

```
conf.type="plain"
```

Hay varios tipos de comparación de curvas. 

Curva de supervivencia. El eje X es el tiempo. El eje Y es la supervivencia. Todos empiezan en el 100%.

El test de logrank es el que se utiliza para comparar curvas de supervivencia.

```
ggsurvplot(fit, pval = T, pval.method = T)
```

Para evaluar las diferencias entre pares se utiliza:

```
pairwise_survdiff(fit)
```

Regresión de Cox. Vivo -> OHE

Modelo de Riesgos Competitivos. Vivo -> OHE | Vivo -> Trasplante.

Modelo multiestado. Vivo -> OHE | Vivo -> OHE | OHE -> Trasplante | OHE -> Muerte.

Lo que vemos son hazard ratios.

- < 1 -> Protectora contra el evento.
- > 1 -> Factor de riesgo.
- 1 -> SIn efecto

El `exp(coef)` es el Hazard Ratio. Primero vemos el p-valor y luego vemos si es factor de riesgo o no. Si es p-valor > 0.05, da igual lo que salga de HR.

```text
(HR: 1.583, CI95%: 1.21-2.07, p < 0.001)
```

La función `step()` utiliza el [[Criterio de información de Akaike (AIC)]] para seleccionar el mejor modelo, con las mejores variables.