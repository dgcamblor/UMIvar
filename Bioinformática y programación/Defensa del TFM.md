**¿Por qué evaluaste los programas en los dos tipos de datos diferentes?**

La cobertura (el número de lecturas cubriendo una determinada posición) que se genera en uno y otro tipo de datos es diferente. La cobertura en datos de secuenciación es mayor, lo que implica que podemos detectar variantes a frecuencias menores. Ahora bien, coberturas mayores se relacionan con grupos de duplicados de mayor tamaño, es decir, los UMIs proporcionan más información en este tipo de datos. En cambio, en los datos de secuenciación de exomas, las grupos de duplicados son mucho más reducidos. En concreto, el grupo de datos con el que trabajaba tenía como máximo dos lecturas duplicadas por cada UMI. Quería ver si seguía existiendo beneficio a pesar de tener información limitada de UMIs, que fue lo que vi, especialmente porque en mi unidad trabajamos más con este tipo de datos.

**¿Por qué no probaste los métodos de corrección de k-meros?**

Al ser datos simulados, no confiaría en los resultados de una aprocimación como esa, especialmente al depender tan fundamentalmente de la secuencia.

**¿Por qué los métodos de corrección de k-meros no se podrían utilizar?**

Estos métodos hacen un conteo de k-meros en las lecturas. Se asigna una mayor confianza a los k-meros más frecuentes, mientras que los más infrecuentes son considerados como ruido tecnológico. De nuevo, si tenemosen cuenta la distribución de frecuencias de las variantes verdaderas en biopsias líquidas, la corrección haría que se eliminaran muchas variantes verdaderas.

**¿No se podrían utilizar métodos basados en la calidad de las bases?**

Por supuesto, la calidad de las bases es muy importante para la corrección del ruido tecnológico. Sin embargo, no todo el ruido tecnológico concurre con una reducción de la calidad. Solamente los errores de secuenciación, y no siempre. Por ejemplo, un error de amplificación llevará a un cambio que se leerá por el secuenciador de forma indistinta de cómo leería una base verdadera.

**¿Por qué la estrategia de llamada de variantes sobre lecturas sin deduplicar es más efectiva?**

Al hacer la llamada de variantes sobre lecturas pueden emplearse la información de las moléculas individuales de cara a evaluar qué tanta confianza tenemos en que una variante sea cierta o no. Por ejemplo:
- Las calidades de cada duplicado.
- El número de moléculas que se encuentran apoyando una variante (cuantas más mejor).