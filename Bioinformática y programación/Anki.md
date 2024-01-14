
En Anki, el funcionamiento está basado en **Notas** y **Tarjetas**:

- Notas. Es la información que se quiere aprender, y consiste en un conjunto de campos.
- Tarjetas. Las tarjetas son la forma en que la información de las Notas es presentada al usuario. Están formadas por un anverso y un reverso, codificados en HTML, y un estilo, codificado en CSS.

## Tarjetas

Por ejemplo, la tarjeta "Completo" para la nota "Palabra". Hace uso de los campos de reemplazamiento condicional ([Card Generation - Anki Manual](https://docs.ankiweb.net/templates/generation.html?highlight=condit#conditional-replacement)).

Plantilla del anverso:

```html
{{Palabra/s}}
```

Plantilla del reverso:

```html
{{FrontSide}}

<hr id=answer>

{{Español}}

<br><br>

{{#Descripción}}
<div class="cont">
<h3>Descripción</h3>
{{Descripción}} 
</div>
<br>
{{/Descripción}}

{{#URL Diccionario}} 
<div class="cont">
<h3>URL Diccionario</h3>
{{URL Diccionario}} 
</div>
<br>
{{/URL Diccionario}}

{{#Adicional}} 
<div class="cont">
<h3>Adicional</h3>
{{Adicional}} 
</div>
<br>
{{/Adicional}} 
```

Estilo (CSS):

```css
.card {
    font-family: arial;
    font-size: 20px;
    text-align: center;
    color: black;
    background-color: white;
}

div.cont{
		border: 2px solid black;
		border-radius: 25px;
		background-color: rgba(255, 255, 255, 0.1);
		padding: 0px 20px 20px 20px;
}
```