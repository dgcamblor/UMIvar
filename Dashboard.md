---
cssclasses:
  - dashboard
banner: "![[guillermo-ferla-kEEl9csCutg-unsplash.jpg]]"
banner_x: 0.5
banner_y: 0.54
obsidianUIMode: preview
---

## Dashboards

- 🖱️Bioinformática y programación
	- [[Pipeline list]]
- 🎓Doctorado
	- [[PhD dashboard]]

## Por escribir

```dataview
TABLE 
from #write 
```

## Información de la bóveda

- 🗄️ Recent file updates
 `$=dv.list(dv.pages('').sort(f=>f.file.mtime.ts,"desc").limit(4).file.link)`
- 🔖 Tagged:  favorite 
 `$=dv.list(dv.pages('#favorite').sort(f=>f.file.name,"desc").limit(4).file.link)`
- 〽️ Stats
	-  Número de notas: `$=dv.pages().length`
	-  Número de notas diarias: `$=dv.pages('"Notas diarias"').length`