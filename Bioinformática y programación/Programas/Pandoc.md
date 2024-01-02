---
tags:
  - program
---


Pandoc is a document conversor, capable of converting several document types like Markdown, HTML, LaTeX and Word.

Install with:

```bash
sudo apt install pandoc
```

Converting from Markdown to LaTeX:

```bash
pandoc document.md -f markdown -t latex -s -o document.tex
```