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

Some key parameters are:

- `-s`: To create a standalone document. This is useful, for example, when converting a Markdown file to HTML or LaTeX, as it ensures that the output file is a complete, self-standing document, including any necessary headers and footers.

- `--citeproc`: To add citations. This requires a bibliography file in BibTeX format, and the Markdown file must be written in a format that supports citations, such as Pandoc’s Markdown.

- `--bibliography`: To add a bibliography file.

- `--lua-filter`: To add a Lua filter. This is a Lua script that modifies the document before conversion.

- `--metadata`: To add metadata to the document. This is a key-value pair, such as `--metadata title="My Document"`.