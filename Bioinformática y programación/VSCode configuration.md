---
tags:
  - cheatsheet
---
## Shortcuts

CTRL + SHIFT + P -> Open Keyboard Shortcuts

```json
// Place your key bindings in this file to override the defaults
[
{
"key": "ctrl+alt+i",
"command": "editor.action.insertSnippet",
"args": {
"snippet": "```{r}\n$0\n```"
},
"when": "editorLangId == markdown || editorLangId == rmd"
},
{ // Insert %>% in R
"key": "ctrl+shift+m",
"command": "editor.action.insertSnippet",
"args": {
"snippet": "%>%"
},
"when": "editorLangId == rmd || editorLangId == r"
},
]
```

## Language specific settings

CTRL + SHIFT + P -> Language specific settings

## Snippets

## TODOs

Adding `#!` tracking as a TODO. Add space after `!`.

```
"! ": {
"type": "none"
},
```