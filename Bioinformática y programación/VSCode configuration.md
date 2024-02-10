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

[Snippets in Visual Studio Code](https://code.visualstudio.com/docs/editor/userdefinedsnippets)

> To create or edit your own snippets, select **Configure User Snippets** under **File** > **Preferences**, and then select the language (by [language identifier](https://code.visualstudio.com/docs/languages/identifiers)) for which the snippets should appear, or the **New Global Snippets file** option if they should appear for all languages. VS Code manages the creation and refreshing of the underlying snippets file(s) for you.

## TODOs

Adding `#!` tracking as a TODO. Add space after `!`.

```
"! ": {
"type": "none"
},
```