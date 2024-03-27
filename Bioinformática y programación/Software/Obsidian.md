[Documentation](https://help.obsidian.md/Home)
[CSS variables](https://docs.obsidian.md/Reference/CSS+variables/CSS+variables)
[Developer documentation](https://help.obsidian.md/Contributing+to+Obsidian/Developers)

## Tags

**YAML tags vs. normal tags.** YAML tags should be general page metainformation while normal tags aim for paragraph searchability. 

## Internal links

```text
[[<Page>#^<paragraph>]]
```

You can embed a reference with:
```text
![[<Page>]]
```

## External links

```
[<text>](<link>)
```

Opening local apps inside Obsidian:
```
notion://<link>
```

## Callouts

https://help.obsidian.md/Editing+and+formatting/Callouts#Supported+types
Callouts are written with the syntax:

```
> [!<Callout>]
```

The callout types are:

- Note
- Abstract
- Info
- Todo
- Tip
- Success
- Warning
- Failure
- Danger
- Bug
- Example
- Quote

## Synced block

Synced blocks are not implemented in the sense of two-way edition. However, blocks can be embedded adding a `!` before the page mention.

```
![[^]]
```
