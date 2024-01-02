---
tags:
  - command
  - cheatsheet
---


A screen can be initiated with:

```bash
screen -S <name>
```

To detach from a screen, press `Ctrl + A` and then `Ctrl + D`. To scroll in the current screen, press `Ctrl + A`, then `Esc`.

To reattach to a screen, use:

```bash
screen -r <name>
```

To list all screens, use:

```bash
screen -ls
```

Finally, to kill a screen session:

```
screen -S <name> -X quit
```

- `-X`: Execute `<cmd>` as a screen command in the specified session.

Other parameters are:

- `-L`: Save a log (`screenlog.0`) in the current directory.