---
tags:
  - command
  - cheatsheet
---

## df

```
df -h .
```

## du

20 heaviest files:

```
sudo du -a / 2>/dev/null | sort -n -r | head -n 20
```

