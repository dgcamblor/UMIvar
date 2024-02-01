---
tags:
  - language
  - cheatsheet
---
## Global configuration of the script

Adding the shebang `#!/bin/bash` instead of `#!/bin/sh`  allows to use bash, which in fact is a superset of sh with extended features.

```bash
#!/bin/bash
```

Use of aliases inside a script:

```bash
shopt -s expand_aliases
```

Sourcing a configuration file:

```bash
source config.conf
```

Exiting the script when any command fails:

```bash
set -e
```

Activate a [[Conda]] environment for the script:

```bash
eval "$(conda shell.bash hook)"; conda activate ENV
```

## Variable assignment

To assign the result of a Linux command to a variable, you can use command substitution or the `$(...)` syntax. 

```bash
variable=$(command)
```

## Block comment

Block comments are not *per se* implemented in bash, but workarounds can be used.

```bash
: <<'END'
bla bla
blurfl
END
```

