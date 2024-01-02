---
tags:
  - command
  - cheatsheet
---

To tar and compress a directory:

```bash
tar -czvf name-of-archive.tar.gz /path/to/directory-or-file
```

- `-c`: Create an archive.
- `-z`: Compress the archive with gzip.
- `-v`: Display progress in the terminal while creating the archive, also known as “verbose” mode. The v is always optional in these commands, but it’s helpful.
- `-f`: Allows you to specify the filename of the archive.

To extract a tar.gz compressed file:

```bash
tar -xzvf name-of-archive.tar.gz
```

- `-x`: Extract the archive.

## tar all directories in current directory

To tar all the directories in the current directory and remove them:

```
for folder in */; do folder_name=$(basename "$folder"); sudo tar -czvf "${folder_name}.tar.gz" "$folder" && yes | sudo rm -r "$folder"; done
```