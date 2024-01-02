---
tags:
  - program
  - cheatsheet
---

## Running a container

### docker run

Start a new container from an image.

```bash
docker run <image>
```

- `--name <name>`: name of the container.
- `-d`: run the container in background.
- `-it`: run the container in interactive mode.

## Manage containers

### docker ps

List running containers.

```bash
docker ps
```

### docker system

Delete everything.

```bash
docker system prune -a --volumes
```

## Image construction

`ENTRYPOINT` -> The executable program that will be the entry point to the image.