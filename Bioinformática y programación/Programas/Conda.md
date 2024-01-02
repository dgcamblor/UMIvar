---
tags:
  - program
---

Conda can create isolated environments in which install specific packages and dependencies, independently of each other.

## Installation

The easiest way to install conda is using miniconda.

```bash
# Install
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

## Configuration

Adding channels:

```bash
# Add channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels anaconda
conda config --add channels r
```

To avoid the base conda environment from activating at startup:

```bash
conda config --set auto_activate_base false
```

## Creation

```
conda create --name myenv python=3.7 numpy pandas
conda activate myenv
```

## Installing packages

```
conda install matplotlib
```

Packages are installed under the `~/miniconda3/envs/env_name` directory.


## Using pip

When using `pip` to install packages, packages are installed to the corresponding environment.