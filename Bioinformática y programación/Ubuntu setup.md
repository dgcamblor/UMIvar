## nvidia drivers

```bash
sudo apt install nvidia-driver-525
```

```bash
sudo snap install code --classic
sudo 
```

## Dual boot problems

### Customizing grub (changing boot order, etc.)

```bash
sudo add-apt-repository ppa:danielrichter2007/grub-customizer  # If Ubuntu 22.04
sudo apt install grub-customizer
```

### Fixing wrong time in Windows

```bash
sudo timedatectl set-local-rtc 1  # 0 to revert
```

### ACPI BIOS Error / AE_ALREADY_EXISTS

The notice is benign, is not the cause of a hangup, All it is is a mismatch between the bios and the kernel where the bios does not live up to its specifications.

You can get rid of those ACPI Error messages by

Open `/etc/default/grub` in an editor with root access. In your case I believe Ubuntu uses gedit as it's text editor.

```
sudo gedit /etc/default/grub
```

The line with `GRUB_CMDLINE_LINUX_DEFAULT`, add the `loglevel=3` part. The original looks like

```
GRUB_CMDLINE_LINUX_DEFAULT='quiet splash'
```

Change it to this:

```
GRUB_CMDLINE_LINUX_DEFAULT='quiet splash loglevel=3'
```

Then save the changes and close it, now open a terminal and run:

```
sudo update-grub
```

## Google Chrome

Deleting firefox with:

```
sudo apt purge firefox
sudo snap remove firefox
```

Install Google Chrome at https://www.google.com/chrome/.

## onedrive

Instructions for Ubuntu 22.04. For other versions, see [here](https://github.com/abraunegg/onedrive/blob/master/docs/ubuntu-package-install.md#known-issues-with-installing-from-the-above-packages).

```bash
wget -qO - https://download.opensuse.org/repositories/home:/npreining:/debian-ubuntu-onedrive/xUbuntu_22.04/Release.key | gpg --dearmor | sudo tee /usr/share/keyrings/obs-onedrive.gpg > /dev/null
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/obs-onedrive.gpg] https://download.opensuse.org/repositories/home:/npreining:/debian-ubuntu-onedrive/xUbuntu_22.04/ ./" | sudo tee /etc/apt/sources.list.d/onedrive.list
sudo apt-get update
sudo apt install --no-install-recommends --no-install-suggests onedrive
```

Configure onedrive:

```bash
onedrive  # Sign in
onedrive --resync --synchronize
```

TODO: Webhooks

## General dependencies

```
sudo apt install build-essential 
sudo apt install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libxml2-dev  # tidyverse dependencies
sudo apt install cmake
```

## General tools and software

```bash
sudo apt install vim 
sudo apt install pigz  # Gzip parallel implementation
sudo apt install alacarte  # Edit applications in the menu
sudo snap install discord 
sudo snap install thunderbird 
sudo snap install bitwarden 
sudo snap install kuro-desktop
```

```bash
sudo snap install spotify  #?
```

NOTE: SpotX-Bash does not work with the snap version of Spotify. Use the .deb version instead.

```bash
curl -sS https://download.spotify.com/debian/pubkey_7A3A762FAFD4A51F.gpg | sudo gpg --dearmor --yes -o /etc/apt/trusted.gpg.d/spotify.gpg
echo "deb http://repository.spotify.com stable non-free" | sudo tee /etc/apt/sources.list.d/spotify.list
sudo apt-get update && sudo apt-get install spotify-client
```

## Bioinformatic tools

```bash
sudo apt install bedtools
```

## R and RStudio

### R installation

Setting up the R installation in Ubuntu 22.04 (instructions from https://cran.r-project.org/)

```bash
# update indices
sudo apt update -qq
# install two helper packages we need
sudo apt install --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
```

Finally installing R:

```bash
# Base R
sudo apt install --no-install-recommends r-base

# For compillation of packages
sudo apt-get install r-base-dev
```

### RStudio installation

Installing RStudio at https://posit.co/download/rstudio-desktop/

### R packages

See [[Basic R package installation]].
TODO: Continue R package installation

## Python and conda

Installing conda with miniconda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh 
```

Add conda to path:

```bash
vim ~/.bashrc
```

And adding the following line:

```bash
export PATH=~/miniconda3/bin:$PATH
```

Adding conda channels:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels anaconda
```

![[2023-10-13#On conda's default solver]]

## Zotero

Follow steps at: [installation - Zotero Documentation](https://www.zotero.org/support/installation)

```bash
sudo mkdir /opt/zotero
```

Zotero configuration: [[Zotero configuration]]
## Website software

- Google Chrome
- Obsidian