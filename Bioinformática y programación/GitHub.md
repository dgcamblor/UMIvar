Esentially, there are two ways of authenticating and cloning repositories in GitHub:

- HTTPS. To clone a repository using HTTPS, you can use the following command:
  
```bash
git clone https://github.com/ACCOUNT/REPO
```

- SSH. To clone a repository using SSH, you can use the following command:

```bash
git clone git@github.com:ACCOUNT/REPO.git
```

or

```bash
git clone ssh://github.com/ACCOUNT/REPO
```

You can check the specific remote that is being used (HTTPS/SSH) with:

```bash
git remote -v
origin	git@github.com:dgcamblor/Notas-PKMS.git (fetch)
origin	git@github.com:dgcamblor/Notas-PKMS.git (push)
```