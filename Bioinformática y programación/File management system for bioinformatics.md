%%TODO
- [ ] Revisar [How Do You Manage Your Files & Directories For Your Projects ?](https://www.biostars.org/p/821/)
%%

#write 

We want to follow the rule: 

> One project, one folder

The folder itself should be 

## INCLIVA structure

- `analysis` -> Analysis results.
- `rawdata` -> Symbolic link to FASTQ files (stored in `/nfs`).
	- `_1`
	- `_2`
	- `_UMI`
- `scripts` -> `inc_pipe` and `.conf` file
- `Sheila` -> Folder for Sarek scripts (germline variant detection enhancing).
- `tmp` -> Temporal files generated along the `inc_pipe`.