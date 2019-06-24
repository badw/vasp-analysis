# vasp-analysis
This is mostly an area for me (BADW) to learn python - but also make something useful out of it.

### setup

simply run:
``` 
git clone git@github.com:badw/vasp-analysis.git

cd ./vasp-analysis

pip3 install --user .
```
### parabo
this should hopefully give you:

a) a table of cell parameters:

* a,b,c
* alpha,beta,gamma
* volume

b) a list of bond lengths (in Å)

simply run:
```
param-bonds (-f POSCARFILE)
```
options for:

* changing the cutoff distance for the bond lengths with `-c` (current cutoff is 4 Angstroms)
* changing the accuracy to the next decimal place with `-d`
* output a `.csv` with `-x` file to easily import into Microsoft Excel
* output the data into a latex table format with `-l` (makes a parameters.tex file)
* 
