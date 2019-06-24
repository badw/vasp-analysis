# vasp-analysis
This is mostly an area for me (BADW) to learn python - but also make something useful out of it.

### setup

simply run:
``` 
git clone git@github.com:badw/vasp-analysis.git

cd ./vasp-analysis

pip3 install --user .
```
### param-bonds
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
