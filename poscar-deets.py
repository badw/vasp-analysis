#!/usr/bin/env python3
import sys
import argparse
import itertools
from itertools import combinations
from pymatgen.core.structure import Structure
from prettytable import PrettyTable
import numpy as np
import pandas as pd
from tabulate import tabulate

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', type=str, default='POSCAR', help='path to input file')
parser.add_argument('-c', '--cutoff', default=3, type=float, help='nearest neighbour radius cutoff (default 3Å)')
args = parser.parse_args()

poscar = Structure.from_file(args.file)
posdict = poscar.as_dict()
possites = posdict['sites'] 
atoms = len(possites)


t1 = PrettyTable(['parameter', 'value'])

a = posdict['lattice']['a']
b = posdict['lattice']['a']
c = posdict['lattice']['a']
alpha = posdict['lattice']['alpha']
beta = posdict['lattice']['beta']
gamma = posdict['lattice']['gamma']
volume = posdict['lattice']['volume']

if a == b and b == c:
    t1.add_row(['a = b = c / Å', round(a,2)])
elif a == b and b != c:
    t1.add_row(['a = b / Å', round(a,2)])
    t1.add_row(['c / Å', round(c,2)])
elif a == c and b != a:
    t1.add_row(['a = c / Å', round(a,2)])
    t1.add_row(['b / Å', round(b,2)])
elif a != b and b == c:
    t1.add_row(['a  / Å', round(a,2)])
    t1.add_row(['b = c / Å', round(b,2)])
elif a != b and b != c:
    t1.add_row(['a / Å', round(a,2)])
    t1.add_row(['b / Å', round(b,2)])
    t1.add_row(['c / Å', round(c,2)])
    
if alpha == beta and beta == gamma:
    t1.add_row(['α = β = γ / º', round(alpha,2)])
elif alpha == beta and beta != gamma:
    t1.add_row(['α = β/ º', round(alpha,2)])
    t1.add_row(['γ / º', round(gamma,2)])
elif alpha == gamma and beta != gamma:
    t1.add_row(['α = γ / º', round(alpha,2)])
    t1.add_row(['β / º', round(beta,2)])
elif a != b and b == c:
    t1.add_row(['α / º', round(alpha,2)])
    t1.add_row(['β = γ/ º', round(beta,2)])
elif alpha != beta and beta != gamma:
    t1.add_row(['α / º', round(alpha,2)])
    t1.add_row(['β / º', round(beta,2)])
    t1.add_row(['γ / º', round(gamma,2)])
    
t1.add_row(['Vol. /Å^3', round(volume,2)])

'''Bond lengths using a pandas dataframe'''

p = combinations(range(atoms),2)

df = pd.DataFrame(columns=['bond', 'length / Å'])

for elements in p:
    i = elements[0] 
    j = elements[1]
    dist = poscar.get_distance(i,j).round(2)
    if dist <= args.cutoff:
        sitei = possites[i]['label']
        sitej = possites[j]['label']
        if sitei != sitej:
            df = df.append({'bond' : '{}-{}'.format(sitei,sitej), 'length / Å' : dist},ignore_index=True)



df = df.drop_duplicates()     
df = df.set_index('bond')
print(t1)   
print(tabulate(df,headers='keys',tablefmt='psql'))
