import sys
import argparse

import numpy as np
import pandas as pd
from itertools import combinations
from pymatgen.core.structure import Structure
from prettytable import PrettyTable
from tabulate import tabulate

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, default='POSCAR',
                        help='path to input file')
    parser.add_argument('-c', '--cutoff', default=4, type=float,
                        help='nearest neighbour radius cutoff (default 3Å)')
    parser.add_argument('-x', '--excel', default=False, action='store_true',
                        help='output data to a .csv file to import to Microsoft Excel')
    args = parser.parse_args()
    
    poscar = Structure.from_file(args.file)
    pdict = poscar.as_dict()
    psites = pdict['sites'] 
    atoms = len(psites)
    ''' Table Structural Parameters ''' 
    t1 = PrettyTable(['parameter', 'value'])
    
    a = pdict['lattice']['a']
    b = pdict['lattice']['a']
    c = pdict['lattice']['a']
    alpha = pdict['lattice']['alpha']
    beta = pdict['lattice']['beta']
    gamma = pdict['lattice']['gamma']
    volume = pdict['lattice']['volume']
    
    if a == b and b == c:
        t1.add_row(['a = b = c / Angst', round(a,2)])
    elif a == b and b != c:
        t1.add_row(['a = b / Angst', round(a,2)])
        t1.add_row(['c / Angst', round(c,2)])
    elif a == c and b != a:
        t1.add_row(['a = c / Angst', round(a,2)])
        t1.add_row(['b / Angst', round(b,2)])
    elif a != b and b == c:
        t1.add_row(['a  / Angst', round(a,2)])
        t1.add_row(['b = c / Angst', round(b,2)])
    elif a != b and b != c:
        t1.add_row(['a / Angst', round(a,2)])
        t1.add_row(['b / Angst', round(b,2)])
        t1.add_row(['c / Angst', round(c,2)])
        
    if alpha == beta and beta == gamma:
        t1.add_row(['α = β = γ / deg', round(alpha,2)])
    elif alpha == beta and beta != gamma:
        t1.add_row(['α = β / deg', round(alpha,2)])
        t1.add_row(['γ / deg', round(gamma,2)])
    elif alpha == gamma and beta != gamma:
        t1.add_row(['α = γ / deg', round(alpha,2)])
        t1.add_row(['β / deg', round(beta,2)])
    elif a != b and b == c:
        t1.add_row(['α / deg', round(alpha,2)])
        t1.add_row(['β = γ/ deg', round(beta,2)])
    elif alpha != beta and beta != gamma:
        t1.add_row(['α / deg', round(alpha,2)])
        t1.add_row(['β / deg', round(beta,2)])
        t1.add_row(['γ / deg', round(gamma,2)])
        
    t1.add_row(['Vol. / Angst^3', round(volume,2)])
    
    '''Table of Bond lengths'''
    
    p = combinations(range(atoms),2)
    
    df = pd.DataFrame(columns=['bond', 'length / Angst'])
    
    for elements in p:
        i = elements[0] 
        j = elements[1]
        dist = poscar.get_distance(i,j).round(2)
        if dist <= args.cutoff:
            sitei = psites[i]['label']
            sitej = psites[j]['label']
            if sitei != sitej:
                df = df.append({'bond' : '{}-{}'.format(sitei,sitej), 'length / Angst' : dist},ignore_index=True)
    
    df = df.drop_duplicates()     
    df = df.set_index('bond')
    print(t1)   
    print(tabulate(df,headers='keys',tablefmt='psql', numalign="center"))

    if args.excel == True:
        df.to_csv('bond-lengths.csv',index=True, header=True)
        print("output saved to bond-lengths.csv")
if __name__ == "__main__":
    main()
