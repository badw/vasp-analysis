import sys
import argparse

import numpy as np
import pandas as pd
from itertools import combinations
from pymatgen.core.structure import Structure
from tabulate import tabulate

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, default='POSCAR',
                        help='path to input file')
    parser.add_argument('-c', '--cutoff', default=4, type=float,
                        help='nearest neighbour radius cutoff (default=3Å)')
    parser.add_argument('-d', '--decimal', default=2, type=int,
                        help='how many decimal places to round the data to (default=2) ')
    parser.add_argument('-x', '--excel', default=False, action='store_true',
                        help='output data in .csv format to easily paste into Microsofrt Excel')
    parser.add_argument('-l', '--latex', default=False, action='store_true',
                        help='output data to a LaTEX table to put into a .tex file')
    args = parser.parse_args()
    
    poscar = Structure.from_file(args.file)
    pdict = poscar.as_dict()
    psites = pdict['sites'] 
    atoms = len(psites)
  
    ''' Table Structural Parameters ''' 
    
    df1 = pd.DataFrame(columns=['parameter', 'value']) 

    a = pdict['lattice']['a']
    b = pdict['lattice']['b']
    c = pdict['lattice']['c']
    alpha = pdict['lattice']['alpha']
    beta = pdict['lattice']['beta']
    gamma = pdict['lattice']['gamma']
    volume = pdict['lattice']['volume']
    
    if a == b and b == c:
        df1 = df1.append({'parameter' : 'a = b = c / Å' , 'value' : round(a,args.decimal)},ignore_index=True)
    elif a == b and b != c:
        df1 = df1.append({'parameter' : 'a = b / Å' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'c / Å' , 'value' : round(c,args.decimal)},ignore_index=True)
    elif a == c and b != a:
        df1 = df1.append({'parameter' : 'a = c / Å' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'b / Å' , 'value' : round(b,args.decimal)},ignore_index=True)
    elif a != b and b == c:
        df1 = df1.append({'parameter' : 'a / Å' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'b = c/ Å' , 'value' : round(b,args.decimal)},ignore_index=True)
    elif a != b and b != c:
        df1 = df1.append({'parameter' : 'a / Å' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'b / Å' , 'value' : round(b,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'c / Å' , 'value' : round(c,args.decimal)},ignore_index=True)
        
    if alpha == beta and beta == gamma:
        df1 = df1.append({'parameter' : 'alpha = beta = gamma  / o' , 'value' : round(alpha,args.decimal)},ignore_index=True)
    elif alpha == beta and beta != gamma:
        df1 = df1.append({'parameter' : 'alpha = beta / o' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'gamma / o' , 'value' : round(gamma,args.decimal)},ignore_index=True)
    elif alpha == gamma and beta != gamma:
        df1 = df1.append({'parameter' : 'alpha = gamma / o' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'beta / o' , 'value' : round(beta,args.decimal)},ignore_index=True)
    elif a != b and b == c:
        df1 = df1.append({'parameter' : 'alpha / o' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'beta = gamma / o' , 'value' : round(beta,args.decimal)},ignore_index=True)
    elif alpha != beta and beta != gamma:
        df1 = df1.append({'parameter' : 'alpha / o' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'beta / o' , 'value' : round(beta,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'gamma / o' , 'value' : round(gamma,args.decimal)},ignore_index=True)
       
    df1 = df1.append({'parameter' : 'Vol. / Å^3', 'value' : round(volume,args.decimal)},ignore_index=True)
    
    '''Table of Bond lengths'''
    
    df2 = pd.DataFrame(columns=['bond', 'length / Å'])
    
    p = combinations(range(atoms),2)
    
    for elements in p:
        i = elements[0] 
        j = elements[1]
        dist = poscar.get_distance(i,j).round(args.decimal)
        if dist <= args.cutoff:
            sitei = psites[i]['label']
            sitej = psites[j]['label']
            if sitei != sitej:
                df2 = df2.append({'bond' : '{}-{}'.format(sitei,sitej), 'length / Å' : dist},ignore_index=True)
    
    df1 = df1.drop_duplicates()
    df1 = df1.set_index('parameter')
    df2 = df2.drop_duplicates()     
    df2 = df2.set_index('bond')

    if args.excel == True:
        print(df1.to_csv(index=True,header=True))
        print(df2.to_csv(index=True,header=True))
    elif args.latex == True:
        print(df1.to_latex(index=True))
        print(df2.to_latex(index=True))
    else:
        print(tabulate(df1,headers='keys',tablefmt='psql', numalign="center"))
        print(tabulate(df2,headers='keys',tablefmt='psql', numalign="center"))

if __name__ == "__main__":
    main()
