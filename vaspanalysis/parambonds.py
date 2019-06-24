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
                        help='nearest neighbour radius cutoff (default 3 Angstroms)')
    parser.add_argument('-d', '--decimal', default=2, type=int,
                        help='how many decimal places to round the data to (default = 2) ')
    parser.add_argument('-x', '--excel', default=False, action='store_true',
                        help='output data to a .csv file to import to Microsoft Excel')
    args = parser.parse_args()
    
    poscar = Structure.from_file(args.file)
    pdict = poscar.as_dict()
    psites = pdict['sites'] 
    atoms = len(psites)
  
  ''' Table Structural Parameters ''' 
    
    df1 = pd.DataFrame(columns=['parameter', 'value']) 

    a = pdict['lattice']['a']
    b = pdict['lattice']['a']
    c = pdict['lattice']['a']
    alpha = pdict['lattice']['alpha']
    beta = pdict['lattice']['beta']
    gamma = pdict['lattice']['gamma']
    volume = pdict['lattice']['volume']
    
    if a == b and b == c:
        df1 = df1.append({'parameter' : 'a = b = c / Angst' , 'value' : round(a,args.decimal)},ignore_index=True)
    elif a == b and b != c:
        df1 = df1.append({'parameter' : 'a = b / Angst' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'c / Angst' , 'value' : round(c,args.decimal)},ignore_index=True)
    elif a == c and b != a:
        df1 = df1.append({'parameter' : 'a = c / Angst' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'b / Angst' , 'value' : round(b,args.decimal)},ignore_index=True)
    elif a != b and b == c:
        df1 = df1.append({'parameter' : 'a / Angst' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'b = c/ Angst' , 'value' : round(b,args.decimal)},ignore_index=True)
    elif a != b and b != c:
        df1 = df1.append({'parameter' : 'a / Angst' , 'value' : round(a,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'b / Angst' , 'value' : round(b,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'c / Angst' , 'value' : round(c,args.decimal)},ignore_index=True)
        
    if alpha == beta and beta == gamma:
        df1 = df1.append({'parameter' : 'alph = bet = gam  / deg' , 'value' : round(alpha,args.decimal)},ignore_index=True)
    elif alpha == beta and beta != gamma:
        df1 = df1.append({'parameter' : 'alph = bet / deg' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'gam / deg' , 'value' : round(gamma,args.decimal)},ignore_index=True)
    elif alpha == gamma and beta != gamma:
        df1 = df1.append({'parameter' : 'alph = gam / deg' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'beta / deg' , 'value' : round(beta,args.decimal)},ignore_index=True)
    elif a != b and b == c:
        df1 = df1.append({'parameter' : 'alph / deg' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'beta = gam / deg' , 'value' : round(beta,args.decimal)},ignore_index=True)
    elif alpha != beta and beta != gamma:
        df1 = df1.append({'parameter' : 'alph / deg' , 'value' : round(alpha,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'beta / deg' , 'value' : round(beta,args.decimal)},ignore_index=True)
        df1 = df1.append({'parameter' : 'gam / deg' , 'value' : round(gamma,args.decimal)},ignore_index=True)
       
    df1 = df1.append({'parameter' : 'Vol. / Angst^3', 'value' : round(volume,args.decimal)},ignore_index=True)
    
    '''Table of Bond lengths'''
    
    df2 = pd.DataFrame(columns=['bond', 'length / Angst'])
    
    p = combinations(range(atoms),2)
    
    for elements in p:
        i = elements[0] 
        j = elements[1]
        dist = poscar.get_distance(i,j).round(args.decimal)
        if dist <= args.cutoff:
            sitei = psites[i]['label']
            sitej = psites[j]['label']
            if sitei != sitej:
                df2 = df2.append({'bond' : '{}-{}'.format(sitei,sitej), 'length / Angst' : dist},ignore_index=True)
    
    df1 = df1.drop_duplicates()
    df1 = df1.set_index('parameter')
    df2 = df2.drop_duplicates()     
    df2 = df2.set_index('bond')
    print(tabulate(df1,headers='keys',tablefmt='psql', numalign="center"))
    print(tabulate(df2,headers='keys',tablefmt='psql', numalign="center"))

    if args.excel == True:
        df1.to_csv('cell-parameters.csv',index=True,header=True)
        df2.to_csv('bond-lengths.csv',index=True, header=True)
        print("output saved to cell-parameters.csv and bond-lengths.csv")

if __name__ == "__main__":
    main()
