
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sparse
import sys 
import os
from subprocess import check_call

def usage():
    print('Usage: python convert.py [h5ad_file]')

def convert():
    filename =  os.path.splitext(sys.argv[1])[0]
    if not os.path.lexists(filename):
        os.makedirs(filename)
    #Create dir
    h5ad_file = sc.read(filename= sys.argv[1])
    h5ad_file.obs.to_csv(f'./{filename}/{filename}_metadata.txt.gz', compression='gzip',sep='\t', index=True)
    
    #write metadata
    if 'counts' in h5ad_file.layers.keys():
        h5ad_file.X = h5ad_file.layers['counts']
    sio.mmwrite(f'./{filename}/{filename}.mtx',sparse.csr_matrix(h5ad_file.X.T))
    #write sparce matrix
    h5ad_file.var.to_csv(f'./{filename}/{filename}_features.txt.gz',compression='gzip',sep='\t')
    #write features 
    check_call(['gzip', f'./{filename}/{filename}.mtx'])
    #gzip files
    print('Converted finish')
if __name__ == '__main__':
    try:
        convert()
    except IndexError:
        usage()
