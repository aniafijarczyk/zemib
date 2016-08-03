'''
Created on Jul 10, 2015

@author: ania.fijarczyk

This script calculates rough estimate of FMR: Fraction of Mapped Reads. This is a fraction of all mapped reads 
which map to the particular reference in an individual. It estimates how coverage is distributed 
throughout references. It reads sam file which is generated with --no-unal function 
(i.e. all reads in sam are considered mapped).

It also outputs number of reads covering reference. It sums all reads mapped to a given 
reference (SN). To obtain a rough estimate of coverage per bp, you have to multiply this 
value by the length of reads and divide by the length of reference (LN).

USAGE: python sam_coverage.py
By default it reads all *.sam files

OUTPUT:
sam_fmr.out and sam_coverage.out
 
'''
from collections import defaultdict


def formatHeader(lista):
    a = [ele.split(':')[1] for ele in lista if 'SN' in ele]
    b = [ele.split(':')[1] for ele in lista if 'LN' in ele]
    new_list = a+b
    return new_list

def readMe(filename):
    fh = open(filename,'r')
    linie = fh.readlines()
    
    # Create dictionary with references and their lengths extracted from sam file
    
    len_list = [formatHeader(ele.split()) for ele in linie if '@SQ' in ele]
    L = {ref:int(length) for ref, length in len_list}
    
    
    # Create dictionary with reference and number of mapped reads
    
    k = [ele.split() for ele in linie if ele.startswith('@') == False]
    suma_readow = len(k)
    d = defaultdict(list) 
    for a in k:
        d[a[2]].append(a[1])
    
    # Dictionary with FMR
    
    D = {mip: len(d[mip])/float(suma_readow) for mip in d.keys()}
    
    # Create dictionary with FMR, number of reads, and number of reads per reference length
    
    Z = {ele: (0,0,0) for ele in L.keys()} # dictionary with all refs and 0 values
    
    for ref_seq in D.keys():
        Z[ref_seq] = D[ref_seq], len(d[ref_seq]), len(d[ref_seq])/float(L[ref_seq])
    
    return Z


if __name__ == '__main__':
    import glob
       
    #plik = "test"
    #r = readMe(plik)
    #print r
    
    flist = glob.glob('*.sam')
    
    list_fmr = defaultdict(list)
    list_reads = defaultdict(list)
    list_coverage = defaultdict(list)
    
    
    # list of fmr
    i = 0
    for plik in flist:
        r = readMe(plik)
        f = {ele : r[ele][0] for ele in r.keys()}
        re = {ele : r[ele][1] for ele in r.keys()}
        cov = {ele : r[ele][2] for ele in r.keys()}
            
        for refer in re.keys():
            p = (i, plik, re[refer])
            q = (i, plik, f[refer])
            list_reads[refer].append(p)
            list_fmr[refer].append(q)    
        i+=1
    
    #print list_reads
    '''    
    # ref names
    one_file = flist[0]
    refs_names = ['header'] + readMe(one_file).keys()
    
    fmrs = zip(*[refs_names] + list_fmr)
    reads = zip(*[refs_names] + list_reads)
    coverage = zip(*[refs_names] + list_coverage)

    #print fmrs
    '''
    
    # Write down FMR
    head1 = list_fmr.values()[0]
    head_sorted1 = sorted(head1, key=lambda x: int(x[0]))
    header1 = [k[1] for k in head_sorted1]
    
    wh1 = open("sam_fmr.out",'w')
    h1 = 'ref' + '\t' + '\t'.join(header1)+'\n'
    wh1.write(h1)

    for ele in list_fmr.keys():
        lista = list_fmr[ele]
        lista_sorted = sorted(lista, key=lambda x: int(x[0]))
        lista_format = [str(j[2]) for j in lista_sorted]
        w = [ele] + lista_format
        W = '\t'.join(w)+'\n'
        wh1.write(W)
        
    wh1.flush()
    wh1.close()
    
    
    # Write down number of reads in reference
    head = list_reads.values()[0]
    head_sorted = sorted(head, key=lambda x: int(x[0]))
    header = [k[1] for k in head_sorted]
    
    wh3 = open("sam_reads.out",'w')
    h = 'ref' + '\t' + '\t'.join(header)+'\n'
    wh3.write(h)

    for ele in list_reads.keys():
        lista = list_reads[ele]
        lista_sorted = sorted(lista, key=lambda x: int(x[0]))
        lista_format = [str(j[2]) for j in lista_sorted]
        w = [ele] + lista_format
        W = '\t'.join(w)+'\n'
        wh3.write(W)
    wh3.flush()
    wh3.close()
        
        
    #print H
