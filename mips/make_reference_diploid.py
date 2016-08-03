'''
Created on Jun 17, 2014

@author: ania.fijarczyk

USAGE:
python make_reference_diploid.py
by default it reads all <name>.fasta files (sequential) in the directory 

OUTPUT:
same fasta files with consensus seuence "chr"<name>.fa

The script reads fasta files of the type <name>.fasta
For each file it reads the fasta alignment (can be a single sequence) and creates a consensus sequence
Generally it is meant to create a consensus sequence from the alignment of diploid sequences with SNPs encoded with IUPAC code
Individual ID's are not considered, all sequences are treated equally
Each base of each sequence is resolved into a couple of bases, i.e. Y -> T, C; T -> T, T
For the purpose of mips design it adds "chr" prefix to each sequence name


Steps:
1. For each position in each sequence, a base is resolved into a couple of bases (e.g. i.e. [Y] -> [T, C]; [T] -> [T, T])
2. Each variant, including 'N' (that includes 'n') is counted
3. If proportion of missing data exceeds or is equal to some value X=0.66 (i.e. proportion of 'N' in the sequence), the consensus base is 'N'
   If proportion of missing data is lower than X=0.66, the consensus base is the base with the highest count

   If more variants have the same count, one is chosen randomly
   Positions '-' are ignored by this script


Generally it is meant to create a consensus sequence from the alignment of diploid sequences with SNPs encoded with IUPAC code
If proportion of misssing data >= 0.66 (can be adjusted), then consensus has 'N'

'''


import random

diploid = {'K':('G','T'),'M':('A','C'),'R':('A','G'),'Y':('C','T'),'S':('C','G'), \
        'W':('A','T'),'V':('A','C','G'),'H':('A','C','T'),'D':('A','G','T'),'B':('C','G','T')}

ambig = ['M','R','W','S','Y','K','V','H','D','B']


def uniq(inlist):    #wybor unikatowych elementow z listy redundantnych elementow
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

def cleaning(lista):
    clean_lista = [ele for ele in lista if (ele != '-') and (ele != 'n') and (ele != 'N')]
    return clean_lista

def resolve_ambigs(lista):
    new_lista = []
    for base in lista:
        if base in ambig:
            nts = diploid[base]
            new_lista+=nts
        else:
            b = base,base
            new_lista+=b
    return new_lista

def readMe(fname):
    
    fh = open(fname, 'r')
    linie = fh.readlines()
    linie2 = [a for a in linie if a != '\n']
    seq=[ele.split()[0].upper() for ele in linie2 if '>' not in ele]

        
    # dictionary with all variants per position
    P={}
    j=0
    all_sample_variants = zip(*seq)
    for position in all_sample_variants:
        resolved = resolve_ambigs(position)
        P[j] = resolved
        j+=1
        
    return P

class POS:
    def __init__(self, kolumna):
        self.lista = kolumna
        self.length = len(kolumna)
        self.called = cleaning(kolumna)
        self.nr = uniq(cleaning(kolumna))
        self.counts = {}
        self.consensus = ''
        
        self.count(self.lista)
        self.choose_consensus(self.counts,self.length)
        
    def count(self, lista):
        
        upper_lista = [ele.upper() for ele in lista]

        new_dic = {}
        new_dic['A'] = upper_lista.count('A')
        new_dic['C'] = upper_lista.count('C')
        new_dic['T'] = upper_lista.count('T')
        new_dic['G'] = upper_lista.count('G')
        new_dic['N'] = upper_lista.count('N')
        
        self.counts = new_dic
        
    def choose_consensus(self, leks, length):
   
        new_leks = {i:leks[i] for i in leks if i != 'N'}
        base_counts = new_leks.values()
        major_base = max(base_counts)
        n_counts = leks['N']
        t = float(n_counts)/float(length)
        if t >= 0.66:
            N = 'N'
        elif major_base == 0:
            N = 'N'
        elif major_base > 0:
            n = [ele for ele in new_leks.keys() if new_leks[ele] == major_base]
            N = random.sample(n,1)[0]
        
        self.consensus = N
        
if __name__ == '__main__':
    import glob
    flist = glob.glob('*.fasta')
    #flist = ('adar.fasta','agl.fasta','appl.fasta')
    #flist = ['test.fasta']
    
    #petla po plikach
    for fname in flist:
        
        r = readMe(fname)
        nazwa = 'chr'+fname.split('.')[0]

        C = []
        for position in r.keys():
            
            pos = POS(r[position])
            c = pos.consensus
            C.append(c)
    
        seq = ''.join(C)
        #writing fasta file with reference
        wh =  open(nazwa+'.fa','w')
        wh.write('>'+nazwa+'\n')
        wh.write(seq+'\n')
    
    #print C
    
