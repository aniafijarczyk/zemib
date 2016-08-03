'''
Created on May 23, 2014

@author: ania.fijarczyk

USAGE:
python mip_dbsnp.py
By default it reads all <name>.fasta and <chrfasta>.fa files in the directory

OUTPUT:
dbSNP.rod

Script reads fasta files with alignments (<file>.fasta) and fasta files with references (<chrfile>.fa) and creates simple dbSNP.rod file necessary for mip design
SNPs in positions where one variant (singleton) is counterbalanced by minimum 5 other variants (called bases), are IGNORED  
If sum of all variants is less than 6, the SNP is kept


1. Reads references; in case there is an ambiguous base in the reference, one nucleotide is chosen randomly (e.g. for Y, there will be C or T)
2. In each position of the alignment
	a. If there are any SNPs, they are resolved in two bases, non-SNPs are kept single
	b. Any missing variants (-, n or N) are removed
	c. If sum of all alleles < 6:
		All alleles are kept, and in case there are multiple variants, SNP is generated in this position
	   If sum of all alleles >= 6:
		If there are singletons, they are ignored, and the position is invariant, otherwise (if they are not singletons) SNP is generated in this position

Additionally it produces default bed file, where each region referrs to one reference, starts at 0 and end at the last base pair

'''

# filter for singletons, singletons are removed if there are 6 samples or more


import random

#Molecule Type: Sample used to find this variant
#    Genomic - variant discovered using a genomic template
#    cDNA - variant discovered using a cDNA template
#    Unknown - sample type not known 
molType = 'cDNA'

#Class: Describes the observed alleles
#Class = ['unknown', 'single', 'in-del', 'het', 'microsatellite', 'named', 'mixed', 'mnp', 'insertion', 'deletion']
Class = 'single'

#Validation:  Method used to validate the variant (each variant may be validated by more than one method)
#valid =['unknown', 'by-cluster', 'by-frequency', 'by-submitter', 'by-2hit-2allele', 'by-hapmap']
valid = 'unknown'

#avHet: Average heterozygosity from all observations (float)
avHet = '0'

#avHetSE: Standard Error for the average heterozygosity (float)
avHetSE = '0'

#Function: Predicted functional role (each variant may have more than one functional role)
#'unknown', 'coding-synon', 'intron', 'cds-reference', 'near-gene-3', 'near-gene-5', 'nonsense', 'missense', 'frameshift', 'untranslated-3', 'untranslated-5', 'splice-3', 'splice-5'
func = 'untranslated-3'

# ??? locType: Type of mapping inferred from size on reference; may not agree with class
#'range', 'exact' (when reference is known), 'between' (when reference is '-'), 'rangeInsertion', 'rangeSubstitution', 'rangeDeletion'
#locType = 'between'

#Weight: Alignment quality assigned by dbSNP
#1 = unique mapping, 2 = non-unique, 3 = many matches, 10 = excluded from the data set
weight = '1'

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
            new_lista.append(base)
    return new_lista

def ref_ambig(seq):
    new_seq = []
    for base in seq:
        if base in ambig:
            new_base = random.sample(diploid[base],1)[0]
        else: new_base = base
        new_seq.append(new_base)
        
    ref_seq =''.join(new_seq)
    return ref_seq

def readMe(filename):
    fh=open(filename,'r')
    linie=fh.readlines()
    name=[ele.split()[0] for ele in linie if '>' in ele]
    seq=[ele.split()[0].upper() for ele in linie if '>' not in ele]
    combi=zip(name,seq)
    D={}
    for a,b in combi:
        D[a[1:]]=b
    
    return D

def filter_singletons(lista):
    if lista == []: pass
    nowa_lista = []
    singletons = []
    n = len(lista)

    if n < 6:
        nowa_lista = lista
    elif n >= 6:
        for base in uniq(lista):
            c = float(lista.count(base))/n
            if c == 1./n:
                singletons.append(base)
        nowa_lista = [base for base in lista if base not in singletons]
        
    return nowa_lista

def get_column(leks_snip, leks_ref):
    
    #dictionary with reference nucleotides
    R={}
    i=0 # BED files have 0 base start
    ref_seq = leks_ref.values()[0]
    ref_seq_no_ambig = ref_ambig(ref_seq)
    for position in ref_seq_no_ambig:
        R[i]=position
        i+=1
    
    #dictionary with only variants per position
    P={}
    j=0 # BED files have 0 base start 
    all_sample_variants = zip(*leks_snip.values())
    for position in all_sample_variants:
        resolved = resolve_ambigs(position)
        selected = cleaning(resolved)
        filtered = filter_singletons(selected)
        if len(uniq(filtered))>1:
                P[j] = R[j], '/'.join(uniq(filtered))
        j+=1
    
    return P

def write_dbFile(lista):
    fh=open('dbSNP.rod','w')
    header = 'bin\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\trefNCBI\trefUCSC\tobserved\tmolType\tclass\tvalid\tavHet\tavHetSE\tfunc\tlocType\tweight\n'
    fh.write(header)
    
    for linijka in lista:
        p='\t'.join(linijka)+'\n'
        fh.write(p)
        
    fh.flush()
    fh.close()
    
def write_Bed(name_and_seq):
    wh=open('reference.bed','w')
    for name,seq in name_and_seq:
        start = '0'
        stop = str(len(seq) - 1)
        wh.write(name+'\t'+start+'\t'+stop+'\n')
        
    wh.flush()
    wh.close() 
    
    
class SNP:
    def __init__(self, D_key, D_val):
        pos = D_key
        lista = D_val
        
        self.all_nt = lista
        self.ref = lista[0]
        self.snips = lista[1]
        self.start = str(pos)
        self.stop = str(int(pos) + 1)
             
if __name__ == '__main__':
    import glob,os
    
    flist = glob.glob("*.fasta")  
    #flist = ("adar.fasta","agl.fasta")
    
    #bin?
    binCalc = 500
    
    all_W = []
    seqs_for_bed = []
    Contigs =[]
    path_to_fasta =[]
    
    for fname in flist:
        
        #nazwa referencji
        contig_name = 'chr'+fname.split(".")[0]
        Contigs.append(contig_name)
        
        # reading in alignments and references
        D_snip = readMe(fname)
        D_ref = readMe(glob.glob(contig_name+'.fa')[0])
        snip_column = get_column(D_snip,D_ref)
        
        #print snip_column
        
        # for BED file
        r_seq = D_ref.values()[0]
        pairs = contig_name, r_seq
        seqs_for_bed.append(pairs)
        
        # filtr na puste leksykony!!!
        if snip_column == {}:
            continue
        
        # klasa SNP dla kazdej pozycji ze snipem
        W=[]
        j=1
        for variant_list in sorted(snip_column.keys()):
            
            db = SNP(variant_list, snip_column[variant_list])
            lista = db.all_nt
            reference = db.ref
            snips = db.snips
            chromStart = db.start
            chromStop = db.stop
            
            snp_name = 'snp_'+contig_name+'%05d' %j
            
            if reference == '-': locType = 'between'
            else: locType = 'exact'
                
            wers = str(binCalc), contig_name, chromStart, chromStop, snp_name, '0', '+', reference, reference, snips, molType, Class, valid, avHet, avHetSE, func, locType, weight
            W.append(wers)
            
            j+=1
        binCalc +=1
        
        all_W+=W
        
    # printing dbSNP.rod file
    w = write_dbFile(all_W)
    
    # printing BED file
    w_bed = write_Bed(seqs_for_bed)


