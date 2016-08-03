'''
Created on Jul 11, 2014

@author: ania.fijarczyk

Script reads *global.gff3 file with coordinates of model exons on the reference,
merges coordinates of the same exons together, 
and if exons are overlapping on the reference, extracts coordinates of regions inbetween exon boundaries 
outputs table with start and stop positions for each exon plus some additional information

example:
ref:      ***************************
exon 1:   xxxxxxxxxxxxxxxxx
exon 2:            xxxxxxxxxxxxxxxxxx
                      
output:   xxxxxxxxx
                   xxxxxxxx
                           xxxxxxxxxx

USAGE:  python exon_intersection.py

OUTPUT: exon_intersection.out

'''

from collections import defaultdict
from sets import Set
import itertools

def uniq(inlist):    #wybor unikatowych elementow z listy redundantnych elementow
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

def merge_lists(big_list):
    new_list = []
    for ele in big_list:
        new_list+=ele
    return new_list

def format_list(lista):
    h = [item.split('=') for item in lista]
    n = {a:b for a,b in h} 
    return n
        

def readMe(filename):
    fh = open(filename,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    ksel = [(ele[0],ele[3],ele[4],format_list(ele[8].split(';'))) for ele in k]
    
    exon_lista = [(ele[3]['ExonName'],ele[3]['GeneName']) for ele in ksel]
    exon_dic = {a:b for a,b in exon_lista}
    
    ref_dic = defaultdict(list)
    for ele in ksel:
        ref_dic[ele[0]].append((ele[1],ele[2],ele[3]['ExonName'],ele[3]['GeneName']))

    return ref_dic

def sort_columns(lista):
    My_table=sorted(lista, key=lambda a: (a[0],a[1]))
    return My_table


class EXON:
    
    def __init__(self,lista):
        self.info = lista
        self.merged = {}
        self.ranges = {}
        self.pairs = []
        self.uniq_list = []
        self.singles_list = []
        self.singles_exons = []
        self.intersection = []
        
        self.merge_exons(self.info)
        self.get_pairs(self.merged)
        self.group_overlapping_exons(self.uniq_list, self.pairs, self.merged, self.ranges)
        
        
        
    def merge_exons(self, hit_list):
        d = defaultdict(list)
        for start,stop,exon,gene in hit_list:
            d[exon].append((start,stop))
        r = {}
        for ex in d.keys():
            r[ex] = min([int(ele[0]) for ele in d[ex]]),max([int(ele[1]) for ele in d[ex]])
        self.merged = r
        
    def get_pairs(self, exon_dict):
        exon_range = {}
        overlapping_pairs = []

        for exon in exon_dict.keys():
            exon_range[exon] = range(exon_dict[exon][0],exon_dict[exon][1]+1)
        exon_ids = exon_range.keys()
        pary = itertools.combinations(exon_ids, 2)
        p = list(pary)
        
        for dwojka in p:
            first_set = Set(exon_range[dwojka[0]])
            second_set = Set(exon_range[dwojka[1]])
            common_part = first_set & second_set
            if len(common_part):
                overlapping_pairs.append(dwojka)
        u = []        
        for ele in overlapping_pairs:
            u+=ele
        uniq_overlapping_exons = uniq(u)
        singles = [ele for ele in exon_ids if ele not in uniq_overlapping_exons]
        
        s = []
        for name in singles:
            s.append((name,exon_dict[name][0],exon_dict[name][1]))
                
        self.ranges = exon_range
        self.pairs = overlapping_pairs
        self.uniq_list = uniq_overlapping_exons
        self.singles_list = singles
        self.singles_exons = s
        
    def group_overlapping_exons(self, uniq_overlapping_list, pairs_of_overlapping_exons, merged, ranges):
        N = {}
        for exon in uniq_overlapping_list:
            n = merge_lists([ele for ele in pairs_of_overlapping_exons if (ele[0] == exon) or (ele[1] == exon)])
            N[exon] = sorted(uniq(n))
            
        X = {}    
        for klucz in N.keys():
            m = merge_lists([ele for ele in N.values() if klucz in ele])
            X[klucz] = sorted(uniq(m))
          
        N_uniq = uniq(X.values())
              
        intersection = []

        if N_uniq == []:
            self.intersection = []
        else:
            
            for grupa in N_uniq:
               
                overlap_pairs = []
                # put exons in pairs
                p = itertools.combinations(grupa, 2)
                pary = list(p)
                
                # take pairs that overlap and their intersection
                for dwojka in pary:
                    first_set = Set(ranges[dwojka[0]])
                    second_set = Set(ranges[dwojka[1]])
                    common_part = first_set & second_set
                    if len(common_part):
                        
                        overlap_pairs.append(('|'.join(dwojka),list(common_part)))
                
                #remove redundant exon pairs
                op = []
                z = {}
                z_filtr = {}
                for a,b in overlap_pairs:
                    z[a] = min(b), max(b)
                zakresy_uniq = uniq(z.values())
                for ele in zakresy_uniq:
                    nazwy = [n for n in z.keys() if z[n] == ele]
                    nazwy_razem = '|'.join(nazwy)
                    z_filtr[nazwy_razem] = ele
                for exon in z_filtr.keys():
                    start = z_filtr[exon][0]
                    stop = z_filtr[exon][1] + 1
                    op.append((exon,range(start,stop)))
                
                argu = 000
                #self.intersection = op
                exon_ranges = merge_lists([ele[1] for ele in op])
                if len(uniq(exon_ranges)) == len(exon_ranges):
                    
                    i = [(ele[0],min(ele[1]),max(ele[1])) for ele in op]
                    intersection.append(i)
                    #self.intersection = op
                    argu = 111
                    
                else:
                    while argu == 000:
                        # put in pairs that have intersection
                        new_op = []
                    
                        p = itertools.combinations(op, 2)
                        pary = list(p)
                
                        # take pairs that overlap and their intersection
                        for dwojka in pary:
                            first_set = Set(dwojka[0][1])
                            second_set = Set(dwojka[1][1])
                            common_part = first_set & second_set
                            if len(common_part):
                                new_op.append(('|'.join([dwojka[0][0],dwojka[1][0]]),list(common_part)))
                    
                        #self.intersection = new_op
                    
                    
                        # remove redundant pairs
                        red_op = []
                        z = {}
                        z_filtr = {}
                        for a,b in new_op:
                            z[a] = min(b), max(b)
                        zakresy_uniq = uniq(z.values())
                        for ele in zakresy_uniq:
                            nazwy = [n for n in z.keys() if z[n] == ele]
                            nazwy_razem = '|'.join(nazwy)
                            z_filtr[nazwy_razem] = ele
                        for exon in z_filtr.keys():
                            start = z_filtr[exon][0]
                            stop = z_filtr[exon][1] + 1
                            red_op.append((exon,range(start,stop)))
                    
                    
                    
                        red = [('|'.join(uniq(ele[0].split('|'))), ele[1]) for ele in red_op]
                        
                        exon_ranges = merge_lists([ele[1] for ele in red])
                        op = red
                        i = [(ele[0],min(ele[1]),max(ele[1])) for ele in op]
                        if len(uniq(exon_ranges)) < len(exon_ranges):
                            argu = 000
                        else:
                            argu = 111
                        
                           
                    intersection.append(i)
                    
        #i = [(ele[0], min(ele[1]), max(ele[1])) for ele in intersection]
        inter = merge_lists(intersection)    
        self.intersection = uniq(inter)


     
if __name__ == '__main__':
    
    fname = 'exons_alignment_by_blast_out_global.gff3'
    #fname = 'sample_tes1.gff3'
    
    ref_dict = readMe(fname)
    R = []
    I = []
    for ref in ref_dict.keys():
        Sing = []
        Over = []
        exon = EXON(ref_dict[ref])
        
        gene_name = ref_dict[ref][0][3]
        
        Sing.append(exon.singles_exons)
        Over.append(exon.intersection)
        all_exons = merge_lists(Sing + Over)
        I.append(all_exons)
        for ele in all_exons:
            r = ref, ele[1], ele[2], ele[0], gene_name
            R.append(r)
      
    Rsorted = sort_columns(R)
    
    W = [(ele[0], str(ele[1]), str(ele[2]), '_'.join(ele[3:])) for ele in Rsorted]
            
    wh = open('exon_intersection.out','w')
    for hit in W:
        wh.write('\t'.join(hit)+'\n')
        
    wh.flush()
    wh.close()
    
    #print Rsorted
    
