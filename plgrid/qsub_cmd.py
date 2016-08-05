'''
Created on Dec 29, 2014

@author: ania.fijarczyk
'''

import os,glob,sys

sciezka = os.getcwd()
lista = glob.glob('qsub_*.sh')
katalog = sciezka.split('/')[-1]

if __name__ == '__main__':
    
    wh = open('kolejkuj_'+katalog+'.sh','w')
    for skrypt in sorted(lista):
        wh.write('qsub -A newt001 '+sciezka+'/'+skrypt+'\n')
	#wh.write('sleep 1\n')
    wh.flush()
    wh.close()
