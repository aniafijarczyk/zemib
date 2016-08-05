'''
Created on Dec 29, 2014

@author: ania.fijarczyk
'''

import os,glob,sys,random

from optparse import OptionParser


random_list = random.sample(range(10000,1000000),2000)

sciezka = os.getcwd()

ABC_summary_prog = 'ABC_summary_one_locus.sh'
ABCtoolbox_prog = 'ABCtoolbox'
ArpToFatsa_prog = 'ArpToFasta_one_fast.py'
est_file = 'exampleDNA_1.est'
input_file = 'exampleDNA.input'
par_file = 'exampleDNA.par'
fastsimcoal_prog = 'fastsimcoal21'
stats_file = 'Lm_Lvg_SS.obs'
stats_prog = 'stats_counter_one_faster.py'
mstatspop_prog = 'mstatspop'

lista = glob.glob('aniaf*')

def get_queue(que_int):

    if que_int == 0:
        return "plgrid"
    elif que_int == 1:
        return "plgrid-long"
    elif que_int == 2:
        return "l_short"
    elif que_int == 3:
        return "l_long"
    elif que_int == 4:
        return "l_prio"
    elif que_int == 5:
        return "l_test"
    else:
        return "plgrid"

    
def main(args=[]):
    
   
    usage = "usage: %prog [options] arg \nProgram generates qsub shells for ABC simulations (for abc_gwh workflow)"
    parser = OptionParser(usage, version='%prog version 1.0')
    parser.add_option("-q", "--queue_type", dest="QUEUE_TYPE", help="queue type", default= 0)
    
    (options, arg) = parser.parse_args(args)
    
####################### entering program ###########################
    
    QUEUE_TYPE = int(options.QUEUE_TYPE)
    
    i = 0
    for plik in lista:

        plgrid_folder_name = plik+'_plgrid'
    
        wh = open('qsub_'+plik+'.sh','w')
        W = []
        W.append('#!/bin/sh\n')
        W.append('#PBS -q '+get_queue(QUEUE_TYPE)+'\n\n')
  	
        W.append('cd $PLG_USER_SCRATCH_SHARED\n')
        W.append('mkdir -p '+plgrid_folder_name+'\n')
        W.append('cp -rp '+sciezka+'/'+plik+'/* ./'+plgrid_folder_name+'\n')
        W.append('cd '+plgrid_folder_name+'\n')
    
        W.append('chmod +x ./'+ABC_summary_prog+'\n')
        W.append('chmod +x ./'+ABCtoolbox_prog+'\n')
        W.append('chmod +x ./'+fastsimcoal_prog+'\n')
        W.append('chmod +x ./'+mstatspop_prog+'\n')    
    
        W.append('./'+ABCtoolbox_prog+' ./'+input_file+' addToSeed='+str(random_list[i])+' verbose\n')
	W.append('cd ../\n')
    
        W.append('mkdir -p $PLG_USER_STORAGE/ABC_results\n')
        W.append('mkdir $PLG_USER_STORAGE/ABC_results/'+plik+'_sim\n')
        W.append('cp -rp ./'+plgrid_folder_name+'/example_output_sampling1.txt $PLG_USER_STORAGE/ABC_results/'+plik+'_sim\n')
        #W.append('cp -rp ./'+plgrid_folder_name+'/seed.txt $PLG_USER_STORAGE/ABC_results/'+plik+'_sim\n')
        W.append('rm -rf $PLG_USER_SCRATCH_SHARED/'+plgrid_folder_name+'\n')
 	
	i+=1
   
        w = ''.join(W)
        wh.write(w)
        wh.flush()
        wh.close()
        #print w 
  
  
    
if __name__ == '__main__':
    main(sys.argv[1:])

    
    
    
