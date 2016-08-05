'''
Created on Dec 30, 2014

@author: ania.fijarczyk
'''


import os,sys,glob

lista = glob.glob('M*')

def make_example_file(job,sim):
    
    p = '\
//Inputfile for the program ABCsampler\n\
//-----------------------------------------------------------------------\n\
samplerType standard \n\
//-----------------------------------------------------------------------\n\
estName exampleDNA_1.est \n\
obsName Lm_Lvg_SS.obs\n\
simDataName exampleDNA-temp_1_1.arp\n\
outName example_output\n\
separateOutputFiles 0\n\
addToSeed 0\n\
nbSims '+str(sim)+'\n\
writeHeader 1\n\
simulationProgram fastsimcoal21\n\
simInputName exampleDNA.par \n\
simParam -i#exampleDNA-temp.par#-n#1#-q\n\
launchAfterSim ABC_summary_one_locus.sh\n\
sumStatFile Lm_Lv_output.txt\n'
    
    wh = open('exampleDNA.input','w')
    wh.write(p)
    wh.flush()
    wh.close()


def multiply(n,sims):
    job = int(n)
    x=1
    for model in lista:
        i=1
        while i <= job:
            cmd = 'cp -r '+model+' aniaf_'+model+'_'+str(i)
            os.system(cmd)
            
            # create example file and move it to the new folder
            
            make_example_file(x,sims)
            cmd1 = 'mv exampleDNA.input aniaf_'+model+'_'+str(i)
            os.system(cmd1)
            
            i+=1
            x+=1
            
    

if __name__ == '__main__':
    njobs = sys.argv[1]
    nsims = sys.argv[2]
    #njobs = '2'
    #nsims = '400'
    multiply(njobs,nsims)
