#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt


PROG = 'QDMLL'
PATH = "/net/granat/users/meinecke/QDMLL/GA"

### define hyperparameters
n_pop = 100
n_max_gen = 500
n_elits = 4
n_challengers = 8


#PG_bounds = [0.005, 0.3]
JG_bounds = [1e8,3e9]
U_bounds = [4.0,10.0]
rL_bounds = [0.03,0.95]
rR_bounds = [0.03,0.95]
lA_bounds = [0.01,0.3]
aL_bounds = [1,15]
inh_bounds = [0.010,0.080]
#ai_bounds = [1.5,10.0]

input_bounds = np.array([JG_bounds, U_bounds, rL_bounds, rR_bounds, lA_bounds, aL_bounds, inh_bounds])
n_genes = input_bounds.shape[0]


p_cross = 0.8
p_mut = 1.0/(1.0*n_genes)
p_mut_uni = 0.5
p_mut_ch = 1.0/n_genes
k_tourn = 3

verb = True

### random generator
rng = np.random.default_rng(817)

### set up population and scores matrices to keep track of the evolution
pop_evoMat = np.zeros((n_max_gen,n_pop,n_genes))
scores = np.zeros((n_max_gen,n_pop))



### read utils for output files generated by the c++ program
def get_float_from_datamatrix(data,floatname):
  for k in range(data.shape[0]):
    #print(data[k][0])
    if data[k][0] == floatname:
      return float(data[k][1])
  return np.NAN

def get_float_from_file(filename,floatname):
  data = np.genfromtxt(filename, dtype=None, delimiter=r': ', encoding='UTF-8')
  return get_float_from_datamatrix(data,floatname)



### calc scores via condor parallelization
import subprocess
import os
import time
def calc_scores(pop):
    
    ### condor paras
    MEM = 0
    #PROG = 'QDMLL'
    #PATH = "/net/granat/users/meinecke/QDMLL/GA"
    if( not os.path.exists(PATH)  ):
      print(" Error: directory " + PATH + " does not exists" )
    

    argfile = open('GA_pyargfile','w')
    
    for k in range(pop.shape[0]):
        argfile.write('-fitness -outTime 10000 -intTime 60000')
        argfile.write(' -noNoise')
        #argfile.write(" -loadHist -histFile " + PATH + "/Xhist.bin" )
        argfile.write(' -filename individual_' + "{:03d}".format(k+1))
        argfile.write(' -J_G ' + "{:1.6f}".format(pop[k,0]))
        argfile.write(' -U ' + "{:1.6f}".format(pop[k,1]))
        argfile.write(' -r_L ' + "{:1.6f}".format(pop[k,2]))
        argfile.write(' -r_R ' + "{:1.6f}".format(pop[k,3]))
        argfile.write(' -l_A ' + "{:1.6f}".format(pop[k,4]))
        argfile.write(' -a_L ' + "{:1.6f}".format(int(pop[k,5])))
        argfile.write(' -QD_inh ' + "{:1.6f}".format(pop[k,6]))
        argfile.write("\n")
    argfile.close()
        
    ### condor submission - do not collect output
    substr = "qsub -mem " + str(MEM) + " -m n -speed 4 -w " + PATH + " -argfile GA_pyargfile " + PROG
    #print("Submit string: ")
    #print(substr)
    process = subprocess.Popen(substr, shell=True, stdout=subprocess.PIPE)
    #process = subprocess.Popen(substr, shell=True)
    process.wait()
    
    ### check every few second whether jobs have finished
    jobs_completed = False
    while not jobs_completed:
      time.sleep(5.0)
      process = subprocess.Popen('qstat | grep meinecke', shell=True, stdout=subprocess.PIPE)
      process.wait()
      output = str(process.stdout.read())
      out = output.strip("b'")
      #print(out)
      if out == '':
        jobs_completed = True
    
    ### delete condor job files
    subprocess.Popen('rm ' + PATH +'/job*', shell=True, stdout=subprocess.PIPE)
    
    ### read fitness scores from output files
    scrs = np.zeros(pop.shape[0])
    for k in range(pop.shape[0]): 
      scrs[k] = get_float_from_file(PATH + "/individual_" + "{:03d}".format(k+1), 'fitness')
    
    #plt.plot(scrs)
    #plt.show()
    
    return scrs
  
  
  
### genetic algorithm subroutines

def select_tournament(pop, scores):
    rinds = rng.integers(n_pop,size=k_tourn)
    t_scores = scores[rinds]
    rind_winner = np.argmax(t_scores)
    
    return pop[rinds[rind_winner]]


def crossover(p1,p2,p_cross):
    if rng.random() < p_cross:
        cross_ind = rng.integers(1,n_genes)

        tmp = p1[cross_ind:]
        p1[cross_ind:] = p2[cross_ind:]
        p2[cross_ind:] = tmp

def mutate_uniform(cr,p_mut):
    rn = rng.random(size=n_genes)
    for n in range(n_genes):
        if rn[n] < p_mut:
            cr[n] = rng.uniform(low=input_bounds[n,0], high=input_bounds[n,1])
            
def mutate_centered_normal(cr,p_mut,rel_scale=0.2):
    rn = rng.random(size=n_genes)
    for n in range(n_genes):
        if rn[n] < p_mut:
            gmin = input_bounds[n,0]
            gmax = input_bounds[n,1]
            #cr[n] = np.clip(rng.normal(loc=cr[n], scale=(gmax-gmin)*rel_scale), a_min = gmin, a_max = gmax)
            candidate = rng.normal(loc=cr[n], scale=(gmax-gmin)*rel_scale)
            while (candidate < gmin or candidate > gmax):
              candidate = rng.normal(loc=cr[n], scale=(gmax-gmin)*rel_scale)
            cr[n] = candidate
            

def mutate(cr, p_mut=0.1, p_mut_uni=0.5, rel_scale=0.2):
        if(rng.random() < p_mut_uni):
          mutate_uniform(pop_evoMat[k+1,m], p_mut)
        else:
          mutate_centered_normal(pop_evoMat[k+1,m], p_mut, rel_scale)


def challenger(pop_elits, p_mut_ch=0.2, rel_scale=0.2):
    #pick elite
    cr = pop_elits[rng.integers(n_elits)]
    rn = rng.random(size=n_genes)
    for n in range(n_genes):
        if rn[n] < p_mut_ch:
          gmin = input_bounds[n,0]
          gmax = input_bounds[n,1]
          #cr[n] = np.clip(rng.normal(loc=cr[n], scale=(gmax-gmin)*rel_scale), a_min = gmin, a_max = gmax)
          candidate = rng.normal(loc=cr[n], scale=(gmax-gmin)*rel_scale)
          while (candidate < gmin or candidate > gmax):
            candidate = rng.normal(loc=cr[n], scale=(gmax-gmin)*rel_scale)
          cr[n] = candidate
    return cr



### initialize first generation
pop_evoMat[0] = rng.uniform(low=input_bounds.T[0], high=input_bounds.T[1], size=(n_pop,n_genes))


### evolve
for k in range(n_max_gen-1):
    
    ### calculate scores
    scores[k] = calc_scores(pop_evoMat[k])
    if verb:
        print('generation: ' + str(k))
        print('mean score: ' + str(np.mean(scores[k])))
        bestInd = np.argmax(scores[k])
        print('best score: ' + str(scores[k,bestInd]))
        print('best individual: ' + str( pop_evoMat[k,bestInd]))
    
    ### sort
    sorted_ind = np.flip(np.argsort(scores[k]))
    scores[k] = scores[k,sorted_ind]
    pop_evoMat[k] = pop_evoMat[k,sorted_ind]
    #print(scores[k])
    #plt.plot(scores[k])
    #plt.show()
    
    ### elitism
    for m in range(n_elits):
        pop_evoMat[k+1,m] = pop_evoMat[k,m]
    
    ### challengers
    rel_scale_ch = 0.1*(n_max_gen-k)/n_max_gen + 0.01*k/n_max_gen
    for m in range(n_elits, n_elits + n_challengers):
        pop_evoMat[k+1,m] = challenger(pop_evoMat[k,:n_elits], p_mut_ch, rel_scale = rel_scale_ch)
    
    ### selection
    for m in range(n_elits + n_challengers,n_pop):
        pop_evoMat[k+1,m] = select_tournament(pop_evoMat[k],scores[k])
    
    ### crossover
    for m in range(n_elits + n_challengers,n_pop-1,2):
        crossover(pop_evoMat[k+1,m],pop_evoMat[k+1,m+1],p_cross)
    
    ### mutation
    rel_scale_mut = 1.0*(n_max_gen-k)/n_max_gen + 0.1*k/n_max_gen
    for m in range(n_elits + n_challengers,n_pop):
        mutate(pop_evoMat[k+1,m], p_mut, p_mut_uni = p_mut_uni, rel_scale = rel_scale_mut)
    
    
    np.savetxt(PATH + '/pop_k'+ "{:03d}".format(k+1),pop_evoMat[k])
    np.savetxt(PATH + '/scores_k'+ "{:03d}".format(k+1),scores[k])
    
    
### evaluate last generation
scores[-1] = calc_scores(pop_evoMat[-1])
sorted_ind = np.flip(np.argsort(scores[-1]))
scores[-1] = scores[-1,sorted_ind]
pop_evoMat[-1] = pop_evoMat[-1,sorted_ind]

np.savetxt(PATH + '/pop_k'+ "{:03d}".format(n_max_gen),pop_evoMat[n_max_gen-1])
np.savetxt(PATH + '/scores_k'+ "{:03d}".format(n_max_gen),scores[n_max_gen-1])


print('final best score: ' + str(scores[-1,0]))
print('final best individual: ' + str(pop_evoMat[-1,0]))


### visualize results
mean = np.mean(scores, axis=1)
std = np.std(scores, axis=1)

### plot results
plt.plot(np.arange(n_max_gen),mean)
plt.fill_between(np.arange(n_max_gen), mean-std, mean+std, alpha=0.3)
plt.plot(np.amax(scores, axis=1))
plt.show()

plt.imshow(scores, aspect='auto')
plt.colorbar()
plt.show()
