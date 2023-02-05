import numpy as np
import itertools as it

import allel
import pandas as pd
import os
from datetime import datetime

import collections
def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


from tools.sfs_utilities import(
	single_gen_matrix, freq_progression
	)

#############


def get_fixedtally(tally_array, total_prop= 1):
    '''
    get total propotion of fixed alleles.
    '''
    
    survived= []
    for idx in range(len(tally_array)):

        survived.append(tally_array[idx] * total_prop)
        
        total_prop= total_prop * (1 - tally_array[idx])
    
    survived= sum(survived)
    return total_prop, survived

#############
fig_dir= 'Figures'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'

#############
#############
model_dict= {
    'test': {
        'N': 1000,
        'T': 5000,
        'scale': 1
    },
    'schweinfurthii': {
        'N':10528, 
        'T':17925,
        'scale': 1
    },
    'troglodytes': {
        'N': 72001, 
        'T': 17925,
        'scale': 1 / 7
        },
    'ellioti': {
        'N': 6033,
        'T': 33402,
        'scale': 1 
    },
    'verus': {
        'N': 5710,
        'T': 33402,
        'scale': 1 
    }
}

############
############
scale_sim= True

if scale_sim:
    fig_dir= fig_dir + 'scaled_'

muGlobal= 1.08e-8

seqL= 1e6

nbins= 50

Nreps= 20
sample_range= [2,1000]
Nsteps= 50
range_use= np.linspace(*sample_range,Nsteps,dtype= int)

##########
##########

instance_dict= {}

for pop_select in model_dict.keys():
	print('going on: {}'.format(pop_select))
	Ne= model_dict[pop_select]["N"]

	if scale_sim:
	    scale_gen= model_dict[pop_select]['scale']
	    fig_dir= fig_dir + 'scaled_'
	    muNuc= muGlobal / scale_gen
	    Ne= int(Ne * scale_gen)

	else: muNuc= muGlobal

	mu= muNuc * seqL
	Theta= 4 * Ne * mu
	rate_Pmut= Theta / 2

	branch_len= model_dict[pop_select]['T']
	MRCA= branch_len / Ne / 2

	Pexp= 2 * Ne * MRCA * Theta / 2 ## Multiply the standard poisson rate by 2Ne to do all the pop.

	muts= np.random.poisson(Pexp,1)[0]
	bins= np.linspace(0,MRCA,nbins)
	bins= np.round(bins,4)
	bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]
	t_list= [sum(x)/2 for x in bins]

	gens= [x * branch_len / MRCA for x in t_list]
	gens= np.array(gens,dtype= int)

	precision= 2 * Ne
	sel_coeff= 0
	remove_tails= True
	multiplier= precision / Ne

	freqi= 1

	freq_matrix= single_gen_matrix(Ne= Ne,precision= precision,s=sel_coeff)

	freq_prog= [freq_progression(freqi,n_gens= x, freq_matrix= freq_matrix,remove_tails= remove_tails) for x in gens]
	fixed_tally= [x[1] for x in freq_prog]
	freq_prog= [x[0] for x in freq_prog]


	########## plot SFS per time bin
	surface= 'N'
	surf= np.linspace(0,1,freq_matrix.shape[0])
	surf_dict= {
	    'freq': surf,
	    'N': [int(x * Ne) for x in surf]
	}

	plt.figure(figsize=(10,10))

	for x in list(range(len(gens))):
		plt.plot(surf_dict[surface],freq_prog[x][0])

	plt.xlabel('freq',fontsize=20)
	plt.ylabel('probability',fontsize=20)
	plt.xticks(fontsize= 15)
	plt.yticks(fontsize= 15)
	plt.title(pop_select,fontsize=20)

	plt.legend()
	plt.savefig(fig_dir + '{}_SFSbins_Nsamp{}.png'.format(pop_select,sample_range[1]),bbox_inches='tight')
	plt.close()

	##########

	tallies= [get_fixedtally(x) for x in fixed_tally]

	mut_bin= int(muts / (nbins - 1))
	nseg_tal= [x[0] * mut_bin for x in tallies]
	nseg_tal= np.array(nseg_tal,dtype= int)


	########## Overall SFS

	stand_num= np.array(nseg_tal) / sum(nseg_tal)
	seg_total_matrix= np.array([x[0] for x in freq_prog])
	seg_total_matrix= stand_num.reshape(-1,1) * seg_total_matrix
	seg_total_matrix= np.sum(seg_total_matrix,axis= 0)
	seg_total_matrix= seg_total_matrix / np.sum(seg_total_matrix)

	########## plot overall SFS

	plt.figure(figsize=(10,10))
	plt.plot(surf_dict[surface],seg_total_matrix)

	plt.xlabel('freq',fontsize=20)
	plt.ylabel('probability',fontsize=20)
	plt.xticks(fontsize= 15)
	plt.yticks(fontsize= 15)
	plt.title(pop_select,fontsize=20)

	plt.legend()
	plt.savefig(fig_dir + '{}_SFStotal_Nsamp{}.png'.format(pop_select,sample_range[1]),bbox_inches='tight')
	plt.close()

	######

	freq_surface= np.linspace(0,1,freq_matrix.shape[0])

	## sampling parameters

	ksamp_dict= {
	    x: [] for x in range_use
	}

	for si in range_use:
	    
	    for proxy in range(Nreps):
	        
	        ###
	        ###
	        count_bin= 0
	        for idx in range(nbins-1):
	            probs= np.random.choice(freq_surface,nseg_tal[idx], p= freq_prog[idx][0])
	            binary= [np.random.choice([0,1],si,p= [1-x,x]) for x in probs]
	            binary= np.array(binary).T
	            
	            if not len(binary):
	            	continue

	            binary= np.sum(binary,axis= 0)
	            binary[binary > 0]= 1
	            
	            binary= sum(binary)
	            count_bin += binary
	        
	        ksamp_dict[si].append(count_bin)

	###############
	###############

	stats_dict= {
	    z: {
	        'mean': np.mean(g),
	        'std': np.std(g)
	    } for z,g in ksamp_dict.items()
	}

	samp_order= sorted(ksamp_dict.keys())


	plt.figure(figsize=(10,10))

	global_y= [stats_dict[x]['mean'] for x in samp_order]

	global_error= [stats_dict[x]['std'] for x in samp_order]

	plt.errorbar(samp_order,global_y,yerr=global_error)

	plt.xlabel('samp size',fontsize=20)
	plt.ylabel('PA seg #',fontsize=20)
	plt.xticks(fontsize= 15)
	plt.yticks(fontsize= 15)
	plt.title(pop_select,fontsize=20)

	plt.legend()
	plt.savefig(fig_dir + '{}_sampN{}.png'.format(pop_select,sample_range[1]),bbox_inches='tight')
	plt.close()

	instance_dict[pop_select]= {
		'mean': global_y,
		'std': global_error
	}


plt.figure(figsize=(10,10))
for pop_select in model_dict.keys():
	if pop_select == 'test':
		continue


	global_y= instance_dict[pop_select]['mean']

	global_error= instance_dict[pop_select]['std']

	plt.errorbar(samp_order,global_y,yerr=global_error,label= pop_select)

plt.xlabel('samp size',fontsize=20)
plt.ylabel('PA seg #',fontsize=20)
plt.xticks(fontsize= 15)
plt.yticks(fontsize= 15)
plt.title(pop_select,fontsize=20)

plt.legend()
plt.savefig(fig_dir + 'combined_sampN.png',bbox_inches='tight')
plt.close()

instance_dict[pop_select]= {
	'mean': global_y,
	'std': global_error
}



