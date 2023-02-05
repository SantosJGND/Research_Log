from scipy.stats import binom
import numpy as np
import itertools as it
import math

def nCr(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)


def get_grid(P= 20,range_height= [0,1],range_length= [0,1]):

    traces= [x for x in it.product(range(P),range(P))]

    i_coords, j_coords = np.meshgrid(np.linspace(range_height[0],range_height[1],P),
                          np.linspace(range_length[0],range_length[1],P),indexing= 'ij')

    traces= [x for x in it.product(range(P),range(P))]

    i_coords, j_coords = np.meshgrid(np.linspace(range_height[0],range_height[1],P),
                          np.linspace(range_length[0],range_length[1],P),indexing= 'ij')

    background= np.array([i_coords, j_coords])

    background= [background[:,c[0],c[1]] for c in traces]
    background=np.array(background)
    
    return background,i_coords, j_coords


def single_gen_matrix(Ne= 1092, ploidy= 2,precision= 1092,s=0):
    '''
    define transition probabilities for alleles at any given frequency from 1 to Ne.
    '''
    pop_size= Ne * ploidy
    space= np.linspace(0,pop_size,precision,dtype=int)
    t= np.tile(np.array(space),(precision,1))
    
    probs= [x / (pop_size) for x in space]
    
    ## 
    sadjust= [x * (1+s) for x in probs]
    scomp= [(1-x) for x in probs]
    
    new_probs= [sadjust[x] / (scomp[x]+sadjust[x]) for x in range(precision)]
    new_probs= np.nan_to_num(new_probs)
    
    probs= new_probs
    
    freq_matrix= binom.pmf(t.T,pop_size,probs)
    
    return freq_matrix


def freq_progression(fr= 0,n_gens= 20,freq_matrix= {},remove_tails= False):
    '''frequency distribution after n generations given initial freq'''
    fixed_tally= []
    if isinstance(fr,int):
        freq_ar= [0] * freq_matrix.shape[0]
        freq_ar[fr]= 1
    else:
        freq_ar= fr
    
    for idx in range(n_gens):
        
        freq_ar= np.array(freq_ar) @ freq_matrix.T
        freq_ar= freq_ar.reshape(1,-1)
        
        prop_fixed= sum([freq_ar[0,0],freq_ar[0,-1]])
        prop_fixed= prop_fixed / np.sum(freq_ar,axis= 1)
        fixed_tally.append(prop_fixed[0])
        
        if remove_tails:
            freq_ar[0,0]= 0
            freq_ar[0,-1]= 0
        
        #print(freq_ar.shape)
            
        freq_ar= freq_ar.T / np.sum(freq_ar,axis= 1)
        
        freq_ar= freq_ar.T
    
    return freq_ar, fixed_tally


###############
###############

from scipy.stats import binom

def single_gen_matrix_v2(Ne1= 1092,Ne2= 1092,prec1= 1092,prec2= 1092,ploidy= 2,s=0):
    '''
    define transition probabilities for alleles at any given frequency from 1 to Ne.
    '''
    pop_size1= Ne1 * ploidy
    pop_size2= Ne2 * ploidy
    
    space= np.linspace(0,pop_size2,prec2,dtype=int)
    t= np.tile(np.array(space),(prec1,1))
    
    probs= np.linspace(0,pop_size1,prec1,dtype=int)
    probs= [x / (pop_size1) for x in probs]
    
    ## 
    sadjust= [x * (1+s) for x in probs]
    scomp= [(1-x) for x in probs]
    
    new_probs= [sadjust[x] / (scomp[x]+sadjust[x]) for x in range(len(probs))]
    new_probs= np.nan_to_num(new_probs)
    
    freq_matrix= binom.pmf(t.T,pop_size2,new_probs)
    
    return freq_matrix



def freq_progr_func(theta_dict,Ne= 1000,Ne0=1000,T= 2000,ploidy= 2, s= 0,remove_tails= True, fr= 1,
    return_spec= False, fixed_tally= 1):
    """
    model frequency evolution using Marvkov model predicated on demographic data.
    """
    theta_func= theta_dict['func']
    theta_args= theta_dict['kargs']

    Ne_t= np.array(range(T)) + 1
    Ne_t= theta_func(Ne_t,Ne=Ne,**theta_args)
    Ne_t= np.array(Ne_t,dtype= int)
    Ne_process= int(Ne0)
    
    #####
    trans_matrix= single_gen_matrix_v2(Ne1= Ne,prec1=Ne,
                                       Ne2=Ne,prec2=Ne,ploidy= ploidy,s=s)
    
    #####
    if isinstance(fr,int):
        freq_ar= [0] * trans_matrix.shape[0]
        freq_ar[fr]= 1
        freq_ar= np.array(freq_ar).reshape(-1,1).T
    else:
        freq_ar= fr
    
    if isinstance(fixed_tally,int)==1:
        fixed_tally= np.ones((1,1))
    
    array_spec= []

    for idx in range(T):

        Ne_now= Ne_t[idx]
        
        if Ne_now != Ne_process:
            trans_matrix= single_gen_matrix_v2(Ne1= Ne_process,prec1=Ne_process,
                                              Ne2=Ne_now,prec2=Ne_now,ploidy= ploidy,s=s)
            
            Ne_process= Ne_now
        
        if freq_ar.shape[1] != trans_matrix.shape[1]:
            old_shape= freq_ar.shape[1]
            trans_matrix= single_gen_matrix_v2(Ne1= Ne_process,prec1=old_shape,
                                              Ne2=Ne_now,prec2=Ne_now,ploidy= ploidy,s=s)
        
        freq_ar= freq_ar @ trans_matrix.T
               
        ### proportion of fixed
        prop_fixed= freq_ar[:,[0,-1]]
        prop_fixed= np.sum(prop_fixed,axis= 1)
        prop_fixed= prop_fixed.T / np.sum(freq_ar,axis= 1)
        ## update segregating tally
        fixed_tally= fixed_tally * (1-prop_fixed.reshape(-1,1))
        new_tally= np.ones((1,1))
        fixed_tally= np.concatenate((fixed_tally,new_tally))
        
        if remove_tails:
            freq_ar[:,0]= 0
            freq_ar[:,-1]= 0
        
        freq_ar= freq_ar.T / np.sum(freq_ar,axis= 1)

        freq_ar= freq_ar.T
        
        new_freq_ar= [0] * trans_matrix.shape[0]
        new_freq_ar[1]= 1
        new_freq_ar= np.array(new_freq_ar).reshape(-1,1).T
        freq_ar= np.concatenate((freq_ar,new_freq_ar))
        
        if return_spec:
            array_spec.append(freq_ar)
       
    #
    
    return freq_ar, fixed_tally, array_spec



def get_fixedtally_v2(tally_array, total_prop= 1):
    '''
    get total propotion of fixed alleles.
    '''
    
    fixed= []
    for idx in range(len(tally_array)):

        fixed.append(tally_array[idx] * total_prop)
        
        total_prop= total_prop * (1 - tally_array[idx])
    
    fixed_array= [sum(fixed[x+1:]) for x in range(len(fixed)-1)] + [fixed[-1]]
    total_fixed= sum(fixed)
    return total_fixed, fixed_array

