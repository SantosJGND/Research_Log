import numpy as np
#######
#######
def lineAssign(new_array,mut_idx,nmuts = 192):
    new_array= [0]*nmuts
    
    for idx in range(len(new_array)):
        nmut= mut_idx[idx]
        new_array[nmut] += new_array[idx]
        
    return new_array


####### v1
def count_popKmers_bench(Window, mut_matrix, pop_dict, single= True, frequency_range= [0,1],row=32,col=3,segregating= False,
                        scale= 1, prop_gen_used= 1, return_private= False,PA= {},pop_tag= '_ss', counted= False,
                        mut_list= [],nmuts= 192):
    '''
    Extract population mutation counts from _ind x kmer_ mutation matrix. 
    '''
    pop_counts_matrix= {}
    pop_counts_lin= {}
    times_dict= {}
    pop_list= list(pop_dict.keys())
    
    ncol= Window.shape[1]
    if not counted:
        ncol= mut_matrix.shape[0]

    array_counts= np.zeros((Window.shape[0],ncol))

    for pop in pop_list:
        t0= time.time()
        pop_ori= pop
        if pop_tag in pop:
            pop_ori= pop[len(pop_tag):].split('.')[0]
        ####
        #### filtering part
        klist= pop_dict[pop]
        pop_gen= Window[klist,:]

        if PA: 
            shared= [x for x in range(pop_gen.shape[1]) if PA[pop_ori][x] == 0]
            pop_gen[:,shared] = 0
        else:
            freqs= np.sum(pop_gen,axis= 0) / pop_gen.shape[0]
            ## discount alleles outside freq range.
            in_out= (freqs <= frequency_range[0]) | (freqs >= frequency_range[1])

        if single: 
            pop_gen= np.sum(pop_gen,axis= 0) > 0
            pop_gen= np.array(pop_gen,dtype= int).reshape(1,len(pop_gen))

        pop_seg_ori= np.sum(pop_gen,axis= 0) > 0
        pop_seg_ori= np.array(pop_seg_ori,dtype= int).reshape(1,len(pop_seg_ori))
        pop_seg_ori= pop_seg_ori * scale * prop_gen_used
        
        if not PA:
            pop_gen[:,in_out]= 0
        t1= time.time()
        ####
        #### benchmarking

        pop_collapsed_mat= geno_muts_v2(pop_gen, mut_matrix)

        t2= time.time()
        
        pop_coll_II= [lineAssign(list(x),mut_idx,nmuts = 192) for x in pop_gen]
        
        t3= time.time()
        ###
        ###
        #pop_counts_matrix= pop_collapsed_mat
        #pop_counts_lin= pop_coll_II
        
        times_dict= {
            'size': len(klist),
            'filter': t1-t0,
            'mat': t2-t1,
            'lin': t3-t2
        }
        
        ######
    
    
    return pop_counts_matrix, pop_counts_lin, times_dict

############################