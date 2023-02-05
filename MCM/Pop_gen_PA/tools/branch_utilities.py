
import numpy as np 


###################
################### utility co-factors

def get_edge_dict(tree_summ, directed= True):
    '''
    get dictionary of edges.
    '''
    edge_dict= {}
    for ed in tree_summ['edges']:
        
        for idx in range(1-int(directed)+1):
            
            parent= ed[idx]
            daughter= ed[1-idx]
            if parent in edge_dict.keys():
                edge_dict[parent].append(daughter)
            else:
                edge_dict[parent]= [daughter]

    return edge_dict


def node_assign(tree_demo, tree_summ):
    '''
    '''
    edge_dict= get_edge_dict(tree_summ)
    
    parent_dict= {}
    for z,g in edge_dict.items():
        for nd in g:
            parent_dict[nd]= z
    
    
    new_nodes= {}
    
    for node in tree_demo:
        children= list(node['N'].keys())
        parent= parent_dict[children[0]]
        nnode= {z:g for z,g in node.items() if z != 'node'}
        new_nodes[parent]= nnode
    
    return new_nodes

#####

from tools.ABC_utilities import read_demofile,tree_fill_list,get_tree_nodes

def demo_file_branchProcess(demo_file,anc_r= '0'):

    tree, demo_data= read_demofile(demo_file)

    tree_summ= get_tree_nodes({anc_r:tree},nodes=[],edges= [],leaves= [])
    #
    anc_name= anc_r
    anc_size= demo_data['N'][anc_r]

    times_order= [0] + sorted(demo_data['T'].keys(),reverse= True)
    int_sizes= [x for x in demo_data['N'].keys() if x[0] == 'T']
    int_sizes= sorted(int_sizes)

    tree_demo= tree_fill_list(times_order,tree, demo_data, tree_summ, int_sizes= int_sizes,demo_tape= [])
    ###

    pops= tree_summ['leaves']

    tree_demo= [{z:g for z,g in x.items() if len(g)} for x in tree_demo]
    
    return tree, demo_data, tree_summ, tree_demo



####################
#################### Branch Progress and Recursive

from tools.ABC_utilities import (
    sample_dist_beta, return_replica
)

from tools.sfs_utilities import (
    freq_progr_func, get_fixedtally_v2
)


def branch_progress(theta_dict,Ne= 1000,Ne0= 1000,T= 200,muG= 1.08e-8, ploidy= 2,s= 0,
                    seqL= 1e6, scale_sim= False,model_dict= {},fr= 1,fixed_tally= 1):
    '''
    Calculate frequency evolution on a single branch given provided parameters and theta function.
    '''
    muNuc= muG

    if scale_sim:
        scale_gen= model_dict[pop_select]['scale']
        muNuc= muNuc / scale_gen
        Ne= int(Ne * scale_gen)

    ##### get number of mutations
    mu= muNuc * seqL
    Theta= 4 * Ne * mu
    rate_Pmut= Theta / 2

    MRCA= T / Ne / 2

    Pexp= 2 * Ne * MRCA * Theta / 2 ## Multiply the standard poisson rate by 2Ne to do all the pop.

    Poisson_conv= theta_dict['func'](T,Ne= Pexp,**theta_dict['kargs'])
    
    muts= np.random.poisson(Poisson_conv,1)[0]

    #####
    freq_ar, fixed_tally, array_spec= freq_progr_func(theta_dict,fr= fr,Ne= Ne,Ne0= Ne,
                                          T= T,ploidy= ploidy, s= s, return_spec= False,fixed_tally= fixed_tally)
    
    return freq_ar, array_spec, fixed_tally, muts



def traverse_chose(child,Nec,tree_summ, node_dict,theta_dict,freq_ar,array_spec,fr= 1,Ne= 500,T= 200,muG= 1.08e-8, ploidy= 2,s= 0,
                    seqL= 1e6, scale_sim= False,model_dict= {},sample_func= sample_dist_beta,
                    med_samp= True,fixed_tally= 1):
        
    if child in tree_summ['leaves']:
        
        freq_ar_c, array_spec_c, fixed_c, muts_c= branch_progress(theta_dict,fr=fr,Ne= Nec,Ne0=Ne,T= T,muG= muG, 
                                                    ploidy= ploidy,s= s, seqL= seqL, 
                                                    scale_sim= scale_sim,model_dict= model_dict,fixed_tally= fixed_tally)
        
        fixed_c= fixed_c.reshape(-1,1)

        return {
            'node': child,
            'specs':freq_ar_c,
            'T': T,
            'fixed': fixed_c,
            'muts': muts_c
            }


    else:
        next_t= node_dict[child]['T']

        next_t= T - next_t
        return traverse_sfs(node_dict,tree_summ,theta_dict,node= child,fr= fr,Ne= Ne,T= next_t,
                                   muG= muG, ploidy= ploidy,s= s, sample_func= sample_func,
                                    seqL= seqL, scale_sim= scale_sim,model_dict= model_dict,
                                    med_samp= med_samp,fixed_tally= fixed_tally)


def traverse_intermediate(node,tree_summ, node_dict,theta_dict,freq_ar,array_spec,fr= 1,Ne= 500,T= 200,muG= 1.08e-8, ploidy= 2,s= 0,
                    seqL= 1e6, scale_sim= False,model_dict= {},sample_func= sample_dist_beta,
                    med_samp= True,fixed_tally= 1):
    
    
    return {
        child: traverse_chose(child,int(Nec),tree_summ, node_dict,theta_dict,freq_ar,array_spec,
                                    fr= freq_ar,Ne= int(Ne),T= T, seqL= seqL, 
                                   muG= muG, ploidy= ploidy,s= s, sample_func= sample_func,
                                    scale_sim= scale_sim, model_dict= model_dict,
                                    med_samp= med_samp,fixed_tally= fixed_tally) for child,Nec in node_dict[node]["N"].items()
           }




def get_SFS(Ne= 1000, mu= 1e-8,Nsamp= 100,sinks= True):
    '''
    get expected_SFS
    '''
    Theta= 4 * Ne * mu
    freq_exp= [Theta / x for x in range(1,Nsamp+1)]
    if sinks:
        freq_exp= [0] + freq_exp + [0]
    freq_exp= np.array(freq_exp) / np.sum(freq_exp).reshape(-1,1)
    
    return freq_exp



## GGVE p. 59 eq. 2.30-2.31

from math import factorial
from operator import mul
from fractions import Fraction
from functools import reduce

def nCk(n,k):
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def exact_SN(k,n,Theta):
    
    if k==0:
        dis= [Theta + x for x in range(n)]
        dis= reduce((lambda x, y: x * y), dis)
        
        su= factorial(n - 1) / dis
        
        return su
    
    preff= (n-1) / Theta
    weight= 0
    
    for i in range(1,n):
        
        mu= (-1)**(i-1)
        mu= mu * nCk(n-2,i-1)
        
        mu= mu * (Theta / (i + Theta))**(k+1)
        
        weight += mu
    
    Sk= preff * weight

    return Sk


def get_Nsamps(Nsamp,Theta,range_k=[3,1e4],steps= 100):

    kvec= np.linspace(*range_k,steps)
    probs= [exact_SN(k,Nsamp,Theta) for k in kvec]
    probs= np.array(probs) / sum(probs)
    probs[probs < 0] = 0

    Nk= np.random.choice(kvec,size= 1,p= probs)[0]

    return int(Nk)


def traverse_sfs(node_dict,tree_summ,theta_dict,node= '0',fr= 1,Ne= 500,Ne0=500,T= 200,muG= 1.08e-8, ploidy= 2,s= 0,
                seqL= 1e6, scale_sim= False,model_dict= {},sample_func= sample_dist_beta,
                med_samp= True, node_stat_dict= {},fixed_tally= 1):

    if isinstance(T,str):
        T= 1
        ## expected number of mutations segregating in full population. Average coal time = 2Ne.
        Theta= 4 * Ne0 * muG * seqL
        MRCA= 2
        Pexp= 2 * Ne * MRCA * Theta / 2 
        muts= np.random.poisson(Pexp,1)[0]
        #muts=  get_Nsamps(Ne0,Theta= Theta)
        print(Theta,muts,Ne0)
        fixed_tally= np.ones((1,1))
        array_spec= []
        freq_ar= get_SFS(Ne= Ne0,mu= muG,Nsamp= Ne0)

    else:
        freq_ar, array_spec, fixed_tally, muts= branch_progress(theta_dict,fr=fr,Ne= Ne,Ne0= Ne0,T= T,muG= muG, ploidy= ploidy,s= s,
                                                                scale_sim= scale_sim,model_dict= model_dict)

    New_time= node_dict[node]['T']

    return {
        'node': node,
        'T': T,
        'muts': muts,
        'branch': traverse_intermediate(node,tree_summ, node_dict,theta_dict,freq_ar,array_spec,
                                        fr= freq_ar,Ne= int(Ne),T= New_time, seqL= seqL, 
                                       muG= muG, ploidy= ploidy,s= s, sample_func= sample_func,
                                        scale_sim= scale_sim, model_dict= model_dict,
                                        med_samp= med_samp,fixed_tally= fixed_tally)
    }
 
        
#######################################
#######################################
##### SFS_TRAVERSE PROCESSING #########


def get_sum_stats(node_dict,node_info= []):
    '''
    get single list of sum_stats.
    '''
    node_inst= {
        z:g for z,g in node_dict.items() if z != 'branch'
    }
    node_info.append(node_inst)
    
    if 'branch' in node_dict.keys():
        for br in node_dict['branch'].keys():
            get_sum_stats(node_dict['branch'][br], node_info= node_info)
    

    return node_info

def get_setlist(mut_t,T= 0):
    if not T:
        return []
        
    mut_t= sorted(mut_t)
    mut_now= mut_t[0]
    
    timetable= {idx: 0 for idx in range(T)} 
    for idx in range(len(mut_t)):
        mut= mut_t[idx]

        timetable[mut] += 1
    
    timetable= [[z,g] for z,g in timetable.items()]
    
    return timetable


def node_mut_def(node_stats):
    
    for branch in node_stats:
        muts= branch['muts']
        T= branch['T']
        
        mut_t= np.random.choice(list(range(T)),size= muts)
        
        timetable= get_setlist(mut_t,T= T)
        
        timetable= np.array(timetable)
        branch['timetable']= timetable
    
    return node_stats

def get_up(leaf,parent_dict,uplist= []):
    """
    """
    uplist.append(leaf)
    parent= parent_dict[leaf]
    
    if parent in parent_dict.keys():
        get_up(parent,parent_dict,uplist= uplist)
    else:
        uplist.append(parent)
    
    return uplist

def get_track(leaf,tree_summ):

    edge_dict= get_edge_dict(tree_summ)

    parent_dict= {}
    for z,g in edge_dict.items():
        for nd in g:
            parent_dict[nd]= z

    leaf_track= get_up(leaf,parent_dict,uplist= [])
    
    return leaf_track


####################################
####################################




def get_probs(sim_sfs,tree_summ):

    leaf_tracks= {
        l: get_track(l,tree_summ) for l in tree_summ['leaves']
    }


    node_stats= get_sum_stats(sim_sfs,node_info= [])
    node_stats= node_mut_def(node_stats)

    node_stats= {
        z["node"]: {g:f for g,f in z.items() if g != "node"} for z in node_stats
    }
    
    return node_stats, leaf_tracks

def merge_branch_info(node_stats,leaf_tracks,tree_summ):
    '''

    '''
    pick_dict= {}
    pick_array= []
    pseg_array= []
    gen_track= []
    sfs_stack= []
    for leaf in tree_summ['leaves']:

        upsteps= leaf_tracks[leaf]
        sfs_array= node_stats[leaf]['specs']
        fixed_array= node_stats[leaf]['fixed']

        timestamps= []
        upsteps_track= []
        for upst in upsteps:
            ntable= node_stats[upst]['timetable']
            if len(ntable):
                timestamps.append(ntable)
            upsteps_track.extend([upst]*node_stats[upst]['T'])


        timestamps= timestamps[::-1]    

        timestamps= np.concatenate(tuple(timestamps))

        N_t= timestamps[:-1,1]
        fixed_tal= fixed_array[:-1]

        pick= fixed_tal * N_t.reshape(-1,1)
        pick= np.array(pick,dtype= int)

        pick_dict[leaf]= pick
        gen_track.append(upsteps_track)
        pseg_array.append(fixed_tal)
        pick_array.append(N_t)


    gen_track= np.array(gen_track).T[::-1]
    pseg_array= np.concatenate(tuple(pseg_array),axis= 1)
    pseg_array= np.nan_to_num(pseg_array)

    pick_array= np.array(pick_array).T
    
    return gen_track, pseg_array, pick_array


def pop_frequencies(gen_track,pseg_array,pick_array,node_stats,leaves):

    freq_dict= {
        pop: [] for pop in leaves
    }

    for idx in range(gen_track.shape[0]-1):
        vect_comp= list(gen_track[idx])
        comp_dict= {z: [x for x in range(gen_track.shape[1]) if vect_comp[x]==z] for z in list(set(vect_comp))}

        for gp,pp in comp_dict.items():
            prob_vec= [pseg_array[idx,z] for z in pp]
            if not np.sum(prob_vec):
                continue
            prob_vec= np.array(prob_vec) / np.sum(prob_vec)

            nmuts= [pick_array[idx,z] for z in pp]
            nmuts= list(set(nmuts))[0]
            if not nmuts:
                continue

            tt= [np.random.choice([0,1],size= nmuts,p= [1-y,y]) for y in prob_vec]
            tt= np.array(tt).T
            

            for mt in range(nmuts):
                gladly= tt[mt]

                if not sum(gladly):
                    continue

                for pidx in range(gen_track.shape[1]):
                    pop= leaves[pidx]
                    if pidx not in pp:
                        freq_dict[pop].append(0)
                        continue

                    if gladly[pp.index(pidx)] == 0:
                        freq_dict[pop].append(0)
                        continue    
                    
                    freq= node_stats[pop]['specs'][idx]
                    freq= np.random.choice(list(range(len(freq))),size= 1,p= freq)[0] / len(freq)

                    freq_dict[pop].append(freq)
    
    freq_array= np.array(list(freq_dict.values()))
    
    return freq_array
    
    
def hap_sample(freq_array,Sizes=50):
    '''
    '''
    if isinstance(Sizes,int):
        Sizes= [Sizes] * freq_array.shape[0]

    N_pops= freq_array.shape[0]

    pop_names= ['pop{}'.format(x+1) for x in range(len(Sizes))]
    labels= np.repeat(np.array([pop_names]),Sizes)

    data= []

    for k in range(N_pops):

        probs= freq_array[k]
        probs[(probs > 1)]= 1

        m= Sizes[k]
        Haps= [[np.random.choice([1,0],p= [1-probs[x],probs[x]]) for x in range(len(probs))] for acc in range(m)]

        data.extend(Haps)

    data= np.array(data)

    return data, pop_names, labels
