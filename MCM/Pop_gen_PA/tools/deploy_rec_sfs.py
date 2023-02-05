

from tools.sfs_utilities import (
    single_gen_matrix_v2, freq_progr_func, get_fixedtally_v2
)

from tools.ne_util import (
    theta_constant, theta_exp
)

from tools.branch_utilities import (
    demo_file_branchProcess, node_assign, get_edge_dict,
    branch_progress, traverse_sfs,
    get_probs, merge_branch_info, pop_frequencies, hap_sample
)

from tools.ABC_utilities import (
    sample_dist_beta, return_replica
)



###
def deploy_sim(demo_file,
				anc_r= '0',
                Nsamp= 1,
                sizes= 1000,
                med_samp= True,
                rescale_dict= {},
                directed= False,
                M_convert= True,
                sample_func= sample_dist_beta,
                scale_sim= False,
                burnin= 20,
                seqL= 1e6,
                muG= 1.08e-9,
                s= 0,
                ploidy= 2,
                Sizes= 50,
                ):
    '''
    deploy sims.
    '''
    tree, demo_data, tree_summ, tree_demo= demo_file_branchProcess(demo_file)


    replic= [return_replica(x,sample_func=sample_func,func=int,rescale_dict= rescale_dict,med_samp= med_samp) for x in tree_demo]

    node_dict= node_assign(replic, tree_summ)

    theta_dict= {
        'func': theta_constant,
        'kargs': {}
    }

    anc_name= '0'

    anc_size= demo_data['N'][anc_r]
    anc_sample= sample_func(1,*anc_size,func= int,func_args= [],med_samp= med_samp)
    asize= int(anc_sample[0])

    sim_sfs= traverse_sfs(node_dict,tree_summ,theta_dict,node= '0',fr= 1,Ne= asize,Ne0=asize,
                            T= '0',muG= muG, ploidy= ploidy,s= s, seqL= seqL,
                            scale_sim= scale_sim,sample_func= sample_func,
                            med_samp= med_samp)

    node_stats, leaf_tracks= get_probs(sim_sfs,tree_summ)
    gen_track, pseg_array, pick_array= merge_branch_info(node_stats,leaf_tracks,tree_summ)

    leaves= tree_summ['leaves']

    freq_array= pop_frequencies(gen_track, pseg_array, pick_array,node_stats,leaves)

    data, pop_names, labels= hap_sample(freq_array,Sizes= 50)
    
    return data, pop_names, labels, freq_array

