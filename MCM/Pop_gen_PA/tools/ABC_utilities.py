import numpy as np
import itertools as it
import os

import collections
def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


def str_simple(string,tag='(',end= ')'):
    '''
    return index of first comma in the first layer of a nested tuple string.
    '''
    left= 0
    right= 0
    
    for s in range(1,len(string)):
        
        if string[s] == tag:
            left += 1
        if string[s] == end:
            right += 1
        
        if left == right:
            return s + 1


def split_middle_comma(tree_sr):
    '''
    split a 3-element tuple string.
    '''
    line= str(tree_sr)
    
    if line[1] != '(':
        #
        line= line.strip(')').strip('(')
        line= line.split(',')
        return line
    
    indx_split= str_simple(tree_sr)
    
    pos1= tree_sr[:indx_split][1:]
    pos2= tree_sr[indx_split+1:][:-1]
    
    if pos2[0] != '(':
        pos2,pos2b= pos2.split(',')
        return [pos1,pos2,pos2b]
    
    
    pos2= '({})'.format(pos2)
    
    pos2_idx= str_simple(pos2)
    
    pos2b= pos2[:pos2_idx][1:]
    pos2c= pos2[pos2_idx+1:][:-1]
    
    return [pos1,pos2b,pos2c]



def det_parenthesis(node_str,det= '('):
    if node_str[0]== det:
        return tuple_tree_rec(node_str)
    else:
        return node_str

    
def tuple_tree_rec(line):
    '''
    read binary trees in ((node1,node2),node3). 
    Recursive = inner branch number limited by machine recursion limit.
    
    '''
    test= split_middle_comma(line)
    
    return {
        test[0]: det_parenthesis(test[0]),
        test[1]: det_parenthesis(test[1]),
        'node': test[2]
    }


## read dems file

def read_demofile(filename):
    '''
    read demographic history table. Simple.
    first line must describe format in ((n1,n2),n3) format.
    1st column: T,N or M for time, size and migration.
    2nd colomn: ordered times, population names and concatenated (into tuples) population names, respectively.
    Further columns are values. 3 = mean or sd, 2-4 = lower & upper confidence intervals respectively.
    
    return tree from first line, dictionary of values w/ levels 1st then second columns.
    '''
    args= recursively_default_dict()
    
    with open(filename,'r') as fp:
        lines= fp.readlines()
    
    tree= lines[0].strip()
    tree= tuple_tree_rec(tree)
    
    mig_list= []
    for line in lines[1:]:
        line= line.split()
        if line[0] != 'M':
            args[line[0]][line[1]]= np.array(line[2:],dtype= int)
        else: 
            mig_list.append(line[1:])
        
    
    M_dict= recursively_default_dict()
    
    for line in mig_list:
        source= line[0].split('-')[0]
        #source= ''.join(source.split(','))
        sink= line[0].split('-')[1]
        #sink= ''.join(sink.split(','))
        
        M_dict[source][sink]= np.array(line[1:],dtype= float)
        
    args['M']= M_dict
    
    return tree, args


def split_unzip(str_key):
    '''
    return elements of tuple.
    '''
    sp_idx= str_simple(str_key)
    sp1= str_key[1:sp_idx]
    
    sp2= '(' + str_key[sp_idx+1:]
    
    sp2_idx= str_simple(sp2)
    sp2= sp2[1:sp2_idx]
    
    return sp1,sp2
    
def remove_node(str_key):
    '''
    remove 3d position from tuple.
    '''
    
    if str_key[0] != '(':
        return str_key
    
    sp1, sp2= split_unzip(str_key)
    
    return '({})'.format(sp1 + ','+ sp2)


def pops_rec_decide(str_key,pops= []):
    '''
    co-factor to get_pops. If no enclosure, append to listn else launch get_pops within the next level.
    '''
    if str_key[0] != '(':
        pops.append(str_key)
    else:
        get_pops(str_key,pops= pops)

def get_pops(str_key,pops= []):
    '''
    Extract pops from tuples ignoring 3d position. Allows for nested loops.
    return list of pops. returns string if not enclosed in parentheses.
    '''
    pps= remove_node(str_key)
    
    if str_key[0] != '(':
        pops.append(str_key)
        if ',' in str_key:
            pops.extend(str_key.split(','))
    else:
        pps= split_unzip(str_key)

        for p in pps:
            pops_rec_decide(p,pops= pops)
    
    return pops


def make_pop(str_key):
    '''
    make pop name as concatenation of pops in nested tuple.
    '''
    pop_name= get_pops(str_key,pops= [])
    #print('kj: {}'.format(pop_name))
    pop_name= 'p' + ''.join(pop_name)
    return pop_name

def decide_rec_list(key_str,tree_dict,tree_summ,times,demo_data,func,int_sizes= [],demo_tape= []):
    
    if key_str[0] == '(':
        tree_fill_list(times,tree_dict[key_str],demo_data,tree_summ,int_sizes= int_sizes,demo_tape=demo_tape)
    

def tree_fill_list(times,tree_dict, demo_data, tree_summ,int_sizes= [], demo_tape= []):
    '''
    returns list of tree node information.
        - 'T': node index in time array. 
        - 'M': migratory events concerning new populations.
        - 'N': create population events. 
    
    '''
    #print(demo_tape)
    tree_nodes= {
        x: remove_node(x) for x in tree_dict.keys() if x != 'node'
    }
    
    nod_t= int(tree_dict['node'])
    nod_t= times[nod_t]
    
    node_mig= {}
    node_n= {}
    
    gnode_sort= {}
    
    for k,g in tree_nodes.items():
        
        if g[0] != '(':
            node_n[g]= demo_data['N'][g]
            if g in demo_data['M'].keys():
                
                node_mig[make_pop(k)]= {
                    make_pop(k): f for k,f in demo_data['M'][g].items()
                }
        
        if g[0] == '(':
            gnode= tree_dict[k]['node']
            gnode_sort[g]= gnode
    
    if gnode_sort:
        ## only using the keys, but need them sorted by values.
        gnode_sorted= sorted(gnode_sort,key= gnode_sort.get)
        
        ## from root to leaves: T1-N decreases, Node number increases.
        startat= int(tree_dict['node']) - 2
        
        ## gnode times is used to identify key in demos table, where internal branch names are coded T1-TN
        gnode_t= {
            gnode_sorted[x]:int_sizes[startat - x] for x in range(len(gnode_sorted))
        }
        
        gnode_sizes= {
            v:demo_data['N'][g] for v,g in gnode_t.items()
        }
        
        gnode_parent= {}
        
        for gnode in gnode_sorted:
            ## parents are not used. 
            parent= [x[0] for x in tree_summ['edges'] if gnode == x[1]][0]
            parent= make_pop(parent)
            gnode_parent[gnode]= parent
            
            gnode_pops= get_pops(gnode,pops= [])
            
            gnode_pops= ','.join(gnode_pops)
            
            if gnode_pops in demo_data['M'].keys():

                #print(demo_data['M'][gnode_pops].keys())
                node_mig[make_pop(gnode)]= {
                    'p'+''.join(f.split(',')): x for f,x in demo_data['M'][gnode_pops].items()
                }
                #print(node_mig[make_pop(gnode)])
            
        
        gnode_N= {
            g: gnode_sizes[g] for g in gnode_sorted
        }
        
        node_n.update(gnode_N)
    
    demo_tape.append({
        'node': nod_t,
        'T': demo_data['T'][nod_t],
        'N': node_n,
        'M': node_mig,
    
    })
    
    for edge in tree_nodes.keys():
        decide_rec_list(edge,tree_dict,tree_summ,times,demo_data,tree_fill_list,
                        int_sizes= int_sizes,demo_tape=demo_tape)
    
    return demo_tape




def get_tree_nodes(tree,nodes=[],edges= [],leaves= []):
    '''
    break down binary dictionary tree.
    '''
    
    for node in tree.keys():
        if node == 'node':
            continue
        nodes.append(node)
        
        if tree[node]== node:
            leaves.append(node)
            continue
        
        for dt in tree[node].keys():
            if dt != 'node':
                edges.append((remove_node(node),remove_node(dt)))
        
        get_tree_nodes(tree[node],nodes= nodes,edges= edges,leaves= leaves)
    
    summary_dict= {
        'nodes': nodes,
        'edges': edges,
        'leaves': leaves
    }
    
    return summary_dict


from scipy.stats import norm

def sample_dist(nsample,median,ll_cl,up_cl,assume='norm',func= '',func_args= [3],source= 0):
    '''
    determine mean and sd from UP-CL 95%.
    sample using scipy.
    '''
    
    mean_s= (ll_cl+up_cl) / 2
    sd_s= (up_cl - mean_s) / 2
    
    t= norm.rvs(loc=mean_s,scale= sd_s,size= nsample)
    if func:
        t= [func(x,*func_args) for x in t]
    
    return t



from scipy.stats import beta

def sample_dist_beta(nsample,median,ll_cl,up_cl,blur= 500,assume='norm',func= '',func_args= [3],source= 0,
						med_samp= False, rescale= 1):
    '''
    Use beta distribution to add skew.
    '''

    if not source: 
        blur= 50
    else:
        blur= 500
    
    window= up_cl - ll_cl
    sd_s= (window) / 2
    
    rate= (median - ll_cl) / window
    
    t= np.pi / 2
    
    a= np.sin(rate * t) * blur
    b= np.cos(rate * t) * blur

    if med_samp:
        f= [median] * nsample
        f= np.array(f)
    
    else:
        f= beta.rvs(a, b, size=nsample)
    
    if not source or up_cl < 1:
        f= f * window + ll_cl
    
    f= f * rescale

    if func:
        f= [func(x,*func_args) for x in f]
    
    return f



def rep_chose(val,g,tree_demo,assume= 'norm',func= int,sample_func= sample_dist_beta,source= 0,track= ['M'],
	rescale_dict= {},med_samp= False):
    
    branch= source
    rescale= 1
    if val in rescale_dict.keys():
    	rescale= rescale_dict[val]
        
    if source in rescale_dict.keys():
    	rescale= rescale_dict[source]

    if val in track:
        branch= val
    
    if isinstance(g, dict):
        return return_replica(tree_demo[val],assume= assume,source= branch,med_samp=med_samp,rescale_dict= rescale_dict)
    
    if type(g).__module__ == np.__name__:
        g= sample_func(1,*g,func= round,source= source,rescale= rescale,med_samp=med_samp)[0]
        if func: 
            g= func(g)
    
    if isinstance(g,list):
        g= [g[0],sample_func(1,*g[1],func= round,func_args= [],source= source,rescale= rescale,med_samp=med_samp)[0]]

    return g



def return_replica(tree_demo,assume='norm',func= float,sample_func= sample_dist_beta,source= 0,
						rescale_dict= {},med_samp= False):
    '''
    sample topologie from from probab demos dictionary. 
    '''
    return {
        val: rep_chose(val,g,tree_demo,assume= assume,func= func,sample_func= sample_func,source= source,rescale_dict= rescale_dict,med_samp=med_samp) for val,g in tree_demo.items()
    }



def ancestral_initialize(anc_name= 'p1',anc_size= 20000,return_list= True):
    anc_intro= '''
    sim.addSubpop("{}", {});
    c = sim.chromosome;
    catn("Ancestral: " + paste(c.ancestralNucleotides(format="char")[0:20],sep=""));
    catn();\n'''
    
    anc_intro= anc_intro.format(anc_name,str(anc_size))
    
    anc_intro= """1 {\n""" + anc_intro  + """}\n"""
    
    if return_list:
        anc_intro= anc_intro.split('\n')
        anc_intro= [x + '\n' for x in anc_intro]
        
    return anc_intro


def sample_block(gen= 60000,pops= ['p1','p2'],sizes= [500,500]):
    #pops= ','.join(pops)
    sizes= ','.join([str(x) for x in sizes])
    
    sample_simple= """
    g = c();
    pops= 0:{};
    samples= c({});
    for (x in pops) 
        g= c(g, sim.subpopulations[x].sampleIndividuals(samples[x]).genomes);

    g.outputVCF(vcf_file,simplifyNucleotides=T);
    """
    
    sample_simple= sample_simple.format(len(pops)-1,sizes)
    sample_simple= """{} late() """.format(gen) + """{\n""" + sample_simple
    sample_simple= sample_simple.split('\n')
    sample_simple= [x + '\n' for x in sample_simple]
    sample_simple.append( """}\n""")
    return sample_simple




def demos_to_SLiM(batch, template, tree, demo_data, anc_r= 'anc', Nsamp= 5, sizes= 500, burnin= 5e4,size_key= '\t{}.setSubpopulationSize({});\n',
                                                    sample_func= sample_dist_beta, mig_key= '{}.setMigrationRates(c({}), c({}));\n',
                                                    create_key= 'sim.addSubpopSplit("{}", {}, {});\n',sim_scale= 1,
                                                    rescale_dict= {},med_samp= False,M_convert= True,directed= False):

    ###
    ###
    tree_summ= get_tree_nodes({anc_r:tree},nodes=[],edges= [],leaves= [])

    recipe_dir= '/'.join(template.split('/')[:-1]) + '/'
    #
    anc_name= anc_r
    anc_size= demo_data['N'][anc_r]
    
    times_order= [0] + sorted(demo_data['T'].keys(),reverse= True)
    int_sizes= [x for x in demo_data['N'].keys() if x[0] == 'T']
    int_sizes= sorted(int_sizes)

    tree_demo= tree_fill_list(times_order,tree, demo_data, tree_summ, int_sizes= int_sizes,demo_tape= [])
    ###
    ###

    template_file= template
    
    with open(template_file,'r') as fp:
        template= fp.readlines()

    pops= tree_summ['leaves']
    
    sizes= [sizes] * len(pops)

    files= []
    anc_sample= sample_func(Nsamp,*anc_size,func= int,func_args= [],med_samp= med_samp)

    for idx_sim in range(Nsamp):
        temp= list(template)

        ## get data
        asize= int(anc_sample[idx_sim])
        replic= [return_replica(x,sample_func=sample_func,rescale_dict= rescale_dict,med_samp= med_samp) for x in tree_demo]
        replic= {x['node']: x for x in replic}


        ## existing pops, to decide whether to update or create.
        existing_pops= [anc_name]
        existing_sizes= {
        	anc_name: asize
        }
        trail_migs= []

        time_spent= int(burnin * sim_scale)
        ## setup ancestral init. 
        anc_lines= ancestral_initialize(anc_name= make_pop(anc_name),anc_size= asize)

        temp.extend(anc_lines)

        node_order= sorted(replic.keys(),reverse= True)

        for nod in node_order: 
            new_lines= ['\n']

            Tev= replic[nod]['T'] * sim_scale

            new_lines.append(str(time_spent) + ' {\n')
            time_spent += int(Tev)

            deprecate= []

            for br in replic[nod]['N'].keys():            
                pop_name= make_pop(br)
                #print(pop_name)
                parent= [x[0] for x in tree_summ['edges'] if br == x[1]][0]
                #print(parent)
                parent= make_pop(parent)
                pop_size= int(replic[nod]['N'][br] * sim_scale)

                if pop_name in existing_pops:
                    # change size
                    size_change= '\t' + size_key.format(pop_name,pop_size)
                    new_lines.append(size_change)
                    existing_sizes[pop_name]= pop_size

                else:
                    # create population. 
                    create_pop= '\t' + create_key.format(pop_name,pop_size,parent) ###
                    new_lines.append(create_pop)
                    existing_pops.append(pop_name)
                    existing_sizes[pop_name]= pop_size
                    deprecate.append(parent)


            ### remove old populations
            deprecate= list(set(deprecate))

            for old in deprecate:
                size_change= size_key.format(old,0)
                new_lines.append(size_change)

            ### migration
            new_trail=[]
            for trail in trail_migs:
                if trail[0] in existing_pops and trail[1] in existing_pops:
                    new_mig= trail[2]

                    if M_convert:
                        new_mig= new_mig / existing_sizes[trail[0]]
                    
                    mig_change= '\t' + mig_key.format(trail[0],trail[1],new_mig)
                    new_lines.append(mig_change)
                    if not directed:
                        new_mig= trail[2]
                        
                        if M_convert:
                            new_mig= trail[2] / existing_sizes[trail[1]]

                        mig_change= '\t' + mig_key.format(trail[1],trail[0],new_mig)
                        new_lines.append(mig_change)

                    
                else:
                    new_trail.append(trail)

            trail_migs= new_trail

            for mr in replic[nod]['M'].keys():
                # set migration rate.
                for dr,v in replic[nod]['M'][mr].items():

                    if mr in existing_pops and dr in existing_pops:
                        new_mig= v
                        
                        if M_convert:
                            new_mig= new_mig / existing_sizes[mr]

                        mig_change= '\t' + mig_key.format(mr,dr,new_mig)
                        new_lines.append(mig_change)
                        if not directed:
                            
                            if M_convert:
                                v= v / existing_sizes[dr]

                            mig_change= '\t' + mig_key.format(dr,mr,v)
                            new_lines.append(mig_change)

                    else:
                        trail_migs.append((mr,dr,v))

            new_lines.append('}\n\n')

            temp.extend(new_lines)

        ## sample: 
        
        sample_simple= sample_block(time_spent,tree_summ['leaves'],sizes)
        temp.extend(sample_simple)

        #### write new recipe
        recipe_name= '{}{}.slim'.format(batch,idx_sim)
        recipe_name= recipe_dir + recipe_name
        files.append(recipe_name)
        
        with open(recipe_name,'w') as fp:
            fp.write(''.join(temp))
        
    
    pops= tree_summ['leaves']
    
    return pops, files


############################################################################
############################################################################


def demo_to_recipe(demo_file,template,batch= 'test',anc_r= '0',Nsamp= 5,sizes= 500, burnin= 5e4, sim_scale= 1,recipe_dir= 'Recipes/demos_mat/',
					rescale_dict= {},med_samp= False,M_convert= False,directed= False):
    
    tree, demo_data= read_demofile(demo_file)
    
    pops, files= demos_to_SLiM(batch, template,tree, demo_data, anc_r= anc_r, Nsamp= Nsamp, sizes= sizes, burnin= burnin, sim_scale= sim_scale,
    												rescale_dict= rescale_dict,med_samp= med_samp, M_convert= M_convert, directed= directed,
                                                    size_key= '\t{}.setSubpopulationSize({});\n',
                                                    mig_key= '{}.setMigrationRates(c({}), c({}));\n',
                                                    create_key= 'sim.addSubpopSplit("{}", {}, {});\n')
    
    
    return pops, files


