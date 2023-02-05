import numpy as np
import itertools as it

import plotly.graph_objs as go

from sklearn.decomposition import PCA
import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


import Lab_modules.StructE_tools as Ste


###


def Gen_samples(Pops,Sizes,vector_lib,prior_func,prior_kwargs,return_pca= False,n_comp= 100,prior= 'sinusoid',range_diff= [-10,10],steps= 100):
    
    print('...')
    
    Npops= len(Sizes)
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(vector_lib)
    features = pca.transform(vector_lib)
    
    target= [0,1]

    Windows= recursively_default_dict()
    
    threshold= .005
    P= 30

    Fst_labels= []

    Fst_crawl= []

    Fst_windows= recursively_default_dict()
    
    
    if len(Sizes) > Npops:
        print('Size vector longer than N pops. using first {}.'.format(Npops))
        Sizes= Sizes[:Npops]

    for angle in np.linspace(range_diff[0],range_diff[1],steps):
        coords= features[Pops,:]
        vector2= coords[target[0]] - coords[target[1]]
        
        coords= prior_func(coords,target,vector2,angle,**prior_kwargs)
        
        new_freqs= pca.inverse_transform(coords[:Npops])
        
        N_pops= len(Sizes)
        

        data= []

        for k in range(N_pops):

            probs= new_freqs[k,:]

            probs[(probs > 1)]= 1
            probs[(probs < 0)]= 0
            m= Sizes[k]
            Haps= [[np.random.choice([1,0],p= [1-probs[x],probs[x]]) for x in range(vector_lib.shape[1])] for acc in range(m)]

            data.extend(Haps)

        data= np.array(data)

        if return_pca:
                pca2 = PCA(n_components=n_comp, whiten=False,svd_solver='randomized')

                data= pca2.fit_transform(data)

        Windows[int(angle*1000)]= data

        ##### FSTs
        Pairwise= Ste.return_fsts2(new_freqs)
        Pairwise['angle']= [angle] * Pairwise.shape[0]
        Fst_labels.extend(Pairwise.pops)

        Fst_crawl.extend(Pairwise.fst)



        Fst_windows[int(angle*1000)]= Pairwise

        ### store stuff.
        Windows[int(angle*1000)]= data
    
    
    Windows= {1:Windows}
    Fst_windows= {1:Fst_windows}
    print('Done.')
    
    return Windows, Fst_windows


def Check_Path(Npops,vector_lib,prior_func,prior_kwargs,Pops= [], random= True,n_comp= 100,range_diff= [-10,10],steps= 100):
    
    pca = PCA(n_components=100, whiten=False,svd_solver='randomized').fit(vector_lib)
    features = pca.transform(vector_lib)# * pca.explained_variance_ratio_
    
    target= [0,1]
    
    if random:
        Pops= np.random.choice(vector_lib.shape[0],Npops,replace= False)
    
    
    Windows= recursively_default_dict()
    
    Fst_labels= []
    
    Fst_crawl= []
    
    Fst_windows= recursively_default_dict()
    
    
    for angle in np.linspace(range_diff[0],range_diff[1],steps):
        
        coords= features[Pops,:]
        
        vector2= coords[target[0]] - coords[target[1]]
        
        coords, prior= prior_func(coords,target,vector2,angle,passport= True,**prior_kwargs)
        
        new_freqs= pca.inverse_transform(coords)
        new_freqs[new_freqs > 1] = 1
        new_freqs[new_freqs < 0] = 0
        
        Pairwise= Ste.return_fsts2(new_freqs)
        Pairwise['angle']= [angle] * Pairwise.shape[0]
        #
        Fst_labels.extend(Pairwise.pops)

        Fst_crawl.extend(Pairwise.fst)
        

        Fst_windows[int(angle*1000)]= Pairwise

        ### store stuff
    
    Fst_windows= {1:Fst_windows}

    fig_data= [go.Scatter(
    x= [x for x in Fst_windows[1].keys()],
    y= [Fst_windows[1][x].fst[i] for x in Fst_windows[1].keys()],
    mode= 'markers',
    name= '{}'.format([x for x in it.combinations(range(Npops),2)][i])
    ) for i in range(len([x for x in it.combinations(range(Npops),2)]))
    ]

    layout = go.Layout(
        title= 'Fst across sets. prior: {}'.format(prior),
        yaxis=dict(
            title='fsts',
            range= [0,.5]),
        xaxis=dict(
            title='eucledian distance in feature space')
    )

    fig= go.Figure(data=fig_data, layout=layout)
    
    if random:
        return fig, Pops, prior
    else: return fig, prior



def plot_GenFst(Fst_lib,Npops,Chr):
    
    fig_data= [go.Scatter(
    x= [x for x in Fst_lib[1].keys()],
    y= [Fst_lib[1][x].fst[i] for x in Fst_lib[1].keys()],
    mode= 'markers',
    name= '{}'.format([x for x in it.combinations(range(Npops),2)][i])
    ) for i in range(len([x for x in it.combinations(range(Npops),2)]))
    ]

    layout = go.Layout(
        title= 'Pairwise structure',
        yaxis=dict(
            title='fst values',
            range= [0,.5]),
        xaxis=dict(
            title='data sets')
    )

    fig= go.Figure(data=fig_data, layout=layout)

    iplot(fig)

