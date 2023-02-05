import numpy as np
import itertools as it
import pandas as pd

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import scale


from IPython.display import clear_output


def Distance_analysis(SequenceStore,target,refs_lib,DIMr = 'PCA',
                                            Bandwidth_split= 30,
                                            ncomp_local= 4,
                                            clsize= 15):
    Clover= []
    Coordinates= []
    Clusters_coords= []
    PC_var= recursively_default_dict()

    Dist_vectors= []

    Distances= []
    center_distances= []

    Ref_stats= []
    Ref_stats_lib= recursively_default_dict()
    

    for CHR in SequenceStore.keys():
        print('going on CHR: '+ str(CHR))
        for bl in SequenceStore[CHR].keys():

            print('data set: {}'.format(bl))
            ### PCA and MeanShift of information from each window copied from *FM36_Galaxy.py.
            Sequences = SequenceStore[CHR][bl]

            Sequences= np.nan_to_num(Sequences)

            print(Sequences.shape)

            #### Dimensionality reduction

            if DIMr == 'PCA':
                pca = PCA(n_components=ncomp_local, whiten=False,svd_solver='randomized').fit(Sequences)
                data = pca.transform(Sequences)
                PC_var[CHR][bl]= [x for x in pca.explained_variance_]

            if DIMr == 'NMF':
                from sklearn.decomposition import NMF
                data = NMF(n_components=ncomp_local, init='random', random_state=0).fit_transform(Sequences)

            Accurate = []

            params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
            grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

            ######################################
            ####### TEST global Likelihood #######
            ######################################
            Focus_labels = [z for z in it.chain(*refs_lib.values())]

            Who= refs_lib[target]
            Whose_parents= range(sum([len(x) for x in refs_lib.values()]))

            #Refs_local= [x for x in Whose_parents if x not in Who]
            Who_feat= data[Who,:]
            Ref_feat= data[Whose_parents,:]

            #### Normalize by distance between local centroids (to compensate for bias in sampling number).
            #### identify these clusters using MS.
            #### use reference accessions NOT in the cluster identified.
            Dpool= data[[x for x in Whose_parents if x not in Who],:]
            Pdistances= []

            bandwidth = estimate_bandwidth(Dpool, quantile=0.15)
            if bandwidth <= 0:
                bandwidth= .1
            params = {'bandwidth': np.linspace(np.min(Dpool), np.max(Dpool),30)}
            grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

            ## perform MeanShift clustering.
            ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=False, min_bin_freq=25)
            ms.fit(Dpool)
            labels1 = ms.labels_
            label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}

            centers= [np.mean(Dpool[label_select[z],:],axis= 0) for z in label_select.keys()]

            #### Data set of evenly sampled data. ##
            ## We'll generate 50 new observations from each cluster identified locally. ##
            N= 50
            Proxy_data= []
            label_select_labels= [z for z in it.chain(*[[x] * len(label_select[x]) for x in label_select.keys()])]
            Center_store= {}
            Proxy_indexes= {}
            distance_vecs= []

            for lab in label_select.keys():
                if len(label_select[lab]) < 3:
                    continue

                Quanted_set= Dpool[label_select[lab],:]

                if np.max(pairwise_distances(Quanted_set,metric= 'euclidean')) <= 1e-3:
                    Extract= Quanted_set[np.random.choice(Quanted_set.shape[0],N),:]
                else:
                    grid.fit(Quanted_set)
                    kde = grid.best_estimator_
                    Extract= kde.sample(N)

                center= np.mean(Extract,axis= 0)
                Center_store[lab]= center
                Proxy_indexes[lab]= [x for x in range((len(Center_store) - 1) * N, len(Center_store) * N)]

                Proxy_data.extend(Extract)

            Proxy_data= np.array(Proxy_data)
            ##### Get pairwise distances between centroids.

            for pair in it.combinations(label_select.keys(),2):
                coordinates= [np.mean(Dpool[label_select[z],:],axis= 0) for z in pair]
                coordinates= np.array(coordinates)
                iu_control= np.triu_indices(2,1)
                MS_pair_dist= pairwise_distances(coordinates,metric= 'euclidean')
                MS_pair_dist= MS_pair_dist[iu_control][0]
                Pdistances.append(MS_pair_dist)
            ## 

            reference_centroid= np.mean(centers,axis= 0)

            proxy_distances= pairwise_distances(reference_centroid.reshape(1,-1), Proxy_data,metric= 'euclidean')
            distances_to_center= pairwise_distances(reference_centroid.reshape(1,-1), Ref_feat,metric= 'euclidean')[0]
            self_distances= pairwise_distances(reference_centroid.reshape(1,-1), Who_feat, metric= 'euclidean')

            centroid= np.mean(Who_feat,axis= 0)
            distances_pairwise= pairwise_distances(centroid.reshape(1,-1), Ref_feat,metric= 'euclidean')[0]

            Distances.append(distances_pairwise)
            distances_pairwise= [(x - np.mean(proxy_distances)) / np.std(proxy_distances) for x in distances_pairwise]
            Clover.append(distances_pairwise)
            print(np.array(Clover).shape)

            FC_stats= [np.mean(proxy_distances),np.std(proxy_distances), np.mean(self_distances), np.std(self_distances)]
            Coord= [[CHR,bl,x] for x in Who]

            Ref_stats.append(FC_stats)
            Ref_stats_lib[CHR][bl]= FC_stats

            center_distances.append(distances_to_center)

            Coordinates.extend(Coord)
            Clusters_coords.append([CHR,bl])
            
            clear_output()

    return Distances, Clover, Ref_stats_lib, Ref_stats, center_distances, Coordinates, Clusters_coords





def reduc_clust_summ(Clover,Distances,center_distances,Ref_stats,Decoms,dr_names,kclusters= 5):

    data= recursively_default_dict()

    label_store= []

    for decom in range(len(Decoms)):
        dr= Decoms[decom].fit(Clover)
        COMPS= dr.transform(Clover)

        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=kclusters, random_state=0).fit(COMPS)
        labels1 = kmeans.labels_

        label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}

        Cameo = []

        for cramp in sorted(label_select.keys()):
            Clamp = np.mean(Clover[label_select[cramp],:],axis = 0)
            Fry = [x for x in Clamp]
            Cameo.append(Fry)

        Cameo = np.array(Cameo).T
        
        
        COMPS= pd.DataFrame(COMPS,columns= ['pc{}'.format(x+1) for x in range(COMPS.shape[1])])
        COMPS['label']= labels1

        data[dr_names[decom]]= {
            'features': COMPS,
            'KDE':pd.DataFrame(Cameo),
            'stats': recursively_default_dict(),
            'labels_l1': labels1
        }
        
        labels_second_layer= [-1]*Distances.shape[0]

        for Cl in label_select.keys():
            if len(label_select[Cl]) <= len(label_select):
                continue
            ## retrieve distances to centroids selected
            New_comp= Distances[[x for x in label_select[Cl]]]
            ## identify problem rows.
            NAs_row= [sum(np.isnan(New_comp[x])) for x in range(New_comp.shape[0])]
            NAs_row= [x for x in range(len(NAs_row)) if NAs_row[x] == New_comp.shape[1]]
            ## remove problem rows.
            New_comp= New_comp[[x for x in range(New_comp.shape[0]) if x not in NAs_row]]
            ## retrieve distances to global center for centroids selected
            distance_to_center= center_distances[[x for x in label_select[Cl]]]
            distance_to_center= distance_to_center[[x for x in range(New_comp.shape[0]) if x not in NAs_row]]

            ### PCA on Distances
            pca= PCA(n_components= 3,whiten= False).fit(New_comp)
            new_feat= pca.transform(New_comp)
            clock = pca.components_.T*np.sqrt(pca.explained_variance_)

            ## Kmeans in feature space
            from sklearn.cluster import KMeans
            kmeans = KMeans(n_clusters=8, random_state=0).fit(new_feat)
            new_labels = kmeans.labels_
            another_label_select = {y:[x for x in range(len(new_labels)) if new_labels[x] == y] for y in sorted(list(set(new_labels)))}

            ### set up second layer of labels:
            for color in another_label_select.keys():
                for thot in range(len(another_label_select[color])):
                    if another_label_select[color][thot] not in NAs_row:
                        indexed= label_select[Cl][another_label_select[color][thot]]
                        labels_second_layer[indexed]= color

            ### average distances across clustered profiles.
            Cameo = []
            center_means= []
            for cramp in sorted(another_label_select.keys()):
                Fry = np.mean(New_comp[another_label_select[cramp],:],axis = 0)
                Phillip= np.mean(distance_to_center[another_label_select[cramp],:],axis = 0)

                center_means.append(Phillip)
                Cameo.append(Fry)

            ### prepare data to save
            Cameo = np.array(Cameo).T
            center_means= np.array(center_means).T

            new_feat= pd.DataFrame(new_feat,columns= ['pc{}'.format(x+1) for x in range(new_feat.shape[1])])
            new_feat['label']= new_labels

            clock= pd.DataFrame(clock, columns= ['pc{}'.format(x+1) for x in range(clock.shape[1])])
            #clock['id']= [Fam[x] for x in Parent_list]

            data[dr_names[decom]]['stats'][Cl]= {
                'profiles': new_feat,
                'features': clock,
                'averages': pd.DataFrame(Cameo),
                'shapes': pd.DataFrame(center_means)
            }
        data[dr_names[decom]]['labels_l2']= labels_second_layer



    data['Ref_stats']= Ref_stats
    data['Distances']= Distances
    data['centre_dists']= center_distances
    
    return data


