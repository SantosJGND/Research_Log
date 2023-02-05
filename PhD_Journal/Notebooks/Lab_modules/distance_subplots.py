import numpy as  np
from plotly import tools
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot, iplot
from plotly.graph_objs import *
import plotly.figure_factory as ff



def dens_plot(data,data_groups,titles,Dr= 'PCA',Ncols= 2,width= 1000,height= 400):
    
    Distances= data['Distances']
    Centre_distances= data['centre_dists']
    fig_subplots = tools.make_subplots(rows= int(len(titles) / float(Ncols)) + (len(titles) % Ncols > 0), cols=Ncols,
                             subplot_titles=tuple(titles))



    for pop in range(len(titles)):
        name= titles[pop]

        pos1= int(float(pop) / Ncols) + 1
        pos2= pop - (pos1-1)*Ncols + 1

        label= data_groups[pop] + 1
        print(label)

        select_l1= [label]
        selected1= [x for x in range(Distances.shape[0]) if data[Dr]['labels_l1'][x] + 1 in select_l1]    

        x= np.mean(Distances[selected1,:],axis= 0)
        y= np.mean(Centre_distances[selected1,:],axis= 0)

        colorscale = ['#7A4579', '#D56073', 'rgb(236,158,105)', (1, 1, 0.2), (0.98,0.98,0.98)]

        trace1 = go.Scatter(
            x=x, y=y, mode='markers', name='points',
            marker=dict(color='rgb(102,0,0)', size=5, opacity=0.4)
        )
        trace2 = go.Histogram2dContour(
            x=x, y=y, name='density', ncontours=20,
            colorscale='Hot', reversescale=True, showscale=False
        )
        trace3 = go.Histogram(
            x=x, name='x density',
            marker=dict(color='rgb(102,0,0)'),
            yaxis='y2'
        )
        trace4 = go.Histogram(
            y=y, name='y density', marker=dict(color='rgb(102,0,0)'),
            xaxis='x2'
        )

        fig_subplots.append_trace(trace1, pos1, pos2)
        fig_subplots.append_trace(trace2, pos1, pos2)
        #fig_subplots.append_trace(trace3, pos1, pos2)
        #fig_subplots.append_trace(trace4, pos1, pos2)

        fig_subplots['layout']['yaxis' + str(pop + 1)].update(title='reference distances',showgrid=False,
                zeroline=False)

        fig_subplots['layout']['xaxis' + str(pop + 1)].update(title='target distances',showgrid=False,
                zeroline=False)




    fig_subplots['layout'].update(height= height,width= width)

    iplot(fig_subplots)




def trace_plot(data,data_groups,labels_pops,titles,Dr= 'PCA', Z= 0,Ncols= 2,width= 1000,height= 400):

    Distances= data['Distances']
    Centre_distances= data['centre_dists']
    Ref_stats= data['Ref_stats']

    fig_subplots = tools.make_subplots(rows= int(len(titles) / float(Ncols)) + (len(titles) % Ncols > 0), cols=Ncols,
                             subplot_titles=tuple(titles))
    

    for pop in range(len(titles)):
        name= titles[pop]

        pos1= int(float(pop) / Ncols) + 1
        pos2= pop - (pos1-1)*Ncols + 1

        label= data_groups[pop] + 1
        select_l1= [label]
        selected1= [x for x in range(Distances.shape[0]) if data[Dr]['labels_l1'][x] + 1 in select_l1]    

        #names_index = [[f for f in orderCore.ID].index(x) for x in [str(y) for y in df.iloc[:,1]]]

        ###
        scheme = labels_pops
        coords = {y:[x for x in range(len(scheme)) if scheme[x] == y and x ] for y in list(set(scheme))}

        print(coords.keys())
        pop_refs= ['gp: {}'.format(x) for x in coords.keys()]
        #color_here= color_ref

        ###

        ref_names= ['gp: {}'.format(x) for x in coords.keys()]

        #scheme= [orderCore.loc[x,'sNMF_K3'] for x in names_index]
        #scheme= {z:[x for x in range(len(scheme)) if scheme[x] == z] for z in list(set(scheme))}
        if Z== 0:
            x= np.mean(Distances[selected1,:],axis= 0)
            print(len(scheme))
            y= np.mean(Centre_distances[selected1,:],axis= 0)

            max_dists= [max(Centre_distances[x]) for x in selected1]
            min_centre= [min(Distances[x]) for x in selected1]

            Normalized_mean= np.mean(max_dists)
            Normalized_centre= np.mean(min_centre)

        if Z== 1:
            meansVars= Ref_stats[selected1,:]
            trim= [x for x in range(meansVars.shape[0]) if meansVars[x,1] > 0.05]

            max_dists= [max(Centre_distances[x]) for x in selected1]
            min_centre= [min(Distances[x]) for x in selected1]
            meansVars= meansVars[trim,:]
            select_trim= [selected1[x] for x in trim]

            Normalized_mean = [(max_dists[x] - meansVars[x,0]) / meansVars[x,1] for x in range(meansVars.shape[0])]
            Normalized_mean= np.array(Normalized_mean)
            #Normalized_mean= [x for x in Normalized_mean if np.isnan()]
            Normalized_mean[Normalized_mean > 20] = 10
            Normalized_mean[Normalized_mean < -20] = -10
            Normalized_mean= np.mean(Normalized_mean)

            Normalized_centre = [(min_centre[x] - meansVars[x,0]) / meansVars[x,1] for x in range(meansVars.shape[0])]
            Normalized_centre= np.array(Normalized_centre)
            #Normalized_mean= [x for x in Normalized_mean if np.isnan()]
            Normalized_centre[Normalized_centre > 20] = 10
            Normalized_centre[Normalized_centre < -20] = -10
            Normalized_centre= np.mean(Normalized_centre)

            sel_d= Distances[select_trim,:]
            sel_d= [(sel_d[x] - meansVars[x,0]) / meansVars[x,1] for x in range(len(select_trim))]
            x= np.mean(sel_d,axis=0)

            sel_refd=Centre_distances[select_trim,:]
            sel_refd= [(sel_refd[x] - meansVars[x,0]) / meansVars[x,1] for x in range(len(select_trim))]
            y= np.mean(sel_refd,axis=0)

        for i in coords.keys():
            trace1= go.Scatter(
                x= [x[z] for z in coords[i]],
                y= [y[z] for z in coords[i]],
                mode= 'markers',
                name= ref_names[i],
                marker=dict(color=i, size=7, opacity=0.6)
            ) 
            fig_subplots.append_trace(trace1, pos1, pos2)


        trace1= go.Scatter(
            x= [Normalized_centre],
            y= [Normalized_mean],
            mode= 'markers',
            marker=dict(color='black', size=10, opacity=0.6)
        )

        fig_subplots.append_trace(trace1, pos1, pos2)


        fig_subplots['layout']['yaxis' + str(pop + 1)].update(title= 'distances to centre, {}'.format(['raw','scaled'][Z]))

        fig_subplots['layout']['xaxis' + str(pop + 1)].update(title= 'distances to target, {}'.format(['raw','scaled'][Z]))




    fig_subplots['layout'].update(title= 'distances label {}; Target mean: {}'.format(select_l1[0],round(Normalized_mean,2)),
                                  height= height,width= width)

    iplot(fig_subplots)
