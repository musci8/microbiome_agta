import graph_tool as gt
import numpy as np
import graph_tool.draw
import graph_tool.clustering
import pandas as pd
import glob
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
from itertools import combinations


edges_all = pd.read_csv('EdgesAllGt.csv')
camps_id = dict(pd.read_csv('CampsId.csv')[['id','subcamp']].values)
subcamps = pd.read_csv('CampsId.csv')['camp'].tolist()



g = gt.Graph(directed=False)
# defining edge prop
weight = g.new_ep('float')
weight2 = g.new_ep('float')
ecolor1 = g.new_ep('string')
ecolor2 = g.new_ep('string')
etype = g.new_ep('string')
eorder1 = g.new_ep('int')
eorder2 = g.new_ep('int')
#defining node prop
node_size     = g.new_vp('int')
label         = g.new_vp('string')
group         = g.new_vp('int')
node_color    = g.new_vp('int')
edge_color    = g.new_ep('vector<string>')
nodeId        = g.new_vp('string')
halo          = g.new_vp('bool')
halo_size     = g.new_vp('float')
halo_color    = g.new_vp('string')

g.add_edge_list(edges_all[['fromid','toid']].values)
n_edges = len(edges_all)

#eprop
weight.a = np.sqrt((edges_all.N.astype(float).apply(lambda x: np.log(1+x)/2)).values)
W_p = [np.percentile(weight.a,25), np.percentile(weight.a,50),np.percentile(weight.a,75)]                
for i in g.get_vertices():
    
    node_color[i] = camps_id[i]
    group[i] = int(subcamps[i])

i = 0
ec1 = ['black' if ii<len(edges_all) else 'white' for ii in range(len(g.get_edges()))]
ec2 = ['white' if ii<len(edges_all) else 'black' for ii in range(len(g.get_edges()))]
eo1 = [2 if ii<len(edges_all) else 1 for ii in range(len(g.get_edges()))]
eo2 = [1 if ii<len(edges_all) else 2 for ii in range(len(g.get_edges()))]
for e in (g.get_edges()): 
    
    ecolor1[g.edge(*e)] = ec1[i]
    ecolor2[g.edge(*e)] = ec2[i]
    eorder1[g.edge(*e)] = eo1[i]
    eorder1[g.edge(*e)] = eo2[i]
    if weight[g.edge(*e)]<W_p[0]:
        edge_color[g.edge(*e)] = '#4D585858'

    if weight[g.edge(*e)]>=W_p[0] and weight[g.edge(*e)]<W_p[1]:
        edge_color[g.edge(*e)] = '#99585858'

    if weight[g.edge(*e)]>=W_p[1] and weight[g.edge(*e)]<W_p[2]:
        edge_color[g.edge(*e)] = '#CC585858'

    if weight[g.edge(*e)]>W_p[2]:
        edge_color[g.edge(*e)] = '#FF585858'

    i+=1

pos=gt.draw.sfdp_layout(g, groups=group, mu=100,C=2,K=3,p=2,gamma=0.01, mu_p=1)


pos_old = list(pos)

fig = plt.figure(figsize=(10,10))

gt.draw.graph_draw(g, pos=pos,
                    vprops={'fill_color' : '#377eb8',
                            'color' : '#377eb899',
                            'pen_width' : 4,
                            'size' : 20,
                            'anchor' : 0},
                   eprops={
                          'pen_width': weight,
                       'color' : '#5858584D'
                          },
                   output_size=(700,700),
                  output='/home/user/NetworkSocial.svg'
                  )

bacteria_df2 = pd.read_csv('EdgesBacteriaGt.csv')
scale = (edges_all.N.max()/bacteria_df2.N.max())
g = gt.Graph(directed=False)
# defining edge prop
weight = g.new_ep('float')
weight2 = g.new_ep('float')
ecolor1 = g.new_ep('string')
ecolor2 = g.new_ep('string')
etype = g.new_ep('string')
eorder1 = g.new_ep('int')
eorder2 = g.new_ep('int')
label         = g.new_vp('string')
#defining node prop
node_size     = g.new_vp('int')
group         = g.new_vp('int')
node_color    = g.new_vp('int')
nodeId        = g.new_vp('string')
pos2 = g.new_vp('vector<float>')
halo          = g.new_vp('bool')
halo_size     = g.new_vp('float')
halo_color    = g.new_vp('string')

g.add_edge_list(bacteria_df2[['fromid','toid']].values)
n_edges = len(edges_all)

#eprop
weight.a = np.sqrt((1+bacteria_df2.N.astype(float).values)/2)
                 
for i in g.get_vertices():
    
    node_color[i] = camps_id[i]
    pos2[i] = pos_old[i]
    group[i] = int(subcamps[i])

i = 0
ec1 = ['black' if ii<len(edges_all) else 'white' for ii in range(len(g.get_edges()))]
ec2 = ['white' if ii<len(edges_all) else 'black' for ii in range(len(g.get_edges()))]
eo1 = [2 if ii<len(edges_all) else 1 for ii in range(len(g.get_edges()))]
eo2 = [1 if ii<len(edges_all) else 2 for ii in range(len(g.get_edges()))]
for e in g.get_edges(): 
    
    ecolor1[g.edge(*e)] = ec1[i]
    ecolor2[g.edge(*e)] = ec2[i]
    eorder1[g.edge(*e)] = eo1[i]
    eorder1[g.edge(*e)] = eo2[i]

    i+=1


fig = plt.figure(figsize=(10,10))
gt.draw.graph_draw(g, pos=pos2,
                    vprops={'fill_color' : '#984ea3',
                            'size' : 20,
                            'color' : '#984ea399',
                            'pen_width' : 4,
                            'anchor' : 0},
                   eprops={
                          'pen_width': weight,
                       'color' : '#5858584D'
                          },
                   output_size=(700,700),
                  output='/home/user/NetworkBacteria.svg'
                  )
