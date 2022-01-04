#!/usr/bin/env python3
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import os
import subprocess
import xarray as xr

import troute.nhd_io as nhd_io
import troute.nhd_network as nhd_network
import troute.nhd_network_utilities_v02 as nnu
import troute.routing.diffusive_utils as diff_utils
import troute.routing.diffusive_utils_mainstem as diff_utils_mainstem


# In[2]:


# load MC simulation results, which will be used for initial and boundary conditions
#chrtout_folder = pathlib.Path("/glade/work/adamw/forcing/florence_event")
chrtout_folder = pathlib.Path("/home/APD/inland_hydraulics/florence_testcase/glade/work/jamesmcc/domains/private/florence_v21/FORCING_AnA_channel-only")
chrtout_files = sorted(chrtout_folder.glob("*CHRTOUT_DOMAIN1"))


# In[3]:


with xr.open_mfdataset(
    chrtout_files, 
    combine="nested",
    concat_dim="time"
) as ds:
    
    #  lateral inflows (WRF - generated)
    q_lateral = pd.DataFrame(
        ds['qBucket'].values.T + ds['qSfcLatRunoff'].values.T,
        index=ds['feature_id'].values,
        columns=ds.time.values,
    )
    
    # wrf simulated flow (WRF - generated)
    q_wrf = pd.DataFrame(
        ds['streamflow'].values.T,
        index=ds['feature_id'].values,
        columns=ds.time.values,
    )


# In[4]:


q_wrf


# In[ ]:





# In[5]:


# Hand-selected segments along the Cape Fear River between gages at Tar Heel and Kelly 
diffusive_domain = [8831712,8831726,8831730,8831734,8831736,8831748,8831754,8831752,8831780,8831794,8831792,8831798,8831808,8831810,8831806,
                    8831818,8831820,8831826,8831830,8831832,8831840,8831846,8831848,8831852,8831858,8834892,8834894,8834898,8834904,8834908,
                    8834912,8834918,8834924]

diffusive_domain = [8831712,8831726,8831730] # only headsegments of reaches on a mainstem.


# In[ ]:





# In[6]:


#args = nhd_io.read_custom_input_new("/glade/u/home/adamw/projects/t-route/test/jobs/V03/Cape_Fear_test_MC_V3.yaml")
args = nhd_io.read_custom_input_new("/home/dongha.kim/github/t-route/test/input/yaml/Cape_Fear_test_MC_V3.yaml")
supernetwork_parameters = args[1]


# In[ ]:





# In[7]:


connections, param_df, _, _ = nnu.build_connections(
        supernetwork_parameters
    )

rconn = nhd_network.reverse_network(connections)


            
# In[8]:


# identify the tributary segments to the diffusive domain
trib_segs = []
for s in diffusive_domain:
    us_list = rconn[s]
    for u in us_list:
        if u not in diffusive_domain:
            trib_segs.append(u)

connections_trim = {k: connections[k] for k in (diffusive_domain + trib_segs)}
connections_trim[8831730]=[]  ## temporary fix to make proper TW link.
network_break_segments = set()
independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections_trim,
        network_break_segments,
    )

tw = list(nhd_network.headwaters(rconn))[0]
reach_list=reaches_bytw[tw]
rconn=independent_networks[tw]
connections_trim= nhd_network.reverse_network(rconn)
mainstem_headseg_list= diffusive_domain


# In[ ]:





# In[ ]:





# In[ ]:





# In[10]:


'''
# build initial conditions array
initial_conditions = q_wrf.loc[diffusive_domain].iloc[:,0].to_numpy()
# test: q_wrf.loc[8831712].iloc[0]
# test: q_lateral.loc[8831712].iloc[:]

# mask channel geometry data
geo_cols = param_df.columns
geo_index = param_df.loc[diffusive_domain].index
geo_data = param_df.loc[diffusive_domain].to_numpy()

# lateral inflow data
qlat_data = q_lateral.loc[diffusive_domain].to_numpy()
'''


# In[9]:


diffusive_trib_domain=list(set(diffusive_domain)|set(trib_segs))

#q_wrf.loc[8831712].iloc[0]
#q_lateral.loc[8831712].iloc[:]

# mask channel geometry data
param_df_sub= param_df.loc[
                        diffusive_trib_domain,
                        ["bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
geo_cols = param_df_sub.columns.values
geo_index = param_df_sub.index.values
geo_data = param_df_sub.values

# build initial conditions array
q0 = q_wrf.loc[geo_index].iloc[:,0].to_numpy()
q0_2 = [[val] for val in q0]
initial_conditions = np.asarray(q0_2)

# lateral inflow data
qlat_data = q_lateral.loc[geo_index].to_numpy()
# MC data
q_wrf_data= q_wrf.loc[geo_index].to_numpy()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[11]:


import yaml
custom_input_file='../../../../test/input/yaml/Florence_Benchmark_diff.yaml'
with open(custom_input_file) as custom_file:
    data = yaml.load(custom_file, Loader=yaml.SafeLoader)
diffusive_parameters = data.get("diffusive_parameters", {})  
upstream_results={}
lake_segs=[]
waterbodies_df_sub = pd.DataFrame()
waterbodies_df_sub.values
lake_segs=np.asarray(lake_segs)
wbody_params=np.array(waterbodies_df_sub.values)
qts_subdivisions=12
nsteps=288
dt=300


# In[ ]:





# In[13]:


# build upstream results array
#trib_segs
#wbody_params


# In[14]:


diff_inputs= diff_utils_mainstem.diffusive_mainstem_input_data_v01(
    tw,
    connections_trim,
    rconn,
    reach_list,
    mainstem_headseg_list,
    q_wrf_data,    
    diffusive_parameters,
    geo_cols,
    geo_index,
    geo_data,
    qlat_data,
    initial_conditions,
    upstream_results,
    qts_subdivisions,
    nsteps,
    dt,
    lake_segs,
    wbody_params,
)


# In[ ]:




