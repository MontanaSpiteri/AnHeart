## GMM Unsupervised clustering from GraphSAGE embedded outputs
# PyCave 

import pycave
from pycave.clustering import KMeans
import numpy as np
import torch
import glob
from  pycave.bayes import GaussianMixture
import random
#import matplotlib.pyplot as plt
import pathlib

# set directory
directory = "/graphsage/tmp/subgraph3NB5000/"
pattern = "*_node_meta.csv"  
files = glob.glob(f"{directory}/{pattern}")
len(files)

sample_names = [f_name.split("/")[-1] for f_name in files]
sample_names = [name.split('_node_meta')[0] for name in sample_names]
sample_names

X_list = []
npy_files = [x +"_50_embeddings.npy" for x in sample_names]
npy_files 
len(npy_files)

npy_dir = directory
X_list = [np.load(npy_dir+npy) for npy in npy_files]
[x.shape for x in X_list]
X = np.concatenate(X_list)

print("Individual shapes:", [x.shape for x in X_list])
print("Concatenated shape:", X.shape)

X.shape
sample_size = [x.shape[0] for x in X_list]
sample_size

print(torch.cuda.get_arch_list())
print(torch.cuda.current_device())

k_list = [7,8,9,10,11,12,13]
estimator_gmms = [GaussianMixture(x, trainer_params=dict(accelerator='gpu', devices=1,
                                                       enable_progress_bar = True),
                                init_strategy = "kmeans++", batch_size = 15000000) for x in k_list]
for estimator_gmm in estimator_gmms:
    estimator_gmm.fit(torch.from_numpy(X))
cluster_id_gmms = [estimator_gmm.predict(torch.from_numpy(X)) for estimator_gmm in estimator_gmms]
for k in k_list: 
    list_i = k - min(k_list)
    start_i = 0
    sample_cluster_gmm = list()
    out_dir = npy_dir+"gmm"+str(k)+"_5k_trained/"
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    for s_size, region in zip(sample_size, npy_files):
        sample_cluster_gmm.append( cluster_id_gmms[list_i][start_i:(start_i+s_size)])
        start_i = start_i+s_size

    for c_id, region in zip(sample_cluster_gmm, npy_files):
        f_name = region.split("_50_embeddings")[0]
        np.savetxt(out_dir+f_name+'_noMeans.txt', c_id.numpy().astype(int), fmt="%s")

print(f"Output for {f_name}: {c_id.shape}")
