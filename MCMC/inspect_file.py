import pandas as pd
import numpy as np
import sys

vid = sys.argv[1]

path_to_file = r'/local/data/lava/janine/data'
#path_to_file = r'./Data/'

store = pd.HDFStore(path_to_file + vid + '.h5')
data_vars = store.keys()
store.close()
print(data_vars)

if '/s' in data_vars:
    s = pd.read_hdf(path_to_file + vid + '.h5', key='s')
    meta = pd.read_hdf(path_to_file + vid + '.h5', key='meta')
    K_post = meta.K_post[0]
    K_std = meta.K_std[0]
    tauy_post = meta.tauy_post[0]
    tauy_std = meta.tauy_std[0]
    n_post = meta.n_post[0]
    n_std = meta.n_std[0]
    
    print('K = ' + np.array2string(K_post, precision=1) + ' +- ' + np.array2string(2*K_std, precision=1) +
          ', tau_y = ' + np.array2string(tauy_post, precision=2) + ' +- ' + np.array2string(2*tauy_std, precision=2) +
          ', n = ' + np.array2string(n_post, precision=2) + ' +- ' + np.array2string(2*n_std, precision=2))
    
    print('Samples = ' + str(len(s.K)))