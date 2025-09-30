#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import corner

trace = np.load('chain.npz')
#trace = trace['arr_0'][:,200:,:]
trace = trace['arr_0'][1000:,:,:]
print(trace.shape)

trace = trace.reshape(-1,8)
fig=corner.corner(trace, show_titles=True, nbins=30, title_fmt='.5f', color='b',
              quantiles=[0.025, 0.5, 0.975], verbose=True)
fig.savefig('mcmc.pdf', format='pdf', bbox_inches='tight')
