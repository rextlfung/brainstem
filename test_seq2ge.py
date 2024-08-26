'''
Short script for testing Matteo Cencini's PyPulCeq library
Converts a .seq file to a GE compatible format
'''

#%% Import the libraries
%load_ext autoreload
%autoreload 2
import pypulseq as pps
import pypulceq as ppc

#%% Load in sequence object
fn = '3DEPImultishot_loop.seq'
seq = pps.Sequence()
seq.read(fn)

#%% Convert
ppc.seq2ge('b0',seq,verbose=True)
# %%
