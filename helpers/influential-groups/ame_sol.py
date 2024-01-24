import pickle
import json
import argparse
import numpy as np
from gcm import *
import matplotlib.pyplot as plt

# infection function
beta = lambda n, i, scale, shape: scale * i**shape

# Define the parser
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, help="Experience directory")
parser.add_argument("--conf", type=str, help="Configuration file name")
args = parser.parse_args()

# load configuration file
with open(f"{args.dir}/{args.conf}.json") as filename:
    conf = json.load(filename)

# load dataset
dataset_dir = conf['dataset_dir']
dataset = conf['dataset']

with open(f'{dataset_dir}form_{dataset}.pk', 'rb') as filename:
    data = pickle.load(filename)

gm = data['gm']
pn = data['pn']
nmax = data['nmax']
mmax = data['mmax']
state_meta = get_state_meta(mmax, nmax, gm, pn)

# get invasion threshold and bistability threshold
shape = conf['shape_infection']
scale_c = invasion_threshold_safe(beta, gm, pn, fixed_args=(shape,),
                                  min_param=10**(-14), max_param=1)
shape_c = bistability_threshold_safe(beta, gm, pn, min_params=(10**(-14), 1),
                                     max_params=(1, 7))
print(f"Invasion threshold {scale_c}")
print(f"Bistability threshold {shape_c}")
lower_param = conf['lower_param'] * scale_c
higher_param = conf['higher_param'] * scale_c
print('lower_param', lower_param)
print('higher_param', higher_param)
print('conf_lower_param', conf['lower_param'])
print('conf_higher_param', conf['higher_param'])
print('scale_c', scale_c)

# get the upper branch
# backward direction
param_var = -0.01 * scale_c
Jtol = 0.0001
stable_param, stable_fixed_point, stable_infected_fraction = \
        stable_branch(beta, state_meta, higher_param, param_var,
                      rtol=10**(-12), max_iter=20000,
                      fixed_args=(shape,), Jtol=Jtol, verbose=False)

# get the unstable branch
if 'unstable' in conf and conf['unstable']:
    param_var = abs(stable_param[-2] - stable_param[-1])
    fni = stable_fixed_point[-1][1]
    param_init = stable_param[-1]
    unstable_param, unstable_fixed_point, unstable_infected_fraction = \
            unstable_branch(fni, beta, state_meta, param_init,
                            param_var, fixed_args=(shape,), init_iter=100,
                            h=10**(-2),
                            rtol=10**(-12), Jtol=Jtol,
                            max_iter=10000, verbose=False)

# reverse the direction
param_list_upper = list(reversed(stable_param))
I_list_upper = list(reversed(stable_infected_fraction))
state_list_upper = list(reversed(stable_fixed_point))

# get In for the upper branch
In_list_upper = [np.zeros(nmax + 1) for state in state_list_upper]
for ind, state in enumerate(state_list_upper):
    for n in range(2, nmax + 1):
        In_list_upper[ind][n] = np.sum(state[1][n] * np.arange(nmax + 1)) / n
In_list_upper = np.array(In_list_upper)
param_list_lower = np.linspace(lower_param, scale_c, 2)
I_list_lower = np.zeros(2)

# plot
plt.plot(param_list_upper, I_list_upper)
plt.plot(param_list_upper, In_list_upper[:, nmax],
         label=fr'$n = {nmax}$')
plt.plot(param_list_upper, In_list_upper[:, 2],
         label=fr'$n = 2$')
plt.plot(param_list_upper, In_list_upper[:, 2])
plt.plot(param_list_lower, I_list_lower)
plt.show()

# output results
results = dict()
results['scale_c'] = scale_c
results['param_list_lower'] = param_list_lower
results['param_list_upper'] = param_list_upper
results['I_list_upper'] = I_list_upper
results['I_list_lower'] = I_list_lower
results['In_list_upper'] = In_list_upper
if 'unstable' in conf and conf['unstable']:
    results['param_list_unstable'] = unstable_param
    results['I_list_unstable'] = unstable_infected_fraction


with open(f"{args.dir}/{args.conf}.pk", 'wb') as filename:
    pickle.dump(results, filename)
