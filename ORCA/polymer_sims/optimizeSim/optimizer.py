#### 12/20/24

import numpy as np
import os
import datetime

from le_1D import le_oneD
from le_3D import le_threeD
from trajectory_to_txt import main as ttt
from make_distanceMap import main as dmap

import scipy.stats as stats
from skopt import gp_minimize
from skopt.space import Real

def objective(params):
    # Load in ground truth distance map (ORCA)
    truth = np.loadtxt('Granta519cl27_0hr_MYC_ORCA_aggDistMatrix.tsv', delimiter='\t')

    # Unpack parameters
    B1_cap, B1_rel, prom_cap, prom_rel, B2_cap, B2_rel = params

    name = f"B1_{B1_cap}_{B1_rel}-Promoter_{prom_cap}_{prom_rel}-B2_{B2_cap}_{B2_rel}"
    dir_name = f"simulation_run_{name}_dir"
    os.mkdir(dir_name)

    # Run full simulation
    le_oneD(name, B1_cap, B1_rel, prom_cap, prom_rel, B2_cap, B2_rel)
    le_threeD(dir_name)

    # Make distance map of simulated polymers
    ttt(f"{dir_name}")
    avg_sim_dist = dmap(f"{dir_name}/confs_txt")
    
    # Calculate loss
    coeff = stats.spearmanr(avg_sim_dist.flatten(), truth.flatten())
    
    end_time = datetime.datetime.now()
    with open(f"{dir_name}/{name}_coeff.txt", 'w') as outf:
        outf.write(f"Completed at {end_time}")
        outf.write(f"Spearman corr. coefficient: {coeff}")
    return -coeff[0]

def runner():
    search_space = [
        Real(0.1, 0.5, name='B1_cap'),
        Real(0.01, 0.05, name='B1_rel'),
        Real(0.05, 0.25, name='prom_cap'),
        Real(0.001, 0.01, name='prom_rel'),
        Real(0.1, 0.5, name='B2_cap'),
        Real(0.01, 0.05, name='B2_rel')
    ]

    result = gp_minimize(
    func=objective,  # The function to minimize
    dimensions=search_space,  # Parameter space
    acq_func='EI',            # Acquisition function ('EI', 'PI', 'LCB')
    n_calls=50,               # Number of function evaluations
    n_random_starts=10,       # Number of random initial evaluations
    random_state=42           # For reproducibility
    
    )

    with open('sim_opt_res.txt', 'w') as o:
        o.write(result.x)
        o.write(result.fun)

if __name__ == '__main__':
    runner()
