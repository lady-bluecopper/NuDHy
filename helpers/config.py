# path where the files with the metrics will be stored
met_path = '../out/'
# path where the samples are stored
sample_path = '../out/samples'
# random seed for reproducibility
seed = 0
# folder where the datasets are placed
folder = '../data/'
# max number of workers to use in the parallel computations
max_workers = 2
# number of samples to generate
num_samples = 33

# RECIPROCITY
alpha = 1e-6

# REDI
allow_inters = False

# LAPLACIAN
# max order to consider
order = 20
# number of eigenvalues to return
eig_num = 6
# whether we want the symmetric or non-symmetric Laplacian
symmetric = False
# which eigenvalues we want to compute
# (SM = smallest magnitude, LM = largest magnitude)
which = 'SM'
# tolerance
toll = 1e-6

# PAGERANK
n_iter = 5
damping_factor = 0.85
solver = 'piteration'
# tolerance
tol = 1e-6

# AFFINITY
# lower-bound to the hyperedge sizes to consider
min_size = 2
# upper-bound to the hyperedge sizes to consider
max_size = 15

# GENEPY
# 'head' means we consider the export matrix
# 'tail' means we consider the import matrix
side = 'head'  # could be tail as well
# which index to compute between 'countries' and 'products'
group = 'countries'

# NON-LINEAR CONTAGION
# number of lambda values to test
nb_points = 20
# lower bound to lambda
lower_param = 0.01
# recovery rate
recovery_rate = 1
# initial density of infected nodes
# original paper uses 0.8 for upper branch and 0.02 for lower branch
init_I_list = [0.01]
# history of past states for the quasistationary-state method
nb_history = 50
# decorrelation period for sampling states
burnin_dec_dt = 1
# decorrelation period for measuring states
measure_dec_dt = 10

nl_params = {
    'lyon': {
        'shape_infection': [1, 4],  # exponent of the beta function
        'higher_param': 2,  # highest value of lambda
        'quasi': True,  # run the quasi-state method (only for small graphs)
        'burnin_dt': 10000,  # burn-in period
        'measure_dt': 10000  # sampling period
    },
    'high': {
        'shape_infection': [1, 4],  # exponent of the beta function
        'higher_param': 2,  # highest value of lambda
        'quasi': True,  # run the quasi-state method (only for small graphs)
        'burnin_dt': 10000,  # burn-in period
        'measure_dt': 10000  # sampling period
    },
    'coauth': {
        'shape_infection': [1, 2],  # exponent of the beta function
        'higher_param': 10,  # highest value of lambda
        'quasi': False,  # run the quasi-state method (no for large graphs)
        'burnin_dt': 1000,  # burn-in period
        'measure_dt': 1000  # sampling period
    },
    'email-Eu': {
        'shape_infection': [1, 3],  # exponent of the beta function
        'higher_param': 5,  # highest value of lambda
        'quasi': True,  # run the quasi-state method (no for large graphs)
        'burnin_dt': 10000,  # burn-in period
        'measure_dt': 10000  # sampling period
    },
    'email-Enron': {
        'shape_infection': [1, 3],  # exponent of the beta function
        'higher_param': 5,  # highest value of lambda
        'quasi': False,  # run the quasi-state method (no for large graphs)
        'burnin_dt': 10000,  # burn-in period
        'measure_dt': 10000  # sampling period
    }
}

# frequency thresholds for convergence experiment
fi_threshs = {'eco01100': '1e-3,1e-3',
              'dblp_v9': '1e-3,1e-4',
              'enron': '1,1e-3',
              'qna_math': '1,5e-5',
              'uspto19762016': '1e-3,1',
              'metabolic_iaf1260b': '1e-4,1e-4',
              'metabolic_ijo1366': '1e-4,1e-4',
              'citation_software': '1e-4,1e-4',
              'hs1995_thresh': '5e-2,1e-1',
              'hs2009_thresh': '5e-2,1e-1',
              'hs2019_thresh': '5e-2,1e-1',
              'hs2020_thresh': '5e-2,1e-1'
              }

# names of the samplers
sampler_names = {'samosa.adri.NuDHy_A': 'NuDHy-Degs',
                 'samosa.adri.NuDHy_B': 'NuDHy-Degs-MH',
                 'samosa.adri.NuDHy_C': 'NuDHy-BIOT',
                 'samosa.samplers.NuDHy_Degs': 'NuDHy-Degs',
                 'samosa.samplers.NuDHy_BIOT': 'NuDHy-BIOT',
                 'NuDHy_A': 'NuDHy-Degs',
                 'NuDHy_B': 'NuDHy-Degs-MH',
                 'NuDHy_C': 'NuDHy-BIOT',
                 'NuDHy_Degs': 'NuDHy-Degs',
                 'NuDHy_BIOT': 'NuDHy-BIOT'}

# user-friendly dataset names
dataset_names = {'eco01100': 'ecoli',
                 'dblp_v9': 'DBLP-9',
                 'dblp_v9_rest': 'DBLP-9',
                 'bitcoin_2016': 'Bit-2016',
                 'enron': 'Enron',
                 'qna_math': 'Math',
                 'uspto19762016': 'ORD',
                 'lyon': 'Primary',
                 'high': 'High',
                 'coauth': 'coauth-DBLP',
                 'metabolic_iaf1260b': 'iaf1260b',
                 'metabolic_ijo1366': 'ijo1366',
                 'citation_software': 'Cit-SW',
                 'hs1995_thresh': 'HS1995',
                 'hs2009_thresh': 'HS2009',
                 'hs2019_thresh': 'HS2019',
                 'hs2020_thresh': 'HS2020',
                 }

#  number of swaps to perform to generate the samples
swaps_num = {'eco01100': 79080,
             'metabolic_iaf1260b': 177680,
             'metabolic_ijo1366': 193500,
             'dblp_v9': 9291320,
             'enron': 14819700,
             'qna_math': 5208700,
             'uspto19762016': 53179340,
             'citation_software': 6002840,
             'dblp_v9_rest': 1660900,
             'hs1995_thresh': 12016350,
             'hs2009_thresh': 12724250,
             'hs2019_thresh': 12364000,
             'hs2020_thresh': 12284050,
             'lyon': 57060,
             'high': 363840,
             'email-Enron': 227500,
             'email-Eu': 1714740,
             'coauth': 4308680,
             'congress_bills_S_093': 283040,
             'congress_bills_S_094': 231780,
             'congress_bills_S_095': 220840,
             'congress_bills_S_096': 230700,
             'congress_bills_S_097': 332960,
             'congress_bills_S_098': 431980,
             'congress_bills_S_099': 474940,
             'congress_bills_S_100': 576020,
             'congress_bills_S_101': 611680,
             'congress_bills_S_102': 557000,
             'congress_bills_S_103': 367460,
             'congress_bills_S_104': 208620,
             'congress_bills_S_105': 279480,
             'congress_bills_S_106': 395200,
             'congress_bills_S_107': 349740,
             'congress_bills_S_108': 359820,
             'congress_bills_H_093': 1040100,
             'congress_bills_H_094': 1132420,
             'congress_bills_H_095': 1304020,
             'congress_bills_H_096': 1429320,
             'congress_bills_H_097': 1681000,
             'congress_bills_H_098': 2318980,
             'congress_bills_H_099': 2579480,
             'congress_bills_H_100': 2704920,
             'congress_bills_H_101': 2970040,
             'congress_bills_H_102': 2700440,
             'congress_bills_H_103': 2052860,
             'congress_bills_H_104': 1342800,
             'congress_bills_H_105': 1738700,
             'congress_bills_H_106': 2153440,
             'congress_bills_H_107': 2131980,
             'congress_bills_H_108': 2089200}
