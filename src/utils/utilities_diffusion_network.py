"""
utilities_diffusion_network.py

This module provides utilities for performing network analysis of line-list data

Main Classes:
- DiffusionNetwork: network analysis class
"""

import numpy as np
import pandas as pd
from haversine import haversine_vector, Unit
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow_probability import distributions as tfd
import rasterio as rasterio
from tqdm import tqdm

import matplotlib.pyplot as plt

import networkx as nx

class DiffusionNetwork:
    """
    A class for performing network analysis on line-list epidemiological data.

    Attributes:
        alpha_mean (float): the 'loc' of alpha (temporal decay shape) in the truncated normal distribution
        alpha_scale (float): the 'scale' of alpha in the truncated normal distribution
        beta_mean (float): the 'loc' of beta (spatial decay shape) in the truncated normal distribution
        beta_scale (float): the 'scale' of beta in the truncated normal distribution
        epsilon_mean (float): the 'loc' of epsilon (unobserved source) in the truncated normal distribution 
        epsilon_scale (float) the 'scale' of epsilon in the truncated normal distribution
        kernel (str): Choice of spatial kernel ('Exponential', 'Gaussian', 'Linear')
        input_data_fn (str): path to input data
        parasite (str): parasite name ('P.falciparum', 'P.vivax')
        subset_switch (str): data subset for analysis ('P.falciparum', 'P.vivax', 'P.vivax_relapseReinfection')
        vulnerability_raster_fn (str): file path for vulnerability raster 
        t_mat (np.array): matrix of days between cases
        d_mat (np.array): matrix of distances between cases
        pi (tf.constant): numerical value of pi
        gamma (tf.constant): minimum time (days) between cases to be linked by infection

    Methods:
        initialise_new_run: Load line-list data and process to prepare model fit.
        neg_loglike_fixedEps: Negative log-likelihood function.
        neg_loglike_and_grad: call negative log-likelihood and get gradient.
        run_optimise: Fit network diffusion model.
        transmission_likelihood: Matrix of transmission likelihoods (i.e., un-normalised)
        get_Re: Calculate expected effective preproduction number from transmission network
        get_random_graph: Get a single realisation of a transmission network (unweighted DAG)
        get_weighted_random_graph: give nx.graph of probability network.
        get_index_cases: determine which cases have no obvious sources of infection.
        calculate_regression_weights: calculate sample weights for each case.
        simulate_outbreaks: simulate nsims random graphs and determine outbreak statistics 
        make_outbreak_size_date_fig: Plot outbreak-size distribution vs date.
        view_transmission_network: Not implemented.
        make_dual_fig: plot Re (histogram and time-line). 
        make_Re_fig: plot Re time-line. 
        make_timeline_original: plot time-line of case counts (as 'imported vs indigeneous').
    """
    def __init__(self):
        self.alpha_mean    = 0.003
        self.alpha_scale   = 0.001
        self.beta_mean     = 0.001
        self.beta_scale    = 0.001
        self.epsilon_mean  = 0.001
        self.epsilon_scale = 0.001
        self.kernel = 'Gaussian'
        self.input_data_fn = ''
        self.parasite = ''
        self.subset_switch = ''
        self.vulnerability_raster_fn = ''
        self.t_mat = np.array([0])
        self.d_mat = np.array([0])
        self.pi = tf.constant(np.pi, dtype=tf.float32)
        self.gamma = tf.constant(15, dtype=tf.float32)
        
    def initialise_new_run(self, subset_switch = 'P.falciparum', original_df = None):
        """
        Load data and compute inter-case relationships for network analysis.

        Args:
            subset_switch (str): Which data subset to use. One of ('P.falciparum', 'P.vivax', 
                                'P.vivax_relapseReinfection')
            original_df (pd.DataFrame): [optional] the input data (if self.input_data_fn unspecified)
        """
        self.subset_switch = subset_switch
        if type(original_df)==type(None):
            self.original_df = pd.read_csv(self.input_data_fn, encoding = "ISO-8859-1")
        elif type(original_df)==type(pd.DataFrame({'': []})):
            self.original_df = original_df 
        else: 
            print('Error: original_df argument needs to be a pandas data frame. Object not initialised')
            return
        
        self.original_df = self.original_df[['date_test', 'mp_species','date_onset','imported','LONG', 'LAT', 'case_classification', 'NAME_1', 'NAME_2', 'NAME_3', 'delay_days', 'occu_pt', 'age', 'gender_pt']]
        self.original_df = self.original_df.dropna()
        if subset_switch=='P.falciparum':
            self.parasite = 'P.falciparum'
            self.original_df = self.original_df.loc[self.original_df.mp_species == 'P.falciparum'].reindex()
        elif subset_switch=='P.vivax':
            self.parasite = 'P.vivax'
            self.original_df = self.original_df.loc[self.original_df.mp_species == 'P.vivax'].reindex()
        elif subset_switch=='P.vivax_relapseReinfection':
            #Treats a relapsed p.vivax as a reinfection (i.e. with a source in the data base).
            self.parasite = 'P.vivax'
            self.original_df = self.original_df.loc[(self.original_df.mp_species == 'P.vivax')].reindex()     
            self.original_df.loc[self.original_df.case_classification.str.lower() == 'relapse', 'imported'] = 'N'

        self.original_df = self.original_df.reset_index()
        self.original_df['date_onset']    = pd.to_datetime(self.original_df.date_onset)
        self.original_df['date_test']   = pd.to_datetime(self.original_df.date_test)
        self.original_df['relative_date'] = (self.original_df.date_onset- self.original_df.date_onset.min()).dt.days

        #Ignore rows whos test was > 40 days after onset of symptoms
        self.original_df = self.original_df.loc[(self.original_df.date_test - self.original_df.date_onset).dt.days <= 40]
        self.original_df = self.original_df.reset_index()
        
        #Get matrix of times between cases
        list_of_date_diffs = []
        for i in range(len(self.original_df['relative_date'])):
            series = self.original_df['relative_date'].iloc[i] - self.original_df['relative_date']
            series.loc[series<0] = 0
            list_of_date_diffs.append(list(series))

        self.t_mat = np.c_[list_of_date_diffs]
        self.t_mat = self.t_mat[self.original_df.index[self.original_df['imported']=='N'].tolist()]

        #get matrix of distances between cases
        loc_list = [(self.original_df['LAT'].iloc[i], self.original_df['LONG'].iloc[i]) for i in range(len(self.original_df['LAT']))]
        self.d_mat = haversine_vector(loc_list,loc_list, Unit.KILOMETERS, comb=True)
        self.d_mat = self.d_mat[self.original_df.index[self.original_df['imported']=='N'].tolist()]
        
        #define prior parameter distributions
        self.alpha_dist = tfp.distributions.TruncatedNormal(loc=self.alpha_mean, scale=self.alpha_scale, low=0.000001, high=1)
        self.beta_dist = tfp.distributions.TruncatedNormal(loc=self.beta_mean, scale=self.beta_scale, low=0.0000001, high=1)
        self.epsilon_dist = tfp.distributions.TruncatedNormal(loc=self.epsilon_mean, scale=self.epsilon_scale, low=0.000001, high=1)
        
        T =  tf.constant(self.t_mat, dtype=tf.float32)
        #Number of days 'row' was infected after 'column' became infectious (gamma day lag between infection and infectious)
        self.T_shifted = tf.nn.relu(T - self.gamma) # i.e. rectified linear (max(X, 0))
        T = tf.where(self.T_shifted >0.00001, T, 0)
    
        #Drop cases that are too far apart (exp(-10) \approx 5E-5) in space or time)
        self.D =  tf.constant(self.d_mat, dtype=tf.float32) 
        # self.D = tf.where(self.beta_mean*self.D*self.D < 10, self.D, 100) #TODO: clean here
        self.T_shifted = tf.where(self.T_shifted < 100, self.T_shifted, 0)

        self.mask_tensor = tf.where((self.T_shifted > 0.0001)&(self.beta_mean*self.D*self.D < 10), 1, 0)
        self.mask_indices = tf.where(self.mask_tensor==1)
        
        mask_np = self.mask_tensor.numpy().flatten()

        #Define inital parameter values for optimisation
        self.alpha_start = tf.Variable(tf.constant(0.9*self.alpha_mean*np.ones(self.T_shifted.shape[0]*self.T_shifted.shape[1])[mask_np==1], dtype=tf.float32))
        self.alpha_start_shape = self.alpha_start.shape[0]
        self.beta_start = tf.Variable(tf.constant(self.beta_mean*np.ones(1), dtype=tf.float32))
        return
        
    @tf.function
    def neg_loglike_fixedEps(self, param_vec):
        """
        Calculate negative loglikelihood of cascade given parameter values

        Args:
            param_vec (tf.Variable): parameter values [alpha, beta]

        Returns:
            nll_total (tf.Tensor): value of negative log likelihood
        """
        epsilon = tf.constant(self.epsilon_mean,shape=(self.T_shifted.shape[0],),dtype=tf.float32)
        alpha_sparse = tf.SparseTensor(self.mask_indices, values = param_vec[:self.alpha_start_shape], dense_shape = [self.T_shifted.shape[0],self.T_shifted.shape[1]])
        alpha_orig = tf.sparse.to_dense(alpha_sparse, default_value=0.00000015)
        beta_orig = param_vec[self.alpha_start_shape]

        #Enforce function evaluations within truncated distributions to prevent code crashing during optimisation
        #This is not a robust way to do this (i.e. no incentive to return to optimisation domain).
        alpha = tf.where(alpha_orig > 0.00000015, alpha_orig, 0.00000015)
        beta = tf.where(beta_orig > 0.00000015, beta_orig, 0.00000015)
        alpha = tf.where(alpha_orig < 1, alpha, 0.99999)
        beta = tf.where(beta_orig < 1, beta, 0.99999)

        H0 = epsilon*tf.constant(1.0, shape=(self.T_shifted.shape[0],),dtype=tf.float32) 
        logS0 = epsilon*tf.constant(-1.0,shape=(self.T_shifted.shape[0],), dtype=tf.float32)

        if self.kernel=='Exponential':
            H    = beta*alpha*self.T_shifted*tf.math.exp(-beta*self.D)
            S    = tf.math.exp(-1* alpha*self.T_shifted*self.T_shifted/2)/beta
            logS = -1* alpha*self.T_shifted*self.T_shifted/2 + tf.math.log(beta)
        elif self.kernel=='Gaussian':
            H    = 2*tf.math.sqrt(beta)*alpha*self.T_shifted*tf.math.exp(-beta*self.D*self.D)/tf.math.sqrt(self.pi)
            logS = -1*alpha*self.T_shifted*self.T_shifted/2 + tf.math.log(self.pi)/2 + tf.math.log(0.5) - 1*tf.math.log(beta)/2 
        elif self.kernel=='Linear':
            #Time only case
            H = alpha*self.T_shifted
            S = tf.math.exp(-alpha*self.T_shifted*self.T_shifted/2)
            logS = -alpha*self.T_shifted*self.T_shifted/2
        else:
            H = self.T_shifted
            S = self.T_shifted

        nll_total = -1*tf.reduce_sum(
                    logS0  + tf.reduce_sum(tf.where(self.mask_tensor ==1, logS, 0), axis=1) + tf.reduce_sum(tf.math.log(H0 + tf.reduce_sum(tf.where(self.mask_tensor ==1, H, 0), axis=1)))
                   ) 
        #Add prior likelihoods
        nll_total += -1*tf.reduce_sum(self.alpha_dist.log_prob(param_vec[:self.alpha_start_shape]))          
        nll_total += -1*tf.reduce_sum(self.beta_dist.log_prob(beta))

        return  nll_total
    
    # @tf.function
    # def neg_loglike_and_grad(self, param_vec):
    #     """
    #     Calculate negative loglikelihood and gradient

    #     Args:
    #         param_vec (tf.Variable): parameter values [alpha, beta]

    #     Returns:
    #         val (tf.Tensor): value of negative log likelihood
    #         g (tf.Tensor): gradient
    #     """
    #     with tf.GradientTape() as tape:
    #         tape.reset()
    #         tape.watch(param_vec)
    #         y = self.neg_loglike_fixedEps(param_vec)
    #         g = tape.gradient(y, param_vec,
    #                       unconnected_gradients=tf.UnconnectedGradients.ZERO)
    #     val = self.neg_loglike_fixedEps(param_vec)
    #     return (val, g)

    @tf.function
    def neg_loglike_and_grad(self, param_vec):
        """
        Calculate negative loglikelihood and gradient

        Args:
            param_vec (tf.Variable): parameter values [alpha, beta]

        Returns:
            val (tf.Tensor): value of negative log likelihood
            g (tf.Tensor): gradient
        """
        with tf.GradientTape() as tape:
            tape.watch(param_vec)
            loss = self.neg_loglike_fixedEps(param_vec)
        grad = tape.gradient(loss, param_vec,
                          unconnected_gradients=tf.UnconnectedGradients.ZERO)
        return (loss, grad)


    def run_optimise(self):
        """
        Run optimise function and save results in self.optim_results
        """
        start = tf.concat([self.alpha_start, self.beta_start], 0)
        self.optim_results = tfp.optimizer.lbfgs_minimize(
            self.neg_loglike_and_grad, 
            start, 
            tolerance=1e-8, 
            )
        return 
        
    def transmission_likelihood(self, param_vec):
        """
        Calculate the transmission likelihood matrix

        Args:
            param_vec (tf.Variable): optimal parameter values [alpha, beta]
        """
        epsilon = tf.constant(self.epsilon_mean,shape=(self.T_shifted.shape[0],),dtype=tf.float32)
        epsilon_mat = tf.constant(self.epsilon_mean,shape=(self.T_shifted.shape[0],self.T_shifted.shape[1]),dtype=tf.float32)

        alpha_sparse = tf.SparseTensor(self.mask_indices, values = param_vec[:self.alpha_start_shape], dense_shape = [self.T_shifted.shape[0],self.T_shifted.shape[1]])
        alpha_orig = tf.sparse.to_dense(alpha_sparse, default_value=0.00000015)

        beta_orig = param_vec[self.alpha_start_shape] 

        alpha = tf.where(alpha_orig > 0.00000015, alpha_orig, 0.00000015)
        beta = tf.where(beta_orig > 0.00000015, beta_orig, 0.00000015)
        alpha = tf.where(alpha_orig < 1, alpha, 0.99999)
        beta = tf.where(beta_orig < 1, beta, 0.99999)

        H = 2*tf.math.sqrt(beta)*alpha*self.T_shifted*tf.math.exp(-beta*self.D*self.D)/tf.math.sqrt(self.pi) #+ tf.constant(epsilon_mean,shape=(T_shifted.shape[0],T_shifted.shape[1]),dtype=tf.float32)

        S = tf.math.exp(-1*alpha*self.T_shifted*self.T_shifted/2)*tf.math.sqrt(self.pi)/2/tf.math.sqrt(beta)#*tf.math.exp(-epsilon_mat)

        F = H*S
        return F

    def get_Re(self):
        """
        Append estimate Re value to input data frame
        """
        paddings = [[0,0],[0,1]]
        trans_likelihood = self.transmission_likelihood(self.optim_results.position)
        trans_likelihood = tf.where(trans_likelihood > self.epsilon_mean, trans_likelihood, 0) #Chop small vals
        trans_likelihood = tf.pad(trans_likelihood, paddings, constant_values=self.epsilon_mean)
        self.trans_likelihood_norm = trans_likelihood/tf.reshape(tf.reduce_sum(trans_likelihood, axis=1), (-1, 1))
        self.Re = tf.reduce_sum(self.trans_likelihood_norm, axis=0)
        self.original_df['Re'] = self.Re.numpy()[:(len(self.Re.numpy())-1)]
        return
    
    def get_random_graph(self):
        """
        Sample a transmission network from the probability network

        Returns:
            graph (nx.Graph): graph indicating transmission events
        """
        eps=self.epsilon_mean*2
        chopped_graph = self.trans_likelihood_norm + tf.random.normal(self.trans_likelihood_norm.shape, 0, 0.001*eps, tf.float32)  #Add random noise to distinguish equally likely edges
        chopped_graph = tf.where(chopped_graph > eps,chopped_graph, 0)
        chopped_graph = tf.slice(chopped_graph, [0, 0], [-1, chopped_graph.shape[1] - 1])

        sampled_prob_list = []
        for j in range(chopped_graph.shape[0]):
            if tf.reduce_sum(chopped_graph[j]) > 0.1*eps:
                sampled_prob_list.append(chopped_graph[j][np.random.choice(chopped_graph.shape[1], 1, p=(chopped_graph[j]/tf.reduce_sum(chopped_graph[j])).numpy())[0]].numpy())
            else:
                 sampled_prob_list.append(1.1) #no match

        sampled_prob_list_tensor = tf.constant(sampled_prob_list)
        sampled_prob_list_tensor = tf.reshape(tf.repeat(sampled_prob_list_tensor, chopped_graph.shape[1]), (chopped_graph.shape[0], chopped_graph.shape[1]))
        sampled_prob_list_tensor = tf.where(sampled_prob_list_tensor<0.1*eps, 1.1, sampled_prob_list_tensor) #make zeroes in top row > 1 so there will be no match
        most_likely_transmission_network = tf.cast(tf.equal(chopped_graph, sampled_prob_list_tensor), tf.float32) #Only keep values equal to most likely , unfortunately there's multipel

        #augment to reintroduce imported cases (not infected by anyone so rows are zero
        replace_indices = self.original_df.loc[(self.original_df.imported=='N')].index

        origin_embeddings = tf.zeros((most_likely_transmission_network.shape[1], most_likely_transmission_network.shape[1]),
                                     dtype=tf.dtypes.float32,
                                     name=None)
        indices = tf.constant(replace_indices, dtype=tf.int32)
        most_likely_transmission_network_augmented = tf.tensor_scatter_nd_update(origin_embeddings,
                                                     tf.expand_dims(indices, 1),
                                                     most_likely_transmission_network)
        graph = nx.from_numpy_array(np.array(most_likely_transmission_network_augmented), create_using=nx.DiGraph)
        
        return graph
    
    def get_weighted_random_graph(self):
        """
        Get the weighted graph indicating transmission likelihoods
        """
        #A small parameter to chop improbable edges
        eps = 0.001
        # chopped_graph = tf.nn.relu(self.trans_likelihood_norm - eps) 
        chopped_graph = tf.where(self.trans_likelihood_norm > eps, self.trans_likelihood_norm, 0)
        chopped_graph = tf.slice(chopped_graph, [0, 0], [-1, chopped_graph.shape[1] - 1])

        #augment to reintroduce imported cases (not infected by anyone so rows are zero)
        replace_indices = self.original_df.loc[(self.original_df.imported=='N')].index

        origin_embeddings = tf.zeros((chopped_graph.shape[1], chopped_graph.shape[1]),
                                     dtype=tf.dtypes.float32,
                                     name=None)
        indices = tf.constant(replace_indices, dtype=tf.int32)
        transmission_network_augmented = tf.tensor_scatter_nd_update(origin_embeddings,
                                                     tf.expand_dims(indices, 1),
                                                     chopped_graph)
        graph = nx.from_numpy_array(np.array(transmission_network_augmented), create_using=nx.DiGraph)
        
        return graph

    def get_index_cases(self):
        """
        Determine which cases have no obvious source of infection, group with imported cases.
        Append to input data frame.
        """

        weighted_graph = self.get_weighted_random_graph()
        edges_and_weights = list(weighted_graph.edges(data=True))
        indices_in_time_order = self.original_df.sort_values('relative_date').index.values

        is_index_case = dict()

        for idx in tqdm(indices_in_time_order, desc='Testing for parents'):
            incoming_edges = [(x[1], x[2]) for x in  edges_and_weights if idx==x[0]]
            if len(incoming_edges)==0:
                is_index_case[idx] = 1
            else:
                is_index_case[idx] = 0

        self.original_df['is_index_case'] = 1
        for idx in self.original_df.index.values:
            self.original_df.loc[self.original_df.index==idx, 'is_index_case'] = is_index_case[idx] 
        return
    
    def calculate_regression_weights(self):
        """
        Calculate regression weights and append to dataframe
        """
        raster_object = rasterio.open(self.vulnerability_raster_fn)
        
        weighted_graph = self.get_weighted_random_graph()
        edges_and_weights = list(weighted_graph.edges(data=True))
        indices_in_time_order = self.original_df.sort_values('relative_date').index.values

        re_dict = dict()
        probability_dict = dict()
        is_index_case = dict()
        for idx in tqdm(indices_in_time_order):
            re_dict[idx] = self.original_df.loc[self.original_df.index==idx].Re.values[0]
            coords = [(self.original_df.loc[self.original_df.index==idx].LONG.values[0], self.original_df.loc[self.original_df.index==idx].LAT.values[0])]
            
            incoming_edges = [(x[1], x[2]) for x in  edges_and_weights if idx==x[0]]
            if len(incoming_edges)==0:
                #TODO: get from probability
                probability_dict[idx] = [x[0] for x in raster_object.sample(coords)][0]
                # probability_dict[idx] = [1 for x in raster_object.sample(coords)][0]
                is_index_case[idx] = 1
            else:
                probability_dict[idx] = sum([re_dict[source]*probability_dict[source] * weight['weight'] for source, weight in incoming_edges])
                is_index_case[idx] = 0


        self.original_df['likelihood'] = 0.0
        self.original_df['is_index_case'] = 1
        for idx in self.original_df.index.values:
            self.original_df.loc[self.original_df.index==idx, 'likelihood'] = probability_dict[idx] 
            self.original_df.loc[self.original_df.index==idx, 'is_index_case'] = is_index_case[idx] 
        self.original_df['probability'] = self.original_df['likelihood']/self.original_df.likelihood.sum(axis=0)
        self.original_df['inv_probability'] = 1/self.original_df['probability']
        return 
    
    def simulate_outbreaks(self, nsims=10):
        """
        Simulate transmission pathways to estimate outbreak size statistics

        Args:
            nsims (int): number of networks to simulate
        """
        #Check sorting
        new_column_list = [pd.DataFrame({'date':[], 'size':[], 'simulation_num':[]})]
        for i in tqdm(range(nsims)):
            graph = self.get_random_graph()

            size_of_sub_graph = [len(c) for c in sorted(nx.weakly_connected_components(graph), key=len, reverse=True)]

            parent_node_of_each_component = [min(c) for c in sorted(nx.weakly_connected_components(graph), key=len, reverse=True)]

            relative_date_list = [self.original_df.loc[self.original_df.index==parent_node, 'date_onset'].values[0] for parent_node in parent_node_of_each_component]
            size_of_each_component_sorted = [x[1] for x in sorted(zip(relative_date_list, size_of_sub_graph))]
            relative_date_list_sorted = sorted(relative_date_list)


            new_col_name_size = 'size_' + str(i)
            new_col_name_date = 'date_' + str(i)

            new_column_list.append(pd.DataFrame({'date': relative_date_list_sorted, 'size': size_of_each_component_sorted, 'simulation_num': [i for j in range(len(relative_date_list_sorted))]}))
        
        out_break_size_df = pd.concat(new_column_list, axis=0)

        date_list = []
        mean_list = []
        uci_list = []
        lci_list = []
        for date in sorted(out_break_size_df.date.unique()):
            date_list.append(date)
            mean_list.append(out_break_size_df.loc[out_break_size_df.date==date, 'size'].mean())
            uci_list.append(out_break_size_df.loc[out_break_size_df.date==date, 'size'].quantile(0.975))
            lci_list.append(out_break_size_df.loc[out_break_size_df.date==date, 'size'].quantile(0.025))

        self.date_list = np.array(date_list)
        self.mean_list = np.array(mean_list)
        self.uci_list = np.array(uci_list)
        self.lci_list = np.array(lci_list)

        return 

    def random_discrete_Re_estimate(self):

        eps=self.epsilon_mean*2
        chopped_graph = self.trans_likelihood_norm + tf.random.normal(self.trans_likelihood_norm.shape, 0, 0.001*eps, tf.float32)  #Add random noise to distinguish equally likely edges
        chopped_graph = tf.where(chopped_graph > eps,chopped_graph, 0)
        chopped_graph = tf.slice(chopped_graph, [0, 0], [-1, chopped_graph.shape[1] - 1])

        sampled_prob_list = []
        for j in range(chopped_graph.shape[0]):
            if tf.reduce_sum(chopped_graph[j]) > 0.1*eps:
                sampled_prob_list.append(chopped_graph[j][np.random.choice(chopped_graph.shape[1], 1, p=(chopped_graph[j]/tf.reduce_sum(chopped_graph[j])).numpy())[0]].numpy())
            else:
                 sampled_prob_list.append(1.1) #no match

        sampled_prob_list_tensor = tf.constant(sampled_prob_list)
        sampled_prob_list_tensor = tf.reshape(tf.repeat(sampled_prob_list_tensor, chopped_graph.shape[1]), (chopped_graph.shape[0], chopped_graph.shape[1]))
        sampled_prob_list_tensor = tf.where(sampled_prob_list_tensor<0.1*eps, 1.1, sampled_prob_list_tensor) #make zeroes in top row > 1 so there will be no match
        most_likely_transmission_network = tf.cast(tf.equal(chopped_graph, sampled_prob_list_tensor), tf.float32) #Only keep values equal to most likely , unfortunately there's multipel

        #augment to reintroduce imported cases (not infected by anyone so rows are zero
        replace_indices = self.original_df.loc[(self.original_df.imported=='N')].index

        origin_embeddings = tf.zeros((most_likely_transmission_network.shape[1], most_likely_transmission_network.shape[1]),
                                     dtype=tf.dtypes.float32,
                                     name=None)
        indices = tf.constant(replace_indices, dtype=tf.int32)
        most_likely_transmission_network_augmented = tf.tensor_scatter_nd_update(origin_embeddings,
                                                     tf.expand_dims(indices, 1),
                                                     most_likely_transmission_network)

# pt2, get Re
        # paddings = [[0,0],[0,1]]
        # trans_likelihood = self.transmission_likelihood(self.optim_results.position)
        # trans_likelihood = tf.where(trans_likelihood > self.epsilon_mean, trans_likelihood, 0) #Chop small vals
        # trans_likelihood = tf.pad(trans_likelihood, paddings, constant_values=self.epsilon_mean)
        # trans_likelihood_norm = trans_likelihood/tf.reshape(tf.reduce_sum(trans_likelihood, axis=1), (-1, 1))
        
        Re_vals = tf.reduce_sum(most_likely_transmission_network_augmented, axis=0)
        Re_array = np.array(Re_vals)

        return Re_array
    
    def make_outbreak_size_date_fig(self):
        """
        Plot outbreak size statistics time-line

        Returns:
            fig (plt.fig): outbreak size figure
        """
        fig = plt.figure()
        fig, ax = plt.subplots()
        ax.plot(self.date_list, self.mean_list)
        ax.fill_between(
            self.date_list, self.lci_list, self.uci_list, color='b', alpha=.15)
        ax.set_ylim(ymin=0)
        ax.set_title('Outbreak (community transmission) size over time')
        ax.set_xlabel('Date')
        ax.set_ylabel('Outbreak size')
        return fig

    
    def view_transmission_network(self):
        print('not implemented')
        return
    
    def make_dual_fig(self, upper_lim = 1000):
        """
        Make side-by-side figures of Re histogram and time-line

        Args:
            upper_lim (float): maximum Re value

        Returns:
            fig (plt.fig): time-line figure
        """
        fig = plt.figure()
        fig, axs = plt.subplots(1, 2, figsize = (12,4))
        Re_vals = self.Re.numpy()[:len(self.Re)-1]
        axs[0].hist(Re_vals[Re_vals<upper_lim], bins=20)
        axs[0].set_xlabel('Re')
        axs[0].set_ylabel('Count')
        plt.title('Re values ' + self.subset_switch)
        new_df = self.original_df.sort_values(by='Re')
        new_df = new_df.loc[new_df.Re < upper_lim]
        new_df.plot.scatter(x='date_onset', y='Re', ax=axs[1],c='red')
        plt.xticks(rotation=45);
        plt.suptitle('Re values ' +self.subset_switch);
        return fig
    
    def make_Re_fig(self, upper_lim = 1000):
        """
        Make time-line figure of Re values

        Args:
            upper_lim (float): maximum Re value

        Returns:
            fig (plt.fig): time-line figure
        """
        fig = plt.figure()
        new_df = self.original_df.sort_values(by='Re')
        new_df = new_df.loc[new_df.Re < upper_lim]
        # original_df.plot.scatter(x='date_onset', y='Re', ax=axs[1])
        new_df.plot.scatter(x='date_onset', y='Re',c='red', s=8)
        plt.xticks(rotation=30);
        plt.suptitle('Re values ' + self.subset_switch);
        return fig
    
    def make_timeline_original(self):
        """
        Make time-line figure of case counts 

        Returns: 
            fig (plt.fig): time-line figure
        TODO: make histograms
        """
        fig = plt.figure()
        full = self.original_df.date_onset.dt.to_period('M')
        not_imported = self.original_df.loc[self.original_df.imported=='N'].date_onset.dt.to_period('M')
        imported = self.original_df.loc[self.original_df.imported!='N'].date_onset.dt.to_period('M')

        full.value_counts().sort_index().plot(label='All cases')
        not_imported.value_counts().sort_index().plot(label='Indigenous')
        imported.value_counts().sort_index().plot(label='Imported')
        plt.title('Monthly ' + self.parasite +' cases')
        plt.legend()
        return fig


