###################################################################
### MFPCA simulations
### Random censoring 60%
### simulations for MFPCA-DeepSurv
### Intervals (6, 7.5] and (4,5.5]
################################################################
PATH = "/Users/arnaugarcia/Desktop/TFM_super/ML"
import torch
import torch.optim as optim

import sys
sys.path.append(PATH)
from Models.MFPCA_DeepSurv.DeepSurv import DeepSurv
from Models.MFPCA_DeepSurv.functions import (get_numpy, get_numpy2, get_numpy3, sortByTime, surv_loss, BreslowEstimator)
from Models.metrics import (AUC, Brier, Brier2, IBS)
from Simulation.data_simulation_base import simulate_JM_base

# Hyperparameter tunning
import torch.nn as nn
from sklearn.model_selection import GridSearchCV
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split

## Auxiliar functions hyperparameter tuning:
## Selecting hyperparameters
class DeepSurvEstimator(BaseEstimator):
    def __init__(self, n_features=10, hidden_size=64, dropout=0.2, lr=0.001, epochs=12, batch_size=32):
        self.n_features = n_features
        self.hidden_size = hidden_size
        self.dropout = dropout
        self.lr = lr
        self.epochs = epochs
        self.batch_size = batch_size
        self.model = None

    def fit(self, X, y):
        self.model = DeepSurv(self.n_features, self.hidden_size, self.dropout)
        self.model.train()

        optimizer = optim.Adam(self.model.parameters(), lr=self.lr)

        e_train = y[:, 0]
        t_train = y[:, 1]

        for epoch in range(self.epochs):
            perm = torch.randperm(X.shape[0])
            for i in range(0, X.shape[0], self.batch_size):
                idx = perm[i:i+self.batch_size]
                if len(idx) <= 1:
                    continue
                x_batch = X[idx]
                e_batch = e_train[idx]
                optimizer.zero_grad()
                risk = self.model(x_batch)
                loss = surv_loss(risk, e_batch)
                loss.backward()
                optimizer.step()
        return self

    def score(self, X, y):
        self.model.eval()
        with torch.no_grad():
            e_val = y[:, 0]
            risk = self.model(X)
            loss = surv_loss(risk, e_val).item()
        return -loss  # Negative for GridSearchCV to maximize


import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None

import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
rpy2.robjects.numpy2ri.activate()
ro.r.source("/Users/arnaugarcia/Desktop/TFM_super/ML/Models/MFPCA_DeepSurv/MFPCA.r")
mfpca_train = ro.globalenv['mfpca_train']
mfpca_test = ro.globalenv['mfpca_test']

# Display DataFrame in a pretty table format
#from tabulate import tabulate
#print(tabulate(train_data, headers='keys', tablefmt='pretty'))
#print(tabulate(train_data_unique_id, headers='keys', tablefmt='pretty'))

#Load the ADMIN CENS 30% datasets from .RData file
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
rdata_file = '/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/datasets_scenarios/datasets_randCens_60_MFPCA_06may2025.RData'
ro.r['load'](rdata_file)
r_objects = ro.r.ls()
#Now we have the data sets load in the Python environment and we can work with them
list_DF_randCens = ro.globalenv['list_DF_randCens_60']
list_DF_test_randCens = ro.globalenv['list_DF_test_randCens_60']

#The number of data sets:
num_dat = 100

# Creating arrays to save relevant results
censoring_train = np.zeros(num_dat)
censoring_test = np.zeros(num_dat)
IBS_train = np.zeros(num_dat)
IBS_test = np.zeros(num_dat)
BS_train = np.zeros(num_dat)
BS_test = np.zeros(num_dat)
IBS_train2 = np.zeros(num_dat)
IBS_test2 = np.zeros(num_dat)
BS_train2 = np.zeros(num_dat)
BS_test2 = np.zeros(num_dat)
dropout_select = np.zeros(num_dat)
hidden_size_select = np.zeros(num_dat)
epochs_select = np.zeros(num_dat)

#Stuff to compute time of computation:
import pickle
import time
start = time.time()

for i in range(1, num_dat+1):
    if i == 67:
        continue 

    train_data = list_DF_randCens.rx2(i)
    train_data = pandas2ri.rpy2py(train_data)  

    test_data = list_DF_test_randCens.rx2(i)
    test_data = pandas2ri.rpy2py(test_data)

    #Convert sex and treatment variables (baseline covariates) into [0,1]'s
    train_data['sex'] = train_data['sex'].replace({2: 0})
    test_data['sex'] = test_data['sex'].replace({2: 0})
    train_data['treat'] = train_data['treat'].replace({2: 0})
    test_data['treat'] = test_data['treat'].replace({2: 0})

    #Change the names to fit the names used in the python code?? sex->X1, y1 ->Y1, etc.
    train_data = train_data.rename(columns={'time': 'obstime', 'Time': 'time', 
                                        'y1':'Y1', 'y2':'Y2', 'y3':'Y3',
                                        'y4':'Y4', 'sex':'X1',
                                        'treat':'X2'})
    test_data = test_data.rename(columns={'time': 'obstime', 'Time': 'time', 
                                        'y1':'Y1', 'y2':'Y2', 'y3':'Y3',
                                        'y4':'Y4', 'sex':'X1',
                                        'treat':'X2'})

    #create the "predtime" variable (which is identical to obstime)
    train_data['predtime'] = train_data['obstime']
    test_data['predtime'] = test_data['obstime']

    ## Datasets with unique id
    train_data_unique_id = train_data.drop_duplicates(subset='id')
    test_data_unique_id = test_data.drop_duplicates(subset = 'id')

    #We compute and save % of censoring:
    perc_cens_train = train_data_unique_id['event'].value_counts()[0]/train_data_unique_id.shape[0]
    perc_cens_test = test_data_unique_id['event'].value_counts()[0]/test_data_unique_id.shape[0]
    censoring_train[i-1] = perc_cens_train
    censoring_test[i-1] = perc_cens_test

    minmax_scaler = MinMaxScaler(feature_range=(-1,1))
    train_data.loc[:,["Y1","Y2","Y3","Y4"]] = minmax_scaler.fit_transform(train_data.loc[:,["Y1","Y2","Y3","Y4"]])
    test_data.loc[:,["Y1","Y2","Y3","Y4"]] = minmax_scaler.transform(test_data.loc[:,["Y1","Y2","Y3","Y4"]])
    
    x_long_train, x_base_train, e_train, t_train, obs_time = get_numpy3(train_data)

    e_train = train_data_unique_id['event'].to_numpy()
    t_train = train_data_unique_id['time'].to_numpy()
    
    mfpca_out = mfpca_train(x_long_train, obs_time)
    mfpca_scores = np.array(mfpca_out[0])

    Cms = mfpca_out[1]
    psis = mfpca_out[2]
    x_train = np.concatenate((mfpca_scores, x_base_train), axis=1)
    x_train = torch.FloatTensor(x_train)
    e_train_t = torch.FloatTensor(e_train)
    t_train_t = torch.FloatTensor(t_train)

    y_train = list(zip(e_train.flatten().tolist(),t_train.flatten().tolist()))
    y_train = np.array(y_train,dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

    #Hyperparameter tunning
    y_train2 = torch.stack((e_train_t, t_train_t), dim=1)

    # Define model wrapper
    estimator = DeepSurvEstimator(n_features=x_train.shape[1])

    # Define hyperparameter grid
    param_grid = {
        'hidden_size': [16, 32, 64],
        'dropout': [0.2, 0.4],
        'epochs': [12,14]
    }

    # GridSearchCV
    gs = GridSearchCV(estimator, param_grid=param_grid, cv=3, verbose=0)
    gs.fit(x_train, y_train2)

    dropout_select[i-1] = gs.best_params_["dropout"]
    hidden_size_select[i-1] = gs.best_params_["hidden_size"]
    epochs_select[i-1] = gs.best_params_["epochs"]
        
    ## Train model
    torch.manual_seed(0)
    model = DeepSurv(x_train.shape[1], hidden_size = gs.best_params_["hidden_size"], 
                     dropout = gs.best_params_["dropout"])
    model = model.train()

    optimizer = optim.Adam(model.parameters())
        
    n_epoch = gs.best_params_["epochs"]
    batch_size = 32
        
    loss_values = []
    try:
        for epoch in range(n_epoch):
            running_loss = 0
            permutation = torch.randperm(x_long_train.shape[0])
            for batch in range(0, x_long_train.shape[0], batch_size):
                optimizer.zero_grad()
                    
                indices = permutation[batch:batch+batch_size]
                batch_x, batch_e, batch_t = \
                    x_train[indices,:], e_train_t[indices], t_train_t[indices]
                batch_x, batch_e, batch_t = sortByTime(
                    batch_x, batch_e, batch_t)
                if len(indices)>1: #drop if last batch size is 1
                    risk = model.forward(batch_x)
                    loss = surv_loss(risk, batch_e)
                    loss.backward()
                    optimizer.step()
                    running_loss += loss
            loss_values.append(running_loss.tolist())
        plt.plot(loss_values)

    except Exception as e:
        print(f"Error in iteration {i}: {e}")
        continue

    S_func = BreslowEstimator().fit(model(x_train).detach().numpy(), e_train, t_train)

    ### INTERVAL (4,5.5]
    LT = 4
    pred_windows = 1.5
    pred_times = LT + pred_windows

    # Only keep subjects with survival time > landmark time
    tmp_data_test = test_data.loc[test_data["time"]>LT,:]
    tmp_id_test = np.unique(tmp_data_test["id"].values)
    #tmp_all = train_data_all.loc[train_data_all["id"].isin(tmp_id),:]

    # Only keep longitudinal observations <= landmark time
    tmp_data_test = tmp_data_test.loc[tmp_data_test["obstime"]<=LT,:]
    test_tmp_data_unique_id = tmp_data_test.drop_duplicates(subset='id')

    x_long_tmp_test, x_base_tmp_test, e_tmp_test, t_tmp_test, obs_time_test = get_numpy3(tmp_data_test, max_len=len(obs_time))

    e_tmp_test = test_tmp_data_unique_id['event'].to_numpy()
    t_tmp_test = test_tmp_data_unique_id['time'].to_numpy()

    y_tmp_test = list(zip(e_tmp_test.flatten().tolist(),t_tmp_test.flatten().tolist()))
    y_tmp_test = np.array(y_tmp_test,dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

    # MFPCA on testing data before LT
    mfpca_out_test = mfpca_test(x_long_train, x_long_tmp_test, obs_time_test, Cms, psis)
    mfpca_scores_test = np.array(mfpca_out_test[0])

    x_test = np.concatenate((mfpca_scores_test, x_base_tmp_test), axis=1)
    x_test = torch.FloatTensor(x_test)

    # Survival prediction        
    model_test = model.eval()
    risk_test = model_test(x_test)
    risk_test = risk_test.view(-1).detach().numpy()
    S_hat_test = S_func.get_survival_function(risk_test)
    surv_pred_test = np.array([si(pred_times) for si in S_hat_test])

    #We compute IBS (and BS) for test data and save it 
    bs_test = Brier2(surv_pred_test, e_tmp_test, t_tmp_test, e_train, t_train, LT, pred_windows)
    ibs_test = IBS(surv_pred_test, e_tmp_test, t_tmp_test, e_train, t_train, LT, pred_windows)
    BS_test[i-1] = bs_test
    IBS_test[i-1] = ibs_test

    ### Doing the same with the train data
    # Only keep subjects with survival time > landmark time
    tmp_train = train_data.loc[train_data["time"]>LT,:]
    tmp_id_train = np.unique(tmp_train["id"].values)
    #tmp_all = train_data_all.loc[train_data_all["id"].isin(tmp_id),:]

    # Only keep longitudinal observations <= landmark time
    tmp_train = tmp_train.loc[tmp_train["obstime"]<=LT,:]
    train_tmp_data_unique_id = tmp_train.drop_duplicates(subset='id')

    x_long_tmp_train, x_base_tmp_train, e_tmp_train, t_tmp_train, obs_time_train = get_numpy3(tmp_train, max_len=len(obs_time))

    e_tmp_train = train_tmp_data_unique_id['event'].to_numpy()
    t_tmp_train = train_tmp_data_unique_id['time'].to_numpy()

    y_tmp_train = list(zip(e_tmp_train.flatten().tolist(),t_tmp_train.flatten().tolist()))
    y_tmp_train = np.array(y_tmp_train,dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

    # MFPCA on training data before LT
    mfpca_out_tmp_train = mfpca_test(x_long_train, x_long_tmp_train, obs_time_train, Cms, psis)
    mfpca_scores_tmp_train = np.array(mfpca_out_tmp_train[0])

    x_test_train = np.concatenate((mfpca_scores_tmp_train, x_base_tmp_train), axis=1)
    x_test_train = torch.FloatTensor(x_test_train)

    # Survival prediction        
    model_train = model.eval()
    risk_train = model_train(x_test_train)
    risk_train = risk_train.view(-1).detach().numpy()
    S_hat_train = S_func.get_survival_function(risk_train)
    surv_pred_train = np.array([si(pred_times) for si in S_hat_train])

    bs_train = Brier2(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)
    ibs_train = IBS(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)
    BS_train[i-1] = bs_train
    IBS_train[i-1] = ibs_train

    ### INTERVAL (6,7.5]

    LT = 6
    pred_windows = 1.5
    pred_times = LT + pred_windows

    # Only keep subjects with survival time > landmark time
    tmp_data_test = test_data.loc[test_data["time"]>LT,:]
    tmp_id_test = np.unique(tmp_data_test["id"].values)
    #tmp_all = train_data_all.loc[train_data_all["id"].isin(tmp_id),:]

    # Only keep longitudinal observations <= landmark time
    tmp_data_test = tmp_data_test.loc[tmp_data_test["obstime"]<=LT,:]
    test_tmp_data_unique_id = tmp_data_test.drop_duplicates(subset='id')

    x_long_tmp_test, x_base_tmp_test, e_tmp_test, t_tmp_test, obs_time_test = get_numpy3(tmp_data_test, max_len=len(obs_time))

    e_tmp_test = test_tmp_data_unique_id['event'].to_numpy()
    t_tmp_test = test_tmp_data_unique_id['time'].to_numpy()

    y_tmp_test = list(zip(e_tmp_test.flatten().tolist(),t_tmp_test.flatten().tolist()))
    y_tmp_test = np.array(y_tmp_test,dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

    # MFPCA on testing data before LT
    mfpca_out_test = mfpca_test(x_long_train, x_long_tmp_test, obs_time_test, Cms, psis)
    mfpca_scores_test = np.array(mfpca_out_test[0])

    x_test = np.concatenate((mfpca_scores_test, x_base_tmp_test), axis=1)
    x_test = torch.FloatTensor(x_test)

    # Survival prediction        
    model_test = model.eval()
    risk_test = model_test(x_test)
    risk_test = risk_test.view(-1).detach().numpy()
    S_hat_test = S_func.get_survival_function(risk_test)
    surv_pred_test = np.array([si(pred_times) for si in S_hat_test])

    #We compute IBS (and BS) for test data and save it 
    bs_test2 = Brier2(surv_pred_test, e_tmp_test, t_tmp_test, e_train, t_train, LT, pred_windows)
    ibs_test2 = IBS(surv_pred_test, e_tmp_test, t_tmp_test, e_train, t_train, LT, pred_windows)
    BS_test2[i-1] = bs_test2
    IBS_test2[i-1] = ibs_test2

    ### Doing the same with the train data
    # Only keep subjects with survival time > landmark time
    tmp_train = train_data.loc[train_data["time"]>LT,:]
    tmp_id_train = np.unique(tmp_train["id"].values)
    #tmp_all = train_data_all.loc[train_data_all["id"].isin(tmp_id),:]

    # Only keep longitudinal observations <= landmark time
    tmp_train = tmp_train.loc[tmp_train["obstime"]<=LT,:]
    train_tmp_data_unique_id = tmp_train.drop_duplicates(subset='id')

    x_long_tmp_train, x_base_tmp_train, e_tmp_train, t_tmp_train, obs_time_train = get_numpy3(tmp_train, max_len=len(obs_time))

    e_tmp_train = train_tmp_data_unique_id['event'].to_numpy()
    t_tmp_train = train_tmp_data_unique_id['time'].to_numpy()

    y_tmp_train = list(zip(e_tmp_train.flatten().tolist(),t_tmp_train.flatten().tolist()))
    y_tmp_train = np.array(y_tmp_train,dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

    # MFPCA on training data before LT
    mfpca_out_tmp_train = mfpca_test(x_long_train, x_long_tmp_train, obs_time_train, Cms, psis)
    mfpca_scores_tmp_train = np.array(mfpca_out_tmp_train[0])

    x_test_train = np.concatenate((mfpca_scores_tmp_train, x_base_tmp_train), axis=1)
    x_test_train = torch.FloatTensor(x_test_train)

    # Survival prediction        
    model_train = model.eval()
    risk_train = model_train(x_test_train)
    risk_train = risk_train.view(-1).detach().numpy()
    S_hat_train = S_func.get_survival_function(risk_train)
    surv_pred_train = np.array([si(pred_times) for si in S_hat_train])

    bs_train2 = Brier2(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)
    ibs_train2 = IBS(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)
    BS_train2[i-1] = bs_train2
    IBS_train2[i-1] = ibs_train2

    print(i)

#Now we save the results in a .RData file:
results = pd.DataFrame({
    'BS_test': BS_test,
    'IBS_test': IBS_test,
    'BS_train': BS_train,
    'IBS_train': IBS_train,
    'BS_test6_7': BS_test2,
    'IBS_test6_7': IBS_test2,
    'BS_train6_7': BS_train2,
    'IBS_train6_7': IBS_train2,
    'dropouts': dropout_select,
    'hidden_sizes': hidden_size_select,
    'epochs' : epochs_select
})

#in the end we should put:
end = time.time()
print("total time:", (end-start)/60)

# required time: 7mins

# Save results as CSV files
results.to_csv("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_randCen_60_15may2025.csv", index=False)

