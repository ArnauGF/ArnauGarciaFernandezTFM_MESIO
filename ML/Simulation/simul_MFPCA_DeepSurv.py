#Trying to read simulated data frames from an R script and analyizing with Python
PATH = "/Users/arnaugarcia/Desktop/TFM_super/ML"
import torch
import torch.optim as optim

import sys
sys.path.append(PATH)
from Models.MFPCA_DeepSurv.DeepSurv import DeepSurv
from Models.MFPCA_DeepSurv.functions import (get_numpy, get_numpy2, sortByTime, surv_loss, BreslowEstimator)
from Models.metrics import (Brier, Brier2, IBS)
from Simulation.data_simulation_base import simulate_JM_base
#from Simulation.data_simulation_nonPH import simulate_JM_nonPH

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

#We load the simulated datasets (via R) and we analyze them by mens of MFPCA-DeepSurv 

#Load the datasets from .RData file
import rpy2.robjects as ro
rdata_file = '/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/SceI_L5_LLMs_100df.RData'
ro.r['load'](rdata_file)
r_objects = ro.r.ls()
#Now we have the data sets load in the Python environment and we can work with them

#The number of data sets:
num_dat = 100

# Creating arrays to save relevant results
censoring_train = np.zeros(num_dat)
censoring_test = np.zeros(num_dat)
IBS_train = np.zeros(num_dat)
IBS_test = np.zeros(num_dat)
BS_train = np.zeros(num_dat)
BS_test = np.zeros(num_dat)

#Stuff to compute time of computation:
import pickle
import time
start = time.time()

for i in range(1, num_dat+1):
    train_data = ro.globalenv[f'DF_{i}']
    train_data = pandas2ri.rpy2py_dataframe(train_data)
    test_data = ro.globalenv[f'DF_test_{i}']
    test_data = pandas2ri.rpy2py_dataframe(test_data)
    train_data_unique_id = train_data.drop_duplicates(subset='id')
    test_data_unique_id = test_data.drop_duplicates(subset = 'id')

    #We compute and save % of censoring:
    perc_cens_train = train_data_unique_id['event'].value_counts()[0]/train_data_unique_id.shape[0]
    perc_cens_test = test_data_unique_id['event'].value_counts()[0]/test_data_unique_id.shape[0]
    censoring_train[i-1] = perc_cens_train
    censoring_test[i-1] = perc_cens_test

    #Convert sex and treatment variables (baseline covariates) into [0,1]'s
    train_data['sex'] = train_data['sex'].replace({2: 0})
    test_data['sex'] = test_data['sex'].replace({2: 0})
    train_data['treatment'] = train_data['treatment'].replace({2: 0})
    test_data['treatment'] = test_data['treatment'].replace({2: 0})

    #Delete ind1, ind2 variable. 
    train_data = train_data.drop(['ind1', 'ind2'], axis=1)
    test_data = test_data.drop(['ind1', 'ind2'], axis=1)

    #Change the names to fit the names used in the python code?? sex->X1, y1 ->Y1, etc.
    train_data = train_data.rename(columns={'time': 'obstime', 'Time': 'time', 
                                        'y1':'Y1', 'y2':'Y2', 'y3':'Y3',
                                        'y4':'Y4', 'y5':'Y5', 'sex':'X1',
                                        'treatment':'X2'})
    test_data = test_data.rename(columns={'time': 'obstime', 'Time': 'time', 
                                        'y1':'Y1', 'y2':'Y2', 'y3':'Y3',
                                        'y4':'Y4', 'y5':'Y5', 'sex':'X1',
                                        'treatment':'X2'})

    #create the "predtime" variable (which is identical to obstime)
    train_data['predtime'] = train_data['obstime']
    test_data['predtime'] = test_data['obstime']

    minmax_scaler = MinMaxScaler(feature_range=(-1,1))
    train_data.loc[:,["Y1","Y2","Y3","Y4","Y5"]] = minmax_scaler.fit_transform(train_data.loc[:,["Y1","Y2","Y3","Y4","Y5"]])
    test_data.loc[:,["Y1","Y2","Y3","Y4","Y5"]] = minmax_scaler.transform(test_data.loc[:,["Y1","Y2","Y3","Y4","Y5"]])
    
    x_long_train, x_base_train, e_train, t_train, obs_time = get_numpy2(train_data)
    
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
        
    ## Train model
    torch.manual_seed(0)
    model = DeepSurv(x_train.shape[1])
    model = model.train()

    optimizer = optim.Adam(model.parameters())
        
    n_epoch = 12
    batch_size = 32
        
    loss_values = []
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

    S_func = BreslowEstimator().fit(model(x_train).detach().numpy(), e_train, t_train)

    LT = 2.5
    pred_windows = 1
    pred_times = LT + pred_windows

    # Only keep subjects with survival time > landmark time
    tmp_data_test = test_data.loc[test_data["time"]>LT,:]
    tmp_id_test = np.unique(tmp_data_test["id"].values)
    #tmp_all = train_data_all.loc[train_data_all["id"].isin(tmp_id),:]

    # Only keep longitudinal observations <= landmark time
    tmp_data_test = tmp_data_test.loc[tmp_data_test["obstime"]<=LT,:]

    x_long_tmp_test, x_base_tmp_test, e_tmp_test, t_tmp_test, obs_time_test = get_numpy2(tmp_data_test, max_len=len(obs_time))

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

    x_long_tmp_train, x_base_tmp_train, e_tmp_train, t_tmp_train, obs_time_train = get_numpy2(tmp_train, max_len=len(obs_time))

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

#Now we save the results in a .RData file:
results = pd.DataFrame({
    'BS_test': BS_test,
    'IBS_test': IBS_test,
    'BS_train': BS_train,
    'IBS_train': IBS_train
})

#in the end we should put:
end = time.time()
print("total time:", (end-start)/60)

# Save results as CSV files
results.to_csv("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/DeepSurv_sceI_5lmms_100df.csv", index=False)

### We spend a total amount of 2.67649 mins to do it! (50 datasets)
## SCENARIO II: it tooks 5.073203202088674 mins to do it with 100 data sets.
## SCENARIO I: 4.6806795199712115 mins with 100 data sets
