############################################################
## MFPCA-DeepSurv in PBC dataset
############################################################

PATH = "/Users/arnaugarcia/Desktop/TFM_super/ML"

import sys
sys.path.append(PATH)
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sklearn.model_selection import GridSearchCV
from Models.MFPCA_DeepSurv.functions import (get_numpy, get_numpy2, get_numpy3, BreslowEstimator, sortByTime, surv_loss)
from Models.metrics import (AUC, Brier, Brier2, IBS, MSE)
from Simulation.data_simulation_base import simulate_JM_base
#from Simulation.data_simulation_nonPH import simulate_JM_nonPH

import torch
import torch.optim as optim
from Models.MFPCA_DeepSurv.DeepSurv import DeepSurv

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler, StandardScaler
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt


import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
ro.r.source("/Users/arnaugarcia/Desktop/TFM_super/ML/Models/MFPCA_DeepSurv/MFPCA.r")
mfpca_train = ro.globalenv['mfpca_train']
mfpca_test = ro.globalenv['mfpca_test']

#Load the PBC datasets from .RData file
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
rdata_file = '/Users/arnaugarcia/Desktop/TFM_super/ML/PBC/PBC_data.RData'
ro.r['load'](rdata_file)
r_objects = ro.r.ls()

train_data = ro.globalenv[f'DF']
train_data = pandas2ri.rpy2py_dataframe(train_data)
train_data_unique_id = train_data.drop_duplicates(subset='id')


#Convert sex and treatment variables (baseline covariates) into [0,1]'s
train_data['sex'] = train_data['sex'].replace({2: 0})
train_data['drug'] = train_data['drug'].replace({2: 0})


#Change the names to fit the names used in the python code?? sex->X1, y1 ->Y1, etc.
train_data = train_data.rename(columns={'year': 'obstime', 'years': 'time', 
                                    'serBilir':'Y1', 'albumin':'Y2', 'alkaline':'Y3',
                                    'platelets':'Y4', 'sex':'X1', 'drug':'X2', 'status2':'event'})

#create the "predtime" variable (which is identical to obstime)
train_data['predtime'] = train_data['obstime']

minmax_scaler = MinMaxScaler(feature_range=(-1,1))
train_data.loc[:,["Y1","Y2","Y3","Y4"]] = minmax_scaler.fit_transform(train_data.loc[:,["Y1","Y2","Y3","Y4"]])


x_long_train, x_base_train, e_train, t_train, obs_time = get_numpy3(train_data)

e_train = train_data_unique_id['status2'].to_numpy()
t_train = train_data_unique_id['years'].to_numpy()


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
model = DeepSurv(x_train.shape[1], dropout=0.2)
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

###########################
## first time window (5,8]
###########################
LT = 5
pred_windows = 3
pred_times = LT + pred_windows

# Only keep subjects with survival time > landmark time
tmp_train = train_data.loc[train_data["time"]>LT,:]
tmp_id_train = np.unique(tmp_train["id"].values)

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

###########################
## first time window (7,10]
###########################
LT = 7
pred_windows = 3
pred_times = LT + pred_windows

# Only keep subjects with survival time > landmark time
tmp_train = train_data.loc[train_data["time"]>LT,:]
tmp_id_train = np.unique(tmp_train["id"].values)

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





#######################################################################################
## Hyperparam tunning

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split

## DeepSurv model
class DeepSurv(nn.Module):
    def __init__(self, n_features, hidden_size=64, dropout=0.2):
        super().__init__()
        self.model = nn.Sequential(
            nn.Linear(n_features, hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.BatchNorm1d(hidden_size),
            nn.Linear(hidden_size, hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.BatchNorm1d(hidden_size),
            nn.Linear(hidden_size, 1)
        )

    def forward(self, x):
        return self.model(x)

# Dummy Cox loss (replace this!)
def surv_loss(risk_scores, events):
    return torch.mean(risk_scores)  # Replace with proper loss

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

# Assume you have: x_train, e_train, t_train as PyTorch tensors
# e_train: event, t_train: time
y_train = torch.stack((e_train_t, t_train_t), dim=1)

# Define model wrapper
estimator = DeepSurvEstimator(n_features=x_train.shape[1])

# Define hyperparameter grid
param_grid = {
    'hidden_size': [16, 32, 64],
    'dropout': [0.2, 0.4],
    'epochs': [12,14]
}

# GridSearchCV
gs = GridSearchCV(estimator, param_grid=param_grid, cv=3, verbose=2)
gs.fit(x_train, y_train)

print("Best Params:", gs.best_params_)
print("Best Score:", gs.best_score_)