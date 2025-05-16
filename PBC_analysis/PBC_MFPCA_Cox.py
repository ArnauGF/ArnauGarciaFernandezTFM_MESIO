PATH = "/Users/arnaugarcia/Desktop/TFM_super/ML"

import sys
sys.path.append(PATH)
from sksurv.linear_model import CoxPHSurvivalAnalysis
from Models.MFPCA_DeepSurv.functions import (get_numpy, get_numpy2, get_numpy3, BreslowEstimator)
from Models.metrics import (AUC, Brier, Brier2, IBS, MSE)
from Simulation.data_simulation_base import simulate_JM_base
#from Simulation.data_simulation_nonPH import simulate_JM_nonPH

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler, StandardScaler
pd.options.mode.chained_assignment = None

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


y_train = list(zip(e_train.flatten().tolist(),t_train.flatten().tolist()))
y_train = np.array(y_train,dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

## Train model
cox_fit = CoxPHSurvivalAnalysis().fit(x_train, y_train)

S_func = BreslowEstimator().fit(cox_fit.predict(x_train), e_train, t_train)

###########################
## first time window (5,8]
###########################
LT = 5
pred_windows = 3
pred_times = LT + pred_windows

### Doing the same with the train data
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


# Survival prediction        
x_train_tmp = np.concatenate((mfpca_scores_tmp_train, x_base_tmp_train), axis=1)
risk_train = cox_fit.predict(x_train_tmp)
S_hat_train = S_func.get_survival_function(risk_train)
surv_pred_train = np.array([si(pred_times) for si in S_hat_train])

bs_train = Brier2(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)
ibs_train = IBS(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)


###########################
## Second time window (7,10]
###########################
LT = 7
pred_windows = 3
pred_times = LT + pred_windows

### Doing the same with the train data
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


# Survival prediction        
x_train_tmp = np.concatenate((mfpca_scores_tmp_train, x_base_tmp_train), axis=1)
risk_train = cox_fit.predict(x_train_tmp)
S_hat_train = S_func.get_survival_function(risk_train)
surv_pred_train = np.array([si(pred_times) for si in S_hat_train])

bs_train = Brier2(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)
ibs_train = IBS(surv_pred_train, e_tmp_train, t_tmp_train, e_train, t_train, LT, pred_windows)


