>>conda activate tensorflow
>>python


import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn import linear_model, preprocessing

# From your Illumina 27k, 450k or EPIC array data, select the 20318 CpG sites from the file "CpGsites.csv" in the correct order.
cpgs = np.array(pd.read_pickle('/home/PERSONALE/katarzyn.kwiatkowsk2/AltumAge/multi_platform_cpgs.pkl'))
#df = pd.DataFrame(cpgs)
#df.to_csv("/home/gs66/AltumAge/CpGsites.csv", header=False, index=False)

# Load the BMIQCalibration-normalized methylation data.
#data = pd.read_pickle('/home/gs66/AltumAge/example_dependencies/example_data.pkl')
#data = pd.read_csv('home/gs66/Desktop/epicPainNet_ClocksUp/out_clocks_2023/BMIQcalibrated_beta_values_meat2.csv')
data = pd.read_csv('/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/input_file_for_AltumAge_noBMIQ.csv')
#data = pd.read_csv('/home/gs66/Desktop/BRIC_vaccineHBV/out_horvath_october2023/not_BMIQcalibrated_beta_values.csv')
#data = pd.read_csv('home/gs66/Desktop/epicPainNet_ClocksUp/out_clocks_2023/not_BMIQcalibrated_beta_values.csv')
print(data)

real_age = data.age
methylation_data = data[cpgs]
print(methylation_data)

# Load the scaler, which transforms the distribution of beta values of each CpG site to mean = 0 and variance = 1.
scaler = pd.read_pickle('/home/PERSONALE/katarzyn.kwiatkowsk2/AltumAge/scaler.pkl')
# load AltumAge.
AltumAge = tf.keras.models.load_model('/home/PERSONALE/katarzyn.kwiatkowsk2/AltumAge/AltumAge.h5')

# Scale the beta values of each CpG with sklearn robust scaler.
methylation_data_scaled = scaler.transform(methylation_data)

# replace nan's with 0
methylation_data_scaled = np.nan_to_num(methylation_data_scaled, nan=0.0)


# Predict AltumAge.
pred_age_AltumAge = AltumAge.predict(methylation_data_scaled).flatten()

#pred_age_AltumAge = np.append(pred_age_AltumAge, real_age.to_numpy(), axis=1)

#data.to_csv(r'/home/gs66/AltumAge/example_dependencies/input_data.csv', header=True, index=True, sep=',', mode='w')
#np.savetxt("/home/gs66/AltumAge/example_dependencies/output_preds.csv", pred_age_AltumAge, delimiter=",")

np.savetxt("/home/PERSONALE/katarzyn.kwiatkowsk2/EPIMODE/out_script0/epimode_AltumAge_output_preds_noBMIQ.csv", pred_age_AltumAge, delimiter=",")

# all 3 input versions : with BMIQ_meat2 calibration, BMIQ_meat1 and noBMIQ calibration resulted in the same output with AltumAge predictions
