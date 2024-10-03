## Helper functions to load CoMMpass data 
## This code is not meant to run "as-is." Instead, it is meant 
# as a guide for users to write their own helper functions
# Please email me if you need the "data" folder

## Created Nov 27, 2023
## Author: Anish Simhal

import pandas as pd
import numpy as np 
import scipy


def load_MMRF_clinical_data(): 
    """
    Load clinical data 
    
    Return
    patient_list: (list) subject ids 
    patients_mask: (list) subjects with clinical data
    event_duration: (list)
    censorlist: (list)
    """
    km_clinical_data_df = pd.read_csv("data/km_pfs_data_669.csv")
    censorlist = km_clinical_data_df['death_flag_list'].values
    km_clinical_data_df['days_since_last_visit_list'] = km_clinical_data_df['days_since_last_visit_list'].fillna(-1)
    event_duration = km_clinical_data_df['days_since_last_visit_list'].values/365
    
    if np.sum(km_clinical_data_df['mask_list']) != 659: 
        print("data loaded incorrectly") 
        
    print(np.sum(km_clinical_data_df['mask_list']))
    
    missing_subject_list = ['MMRF_2903','MMRF_2905','MMRF_2908','MMRF_2914','MMRF_2926',\
                            'MMRF_2938', 'MMRF_2939', 'MMRF_2941', 'MMRF_2946', 'MMRF_2947']
    
    patient_list = pd.read_csv('data/subject_list_669.csv')
    patient_list = patient_list['0'].values
    patient_list = patient_list[1:]
    number_of_patients = len(patient_list)
    
    
    patients_mask = np.ones(number_of_patients,)
    for n, patient_id in enumerate(patient_list): 
        if patient_id[0:9] in missing_subject_list: 
            patients_mask[n] = 0
    
    patients_mask = patients_mask > 0
    print(len(missing_subject_list))
    if len(missing_subject_list) != 10: 
        print("data not loaded correctly")

    
    for n, event in enumerate(event_duration): 
        if event < 0: 
            patients_mask[n] = 0
            
    patients_mask = patients_mask > 0
    print(sum(patients_mask))
    
    patient_list = patient_list[patients_mask]
    event_duration = event_duration[patients_mask]
    censorlist = censorlist[patients_mask]

    return patient_list, patients_mask, event_duration, censorlist


def load_rna_data(patients_mask): 
    data_all = scipy.io.loadmat('data/data_may4_2022.mat')
    rna_data = data_all['rna_data']
    rna_data = rna_data[:, patients_mask]

    for n in range(0, rna_data.shape[0]): 
        for m in range(0, rna_data.shape[1]): 
            rna_data[n, m] = rna_data[n, m][0][0]


    gene_list = data_all['gene_list']
    gene_list = [x.strip(' ') for x in gene_list]
    gene_list = np.array(gene_list) 
    print("length gene list: ", len(gene_list))
    print("size of rna_data: ", rna_data.shape)
    return rna_data, gene_list