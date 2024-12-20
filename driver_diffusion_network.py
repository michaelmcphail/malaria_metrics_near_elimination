"""
This script performs network analysis on line-list data using a network diffusion model. 
It is intended to be run as the main script to fit the network diffusion model, predict Re,
and (optional) calcualate weights used to adjust for sampling bias in regression.

Usage:
    python3 driver_diffusion_network.py --config='./config.json' [--get_weights]

Arguments:
    config            Path to the configuration file (required).
    --get_weights     Enable weight calculation (optional, default: False).

"""

# import numpy as np
import pandas as pd
import os as os
import argparse
import json as json

import sys


from src.utils import DiffusionNetwork


def load_config(config_path):
    with open(config_path, 'r') as file:
        config = json.load(file) 
    return config

def main(args):
    """
    Main function to fit diffusion network model, predict Re, (optionally) predict
    regression weights, and save outputs as specified by config file 

    args:
        config (str): Path to the config file (*.json).
        get_weights (bool): if True calculate weights for regression.
    """

    save_data = True
    #Gather arguments from config file.
    config = load_config(args.config)

    case_data_base_fn = config['case_data_base_fn']
    file_path_out = config['case_data_withRe_fn']
    pf_vulnerability_raster_fn = config['pf_vulnerability_index_raster_fn']
    pv_vulnerability_raster_fn = config['pv_vulnerability_index_raster_fn']
    country_name = config['country_name']
    country_iso = config['country_iso']

    dirname = os.path.dirname(file_path_out)
    if not os.path.exists(dirname):
        os.makedirs(dirname)


    df_dict = dict()
    for parasite in ['P.falciparum', 'P.vivax']:
        print('Fitting for ' + parasite)
        diff_net_obj = DiffusionNetwork()
        diff_net_obj.input_data_fn = case_data_base_fn
        diff_net_obj.initialise_new_run(subset_switch = parasite)
        print(parasite + ': fitting network model')
        diff_net_obj.run_optimise()
        diff_net_obj.get_Re()
        
        #Get index cases (cases without parents) for eventual use in  a LGCM model.
        #This is a slower step, TODO: check efficiency of algorithm
        print(parasite + ': getting index cases')
        diff_net_obj.get_index_cases()
        
        if args.get_weights:
            print(parasite + ': calculating regression weights')
            if parasite=='P.falciparum':
                diff_net_obj.vulnerability_raster_fn = pf_vulnerability_raster_fn
            else:
                diff_net_obj.vulnerability_raster_fn = pv_vulnerability_raster_fn
            diff_net_obj.calculate_regression_weights()

        df_dict[parasite] = diff_net_obj.original_df.copy(deep=True)

    if save_data:
        combined_df = pd.concat([df_dict[parasite] for parasite in df_dict])
        combined_df['month_int'] = list(pd.DatetimeIndex(combined_df.date_onset).month)
        combined_df['year'] = list(pd.DatetimeIndex(combined_df.date_onset).year)
        combined_df['admin_level'] = 'POINT'
        combined_df['ADMIN0'] = country_name
        combined_df['ISO'] = country_iso
        combined_df = combined_df.rename(columns={"NAME_1": "ADMIN1","NAME_2": "ADMIN2", "NAME_3": "ADMIN3"})

        
        out_name_subset = ['index', 'date_test', 'mp_species',
           'date_onset', 'imported', 'LONG', 'LAT', 'case_classification',
           'relative_date', 'Re', 'is_index_case', 'month_int', 'year',
           'admin_level', 'ADMIN0','ADMIN1', 'ADMIN2','ADMIN3','ISO', 'delay_days', 'occu_pt', 'age', 'gender_pt']
        if args.get_weights:
            out_name_subset += ['inv_probability']
        combined_df[out_name_subset].to_csv(file_path_out, index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to ...")
    parser.add_argument('--config', required=True, help='Path to the configuration file')
    parser.add_argument('--get_weights', action='store_true', help='Whether or not to calculate geospatial regression weights')

    args = parser.parse_args()

    print("Python Path:")
    print("\n".join(sys.path))

    main(args)
