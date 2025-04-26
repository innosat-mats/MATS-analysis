import pandas as pd
import glob
import pdb

csv_files = glob.glob('/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/mats_project/MATS-analysis/Aglaja/output/mats_location/Dec_22/*.csv')

dataframes = [pd.read_csv(file) for file in csv_files]
combined_df = pd.concat(dataframes, ignore_index=True)
combined_df.to_csv('/Users/aglajaroth/Documents/KIT/master/geoecology/master_thesis/mats_project/MATS-analysis/Aglaja/output/mats_location/Dec_22/December_22_mats.csv', index=False)
