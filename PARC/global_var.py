"""
Author Lilian Bosc 
Latest update: 13/12/2022

Initial variables which will be used in the code. These variables can be modified by the code main.py

commented
"""

import os
import json

# Default global variables (model properties)
filename = "modelain2.2" ## model name
model = "mf2k5" ## modflow executable


# ## MODFLOW
environment_file_path = os.getcwd()+"/"+filename+"_MODFLOW/"+filename+".h5"
discritization_file_path = os.getcwd()+"/"+filename+"_MODFLOW/"+filename+".dis"
heads_file_path = os.getcwd()+"/"+filename+"_MODFLOW/"+filename+".hed"
leakage_file_path = "data/"+filename+"_MODFLOW/"+filename+".ccf"
drawdown_file_path = os.getcwd()+"/"+filename+"_MODFLOW/"+filename+".drw"
