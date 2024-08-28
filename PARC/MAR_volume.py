"""
Program to calculate the influence area and change in head due to introduction of MAR in the model domain.
The wells will be modeified at the given cell id locations (prepared with GIS and imported from SWAT model
[the centroid of the subbasins]). The recharge data is imported from swat model as the available surface 
runoff of each time step. A factor to this transient recharge will be defined such that the water table 
remains below 3 meter.

Author: Ranveer Kumar
        PMRF at Department of civil engineering
        IIT (BHU), varanasi 221005
        email: ranveerkumar.rs.civ20@iitbhu.ac.in

Requirements:

1. Csv file containing the cell ids to be modified in the modflow model.
2. Csv file containing sub basin runoff values as exported from swat model.

"""
from data_management_tools import *
import pandas as pd
from SWAT_output import *

data_well = pd.read_csv('Parc_values.csv')
recharge_rate = pd.read_csv('MAR_recharge_rates_stressp.csv').iloc[:,1:].values
rr_array = np.array(recharge_rate).T
# print(len(rr_array))

sub_id = data_well.iloc[:,2].values
cell_id_mar = data_well.iloc[:,1].values


print('Initializing the first Simulation')
initial_head, ini_drw = Initiate_model()
print(' Initiai simulation succesful')
well_id, well_prty, no_ofBC = h5_backup()
print('data backup created')

print("well manipulation started")

stored_volume = []
 # at given cell and recharge condition

for w in range(len(cell_id_mar)):
    Mar_rates = np.array(rr_array[w])
    print(Mar_rates)
    write_well(Mar_rates,cell_id_mar[w])
    print('Well data is written')
    mf_simulator()
    pure_h5(well_id, well_prty, no_ofBC)
    volume_stored, head_change = MAR_storage(initial_head,cell_id_mar[w])
    volume_pumped = np.sum(Mar_rates*90)
    stored_volume.append([sub_id[w], cell_id_mar[w], volume_pumped, volume_stored, head_change ])




columns = ['Sub','CellID','VolumePumpin','Volumestore','head_change']
df = pd.DataFrame(stored_volume, columns=columns)


df.to_csv('after_mar.csv')
