## Modules


from data_management_tools import *
from PARC_tools import *
from scipy.stats import linregress


## get the cells for PARC calculations
def get_marcells():
        
    """
    Active cells with no source and sink
    i.e. Active cells - boundry cells - river cells - pumping wells

    """
    ac =get_active_cell_id()
    bc = get_boundary_cells()
    sc = get_stream_cell_id()
    wc = get_pwell_ids()
    cells1 = [i for i in ac if i not in bc]
    cells2 = [i for i in cells1 if i not in sc]
    cells3 = [i for i in cells2 if i not in wc]
    print(f'Total no of cells for MAR: {len(cells3)}')

    return cells3


def get_cell_top_bot(wellid):
    # row, col, lay = get_ijk(wellid)

    top1d, x =read_arrays('top', layer=1)
    bot1d, x = read_arrays('bot', layer= 1)
    print(top1d[wellid-1])

    return top1d[wellid-1], bot1d[wellid-1]



## variables
start  =0
end = 2
############


cell_id = get_marcells()

if os.path.isfile(os.getcwd()+'/PARC.npy'):

    PARCain = list(np.load('PARC.npy'))
else:

    PARCain = []

# PARCain = []

for cellid in cell_id[(len(PARCain)+1) :]:
    print(f'{len(PARCain)+1}. Cellid: {cellid}\n\n')
    row, col, lay = get_ijk(cellid)
    ctop, cbot = get_cell_top_bot(wellid=cellid)
    Hper = ctop - 3
    tend=[365]
    skin_proper = [0.125,0.225,2500]
    well_proper = [cbot+7, cbot+1, row, col]
    
    value = OWIR(tend, cellid, Hper, 1.0, skin_proper, well_proper)
    PARCain.append([cellid, value[0][1]])
    # print(PARCain)
    np.save('PARC_AIN.npy',np.array(PARCain))
