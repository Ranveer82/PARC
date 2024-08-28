import global_var
from modules import *  
from data_management_tools import *
import numpy as np
from scipy.special import exp1
from scipy.stats import linregress


def permissible_head(Permissible_depth, cell_id):
    row=global_var.max_row
    col=global_var.max_col
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    top = h5_data.root['Arrays/top1']
    top_array = np.reshape(top,(row,col))
    # cell_ids = np.arange(1,row*col+1,1).reshape(row,col)

    r_array = top_array.copy()
    r_array[r_array != -9999999.0] -= Permissible_depth

    r_array1 = np.reshape(r_array, row*col)
    per_head = r_array1[cell_id]
    h5_data.close()
    return per_head

def free_cells():
    # returns cells with no source or sink in the model

    row=global_var.max_row
    col=global_var.max_col
    active_cells = np.array(get_active_cell_id())

    boundary_cells = np.array(get_boundary_cells())
    pumping_well_cells = np.array(get_pwell_ids())

    Stream_cells = np.array(get_stream_cell_id())
    drain_cells = np.array(get_drain_cell_id())

    comb_cells = np.array(set(np.concatenate((boundary_cells[boundary_cells<row*col],pumping_well_cells,Stream_cells,drain_cells), axis=0)))

    empty_cells  = np.setdiff1d(active_cells,comb_cells)

    return empty_cells

def head_diff(cell_id, Permissible_depth, head_array):
    
    # head_array, _ = read_head()
    cell_head = []
    t=0
    while t<len(head_array):
        cell_head.append(head_array[t][cell_id])
        t+=global_var.glo_n_layers

    head_mean = np.mean(np.array(cell_head))
    per_head = permissible_head(Permissible_depth, cell_id)
    delh =per_head-head_mean

    return delh

def possible_recharge(Bound,step,result_array):
    
    if step==0:
        q = Bound[0]
    elif step ==1:
        q = Bound[1]
    else:
        if result_array[-2]>0 and result_array[-1]<0:
            q = (Bound[0]+Bound[1])/2
            Bound[1]=q
        elif result_array[-2]<0 and result_array[-1]>0:
            q = (Bound[0]+Bound[1])/2
            Bound[0]=q
        elif result_array[-2]<0 and result_array[-1]<0:
            print('Dcreasing the lower bound')
            Bound[1] = Bound[0]
            Bound[0]=Bound[0]/2
            q = (Bound[0]+Bound[1])/2
        else:
            print('Incresing the upper bound')
            Bound[0] = Bound[1]
            Bound[1]=Bound[1]*2
            q = (Bound[0]+Bound[1])/2
    return Bound, q

def theis_solution(Q, T, S, r, t):
    """
    Compute drawdown at radial distance r and time t using Theis solution.
    
    Parameters:
        Q : float
            Pumping rate (L^3/T).
        T : float
            Transmissivity of the aquifer (L^2/T).
        S : float
            Storativity of the aquifer (dimensionless).
        r : float or array-like
            Radial distance from the well (L).
        t : float or array-like
            Time since pumping began (T).
    
    Returns:
        float or array-like:
            Drawdown at radial distance r and time t.
    """
    u = r**2 * S / (4 * T * t)
    W_u = exp1(u)  # Exponential integral function
    drawdown = Q / (4 * np.pi * T) * W_u
    return drawdown



def hantush_jacob_solution(Q, T, S, r, t, b):
    """
    Compute drawdown at radial distance r and time t using Hantush-Jacob solution. for leaky aquifers
    
    Parameters:
        Q : float
            Pumping rate (L^3/T).
        T : float
            Transmissivity of the aquifer (L^2/T).
        S : float
            Storativity of the aquifer (dimensionless).
        r : float or array-like
            Radial distance from the well (L).
        t : float or array-like
            Time since pumping began (T).
        b : float
            Ratio of vertical to horizontal hydraulic conductivity.
    
    Returns:
        float or array-like:
            Drawdown at radial distance r and time t.
    """
    u = r**2 * S / (4 * T * t)
    v = r**2 * S / (4 * T * t) * (1 - b**2 / r**2)
    W_u = exp1(u)  # Exponential integral function
    W_v = exp1(v)  # Exponential integral function
    drawdown = Q / (4 * np.pi * T) * W_u + Q * b / (4 * np.pi * T) * W_v
    return drawdown



def cooper_jacob_approximation(Q, T, S, r, t, b):
    """
    Compute drawdown at radial distance r and time t using Cooper-Jacob approximation.
    
    Parameters:
        Q : float
            Pumping rate (L^3/T).
        T : float
            Transmissivity of the aquifer (L^2/T).
        S : float
            Storativity of the aquifer (dimensionless).
        r : float or array-like
            Radial distance from the well (L).
        t : float or array-like
            Time since pumping began (T).
        b : float
            Ratio of vertical to horizontal hydraulic conductivity.
    
    Returns:
        float or array-like:
            Drawdown at radial distance r and time t.
    """
    u = r**2 * S / (4 * T * t)
    W_u = exp1(u)  # Exponential integral function
    drawdown = Q / (4 * np.pi * T) * W_u * (1 + 2 * b**2 / r**2)
    return drawdown

def papadopulos_cooper_solution(Q, T, S, r, t, b):
    """
    Compute drawdown at radial distance r and time t using Papadopulos-Cooper solution.
    
    Parameters:
        Q : float
            Pumping rate (L^3/T).
        T : float
            Transmissivity of the aquifer (L^2/T).
        S : float
            Storativity of the aquifer (dimensionless).
        r : float or array-like
            Radial distance from the well (L).
        t : float or array-like
            Time since pumping began (T).
        b : float
            Ratio of vertical to horizontal hydraulic conductivity.
    
    Returns:
        float or array-like:
            Drawdown at radial distance r and time t.
    """
    u = r**2 * S / (4 * T * t)
    v = r**2 * S / (4 * T * t) * (1 - b**2 / r**2)
    W_u = exp1(u)  # Exponential integral function
    W_v = exp1(v)  # Exponential integral function
    drawdown = Q / (4 * np.pi * T) * W_u + Q * b / (4 * np.pi * T) * W_v * (1 - b**2 / r**2)
    return drawdown

def dupuit_forchheimer_solution(Q, T, b, rw, r):
    """
    Compute drawdown at radial distance r from a fully penetrating well using Dupuit-Forchheimer assumptions.
    
    Parameters:
        Q : float
            Pumping rate (L^3/T).
        T : float
            Transmissivity of the aquifer (L^2/T).
        b : float
            Thickness of the aquifer (L).
        rw : float
            Radius of the well (L).
        r : float or array-like
            Radial distance from the well (L).
    
    Returns:
        float or array-like:
            Drawdown at radial distance r.
    """
    drawdown = Q / (2 * np.pi * T) * np.log(r / rw)
    return drawdown



def change_array(par_name,par_value,layers):
    """
    par_name= paramter to be updated (HANI, HK, VANI, TOP etc)

    par_value =  layer wise parameter values [x1,x2,x3....xn]
    layers = layers in which the par_values to be updated [1,2,3,4...n]

    """
    go_to_the_root()
    dis = get_discretization()
    row = dis['nrow']
    col = dis['ncol']
    # print(str_p)
    for ln in range(len(layers)):
        param_array = np.ones((row*col))*par_value[ln]
        param_id  = par_name+str(int(layers[ln]))
        write_arrays(param_array, param_id)
    
    print(f'parameter {par_name} updated')

    return None

def simulate_head(Q, cellid,well_properties, skin_properties):
    dis = get_discretization()
    str_p = dis['nper']
    # print(str_p)

    well_rech = np.ones((1,str_p))*Q  # 1 for number of wells
    # well_rech
    go_to_the_root()
    write_MNW2(well_rech, [cellid], well_properties, skin_properties)
    print(f'simulated injection with recharge rate: {Q} m^3/min')

    mf_simulator()

    head1d= head_ts(cellid)
    return head1d[:,1]

def linear_fitted_head(h, Q, h_per):
    """
    h =  head array (1d) at least 2 element
    Q = recharge rate araay (1d) at least 2 element
    h_per = the head at which Q is desired 
    
    """
  

    slope, intercept, r_value, p_value, std_err = linregress(h, Q)

    Qmx = h_per*slope + intercept

    return Qmx, r_value, slope

def head_correction_confined(head_array, Q, a = 10  ,rw = 0.127 ,T = 0.006944):
    """
    corrects the simulated head for well with theim method
    a = cell size
    rw = radius of well
    T  = transmissivity of confined aquifer

    """
    re = 0.208*a  # equivalent radius Ac to theim solution
    ln  =np.log(re/rw)
    cf = (Q*ln)/(2*np.pi*T)
    head = head_array + cf
    return  head

def get_st_head(wellid):
    sth1d, x =read_arrays('StartHead', layer=1)
    return sth1d[wellid-1]


def OWIR(t_end_array, wellid, H_per, tolerance=0.01, skin_proper = [0.127,1.7,0.00011574], well_proper= [50.0, 0.0, 88, 90]):
    if H_per <= get_st_head(wellid):
        out =[[t, -999] for t in t_end_array]
        return out
    
    t_end = t_end_array #np.arange(3600,86400+3600,3600)
    ts = pd.read_csv('time_steps.csv').iloc[:,0].values # model time steps
    # wellid = 330703
    # H_per = 117.24
    tol = tolerance # 0.01

    bound = [100.0,10000.0]

    # skin_proper = [0.127,1.7,0.00011574] #[Rw,Rskin,Kskin]
    # well_proper = [30.0, 0.0, 88, 90] #[Ztop,Zbotm,ROW,COL]

    # head for lower bound
    Hl_arr = simulate_head(bound[0], wellid, well_proper, skin_proper)
    # Hl_arr_c = head_correction_confined(Hl_arr, bound[0])

    # head for upper bound
    Hu_arr = simulate_head(bound[1], wellid,well_proper, skin_proper)
    # Hu_arr_c = head_correction_confined(Hu_arr, bound[1])

    parc_output = []
    # parc_heads = []

    for te in t_end:
        te_index = np.argmin(np.abs(ts-te))
        # ts_e = ts[0:te_index]

        Hl = Hl_arr[te_index]
        Hu = Hu_arr[te_index]

        # fitting the linear relation and find out next Q
        h_arr = np.array([Hl,Hu])
        Q_arr = np.array(bound)
        Qmx , r2, slope=  linear_fitted_head(h_arr, Q_arr, H_per)

        # head at Qmax
        hmax_arr = simulate_head(Qmx, wellid, well_proper, skin_proper)
        # hmax_c = head_correction_confined(hmax_arr, Qmx)
        hmax = hmax_arr[te_index]

        # check the Qmax
        if abs(hmax-H_per)<=tol:
            PARC = Qmx
            parc_output.append([te, PARC])
            # parc_heads.append(hmax_arr)
            # print(PARC)
        else:
            print(f' residual : {hmax - H_per}')
            count =0
            while abs(hmax-H_per)>=tol and count<=30:
                if hmax<H_per:
                    h_arr = np.array([hmax, Hu])
                    Q_arr = np.array([Qmx, bound[1]])

                    Qmx , r2, slope=  linear_fitted_head(h_arr, Q_arr, H_per)

                    # Qmx += slope*abs(hmax-H_per)
                    hmax_arr=  simulate_head(Qmx, wellid,well_proper, skin_proper)
                    hmax = hmax_arr[te_index]
                else:
                    h_arr = np.array([Hl, hmax])
                    Q_arr = np.array([bound[0],Qmx])

                    Qmx , r2, slope=  linear_fitted_head(h_arr, Q_arr, H_per)

                    # Qmx -= slope*abs(hmax-H_per)
                    hmax_arr=  simulate_head(Qmx, wellid,well_proper, skin_proper)
                    hmax = hmax_arr[te_index]
                count+=1
                PARC = Qmx
                print(f'iteration {count} residual: {hmax-H_per}')
            parc_output.append([te, PARC])
            # parc_heads.append(hmax_arr)

    return parc_output


def default_model(base_pars, layers):
    for pars in base_pars.keys():
        par_val = np.ones(len(layers))*base_pars[pars]
        change_array(pars, par_val, layers)
        print('++++++++++++++++')
    return None


def PARC(HK,SS, VANI,HANI,scl,locsc):
    layers  =np.arange(1,10) # confined aquifer
    Ztop = locsc+(scl/2)
    Zbot = locsc-(scl/2)

    wellid = 4418
    row = 47
    col = 48
    Hper = 90
    tend = [86400]

    ## manipulation of array data
    HK_arr = np.ones(len(layers))*HK
    SS_arr = np.ones(len(layers))*SS
    VANI_arr = np.ones(len(layers))*VANI
    HANI_arr = np.ones(len(layers))*HANI
    

    change_array('HK', HK_arr, layers)
    change_array('SS', SS_arr, layers)
    change_array('VANI', VANI_arr, layers)
    change_array('HANI', HANI_arr, layers)

    ## well manipulation
    well_properties = [Ztop, Zbot, row, col]  #[Ztop,Zbotm,ROW,COL]
    skin_properties = [0.127,1.7,0.00011574]  #[Rw,Rskin,Kskin]
    
    #simulation
    owir2d = OWIR(tend, wellid, Hper, 0.01, skin_properties, well_properties)


    return owir2d[0][1]

