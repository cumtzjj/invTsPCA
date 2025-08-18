import numpy as np
import respknt
from utils import modeling_ac


def obj_func(x):
    '''
    objective function
    '''

#    mse = []
    
#    proc_paras(proc_file)
#    read_zac(zacfile,raypfile,stdfile)
    ac_out=modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=0.1)
    maxl = int(zac_lent/time_samp)
    if mse_flag == 'mse':
    #    maxl = int(zac_lent/time_samp)
        bias = int(bias_lent/time_samp)
#    print('maxl=',maxl)
        tmp = zacs[:maxl,select_slow] - ac_out[:maxl,select_slow]
        if std_flag:
            tmp[bias:,select_slow] = tmp[bias:,select_slow] / std[bias:maxl,select_slow]
        tmp2 = np.square(tmp)
        tmp3 = np.sqrt(np.mean(tmp2,axis=0))

        mse = np.mean(tmp3)

    tmp2 = 0
    if mse_flag == 'corr':
        for i in select_slow:
            tmp = np.corrcoef(zacs[:maxl,i],ac_out[:maxl,i])
            tmp2 = tmp2 + tmp[0,1]
        tmp3 = tmp2 / len(select_slow)

        mse = 1 - tmp3


    return mse


from sko.tools import set_run_mode
from sko.PSO import PSO

if __name__ == '__main__':
    
    proc_file = 'process_file_QSPA.txt'  # processing parameters, including time_length (time length of teleseismic records), time_sampling, low_pass, and high_pass 
    zacfile = 'zac_bins_QSPA.npy'  # input binned stacked waveforms of ZAC
    raypfile = 'rayp_bins_QSPA.npy' # raypath parameters of input ZAC bins
    stdfile = 'zac_std_QSPA.npy' # STD of input ZAC bins
    zac_lent = 5 # time length of ZAC
    std_flag = True  # calculate MSE via STD, Ture or False
    nl=4  # maxium layers of inverted model 
    mse_flag = 'corr'  # methd to calculate misfit, "mse" or "corr" 


    zacs=np.zeros((3001,8))
    slow=np.zeros((1,8))
    std=np.zeros((3001,8))


    zacs = np.load(zacfile)
    slow = np.load(raypfile)
    std = np.load(stdfile)

    print(slow)

    x = np.loadtxt(proc_file,usecols=1)
    time_len = x[0]
    time_samp = x[1]
    freq_band = x[2:]
    print(freq_band)

    bias_lent = 0
    std[std<0.04]=0.04

    select_slow=[0,1,2,3,4,5,6,7]

    x=np.loadtxt('model_ice2.txt')

    vp_ref = x[:,1]

    vp_lb=vp_ref-0.1  #[4,4,5,6,7]
    vp_ub=vp_ref+0.1 #[6,6,7,7,8]

#    h_lb=np.full(nl-1,0.2) #[3,3,5,5]
    h_lb=np.zeros(nl-1)
    h_ub=np.full(nl-1,4) #[10,10,20,10]

    k_lb=x[:,2]*0.95 #np.full(nl,1.5) #[1.6,1.6,1.6,1.6,1.6]
    k_ub=x[:,2]*1.05 #np.full(nl,3.0) #[1.9,1.9,1.9,1.9,1.9]

    rho_lb=x[:,3]*0.95   #*1000 #np.full(nl,0.8)
    rho_ub=x[:,3]*1.05  #*1000 #np.full(nl,3.0)

    lb = np.concatenate((vp_lb,h_lb,k_lb,rho_lb))
    ub = np.concatenate((vp_ub,h_ub,k_ub,rho_ub))

    print('low_boundary=',lb)
    print('upper_boundary=',ub)
    
    #bounds=[(i,j) for i,j in zip(lb,ub)]
    #print(bounds)
    
 
    mode = 'multiprocessing'
    set_run_mode(obj_func, mode)

    nloop=100
   # num_particles = 20*len(lb)
   # max_iter = 1000

    for iloop in range(nloop):
        #best_position, best_value_history = pso(obj_func, bounds, num_particles, max_iter, min_best_increase=0.0001)
        
        pso = PSO(func=obj_func,n_dim=len(lb),pop=50*len(lb),max_iter=800,lb=lb,ub=ub,w=0.8,c1=2,c2=2)
        pso.run()
        print('best_position is ', pso.gbest_x)
        fname='out/qspa_corr'+'_best_x'+str(iloop)+'.csv'
        np.savetxt(fname,pso.gbest_x,delimiter=',')
        fname='out/qspa_corr'+'_gbest_y'+str(iloop)+'.csv'
        np.savetxt(fname,pso.gbest_y_hist,delimiter=',')





