import numpy as np
from sko.PSO import PSO
from sko.tools import set_run_mode
import respknt
from scipy import signal


def read_zac(zacfile,raypfile,stdfile,filetype='npy'):
    '''
    read data, include ZACs, rayps, stds
    '''
    if filetype == 'npy':
        zacs = np.load(zacfile)
        slow = np.load(raypfile)
        std = np.load(stdfile)

    return zacs,slow,std

def proc_paras(proc_file,filetype='txt'):
    '''
    processing parameters of calculating ZACs
    time lenth : sec
    time sampling : sec
    frequency band : Hz
    '''
    if filetype == 'txt':
        x = np.loadtxt(proc_file,usecols=1)

#    print(x)

#    freq_band = []
    time_len = x[0]
    time_samp = x[1]
    freq_band = x[2:]
#    print(freq_band)
    

    return time_len,time_samp,freq_band

#proc_paras(proc_file)
#read_zac(zacfile,raypfile,stdfile)

#print(slow)


def inv_paras(inv_file,filetype='txt'):
    '''
    inversion parameters
    lb : low boundary
    ub : up boundary
    size : swarm size / len(lb)
    max_iter : max iteration
    '''
    return



def modeling_ac(model,nl,slow,time_len,time_samp,freq_band,model_smooth=1):
    '''
    modeling ZACs using respknt
    '''
    maxlay                       = 50
    vp_in                        = np.zeros(maxlay, dtype=np.float32)
    vs_in                        = np.zeros(maxlay, dtype=np.float32)
    k_in                         = np.zeros(maxlay, dtype=np.float32)
    rho_in                       = np.zeros(maxlay, dtype=np.float32)
    thick_in                     = np.zeros(maxlay, dtype=np.float32)
# parameter set up
    dt                           = time_samp
    tduring                      = time_len

    vp_in[:nl]                   = model[:nl]  
    if len(model) >= 2*nl:
        k_in[:nl]                    = model[2*nl-1:3*nl-1] 
    else:
        k_in[:nl] = 1.73
       
#    print(k_in)
    vs_in[:nl]=vp_in[:nl]/k_in[:nl]

    if len(model) >= 3*nl:
        rho_in[:nl] = model[-nl:]
    else:
        rho_in=1.6612*vp_in-0.4721*vp_in**2+0.0671*vp_in**3-0.0043*vp_in**4+0.000106*vp_in**5

    thick_in[:nl-1]                = model[nl:2*nl-1] 

    if model_smooth:
        thick_in[thick_in<model_smooth]=0

    nyq_bd = freq_band*dt*2 
    b, a = signal.butter(2, nyq_bd, 'bandpass')
    time                         = np.arange(0, tduring+dt, dt)
    lent                         = len(time)

    ac_out = np.zeros((lent,slow.size))

#    print(ac_out.shape)


    for j in range(slow.size):
        if slow.size == 1:
            slowness = slow
        else:
            slowness                 = slow[j]

        wavez,waver                 = respknt.respknt_interface(slowness=slowness, dt_in=dt, tduring=tduring, nl=nl, rho_in=rho_in, thick_in=thick_in, vp_in=vp_in, vs_in=vs_in)
        Ttz                          = wavez[:lent]
        bp=signal.filtfilt(b, a, Ttz)
        ac=signal.correlate(bp,bp,mode='full')
        #print(ac.shape)
        ac_out[:,j]=-ac[lent-1:2*lent-1]/ac.max()

    return ac_out

def obj_func(x):
    '''
    objective function
    '''

#    mse = []
    
#    proc_paras(proc_file)
#    read_zac(zacfile,raypfile,stdfile)
    ac_out=modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=0.1)

    if mse_flag == 'mse':
        maxl = int(zac_lent/time_samp)
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
            tmp = np.corrcoef(zacs[:400,i],ac_out[:400,i])
            tmp2 = tmp2 + tmp[0,1]
        tmp3 = tmp2 / len(select_slow)

        mse = 1 - tmp3


    return mse




if __name__ == '__main__':
    
#    global proc_file = 'process_file.txt'
#    global zacfile = 'model_icezac_bins_noise.npy'
#    global raypfile = 'model_icerayp_bins_noise.npy'
#    global stdfile = 'model_icezac_std_noise.npy'
#    global zac_lent = 5
#    global std_flag = True

#    global nl=2
    proc_file = 'process_file_QSPA.txt'
    zacfile = 'zac_bins_QSPA.npy'
    raypfile = 'rayp_bins_QSPA.npy'
    stdfile = 'zac_std_QSPA.npy'
    zac_lent = 5
    std_flag = True
    nl=4
    mse_flag = 'corr'


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

    rho_lb=x[:,3]*0.95 #np.full(nl,0.8)
    rho_ub=x[:,3]*1.05 #np.full(nl,3.0)

    lb = np.concatenate((vp_lb,h_lb,k_lb,rho_lb))
    ub = np.concatenate((vp_ub,h_ub,k_ub,rho_ub))

    print('low_boundary=',lb)
    print('upper_boundary=',ub)
    
 
    mode = 'multiprocessing'
    set_run_mode(obj_func, mode)

    nloop=10
    kloop=9

    for jloop in range(nloop):
        iloop=kloop*nloop+jloop
        pso = PSO(func=obj_func, n_dim=len(lb), pop=len(lb)*20, max_iter=600, lb=lb, ub=ub, w=0.4, c1=2.0, c2=2.0)
        pso.run()
        print('best_x is ', pso.gbest_x, 'best_y is', pso.gbest_y)
        fname='out/qspa_corr_vel_best_x'+str(iloop)+'.csv'
        np.savetxt(fname,pso.gbest_x,delimiter=',')
        fname='out/qspa_corr_vel_gbest_y'+str(iloop)+'.csv'
        np.savetxt(fname,pso.gbest_y_hist,delimiter=',')





