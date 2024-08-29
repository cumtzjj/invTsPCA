import numpy as np
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


#    for j in select_slow: 
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