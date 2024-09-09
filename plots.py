import matplotlib.pyplot as plt
import numpy as np
from utils import read_zac,modeling_ac,proc_paras

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def plot_wiggles(traces,times,slow,outfile,scale=0.005,xlim=[0,40],stack_flag=False):
    '''
    plot wiggles 
    '''
    fig = plt.figure(figsize=(3,4))

    grid = plt.GridSpec(11,1,hspace=0.1,wspace=0.9)

    if stack_flag:
        main = fig.add_subplot(grid[2:])
        stack = fig.add_subplot(grid[:2],xticklabels=[],yticklabels=[])
    else:
        main = fig.add_subplot(grid[0:])

    for j in range(traces.shape[1]):
        main.plot(times,traces[:,j]*scale+slow[j],color='black')
        main.fill_between(times,traces[:,j]*scale+slow[j],slow[j],where=traces[:,j]>0,color='black')

    if stack_flag:
        line_stack=np.mean(traces,axis=1)
        rayp_stack=np.mean(slow)
        stack.plot(times, line_stack*scale+ rayp_stack ,color='red')
        stack.fill_between(times, line_stack*scale+rayp_stack,rayp_stack,where=line_stack>0,color='black')
        stack.set_xlim(xlim[0],xlim[1])
        stack.set_ylim(rayp_stack-0.005,rayp_stack+0.005)


    main.set_xlabel('Time (s)')
    main.set_ylabel('Slowness (s/km)')

    main.set_xlim(xlim[0],xlim[1])
    main.set_ylim(min(slow)-0.005,max(slow)+0.005)
    plt.tight_layout()
    plt.savefig(outfile+'.pdf',bbox_inches='tight')
    plt.savefig(outfile+'.png',dpi=600,bbox_inches='tight')

    plt.close(fig)

    return

def plot_wiggles2(traces,times,slow,outfile,scale=0.005,xlim=[0,40],stack_flag=False,figsize=(4,5)):
    '''
    plot wiggles 
    '''
    fig = plt.figure()

    if stack_flag:
        line_stack=np.mean(traces,axis=1)
    grid = plt.GridSpec(10,1,hspace=0.1,wspace=0.9)

    if stack_flag:
        main = fig.add_subplot(grid[1:])
        stack = fig.add_subplot(grid[0],xticklabels=[],yticklabels=[])
    else:
        main = fig.add_subplot(grid[0:])

    for j in range(traces.shape[1]):
        main.plot(times,traces[:,j]*scale+slow[j],color='black')
        main.fill_between(times,traces[:,j]*scale+slow[j],slow[j],where=traces[:,j]>0,color='black')

    if stack_flag:
        line_stack=np.mean(traces,axis=1)
        rayp_stack=np.mean(slow)
        stack.plot(times, line_stack*scale+ rayp_stack ,color='red')
        stack.fill_between(times, line_stack*scale+rayp_stack,rayp_stack,where=line_stack>0,color='black')
        stack.set_xlim(xlim[0],xlim[1])
        stack.set_ylim(rayp_stack-0.01,rayp_stack+0.01)


    main.set_xlabel('Time (s)')
    main.set_ylabel('Slowness (s/km)')

    main.set_xlim(xlim[0],xlim[1])
    main.set_ylim(min(slow)-0.01,max(slow)+0.01)

    plt.tight_layout()
    plt.savefig(outfile+'.pdf',bbox_inches='tight')
    plt.savefig(outfile+'.png',dpi=600,bbox_inches='tight')

    plt.close(fig)

    return


def plot_model(fileref,model,nl,model_ref,nloop=100,ylim=(60,0),figsize=(2,2.5),model_smooth=0.1):
    '''
    plot models; vp, vs, rho, vp/vs... 
    '''
#    w_inc=figsize_cm[0]*0.39
#    h_inc=figsize_cm[1]*0.39
#    nloop=100
    b=np.zeros(nloop)
    fig = plt.figure(figsize=figsize)

    # plot reference model
#    if len(model_ref)>0:
#        mod=np.loadtxt(model_ref)
#        if model == 'vp':
#            v=mod[:,1]
#            plt.xlim(min(v)-2,max(v)+2)
#        elif model == 'vpvs':
#            v=mod[:,2]
#        elif model == 'density':
#            v=mod[:,3]

#        h=mod[:,0]

#        v_ref, dep_int = make_model(v,h)
#    l1,=plt.plot(v_int,dep_int,color='Black',label='Model',linewidth=1.0)

    for iloop in range(nloop):
        fname=fileref+'gbest_y'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        b[iloop]=np.min(x)

    inb=np.argsort(b)[::-1]
    print(inb)
    bb=np.sort(b)[::-1]
    print(bb)

    best10penc = False
    plotbest10 = False

    if best10penc:
        bestloop = int(nloop*0.1)
    else:
        bestloop = nloop

    for iloop in inb[-bestloop:]:
#    for iloop in inb:
        fname=fileref+'best_x'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        if model == 'vp':
            v=x[:nl]
            plt.xlabel('Vp (km/s)')
#            plt.xlim(4,10)
        elif model == 'vpvs':
            v=x[2*nl-1:3*nl-1]
            plt.xlabel('Vp/Vs')
            plt.xlim(1,3)
        elif model == 'density':
            v=x[-nl:]
            plt.xlabel(r'Density ($g/cm^3$)')
            plt.xlim(0.5,4)

        h_lay = x[nl:(2*nl-1)]
        h_lay[h_lay<model_smooth]=0
        h = np.concatenate((h_lay, [0]))
        v_int, dep_int = make_model(v,h)

#        lnum=[]

        if iloop == inb[-1]:
            if len(model_ref)>0:
                lnum=[]
                for i in range(len(model_ref)):
                    mod=np.loadtxt(model_ref[i])
                    if model == 'vp':
                        v=mod[:,1]
                        plt.xlim(min(v)-2,max(v)+2)
                    elif model == 'vpvs':
                        v=mod[:,2]
                    elif model == 'density':
                        v=mod[:,3]

                    h=mod[:,0]

                    v_ref, dep_int = make_model(v,h)
                    lnum.append('l'+str(i))
                    lnum[i],=plt.plot(v_ref,dep_int,color='black',label=model_label[i],linewidth=1.0)
            lnum.append('l'+str(i+1))
            lnum[i+1],=plt.plot(v_int,dep_int,color='green',label='Best',linewidth=1.0)
        else:
            l100,=plt.plot(v_int,dep_int,color='gray',label='Invert',linewidth=0.5)

    if plotbest10:
        xx = np.zeros_like(x)
        for iloop in inb[-bestloop:]:
            fname=fileref+'best_x'+str(iloop)+'.csv'
            x=np.loadtxt(fname,delimiter=',')
            xx += x

        x = xx/bestloop

        if model == 'vp':
            v=x[:nl]
            plt.xlabel('Vp (km/s)')
#            plt.xlim(4,10)
        elif model == 'vpvs':
            v=x[2*nl-1:3*nl-1]
            plt.xlabel('Vp/Vs')
            plt.xlim(1,3)
        elif model == 'density':
            v=x[-nl:]
            plt.xlabel(r'Density ($g/cm^3$)')
            plt.xlim(0.5,4)

        h_lay = x[nl:(2*nl-1)]
        h_lay[h_lay<model_smooth]=0
        h = np.concatenate((h_lay, [0]))
        v_int, dep_int = make_model(v,h)

        l101,=plt.plot(v_int,dep_int,color='blue',label='Mean',linewidth=1.0)

    if len(model_ref)>0:
        plt.legend(handles=lnum)
#        plt.legend(handles=[l3])
    else:
        plt.legend(handles=[l3])
    plt.ylabel('Depth (km)')
#    plt.grid()
#    plt.xlabel('Vp (km/s)')
#    plt.xlim(min(v)-1,max(v)+1)
    plt.ylim(ylim)
#    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(fileref+model+'_invert.pdf',bbox_inches='tight')
    plt.savefig(fileref+model+'_invert.png',dpi=600,bbox_inches='tight')

    return

def make_model(y,h,dh=0.01):
    '''
    make model
    '''
    v_step = np.concatenate([(v, v) for v in y])
    dep = np.cumsum(h)
    dep_int = np.arange(0, 100, dh)
    dep = np.concatenate([(d, d) for d in dep])
    dep_step = np.concatenate([[0], dep[:-1]])
    dep_step[-1] = np.max([150, dep_step[-1] * 2.5])  # half space
    v_int = np.interp(dep_int, dep_step, v_step)

    return v_int, dep_int


def plot_waveforms(synfile,stdfile,raypfile,proc_file,fileref,figsize=(4,5),nloop=100,xlim=(0,30),scale=0.3,select=[0,1,2,3,4,5,6,7],model_smooth=0.1):
    '''
    plot waveforms compareison between synthetic and inverted
    '''
    b=np.zeros(nloop)
    less_lay=0.1
    fig = plt.figure(figsize=figsize)

    # plot syn and std
    zacs, slow, std = read_zac(zacfile,raypfile,stdfile)

    print(zacs.size,slow.size,std.size)

    time_len, time_samp, freq_band = proc_paras(proc_file)

    time = np.arange(0, time_len+time_samp, time_samp)

    # plot inverted waveforms
    for iloop in range(nloop):
        fname=fileref+'gbest_y'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        b[iloop]=np.min(x)

    inb=np.argsort(b)[::-1]
    print(inb)
    bb=np.sort(b)[::-1]
    print(bb)

    best10penc = True
    plotbest10 = False

    if best10penc:
        bestloop = int(nloop*0.1)
    else:
        bestloop = nloop


    for iloop in inb[-bestloop:]:
        fname=fileref+'best_x'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        ac_out = modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=model_smooth)

        for j in select: #range(slow.size):

            if slow.size == 1:
                obs = zacs
                mac = ac_out
                ostd = std
            else:
                obs = zacs[:,j]
                mac = ac_out[:,j]
                ostd = std[:,j]

            if iloop == inb[-1]:
                l1, = plt.plot(time[:len(obs)],obs+scale*j,color='black',label='Observation',linewidth=2)
                l2, = plt.plot(time,mac+scale*j,color='dodgerblue',label='Best estimation',linewidth=1)
                l4 = plt.fill_between(time[:len(obs)],obs+ostd+scale*j,obs-ostd+scale*j,color='gray')
                #l5, = plt.plot(time[:len(obs)],obs-ostd+scale*j,linestyle='--',color='black',linewidth=0.3)
       #     else:
                #l3, = plt.plot(time,mac+scale*j,color='red',linewidth=0.5)

    if plotbest10:
        xx = np.zeros_like(x)
        for iloop in inb[-bestloop:]:
            fname=fileref+'best_x'+str(iloop)+'.csv'
            x=np.loadtxt(fname,delimiter=',')
            xx += x

        x = xx/bestloop

        ac_out = modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=model_smooth)

        for j in range(slow.size):

            if slow.size == 1:
                obs = zacs
                mac = ac_out
                ostd = std
            else:
                obs = zacs[:,j]
                mac = ac_out[:,j]
                ostd = std[:,j]

            l6, = plt.plot(time,mac+scale*j,color='blue',label='Mean',linewidth=0.5)




    plt.legend(handles=[l1,l2])
    plt.grid()
    plt.xlabel('Time (s)')
    if slow.size > 1:
        plt.yticks(np.arange(0,slow.size*scale,scale),labels=np.round(slow,3))
    else:
        plt.yticks([0,1],labels=[np.round(slow,3),1])
    plt.ylabel('Slowness (s/km)')
    plt.ylim((select[0]-1)*scale,(select[-1]+1)*scale)
    plt.xlim(xlim)
    plt.tight_layout()
    plt.savefig(fileref+'com_obs_est.pdf',bbox_inches='tight')  
    plt.savefig(fileref+'com_obs_est.png',dpi=600,bbox_inches='tight')


    return

from scipy.stats import norm,zscore

def uncertainty_analysis(synfile,stdfile,raypfile,proc_file,fileref,model_ref,title,figsize=(2,2.5),nloop=100,ylim=(50,0),xlim=(0,10),w=0.3,select=[0,1,2,3,4,5,6,7],model_smooth=0.1,model='vp',vmax=5):
    '''
    cal_misfit between synthetic and inverted
    '''

    fig = plt.figure(figsize=figsize)

# plot reference model
    if len(model_ref)>0:
        mod=np.loadtxt(model_ref)
        if model == 'vp':
            v=mod[:,1]
#        plt.xlim(min(v)-2,max(v)+2)
        elif model == 'vpvs':
            v=mod[:,2]
        elif model == 'density':
            v=mod[:,3]

        h=mod[:,0]
        
        h[h<model_smooth]=0

        v_ref, dep_int = make_model(v,h)
        l1,=plt.plot(v_ref,dep_int,color='black',label='Model',linewidth=1.0)

# plot best model

    zacs, slow, std = read_zac(zacfile,raypfile,stdfile)

    time_len, time_samp, freq_band = proc_paras(proc_file)
    maxl = int(zac_lent/time_samp)
    print(maxl)
    bias_len=int(bias/time_samp)
    print(bias_len)
    if len(select)>1:
        data_mse = np.mean(std[bias_len:maxl,:])
    else:
        data_mse = np.mean(std[bias_len:maxl])

    print(data_mse)

    data_mse = data_mse * w

#    time_len, time_samp, freq_band = proc_paras(proc_file)

    time = np.arange(0, time_len+time_samp, time_samp)

#     bmse=np.zeros(nloop)
#     for iloop in range(nloop):
#         fname=fileref+'best_x'+str(iloop)+'.csv'
#         x=np.loadtxt(fname,delimiter=',')

#         ac_out = modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=model_smooth)

# #        maxl = int(zac_lent/time_samp)
# #        print('maxl=',maxl)
#         if len(select)>1:
#             tmp = zacs[bias_len:maxl,select] - ac_out[bias_len:maxl,select]
#             bmse[iloop] = np.mean(np.sqrt(np.mean(tmp**2,axis=0)))

#         else:
#             tmp = np.squeeze(zacs[bias_len:maxl]) - np.squeeze(ac_out[bias_len:maxl])
#             bmse[iloop] = np.sqrt(np.mean(tmp**2))

#     minindex = np.argmin(bmse)
    
    b=np.zeros(nloop)

    
    for iloop in range(nloop):
        fname=fileref+'gbest_y'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        b[iloop]=np.min(x)
    
    #print(bmse)
    minindex = np.argmin(b)
    print(minindex)

    fname=fileref+'best_x'+str(minindex)+'.csv'
    x=np.loadtxt(fname,delimiter=',')

    if model == 'vp':
        v=x[:nl]
        plt.xlabel('Vp (km/s)')
    elif model == 'vpvs':
        v=x[2*nl-1:3*nl-1]
        plt.xlabel('Vp/Vs')
    elif model == 'density':
        v=x[-nl:]
        plt.xlabel(r'Density ($g/cm^3$)')

    h_lay = x[nl:(2*nl-1)]
    h_lay[h_lay<model_smooth]=0
    h = np.concatenate((h_lay, [0]))
    v_int, dep_int = make_model(v,h)

    l2,=plt.plot(v_int,dep_int,color='dodgerblue',label='Best',linewidth=1.0)

###
### plot ppd and Mean
###
#    zacs, slow, std = read_zac(zacfile,raypfile,stdfile)
#
#    time_len, time_samp, freq_band = proc_paras(proc_file)
#    maxl = int(zac_lent/time_samp)
#    print(maxl)
#    if len(select)>1:
#        data_mse = np.mean(std[:maxl,:])
#    else:
#        data_mse = np.mean(std[:maxl])
#
#    print(data_mse)

#    data_mse = data_mse * w

#    time_len, time_samp, freq_band = proc_paras(proc_file)

#    time = np.arange(0, time_len+time_samp, time_samp)

    vp_good=[]
    

    for iloop in range(nloop):
        fname=fileref+'best_x'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        ac_out = modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=model_smooth)

#        maxl = int(zac_lent/time_samp)
#        print('maxl=',maxl)

        if len(select)>1:
            tmp = zacs[bias_len:maxl,select] - ac_out[bias_len:maxl,select]
            mse = np.mean(np.sqrt(np.mean(tmp**2,axis=0)))
            misfit = np.mean(np.sqrt(np.sum(tmp**2,axis=0)))
        else:
            tmp = np.squeeze(zacs[bias_len:maxl]) - np.squeeze(ac_out[bias_len:maxl])
            mse = np.sqrt(np.mean(tmp**2))

        #print(tmp)
       
         
        #mse = np.mean(np.sqrt(np.mean(tmp**2,axis=0)))

        #print(mse,misfit)
        if mse < data_mse:
            if model == 'vp':
                v=x[:nl]
            elif model == 'vpvs':
                v=x[2*nl-1:3*nl-1]
            elif model == 'density':
                v=x[-nl:]

            h_lay = x[nl:(2*nl-1)]
            h_lay[h_lay<model_smooth]=0
            h = np.concatenate((h_lay, [0]))
            v_int, dep_int = make_model(v,h)
            vp_good.append(v_int)

#    print(vp_good)
#    print(len(vp_good))

    vp_good_array = np.array(vp_good)


#    print(vp_good_array[:,400])
#    print(vp_good_array.shape)

#    print(dep_int)

    vp_mean = np.mean(vp_good_array, axis=0)

#    print(vp_mean)
#    print(vp_mean.shape)

#    vp_mean = vp_mean[:int(ylim[0]/0.01+1)]

    vp_std = np.std(vp_good_array, axis=0)

 #   vp_std = vp_std[:int(ylim[0]/0.01+1)]
#    print(vp_std[300:400])

 #   z_scores = np.abs(zscore(vp_std))

#    threshold = 3

#    vp_std[z_scores < threshold] = np.median(vp_std)

#    print(vp_std[300:400])

    vp_mean_ds = vp_mean[0:int(ylim[0]/0.01+1):10]
    vp_std_ds = vp_std[0:int(ylim[0]/0.01+1):10]
    dep_int_ds = dep_int[0:int(ylim[0]/0.01+1):10]


#    print(vp_std_ds)
#    z_scores = np.abs(zscore(vp_std_ds))
#
#    print(z_scores)

#    threshold = 1

#    vp_std_ds[z_scores > threshold] = np.median(vp_std_ds)

#    print(vp_std_ds)

    plt.plot(vp_mean,dep_int,label='Mean',color='darkorange')


    if model == 'vp':
        x_range, y_range = np.meshgrid(np.linspace(3,10,701),dep_int_ds)
    elif model == 'vpvs':
        x_range, y_range = np.meshgrid(np.linspace(1,3,201),dep_int_ds)
    elif model == 'density':
        x_range, y_range = np.meshgrid(np.linspace(0.5,4,351),dep_int_ds)


    ppd = np.array([norm.pdf(x,mean,std) for x,mean,std in zip(x_range,vp_mean_ds,vp_std_ds)])
    
    plt.pcolor(x_range,y_range,ppd,cmap='YlGn')
    
#    plt.legend()

    plt.ylim(ylim)
    plt.xlim(xlim)

    plt.ylabel('Depth (km)')
    if len(title)>0:
        plt.title(title)

    plt.colorbar()
    plt.tight_layout()
    plt.savefig(fileref+model+str(w*10)+'_uncertainty.pdf',bbox_inches='tight')
    plt.savefig(fileref+model+str(w*10)+'_uncertainty.png',dpi=600,bbox_inches='tight')


    return 



def plot_gbesty(fileref,nloop=100,figsize=(3,2)):
    '''
    plot gbesty
    '''
    fig = plt.figure(figsize=figsize)
    b=np.zeros(nloop)

    for iloop in range(nloop):
        fname=fileref+'gbest_y'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        b[iloop]=np.min(x)

    inb=np.argsort(b)[::-1]
    print(inb)
    bb=np.sort(b)[::-1]
    print(bb)
#plt.plot(inb)
    best10penc = False
    plotbest10 = False

    if best10penc:
        bestloop = int(nloop*0.1)
    else:
        bestloop = nloop

    for iloop in inb[-bestloop:]:
        fname=fileref+'gbest_y'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        plt.plot(x,linewidth=0.3)

#plt.legend()
    plt.ylabel('Misfit')
    plt.xlabel('iteration')
    plt.tight_layout()
    plt.savefig(fileref+'besty_loop.pdf',bbox_inches='tight')
    plt.savefig(fileref+'besty_loop.png',dpi=600,bbox_inches='tight')

    return

def plot_evs_wforms(outfile,synfile,evsfile,raypfile,scale,xlim,figsize):
    """
    plot events number and waveforms
    """

    events_num = np.load(evsfile)
    traces = np.load(synfile)
    slow = np.load(raypfile)
    time_len, time_samp, freq_band = proc_paras(proc_file)

    times = np.arange(0, time_len+time_samp, time_samp)
    maxl = int(zac_lent/time_samp)+1
    times = times[:maxl]

    fig = plt.figure(figsize=figsize)

    grid = plt.GridSpec(1,41,hspace=0.9,wspace=0.05)

    main = fig.add_subplot(grid[10:],yticks=slow,yticklabels=[])
    evs = fig.add_subplot(grid[:8],yticks=slow,yticklabels=np.round(slow,3),xticks=np.arange(0,max(events_num)+5,10))

    for j in range(traces.shape[1]):
        main.plot(times,traces[:,j]*scale+slow[j],color='black')
        main.fill_between(times,traces[:,j]*scale+slow[j],slow[j],where=traces[:,j]>0,color='black')

    evs_num = evs.barh(slow,events_num,color='cyan',height=0.002)
#    ax.bar_label(evs_num,label_type='center')

    main.set_xlabel('Time (s)')
    evs.set_ylabel('Slowness (s/km)')
    evs.set_xlabel('events count')
#    main.yticks(np.arange(0,slow.size*scale,scale))
#    evs.yticks(np.arange(0,slow.size),labels=np.round(slow,3))
    evs.set_xlim(max(events_num)+5,0)
    evs.set_ylim(min(slow)-0.005,max(slow)+0.005)
    main.set_xlim(xlim[0],xlim[1])
    main.set_ylim(min(slow)-0.005,max(slow)+0.005)
    main.grid()
    evs.grid()
    plt.tight_layout()
    plt.savefig(outfile+'evs_wforms.pdf',bbox_inches='tight')
    plt.savefig(outfile+'evs_wforms.png',dpi=600,bbox_inches='tight')

    plt.close(fig)

    return

def get_uncertainty(fileref,zacfile,raypfile,stdfile,proc_file,model):
    zacs, slow, std = read_zac(zacfile,raypfile,stdfile)
#
    time_len, time_samp, freq_band = proc_paras(proc_file)
    maxl = int(zac_lent/time_samp)
    print(maxl)
    if len(select)>1:
        data_mse = np.mean(std[:maxl,:])
    else:
        data_mse = np.mean(std[:maxl])
#
    print(data_mse)

    vlist = []
    for iloop in range(nloop):
        fname=fileref+'best_x'+str(iloop)+'.csv'
        x=np.loadtxt(fname,delimiter=',')

        ac_out = modeling_ac(x,nl,slow,time_len,time_samp,freq_band,model_smooth=model_smooth)

        if len(select)>1:
            tmp = zacs[:maxl,select] - ac_out[:maxl,select]
            mse = np.mean(np.sqrt(np.mean(tmp**2,axis=0)))
        else:
            tmp = np.squeeze(zacs[:maxl]) - np.squeeze(ac_out[:maxl])
            mse = np.sqrt(np.mean(tmp**2))


        print(mse)
        if mse < data_mse:
            if model == 'vp':
                v=x[:nl]
            elif model == 'vpvs':
                v=x[2*nl-1:3*nl-1]
            elif model == 'density':
                v=x[-nl:]
            elif model == 'h':
                v=x[nl:2*nl-1]
            vlist.append(v)
    #print(vlist)

    varray = np.array(vlist)
    
    varray[varray<model_smooth]=0

    #print(varray)

    vmean = np.mean(varray,axis=0)
    vstd = np.std(varray,axis=0)
    vmedian = np.median(varray,axis=0)

    print(vmean,vstd,vmedian)

    threshold = 2

    vfilter = []

    for i in range(varray.shape[0]):
        if np.abs(varray[i,0]-vmean[0])<=threshold*vstd[0]:
            vfilter.append(varray[i,:])

    print(vfilter)

    vfilter_array = np.array(vfilter)

    print(vfilter_array)

    out_mean = np.mean(vfilter_array,axis=0)
    out_std = np.std(vfilter_array,axis=0)

    print(out_mean,out_std)

    ice_h = []
    if model == 'h':
        for i in range(vfilter_array.shape[0]):
            ice_h.append(np.sum(vfilter_array[i,:2]))
    ice_mean = np.mean(np.array(ice_h))
    ice_std = np.std(np.array(ice_h))

    print(ice_mean,ice_std)

    return

from scipy import io
if __name__ == '__main__':
#    fileref='out/bosa_bias_pop_'
    fileref='out/qspa_corr_'
#    fileref='out/model2_5_'
#    fileref='out/BOSA_corr_'
    #model='model'
    nl=4
#    model_ref='model2.txt'
#    model_label=['Model']
#    model_ref='bosa_crust.txt' #['bosa.txt','bosa_crust.txt']
#    model_label=['James et al., 2003','Crust 1.0']
    model_ref=''
    nloop=100
    xlim=(0,5)
    scale=0.5
    select=[0,1,2,3,4,5,6,7]
#    select=[0,1,2,3]
#    select=[0]
    figsize=(4,5)
    zac_lent=5
    std_flag=True
    model_smooth=0.1
    w = 1
    ylim=(5,0)
#    title=r'$M_3$'
    title='QSPA'
    bias=0
    

    proc_file = 'process_file_QSPA.txt'
#    zacfile = 'model2zac_bins_noise.npy'
#    raypfile = 'model2rayp_bins_noise.npy'
#    stdfile = 'model2zac_std_noise.npy'

#    zacfile = 'model0zac_stak_noise.npy'
#    raypfile = 'model0_stack_rayp.npy'
#    stdfile = 'model0zac_all_std_noise.npy'


    zacfile = 'zac_bins_QSPA.npy' 
#    zacfile = 'bosa_zac_bins.npy'
#    raypfile = 'bosa_rayp_bins.npy'
    raypfile = 'rayp_bins_QSPA.npy'    
#    raypfile = 'model0_stack_rayp.npy'
#    stdfile = 'bosa_std_bins.npy'
    stdfile = 'zac_std_QSPA.npy'
#    stdfile = 'model0zac_all_std_noise.npy'
#    proc_file = 'process_file_model0.txt'
#    evsfile = 'events_num_BOSA.npy'
#    evsfile = 'events_num.mat'
#    evs = io.loadmat(evsfile)
#    evss = evs['events_num']
#    evs_num = evss[0]
#    print(evs_num)
#    np.save('events_num_BOSA.npy',arr=evs_num)
#    zacs_tmp=np.zeros((2401,8))
#    std_tmp=np.zeros((2401,8))

#    zacs_tmp[:801,:]=np.load(zacfile)
#    std_tmp[:801,:]=np.load(stdfile)

#    np.save('bosa_zac_bins.npy',arr=zacs_tmp)
#    np.save('bosa_std_bins.npy',arr=std_tmp)

#    zacfile = 'model0zac_stak_noise.npy'
#    raypfile = 'model0_stack_rayp.npy'
#    stdfile = 'model0zac_all_std_noise.npy'


    #plot_model(fileref,model='vp',nl=nl,model_ref=model_ref,nloop=nloop,ylim=(50,0),model_smooth=model_smooth)
    #plot_model(fileref,model='vpvs',nl=nl,model_ref=model_ref,nloop=nloop,ylim=(50,0))
#    plot_model(fileref,model='density',nl=nl,model_ref=model_ref,nloop=nloop,ylim=(5,0))
    #plot_evs_wforms(fileref,zacfile,evsfile,raypfile,scale,xlim,figsize)
    get_uncertainty(fileref,zacfile,raypfile,stdfile,proc_file,model='h')
    plot_waveforms(zacfile,stdfile,raypfile,proc_file,fileref,figsize=figsize,nloop=nloop,xlim=xlim,scale=scale,select=select,model_smooth=model_smooth)
    plot_gbesty(fileref,nloop=nloop)
    uncertainty_analysis(zacfile,stdfile,raypfile,proc_file,fileref,title=title,model_ref=model_ref,nloop=nloop,ylim=ylim,xlim=(3,5),w=w,select=select,model_smooth=model_smooth,model='vp')
    uncertainty_analysis(zacfile,stdfile,raypfile,proc_file,fileref,title=title,model_ref=model_ref,nloop=nloop,ylim=ylim,xlim=(1,3.0),w=w,select=select,model_smooth=model_smooth,model='vpvs')
    uncertainty_analysis(zacfile,stdfile,raypfile,proc_file,fileref,title=title,model_ref=model_ref,nloop=nloop,ylim=ylim,xlim=(0.5,3.0),w=w,select=select,model_smooth=model_smooth,model='density')
#    plot_misfit(zacfile,stdfile,raypfile,proc_file,fileref,figsize=figsize,nloop=nloop,xlim=(0,20),scale=5,select=select,model_smooth=model_smooth)
