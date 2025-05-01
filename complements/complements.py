## COMPLEMENTARY MODULE FOR 'figures.ipynb' PROGRAM
## ------------------------------------------------
##
## The sole pourpose of this program is to serve as a library containing functions that could make
## look the main program more noisy, like plots, small simulation of random walk, etc.
## You can always check this program for more details
##
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import pandas as pd
import scipy.special as sp

## SIMULATION OF THE SIMPLE RANDOM WALK
## ------------------------------------
##
def individual_simple():
    steps = 20;
    random_walker = np.zeros((steps+1));
    for i in range(steps):
        random_walker[i+1] = random_walker[i] + 2*round(np.random.random())-1;

    plt.close();
    plt.rcParams["font.family"] = "monospace";
    figure, axs = plt.subplots(1,figsize=(6,3),facecolor='#f5f5f5',linewidth=1.5,edgecolor='k')
    plt.rcParams['xtick.labelsize'] = 9; plt.rcParams['ytick.labelsize'] = 9;
    axs.plot(random_walker,'-',c='#333E74',lw=2,alpha=0.5);
    axs.plot(random_walker,'o',c='#333E74',ms=10,alpha=0.9);
    axs.set_xticks(range(steps+1)); axs.set_yticks(range(int(min(random_walker)-2),int(max(random_walker)+2)));
    axs.set_xticks(np.arange(steps+1)-0.5,minor=True);
    axs.set_yticks(np.arange(int(min(random_walker)-2),int(max(random_walker)+2))-0.5,minor=True);
    # axs.set_yticks(minor='true');
    axs.set_xlim([-0.5,20.5]); axs.set_ylim([min(random_walker)-1.5,max(random_walker)+1.5]);
    axs.grid(linestyle='--', which='minor', alpha=0.3);
    axs.set_xlabel(r'$\mathbf{n}$',fontsize=12);
    axs.set_ylabel(r'$\mathbf{X}$',fontsize=12);
    plt.suptitle('Individual Simple Random Walk',x=0.53,y=0.95,fontsize=15,fontweight='bold');

    figure.tight_layout();
    plt.savefig('images/individual_simple.png', dpi=300, bbox_inches='tight', pad_inches=0.1);

    return plt.show()


## SIMULATION OF THE RANDOM WALK WITH RESETTING
## --------------------------------------------
##
def individual_resetting(probability):
    if not(probability>=0 and probability<=1):
        probability = 0.1;
    steps = 20;
    random_walker = np.zeros((steps+1));
    resetting_times = [];
    for i in range(steps):
        if (probability>np.random.random()):
            random_walker[i+1] = 0;
            resetting_times.append(i+1);
        else:
            random_walker[i+1] = random_walker[i] + 2*round(np.random.random())-1;

    plt.close();
    plt.rcParams["font.family"] = "monospace";
    figure, axs = plt.subplots(1,figsize=(6,3),facecolor='#f5f5f5',linewidth=1.5,edgecolor='k')
    plt.rcParams['xtick.labelsize'] = 9; plt.rcParams['ytick.labelsize'] = 9;
    axs.fill_between([-0.5,steps+1],[0.5,0.5],color='k',alpha=0.05,edgecolor='none');
    axs.fill_between([-0.5,steps+1],[-0.5,-0.5],color='k',alpha=0.05,edgecolor='none');
    axs.plot(random_walker,'-',c='#333E74',lw=2,alpha=0.5);
    for i in range(len(resetting_times)): axs.plot([resetting_times[i]-1,resetting_times[i]],\
                                                [random_walker[resetting_times[i]-1],random_walker[resetting_times[i]]],\
                                                    '-',c='r',lw=2,alpha=0.3);
    axs.plot(random_walker,'o',c='#333E74',ms=10,alpha=0.9);
    axs.set_xticks(range(steps+1)); axs.set_yticks(range(int(min(random_walker)-2),int(max(random_walker)+2)));
    axs.set_xticks(np.arange(steps+1)-0.5,minor=True);
    axs.set_yticks(np.arange(int(min(random_walker)-2),int(max(random_walker)+2))-0.5,minor=True);
    # axs.set_yticks(minor='true');
    axs.set_xlim([-0.5,20.5]); axs.set_ylim([min(random_walker)-1.5,max(random_walker)+1.5]);
    axs.grid(linestyle='--', which='minor', alpha=0.3);
    axs.set_xlabel(r'$\mathbf{n}$',fontsize=12);
    axs.set_ylabel(r'$\mathbf{X}$',fontsize=12);
    plt.suptitle('Individual Random Walk with Resetting',x=0.53,y=0.95,fontsize=15,fontweight='bold');

    figure.tight_layout();
    plt.savefig('images/individual_resetting.png', dpi=300, bbox_inches='tight', pad_inches=0.1);

    return plt.show()


## SIMULATION OF THE RANDOM WALK WITH MEMORY
## -----------------------------------------
##
def individual_memory(probability):
    if not(probability>=0 and probability<=1):
        probability = 0.1;
    presteps = 40;
    memory = np.zeros((2*presteps+1));
    memory[presteps] = 1;
    prewalk = 0;
    for i in range(presteps):
        if (probability>np.random.random()):
            z = int(np.random.random()*sum(memory))+1;
            ii = 0;
            while (z>0):
                z = z - memory[ii];
                ii = ii + 1;
            prewalk = (ii-1) - presteps;
        else:
            prewalk = prewalk + 2*round(np.random.random())-1;
        
        memory[presteps+prewalk] = memory[presteps+prewalk] + 1;

    steps = 20;
    random_walker = np.zeros((steps+1));
    memory_times = [];
    for i in range(steps):
        if (probability>np.random.random()):
            memory_times.append(i+1);
            z = int(np.random.random()*sum(memory))+1;
            ii = 0;
            while (z>0):
                z = z - memory[ii];
                ii = ii + 1;
            random_walker[i+1] = (ii-1) - presteps;
        else:
            random_walker[i+1] = random_walker[i] + 2*round(np.random.random())-1;
    
        memory[presteps+prewalk] = memory[presteps+prewalk] + 1;

    plt.close();
    plt.rcParams["font.family"] = "monospace";
    figure, axs = plt.subplots(1,figsize=(6,3),facecolor='#f5f5f5',linewidth=1.5,edgecolor='k')
    plt.rcParams['xtick.labelsize'] = 9; plt.rcParams['ytick.labelsize'] = 9;
    axs.plot(random_walker,'-',c='#333E74',lw=2,alpha=0.5);
    for i in range(len(memory_times)): axs.plot([memory_times[i]-1,memory_times[i]],\
                                                [random_walker[memory_times[i]-1],random_walker[memory_times[i]]],\
                                                    '-',c='r',lw=2,alpha=0.3);
    axs.plot(random_walker,'o',c='#333E74',ms=10,alpha=0.9);
    axs.set_xticks(range(steps+1)); axs.set_yticks(range(int(min(random_walker)-2),int(max(random_walker)+2)));
    axs.set_xticks(np.arange(steps+1)-0.5,minor=True);
    axs.set_yticks(np.arange(int(min(random_walker)-2),int(max(random_walker)+2))-0.5,minor=True);
    # axs.set_yticks(minor='true');
    axs.set_xlim([-0.5,20.5]); axs.set_ylim([min(random_walker)-1.5,max(random_walker)+1.5]);
    axs.grid(linestyle='--', which='minor', alpha=0.3);
    axs.set_xlabel(r'$\mathbf{n}$',fontsize=12);
    axs.set_ylabel(r'$\mathbf{X}$',fontsize=12);
    plt.suptitle('Individual Random Walk with Memory',x=0.53,y=0.95,fontsize=15,fontweight='bold');
    figure.tight_layout();
    plt.savefig('images/individual_memory.png', dpi=300, bbox_inches='tight', pad_inches=0.1);

    return plt.show()


def reading_data(name,type):
    short = type[0:3];
    # Reading Files
    distribution = pd.read_csv('data/'+type+'-distribution/distributions/'+name+'_distribution_'+short+'.txt', header=None, delim_whitespace=True, dtype=np.float64)
    variance = pd.read_csv('data/'+type+'-distribution/statistics/'+name+'_variance_'+short+'.txt', header=None, delim_whitespace=True, dtype=np.float64)
    skewness = pd.read_csv('data/'+type+'-distribution/statistics/'+name+'_kurtosis_'+short+'.txt', header=None, delim_whitespace=True, dtype=np.float64)
    kurtosis = pd.read_csv('data/'+type+'-distribution/statistics/'+name+'_kurtosis_'+short+'.txt', header=None, delim_whitespace=True, dtype=np.float64)
    information = pd.read_csv('data/'+type+'-distribution/information/'+name+'_info_'+short+'.txt', header=None,index_col=0, delim_whitespace=True).transpose()
    if (name in ['resetting','memory']):
        cases = round(pd.read_csv('data/'+type+'-distribution/information/'+name+'_cases_'+short+'.txt', header=None, delim_whitespace=True, dtype=np.float64),2)
    else:
        cases = (pd.read_csv('data/'+type+'-distribution/information/'+name+'_cases_'+short+'.txt', header=None, delim_whitespace=True, dtype=np.float64)).astype(int)

    # Adding some details to Pd.Frames
    distribution.rename(index={k: k-information['MaxDistance'][1] for k in distribution.index.values}, inplace=True);
    distribution.rename(columns={k: cases[0][k] for k in distribution.columns.values}, inplace=True);
    variance.rename(columns={k: cases[0][k] for k in variance.columns.values}, inplace=True);
    skewness.rename(columns={k: cases[0][k] for k in skewness.columns.values}, inplace=True);
    kurtosis.rename(columns={k: cases[0][k] for k in kurtosis.columns.values}, inplace=True);

    return information,cases,distribution,variance,skewness,kurtosis


def normalization(distribution,cases):
    for i in cases[0]:
        if not(sum(distribution[i])==1.):
            distribution[i] = distribution[i]/sum(distribution[i]);

    return distribution


def distribution_simple_plotting(name,type,distribution,variance,kurtosis,cases,information,xlim,ylim):
    if (type=='position'):
        text0 = r'${\rm P(x,t)}$';
    else:
        text0 = r'${\rm \mathcal{P}(x,t)}$';
    kwargs = dict(linewidth=1.5,histtype="stepfilled",fc="None", alpha=0.8);
    colors = {0:'#5A5B9f',1:'#009473',2:'#D94f70',3:'#F0C05A',4:'#7BC4C4',5:'#FF6F61'};
    colors_theory = {0:'#393377',1:'#006659',2:'#A92235',3:'#Af9433',4:'#498A91',5:'#Af5642'};
    plt.close();
    fig = plt.figure(figsize=(8,3),facecolor='#f5f5f5',linewidth=1.5,edgecolor='k');
    gs = fig.add_gridspec(2,2,height_ratios=[2,1]);
    ax0 = fig.add_subplot(gs[:,0]); ax1 = fig.add_subplot(gs[0,1]); ax2 = fig.add_subplot(gs[1,1]); axs = [ax0,ax1,ax2];
    axs[1].set_prop_cycle(cycler('color', list(colors.values()))); axs[2].set_prop_cycle(cycler('color', list(colors.values())));
    axs[0].grid(True,linestyle='--', alpha=0.3);
    axs[0].plot([],[],alpha=0.0);
    for i in range(len(cases[0])):
        axs[0].hist(distribution.index.values,weights=distribution[cases[0][i]],bins=2*information['MaxDistance'][1]+1,\
                    edgecolor=colors[i],**kwargs);

        if (type=='position'):
            axs[0].plot(distribution.index.values,np.exp(-(distribution.index.values)**2/(2*cases[0][i]))/np.sqrt(2*np.pi*cases[0][i]),\
            c=colors_theory[i],linestyle=(2,(2,2)),lw=1, label='_nolegend_');
        else:
            axs[0].plot(distribution.index.values,np.sqrt(2)*np.exp(-distribution.index.values**2/(2*cases[0][i]))/np.sqrt(np.pi*cases[0][i])\
            -abs(distribution.index.values)*(1 -sp.erf(abs(distribution.index.values)/np.sqrt(2*cases[0][i])))/cases[0][i],\
            c=colors_theory[i],linestyle=(2,(2,2)),lw=1, label='_nolegend_');
    
    axs[0].set_yscale('log');
    axs[0].set_ylim(1e-4,ylim[0]);
    axs[0].set_xlim([-xlim[0],xlim[0]]);
    axs[0].set_xlabel('Positions');
    axs[0].text(-xlim[0]+0.03*(2*xlim[0]),1e-4+0.6*(ylim[0]-1e-4),text0,fontsize=10);
    leg = axs[0].legend([r'$\bf{Steps:}$']+list(cases[0]),ncol=len(cases[0])+1,frameon=False,\
                        bbox_to_anchor=(1.6, 1.18),fontsize=8);

    if (type=='position'):
        axs[1].plot(variance.index.values,variance.index.values,'-',c='k',alpha=0.9,lw=1);
    else:
        axs[1].plot(variance.index.values,0.5*variance.index.values,'-',c='k',alpha=0.9,lw=1);
    axs[1].grid(True,linestyle='--', alpha=0.3);
    axs[1].plot(variance, '.', markeredgecolor='none',ms=3, alpha=0.8);
    axs[1].set_ylim([0,ylim[1]]);
    axs[1].set_xlim([0,xlim[1]]);
    axs[1].set_xticklabels([]);
    axs[1].text(0+0.03*(xlim[1]),0+0.87*(ylim[1]),r'${\rm M_2}$',fontsize=10);

    axs[2].plot(kurtosis, '.', markeredgecolor='none',ms=3, alpha=0.8);
    axs[2].set_xlim([0,xlim[1]]);
    axs[2].grid(True,linestyle='--', alpha=0.3);
    axs[2].set_ylim([ylim[2],ylim[3]]);
    axs[2].set_xlabel('Time');
    axs[2].text(0+0.03*(xlim[1]),ylim[2]+0.75*(ylim[3]-ylim[2]),r'${\rm M_4}$',fontsize=10);
    plt.suptitle('Simple Random Walk',x=0.5,y=1.1,fontsize=15,fontweight='bold');

    plt.savefig('images/'+name+'_distribution_'+type+'.png',dpi=300,bbox_inches='tight',transparent=False);

    return plt.show();


def distribution_probability_plotting(name,type,distribution,variance,kurtosis,cases,information,xlim,ylim):
    if (type=='position'):
        text0 = r'${\rm P(x,t)}$';
    else:
        text0 = r'${\rm \mathcal{P}(x,t)}$';
    kwargs = dict(linewidth=1.5,histtype="stepfilled",fc="None", alpha=0.8);
    colors = {0:'#5A5B9f',1:'#009473',2:'#D94f70',3:'#F0C05A',4:'#7BC4C4',5:'#FF6F61'};
    colors_theory = {0:'#393377',1:'#006659',2:'#A92235',3:'#Af9433',4:'#498A91',5:'#Af5642'};
    plt.close();
    fig = plt.figure(figsize=(8,3),facecolor='#f5f5f5',linewidth=1.5,edgecolor='k');
    gs = fig.add_gridspec(2,2,height_ratios=[2,1]);
    ax0 = fig.add_subplot(gs[:,0]); ax1 = fig.add_subplot(gs[0,1]); ax2 = fig.add_subplot(gs[1,1]); axs = [ax0,ax1,ax2];
    axs[1].set_prop_cycle(cycler('color', list(colors.values()))); axs[2].set_prop_cycle(cycler('color', list(colors.values())));
    axs[0].grid(True,linestyle='--', alpha=0.3);
    axs[0].plot([],[],alpha=0.0);

    axs[1].grid(True,linestyle='--', alpha=0.3);
    axs[1].plot(variance, '.', markeredgecolor='none',ms=3, alpha=0.8);

    for i in cases[0]:
        ii = np.ceil(i/0.05);
        axs[0].hist(distribution.index.values,weights=distribution[i],bins=2*information['MaxDistance'][1]+1,\
                    edgecolor=colors[ii],**kwargs);

        if (name=='resetting' and int(ii)!=0):
            alpha0 = np.sqrt(2*i/(1-i));
            axs[0].plot(distribution.index.values,alpha0*np.exp(-alpha0*abs(distribution.index.values))/2,\
            c=colors_theory[ii],linestyle=(2,(2,2)),lw=1, label='_nolegend_');
            axs[1].plot([0,xlim[1]],[(1-i)/i]*2,c=colors_theory[ii],linestyle=(2,(2,2)),lw=1, label='_nolegend_');
    
        elif all([name=='memory',type=='position',int(ii)!=0]):
            m2 = (1-i)*(sp.exp1(information['Steps'][1]*i)+np.log(information['Steps'][1]*i)+0.577215)/i;
            axs[0].plot(distribution.index.values,np.exp(-(distribution.index.values)**2/(2*m2))/np.sqrt(2*np.pi*m2),\
            c=colors_theory[ii],linestyle=(2,(2,2)),lw=1, label='_nolegend_')
            axs[1].plot((1-i)*(sp.exp1(variance.index.values[1:len(variance)]*i)+np.log(variance.index.values[1:len(variance)]*i)+0.577215)/i,\
            c=colors_theory[ii],linestyle=(2,(2,2)),lw=1, label='_nolegend_');

    axs[0].set_yscale('log');
    axs[0].set_ylim(1e-4,ylim[0]);
    axs[0].set_xlim([-xlim[0],xlim[0]]);
    axs[0].set_xlabel('Positions');
    axs[0].text(-xlim[0]+0.03*(2*xlim[0]),1e-4+0.6*(ylim[0]-1e-4),text0,fontsize=10);
    leg = axs[0].legend([r'$\bf{\gamma:}$']+list(cases[0]),ncol=len(cases[0])+1,frameon=False,\
                        bbox_to_anchor=(1.67, 1.18),fontsize=8);
    
    axs[1].set_ylim([0,ylim[1]]);
    axs[1].set_xlim([0,xlim[1]]);
    axs[1].set_xticklabels([]);
    axs[1].text(0+0.03*(xlim[1]),0+0.87*(ylim[1]),r'${\rm M_2}$',fontsize=10);

    axs[2].plot(kurtosis, '.', markeredgecolor='none',ms=3, alpha=0.8);
    axs[2].set_xlim([0,xlim[1]]);
    axs[2].grid(True,linestyle='--', alpha=0.3);
    axs[2].set_ylim([ylim[2],ylim[3]]);
    axs[2].set_xlabel('Time');
    axs[2].text(0+0.03*(xlim[1]),ylim[2]+0.75*(ylim[3]-ylim[2]),r'${\rm M_4}$',fontsize=10);
    plt.suptitle('Random Walk with '+name[0].upper()+name[1:len(name)],x=0.5,y=1.1,fontsize=15,fontweight='bold');

    plt.savefig('images/'+name+'_distribution_'+type+'.png',dpi=300,bbox_inches='tight',transparent=False);

    return plt.show();

