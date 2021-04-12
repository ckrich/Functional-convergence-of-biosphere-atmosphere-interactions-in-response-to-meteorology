import numpy as np
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr, GPDC, GPACE, CMIknn, CMIsymb
from tigramite.models import LinearMediation, Prediction

import argparse
import os
import datetime
from datetime import date
from netCDF4 import Dataset as NetCDFFile
import scipy.fftpack
from scipy.signal import argrelextrema
import csv
import numpy.ma as ma

#######################
#commandline argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("Dataset",help = "Which FLUXCOM Dataset was used, long or short?")
parser.add_argument("attribute", help="some characteristic")
#parser.add_argument("variables", help="Which variables shall be used for network construction")
parser.add_argument("tau_min",type=int, help = "Minimum lag that shall be considered for network reconstruction ")
parser.add_argument("tau_max",type=int, help = "Maximum lag that shall be considered for network reconstruction ")
parser.add_argument("mask_type", help= "{y,x,z,xy,xz,yz,xyz}) – Masking mode: Indicators for which variables in the dependence measure I(X; Y | Z) the samples should be masked. If None, ‘y’ is used, which excludes all time slices containing masked samples in Y. ")
parser.add_argument("extent", help = "inter or intra annual(step=365)/seasonal(step=91) networks")
parser.add_argument("minlength", type=int, help = "Which minimum TS length shall be used?")
#parser.add_argument("pc_alpha", type=float, help = "For which pc_alpha values shall PCMCI be calculated")
parser.add_argument("method",help = "which conditional independence test shall be used?")
parser.add_argument("step",type=int, help = "length of TS excerpt?")
parser.add_argument("processing",help = "Do you want the data to be log transformed? Then type log, else leave empty")
#parser.add_argument("hour_mask_use", help= "which masking for the day/hours shall be used?")
parser.add_argument("event_mask_use", help= "any special Event you wanna lock at? Create a mask for it.")
parser.add_argument("quality_value",type=float, help = "minimum quality flag falue")
parser.add_argument("attribute2", help = "information for Networkgenerating process")
#parser.add_argument("aggregate", help = "aggregate halfhourly data to daily data")

#parser.add_argument("MovingWindow", type=bool, help = "Sha")

parsed_args=parser.parse_args()
#########################

########
#function calling PCMCI
def RunPCMCI(con_ind_test,selected_variables,dataframe,pc_alpha,var_names,tau_min,tau_max,alpha_level):
    pcmci= PCMCI(selected_variables = selected_variables,
                  dataframe=dataframe,
                  cond_ind_test=con_ind_test,
                  var_names=var_names,
                      verbosity=1)
    results = pcmci.run_pcmci(tau_min=tau_min,tau_max=tau_max,pc_alpha=pc_alpha, fdr_method='fdr_bh')


    pcmci._print_significant_links(p_matrix=results['p_matrix'],q_matrix = results['q_matrix'],
           val_matrix=results['val_matrix'], alpha_level=alpha_level)

    link_matrix_q = pcmci._return_significant_parents(pq_matrix=results['q_matrix'],
            val_matrix=results['val_matrix'], alpha_level=alpha_level)['link_matrix']
    return results, link_matrix_q

def RunMCI(con_ind_test,selected_variables,dataframe,pc_alpha,var_names,tau_min,tau_max,alpha_level):
    pcmci= PCMCI(selected_variables = selected_variables,
                  dataframe=dataframe,
                  cond_ind_test=con_ind_test,
                  var_names=var_names,
                      verbosity=1)
    results = pcmci.run_mci(tau_min=tau_min,tau_max=tau_max,parents=None)

    pcmci._print_significant_links(p_matrix=results['p_matrix'],
           val_matrix=results['val_matrix'], alpha_level=alpha_level)

    link_matrix_q = pcmci._return_significant_parents(pq_matrix=results['p_matrix'],
            val_matrix=results['val_matrix'], alpha_level=alpha_level)['link_matrix']
    return results, link_matrix_q
##########

####
#create dictionary with cond ind tests
use_mask=True
mask_type=parsed_args.mask_type
con_ind_tests={'parcor': ParCorr(significance='analytic',use_mask=use_mask, mask_type=mask_type, confidence='analytic', verbosity=1),
      'gpdc': GPDC(significance='analytic', gp_params=None,verbosity=1,use_mask=use_mask, confidence='bootstrap', mask_type='y',recycle_residuals=True),
      'gpace' : GPACE(significance='analytic', gp_params=None,verbosity=1,use_mask=use_mask, mask_type='y', ace_version='acepack',recycle_residuals=True),
       'cmiknn' : CMIknn(significance='shuffle_test', knn=50,recycle_residuals=True)
      }
method=parsed_args.method
con_ind_test = con_ind_tests[method]
##################

quality_value=parsed_args.quality_value
attribute=parsed_args.attribute
Dataset=parsed_args.Dataset
attribute2=parsed_args.attribute2 # use _ before : _attribute
processing=parsed_args.processing
step=parsed_args.step
minlength=parsed_args.minlength
event_mask_use=[parsed_args.event_mask_use]

alpha_level=0.01
tau_max=parsed_args.tau_max #int(math.ceil(np.amax(minimizer[2:4,:])))
tau_min=parsed_args.tau_min

if Dataset=="short":
    with open('../Data/FluxtowerLists/Without-9999-short.txt', 'r') as towerfile:
        tower=towerfile.read().split('\n')
elif Dataset=="long":
    with open('../Data/FluxtowerLists/Without-9999-long.txt', 'r') as towerfile:
        #tower=towerfile.read().split('\n')
        tower=os.listdir( "/Net/Groups/BGI/data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/Data/DD/2016_11")

#tower=['DE-Hai.DD.2000.2012.nc']#, 'DE-Tha.DD.1996.2014.nc']

savepath="/Net/Groups/BGI/scratch/ckrich/Projects/EcosystemCharacterisationByCausalNetworks"

for t,fluxtower in enumerate(tower):  # -1 because last line in tower is empty.
    filename='%s%s' % (fluxtower, attribute)
    if Dataset=="long":
        path="/Net/Groups/BGI/data/DataStructureMDI/DATA/Incoming/Fluxnet/berkeley_012016/Data/DD/2016_11"
        if os.path.isfile('%s/%s' %  (path, tower[t])) == False:
            print('not found: %s/%s' %  (path, tower[t]))
            continue
        print("found",fluxtower)
        ts=NetCDFFile('%s/%s' % (path,tower[t]))
        L=len(ts.dimensions["time"])
        if L/minlength<365:
            print("TS too short")
            continue

        print("readingData")
        Rg=np.squeeze(ts.variables["SW_IN_F"][:]*(24*60*60/(10^6))) #conert to MJ per day per m^2
        Rg_q=np.squeeze(ts.variables["SW_IN_F_QC"][:]) #conert to MJ per day per m^2
        Rg_pot=np.squeeze(ts.variables["SW_IN_POT"][:]*(24*60*60/(10^6)))
        T=np.squeeze(ts.variables["TA_F"][:])#[:]
        T_q=np.squeeze(ts.variables["TA_F_QC"][:])#[:]
        GPP=np.squeeze(ts.variables["GPP_NT_CUT_MEAN"][:])
        RECO=np.squeeze(ts.variables["RECO_NT_CUT_MEAN"][:])
        NEE_NT=np.squeeze(ts.variables["NEE_VUT_REF_NIGHT"][:])
        NEE=np.squeeze(ts.variables["NEE_VUT_USTAR50"][:])#[:]
        NEE_q=np.squeeze(ts.variables["NEE_VUT_USTAR50_QC"][:])#[:]
        VPD=np.squeeze(ts.variables["VPD_F"][:])#[:]
        VPD_q=np.squeeze(ts.variables["VPD_F_QC"][:])#[:]
        H=np.squeeze(ts.variables["H_CORR"][:])#[:]
        H_q=np.squeeze(ts.variables["H_F_MDS_QC"][:])#[:]
        LE=np.squeeze(ts.variables["LE_CORR"][:])#[:]
        LE_q=np.squeeze(ts.variables["LE_F_MDS_QC"][:])#[:]
        month=np.squeeze(ts.variables["month"])#
        day=np.squeeze(ts.variables["day"])#;
        year=np.squeeze(ts.variables["year"])#[:];
        P=np.squeeze(ts.variables["P_F"])
        P_q=np.squeeze(ts.variables["P_F_QC"])
        SWC1=np.squeeze(ts.variables["SWC_F_MDS_1"])
        SWC2=np.squeeze(ts.variables["SWC_F_MDS_2"])
        SWC3=np.squeeze(ts.variables["SWC_F_MDS_3"])
        SWC1_q=np.squeeze(ts.variables["SWC_F_MDS_1_QC"])
        SWC2_q=np.squeeze(ts.variables["SWC_F_MDS_2_QC"])
        SWC3_q=np.squeeze(ts.variables["SWC_F_MDS_3_QC"])
        print("finishedreadingData")

    print("Data")
    var_names = ["Rg","T","NEE","VPD","H","LE"]
    Data=np.stack((Rg,T,NEE,VPD,H,LE),axis=1)
    L,W,=Data.shape
    ys=np.unique(year,return_inverse=True)
    if processing=='_meanseason' or processing=='_meanseason_smooth':
        L=365

    if not os.path.exists('%s/Outputs/%s/%s/Data/%s/%s' % (savepath,  var_names, fluxtower, method, step)):
        os.makedirs('%s/Outputs/%s/%s/Data/%s/%s' % (savepath, var_names, fluxtower, method, step))
    if not os.path.exists('%s/Outputs/%s/%s/Plots/%s/%s' % (savepath, var_names, fluxtower, method, step)):
        os.makedirs('%s/Outputs/%s/%s/Plots/%s/%s' % (savepath, var_names, fluxtower, method, step))

    print("creating mask")

    #Create season dataframe
    season_mask_prel = np.ones((L,W,4))
    for i in range(0, L):
        if month[i]==12 or month[i]==1 or month[i]==2:
            season_mask_prel[i,:,0]=False
        elif month[i]==3 or month[i]==4 or month[i]==5:
            season_mask_prel[i,:,1]=False
        elif month[i]==6 or month[i]==7 or month[i]==8:
            season_mask_prel[i,:,2]=False
        elif month[i]==9 or month[i]==10 or month[i]==11:
            season_mask_prel[i,:,3]=False

    season_mask={'': np.zeros((L,W)), 'winter':season_mask_prel[:,:,0], 'spring':season_mask_prel[:,:,1], 'summer':season_mask_prel[:,:,2], 'autumn':season_mask_prel[:,:,3]}

    #####Create an event mask
    event_mask_prel=np.ones((L,W,6))
    for i in range(0,W):
        event_mask_prel[:,i,0]=((month==8)*(day>=1)*(day<=13)==0)
        event_mask_prel[:,i,5]=(((month==7)*(day>=25)+(month==8)*(day<=20))==0)
        event_mask_prel[:,i,1]=((month>=7)*(month<=9)*(year==2003)==0)
        event_mask_prel[:,i,2]=((month>=6)*(month<=10)*(year==2003)==0)
        event_mask_prel[:,i,3]=((month>=6)*(month<=10)*(year==2004)==0)
        event_mask_prel[:,i,4]=((month>=6)*(month<=10)*(year==2002)==0)


    event_mask={'None': np.zeros((L,W)), 'heatwave-2003': event_mask_prel[:,:,0],'25.07-20.08': event_mask_prel[:,:,5], 'Jul-Sep-2003': event_mask_prel[:,:,1], 'Jun-Oct-2003': event_mask_prel[:,:,2],'Jun-Oct-2004': event_mask_prel[:,:,3],'Jun-Oct-2002': event_mask_prel[:,:,4]}

    ##### Yearly mask
    extent=parsed_args.extent
    ys=np.unique(year).astype(int)
    year_mask={'%s'%y: np.repeat(np.reshape(((year==y)==0),(L,1)),W,axis=1) for (i,y) in enumerate(ys)} # add all years to mask
    year_mask["%s-%s"%(min(ys),max(ys))]=np.zeros((L,W)) #add zero mask

    ###### Quality Flag data_mask
    quality_mask=np.ones((L,W))
    for (i,v) in enumerate(var_names):
        quality_mask[:,i]=eval(v+"_q")<0.9
        Data[quality_mask[:,i]==1,i]=-11


### processing Data
    if processing=='_anomalysmooth':
        #year_mask with yearly values is needed, though if extent=inter is choosen, above year_mask is useless. therefore generate an own mask
        ys2=np.unique(year,return_inverse=True)
        year_mask2 = np.ones((L,W,len(ys2[0])))
        dummy=np.ones((365,W,len(ys2[0])))*(-9999)
        meanyear=np.zeros((365,W))
        anomaly=ma.copy(Data)
        #create year mask
        for (i,y) in enumerate(ys2[0].astype(int)):
            year_mask2[np.where(year==y),:,i]=0
        # create year masked arrays
        data_frame={y: pp.DataFrame(Data,mask=year_mask2[:,:,i]) for (i,y) in enumerate(ys2[0].astype(int))}
        # calculate meanyear with masked missing values
        for (i,y) in enumerate(ys2[0].astype(int)):
            dummy[:,:,i]=data_frame[y].values[np.where(data_frame[y].mask[:,0]==0)][0:365]
        meanyear=ma.masked_values(dummy,-9999).mean(axis=2)
        #smooth meanyear
        meanyear_ft= scipy.fftpack.rfft(meanyear,axis=0)
        meanyear_ft[20:365]=0
        meanyear_smoth=scipy.fftpack.irfft(meanyear_ft,axis=0)
        #then subtract season mean from each year [0:365]-> leap year day 366 stays unchanged
        for (i,y) in enumerate(ys2[0].astype(int)):
            anomaly[np.where(ys2[1]==i)[0][0:365]]-=meanyear_smoth
        data=anomaly
        #data_frame={y: pp.DataFrame(data,mask=year_mask[:,:,i]) for (i,y) in enumerate(ys)}  #very important. here ys refers to global year_mask
    #    np.savetxt('%s/Outputs/%s/%s/Data/01_TS_ano_%s.txt' % (savepath,var_names, tower[t], filename), np.concatenate((data,month),1))
    elif processing=="_meanseason":
        #year_mask with yearly values is needed, though if extent=inter is choosen, above year_mask is useless. therefore generate an own mask
        ys2=np.unique(year,return_inverse=True)
        year_mask2 = np.ones((L,W,len(ys2[0])))
        dummy=np.ones((365,W,len(ys2[0])))*(-9999)
        meanyear=np.zeros((365,W))
        anomaly=ma.copy(Data)
        #create year mask
        for (i,y) in enumerate(ys2[0].astype(int)):
            year_mask2[np.where(year==y),:,i]=0
        # create year masked arrays
        data_frame={y: pp.DataFrame(Data,mask=year_mask2[:,:,i]) for (i,y) in enumerate(ys2[0].astype(int))}
        # calculate meanyear with masked missing values
        for (i,y) in enumerate(ys2[0].astype(int)):
            dummy[:,:,i]=data_frame[y].values[np.where(data_frame[y].mask[:,0]==0)][0:365]
        meanyear=ma.masked_values(dummy,-9999).mean(axis=2)
        data=meanyear
    elif processing=="_meanseason_smooth":
        #year_mask with yearly values is needed, though if extent=inter is choosen, above year_mask is useless. therefore generate an own mask
        ys2=np.unique(year,return_inverse=True)
        year_mask2 = np.ones((L,W,len(ys2[0])))
        dummy=np.ones((365,W,len(ys2[0])))*(-9999)
        meanyear=np.zeros((365,W))
        anomaly=ma.copy(Data)
        #create year mask
        for (i,y) in enumerate(ys2[0].astype(int)):
            year_mask2[np.where(year==y),:,i]=0
        # create year masked arrays
        data_frame={y: pp.DataFrame(Data,mask=year_mask2[:,:,i]) for (i,y) in enumerate(ys2[0].astype(int))}
        # calculate meanyear with masked missing values
        for (i,y) in enumerate(ys2[0].astype(int)):
            dummy[:,:,i]=data_frame[y].values[np.where(data_frame[y].mask[:,0]==0)][0:365]
        meanyear=ma.masked_values(dummy,-9999).mean(axis=2)
        #smooth meanyear
        meanyear_ft= scipy.fftpack.rfft(meanyear,axis=0)
        meanyear_ft[20:365]=0
        meanyear_smoth=scipy.fftpack.irfft(meanyear_ft,axis=0)
        data=meanyear_smoth

#Start network construction
    if processing=="_MovingWindow" or processing=="_MovingWindow_anomaly" or processing=="_MovingWindow_normalized":
        print('Start MovingWindow process')
        for pc_alpha in [parsed_args.pc_alpha]: #np.linspace(0.1,1,10):
            alpha_level=0.01
            tau_max=parsed_args.tau_max #int(math.ceil(np.amax(minimizer[2:4,:])))
            tau_min=parsed_args.tau_min
            selected_variables = [0,1,2,3,4,5]

            val_matrizes=np.zeros((W,W,tau_max+1,L-step))
            q_matrizes=np.zeros((W,W,tau_max+1,L-step))
            links_matrizes=np.zeros((W,W,tau_max+1,L-step))
            for m in range(0,L-step,windowstep):
                data_mask = np.ones((L,W))
                data_mask[m:m+step,:]=((window_mask==0)*(hour_mask[h]==0)*(event_mask[e]==0)*(quality_mask==0)==0)
                dataframe = pp.DataFrame(data, mask=data_mask)
                con_ind_test=ParCorr(significance='analytic',use_mask=True, mask_type=mask_type, verbosity=1)
                try:
                    results, links = Run(con_ind_test, selected_variables, dataframe, pc_alpha)
                except:
                    continue

                val_matrizes[:,:,:,m]=results['val_matrix']
                q_matrizes[:,:,:,m]=results['q_matrix']
                links_matrizes[:,:,:,m]=links

            result_gathering={'val_matrizes' : val_matrizes, 'q_matrizes' : q_matrizes, 'links_matrizes' : links_matrizes}
            np.save('../Outputs/%s/%s/Data/%s/%s/01_Results_%s%s_%s_%s_%s_%s_%s_%s%s.npy' % (var_names, tower[t], method, step, filename, processing, mask_type, step, method, pc_alpha,tau_min, tau_max, attribute2), result_gathering)

            print('Produced: 01_Results_%s%s_%s_%s_%s_%s_%s_%s%s.npy' % (filename, processing, mask_type, step, method, pc_alpha, tau_min, tau_max, attribute2))

    else:
        print("Start with %s" % method)
        #con_ind_test=parcor
        #method='parcor'

        #check if Nulldist exist, otherwise generate it.
#        if method='GPDC' or method='GPACE':
#            PATH='../Outputs/%s/%s/Data/%s_nulldists.npz' % (, var_names, tower[t], method)
#            if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
#                print("NulldistFile exists and is readable")
#            else:
#                con_ind_test.generate_and_save_nulldists(sample_sizes=range( Z-2*tau_max-1, Z-2*tau_max+1),
#                null_dist_filename ='../Data/%s_nulldists.npz' % method)
            #####


        if extent=="inter" or processing=="_meanseason_smooth" or processing=="_meanseason":
            year_mask_use = ["%s-%s"%(min(ys),max(ys))]
        elif extent=="intra":
            year_mask_use = ['%s' % y for (i,y) in enumerate(ys)]

        if step==91:
            season_mask_use=["winter", "spring","summer","autumn"] #["summer"] #
        elif step==365:
            season_mask_use=['']

        selected_variables = np.arange(1,W)

        for pc_alpha in [None]: #parsed_args.pc_alpha]: #np.linspace(0.1,1,10):
        #    print('alpha loop')
            for (i,y) in enumerate(year_mask_use):
                print('ys loop %s' % y)
                for (j,s) in enumerate(season_mask_use):
                    print('season loop %s' % s)
                    for (l,e) in enumerate(event_mask_use):
                        print('event loop %s' % e)
                        data_mask=((np.roll(year_mask[y]==0,-31*W))*(season_mask[s]==0)*(event_mask[e]==0)*(quality_mask==0)==0)
                        #data_intermed=data_frame[y].values
                        frame=pp.DataFrame(data,mask=data_mask,missing_flag=-9999)

            #            tp.plot_timeseries(frame.values[np.roll(year_mask[y][:,0]==0,-31*48*W)*(season_mask[s][:,0]==0)*(event_mask[e][:,0]==0),:], figsize=(16,16),  var_names=var_names, use_mask=True, mask=frame.mask[np.roll(year_mask[y][:,0]==0,-31*48*W)*(season_mask[s][:,0]==0)*(event_mask[e][:,0]==0),:],
            #                 grey_masked_samples='data',
                        #         var_units=[r"$\frac{MJ}{m\^2 day}$",r"$\frac{MJ}{m\^2 day}$",r"$^\circ$C", r"$\frac{gC}{day m^2}$", "hPa", r"$\frac{W}{m^2}$", r"$\frac{W}{m^2}$"],
                             #var_units=[r"$\frac{MJ}{m\^2 day}$",r"$\frac{MJ}{m\^2 day}$",r"$^\circ$C", r"$\frac{gC}{day m^2}$", "%", "mm"],
                             #var_units=[r"$\frac{MJ}{m\^2 day}$",r"$^\circ$C", r"$\frac{gC}{day m^2}$", "%","%","%", "mm", "hPa", r"$\frac{W}{m^2}$"],
                            # var_units=[r"$\frac{MJ}{m\^2 day}$",r"$^\circ$C", r"$\frac{gC}{day m^2}$", "%", "mm", "hPa", r"$\frac{W}{m^2}$"],
            #                 var_units=[r"$\frac{MJ}{m\^2 day}$",r"$^\circ$C", r"$\frac{gC}{day m^2}$", "%", r"$\frac{W}{m^2}$", r"$\frac{W}{m^2}$"],
            #                 save_name='../Outputs/%s/%s/Plots/%s/%s/01_TS_%s%s_%s_%s_%s_%s_%s_%s_%s%s.png' % (var_names, tower[t], method, step, filename, processing, y,s,e, method, pc_alpha, tau_min, tau_max, attribute2))
            #            plt.close()

                        try:
                            results, links = RunPCMCI(con_ind_test, selected_variables, frame, pc_alpha, var_names,tau_min,tau_max,alpha_level)
                        except:
                            continue


                        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7,7))
                        fig.subplots_adjust(top=0.88)
                        axes.set_title('%s, %s, %s \n Year %s %s %s, \n %s, PC_alpha %s , tau %s-%s ,%s  ' % (method, tower[t], processing, y,s, e, mask_type, pc_alpha, tau_min, tau_max, attribute2))

                        tp.plot_graph(
                            fig_ax=(fig,axes),
                            val_matrix=results['val_matrix'],
                            link_matrix=links,
                            var_names=var_names,
                            link_colorbar_label='cross-MCI',
                            node_colorbar_label='auto-MCI',
                            arrow_linewidth=50,
                            node_size=40,
                            label_fontsize=15,
                            node_label_size=20,
                        #    curved_radius=0.5,
                            link_label_fontsize=12,
                            show_colorbar=True,
                            )


                        fig.savefig('%s/Outputs/%s/%s/Plots/%s/%s/01_Graph_%s%s_%s_%s_%s_%s_%s_%s_%s_%s_%s%s.pdf' % (savepath,var_names, tower[t], method, step, filename, processing, mask_type,quality_value, y,s, e,method, pc_alpha, tau_min, tau_max, attribute2), bbox_inches='tight')
                        np.save('%s/Outputs/%s/%s/Data/%s/%s/01_Results_%s%s_%s_%s_%s_%s_%s_%s_%s_%s_%s%s.npy' % (savepath,var_names, tower[t], method, step, filename, processing, mask_type,quality_value, y,s, e,method, pc_alpha, tau_min, tau_max, attribute2), results)
                        plt.close()
                        print('Produced: 01_Results__%s%s_%s_%s_%s_%s_%s_%s_%s_%s_%s%s.pdf' % (filename, processing, mask_type,quality_value, y,s, e,method, pc_alpha, tau_min, tau_max, attribute2))



print("finished %s" % tower[t])
print("finished networkConstruction")
