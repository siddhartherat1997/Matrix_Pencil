# -*- coding: utf-8 -*-
"""
 
@author: Prashant Pant, M.Sc
        PhD Candidate
        Centre for Combined Smart Energy Systems (CoSES)
        Research Coordinator-TUM SEED Centre
        Technical University of Munich

"""
from opcua import Client, ua

from time import sleep
from scipy.signal import resample
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import pandas as pd
from math import ceil
from scipy.linalg import  hankel, svd

 
def Pencil(Q,Ts):
    """
    Q a list which contains active power buffer. Ts is sampling time
    
    """
  
    
    N=len(Q)
    
    l=ceil(len(Q)/2) # pencil parameter =range N/2 to 2N/3 #rao limit
    print(l)
    print("N-L",N-l)
    print("L+1",l+1)
    #MAtrix Pencil
    
    #Y=hankel(Q[:l+2],Q[1:N-l])
    Y=hankel(Q[0:N-l],Q[N-l-1:N])
    
    
    u, s, vh = svd(Y,compute_uv=True,full_matrices=True)
    
    #correcting the shape of singular matrix "s"
    smat=np.zeros((u.shape[1],vh.shape[0]),dtype=complex)
    smat[:len(s),:len(s)]=np.diag(s)
    
    #rank=np.linalg.matrix_rank(Y,tol=10e-3) # 
    rank=0
    for i in range (0,len(s)):
        if s[i]/max(s) >10**-5:
            rank+=1
        
    
    
    V = np.conjugate(np.transpose(vh))
    
    
    
    # Filtering the matrices to reduce dimension
    s_new=smat[:,0:rank] # [:,start index, end index] for slicing columns , for slicing rows [star index: end index,:]
    
    v_new=V[:,0:rank]
    
    # creating V1 and V2 from filtered V and taking its complex conjugate
    v1= np.conjugate(np.transpose((np.delete(v_new,(-1), axis=0)))) # Vfilter matrix without last row
    v2=np.conjugate(np.transpose((np.delete(v_new,(0), axis=0)))) # Vfilter matrix without first row
    
    
    Y1=u@ s_new@  v1
    Y2=u@ s_new@  v2
    
    dfil=(np.linalg.pinv(Y1))@Y2
    
    #z=np.linalg.eigvals(dfil)
    z1, _ = np.linalg.eig(dfil)
    z=z1[0:rank]
    
    
    zabs=np.abs(z)
    
    zangle=np.angle(z,deg=False)
    
    zeta=((np.log(zabs))/Ts)[0:rank]
    freq=((zangle)/(2*np.pi*Ts))[0:rank]
    zp = zeta / np.sqrt(zeta**2 + (2 * np.pi * freq)**2)    
    
    
    #Residue calculations
    Z = np.zeros((N, rank), dtype = complex)
    for i in range(N):
        Z[i, :] = pow(z, i)
    epsilon = 1e-10  # Small regularization parameter
    B = np.linalg.inv(Z.T @ Z + epsilon * np.eye(Z.shape[1])) @ Z.T @ Q
 
    mag = 2.0 * np.abs(B)
    
    #Roots calculation for Contineous time equivalent
    roots = []
    for zeta_val, freq_val in zip(zeta, freq):
        real_part = zeta_val * abs(freq_val) * 2 * np.pi
    
        
        imaginary_part = freq_val * 2 * np.pi * np.sqrt(np.clip(1 - zeta_val**2, 0, None))
        roots.append(complex(real_part, imaginary_part))
 
 
# Print the calculated roots
   
    
    
    df = pd.DataFrame({
    "Frequency (Hz)": freq,
    #"Frequency (rad/s)": freq * 2 * np.pi,
    "zeta (%)": zp,
    "Eigen (discrete)": z[:rank], 
    #"Dominance index": abs(z[:rank]),
    "Magnitude": mag})  
    
    df_filtered = df[(abs(df["Frequency (Hz)"]) >= 0.1) & (abs(df["Frequency (Hz)"]) <= 10)]
    
    # Now sort the filtered DataFrame by the "Magnitude" column in descending order
    df_sorted = df_filtered.sort_values(by=['Magnitude'], ascending=False)
    
 
    
    df_filtered_sorted = df_sorted.reset_index(drop=True)
    #print(df_filtered_sorted)
    
  
    return df_filtered_sorted
    
    
# Accessing Server
url="opc.tcp://10.162.231.147:16664" # url of the server to be connected to
 
client=Client(url)
 
client.connect()
 
#print(irradiation)
 
 
##get nodes for sending data
zeta =client.get_node("ns=1;s=zeta")
os_freq =client.get_node("ns=1;s=os_freq")
 
 
q_Dummy=np.zeros(2000)
Q=[q_Dummy]
 
while True:
   
       # reading measurements
    Q_inv1= (client.get_node("ns=1;s=Q_inv1")).get_value()
   
    if np.count_nonzero(Q_inv1)!=2000:    
        pass
 
    else:
       
        Q.append(Q_inv1)
       
        if np.allclose(Q[0],Q[1],rtol=1e-5)==True:
            Q.pop()
           
        else:
            del Q[0]
           
            dfq = pd.DataFrame(Q[0])[0]
            length=len(dfq)
            downsample_factor = length // 200
            
            downsampled_data = dfq.iloc[::downsample_factor].reset_index(drop=True)
            
 
            dum=Pencil(downsampled_data[0:200],1/100)
           
            if dum.empty==True:
                pass
 
            else:
                zeta1 = ua.DataValue(ua.Variant(dum["zeta (%)"][0] , ua.VariantType.Float))
                os = ua.DataValue(ua.Variant(dum["Frequency (Hz)"][0] , ua.VariantType.Float))
                zeta.set_value(zeta1)
                os_freq.set_value(os)
                print(dum)
       
        sleep(0.01)