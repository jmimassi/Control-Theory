import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize



#-----------------------------------      
def LL_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0):
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    
    The function "LL_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation.
    """    
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else :
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+Tlead/Ts)*MV[-1]-(Tlead/Ts)*MV[-2]))
            
    else:
        PV.append(Kp*MV[-1])
        
#-----------------------------------             

#-----------------------------------      
        
def PID_RT(SP,PV,Man,MV_Man,MV_FF,Kc,Ti,Td,alpha, Ts,MVMin,MVMax,MV,MV_P,MV_I,MV_D, E ,ManFF = False, PV_Init = 0,E_init = 0):  
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    Inputs :
    :SP: SP (or SetPoint) vector
    :PV: PV (or Processed Value) vector
    :Man: Man (or Manual controller mode) vector (True or False)
    :MV_Man: MV_Man (or Manual value for MV) vector
    :MV_FF: MV_FF (or FeedForward) vector
    
    
    Parameters : 
    :Kc: controller gain
    :Ti: integral time constant [s]
    :Td: derivative time constant [s]
    :alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s]
    :Ts: sampling period [s]
    :MVMin: Minimum value for MV (used for saturation and anti wind-up)
    :MVMax: Maximum value for MV (used for saturation and anti wind-up)
    :ManFF: Activated FF in manual mode (optional: default value is False)
    :PV_Init: Initial value for PV (optional : default value is 0) : used if PID_RT is ran first in the sequence and no value of PV are available yet.
    
    Outputs: 
    :MV: MV (or Manipulated Value ) vector
    :MV_P: MV_P (or Propotional part of MV) vector
    :MV_I: MV_I (or Integral part of MV) vector
    :MV_D: MV_D (or Derivative part of MV) vector
    :E: (or Controller error) vector
    
    
    
    The function "PID_RT" appends a value to the output vectors "MV", "MV_P", "MV_I", "MV_D".
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up.
    """   
    
    if len(PV) == 0 : 
        E.append(SP[-1]-PV_Init)
    else : 
        E.append(SP[-1]-PV[-1])
    
    #PID 
    ## Proportional action
    
    MV_P.append(Kc*E[-1])
    
    ## Integral action 
    
    if len(MV_I) == 0 : 
        if Ti > 0 :
            MV_I.append(((Kc*Ts)/(Ti))*(E[-1]))
        else : 
            MV_I = 0
    else : 
        if Ti > 0 : 
            MV_I.append(MV_I[-1]+((Kc*Ts)/(Ti))*(E[-1]))
        else : 
            MV_I = 0
    
    ## Derivative action
    
    Tfd = alpha*Td
    
    if len(MV_D) == 0 : 
        if Td > 0 : 
            MV_D.append(((Kc*Td)/(Tfd+Ts))*(E[-1]))
        else : 
            MV_D = 0
    else : 
        if Td > 0 : 
            MV_D.append(((Tfd)/(Tfd+Ts))*MV_D[-1] + ((Kc*Td)/(Tfd+Ts))*(E[-1]-E[-2]))
        else : 
            MV_D = 0
            

    if Man[-1] == True : 
        if ManFF : 
            MV_I[-1] = (MV_Man[-1]-MV_P[-1]-MV_D[-1])
        else : 
            MV_I[-1] = (MV_Man[-1]-MV_P[-1]-MV_D[-1]- MV_FF[-1])
        

    if ((MV_I[-1] + MV_P[-1] ) > MVMax) : 
        MV.append(MVMax - MV_P[-1])
    elif ((MV_I[-1] + MV_P[-1] )  < MVMin ) : 
        MV.append(MVMin - MV_P[-1])
    else :
        MV.append(MV_I[-1] + MV_D[-1] + MV_P[-1] + MV_FF[-1])
    
        

def IMC_Tuning(K, Tlag1, Tlag2,theta,gamma = 0.5) : 
    """
    IMC_Tuning(K,Tlag1,TLag2,theta,gamma=0.5)
    This function computes the optimised IMC PID tuning parameters for a SOPDT Process. 
    
    :Kp: process gain
    :Tlag1: first lag time constant [s]
    :Tlag2: second lag time constant [s]
    :theta: delay [s]
    :gamma: used to compute the closed loop time constant TCLP [s] such as TCLP = gamme*T1p with T1p = main time constant of the process. (range for gamma is [0.2 ... 0.9], default value is 0.5)
    
    
    returns (Kc,Ti,Td) that are the PID controller parameters
    
    """   
    Tc = gamma*Tlag1   ## Calcul ? 
    
    
    KcK = (Tlag1+Tlag2)/(theta + Tc)
    
    Kc = KcK/K
    
    Ti = (Tlag1 + Tlag2)
    
    Td = (Tlag1*Tlag2)/(Tlag1+Tlag2)
    
    return (Kc,Ti,Td)


#-----------------------------------      

#-----------------------------------      

class PID:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 0.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 0.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 1.0

#-----------------------------------      

#-----------------------------------      
        
        
def Margins(P,C,omega, Show = True):
    
    """
    :P: Process as defined by the class "Process".
        Use the following command to define the default process which is simply a unit gain process:
            P = Process({})
        
        A delay, two lead time constants and 2 lag constants can be added.
        
        Use the following commands for a SOPDT process:
            P.parameters['Kp'] = 1.1
            P.parameters['Tlag1'] = 10.0
            P.parameters['Tlag2'] = 2.0
            P.parameters['theta'] = 2.0
        
        Use the following commands for a unit gain Lead-lag process:
            P.parameters['Tlag1'] = 10.0        
            P.parameters['Tlead1'] = 15.0   
            
    :C: Controller as defined by the class "PID".
        Use the following command to define the default process which is simply a unit gain process:
            C = PID({})
        
        Use the following commands for a PID Controller:
            C.parameters['Kc'] = 1.1
            C.parameters['Ti'] = 10.0
            C.parameters['Td'] = 2.0
            C.parameters['alpha'] = 1
        
        Use the following commands for a unit gain PID:
            C.parameters['Ti'] = 10.0        
            C.parameters['Td'] = 15.0 
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ps (P(j omega)) (vector of complex numbers) is returned.
    
    The function "Margins" generates the Bode diagram of the process P*C and prints the margins of the system. 
    """     
    
    s = 1j*omega
    
    Ptheta = np.exp(-P.parameters['theta']*s)
    PGain = P.parameters['Kp']*np.ones_like(Ptheta)
    PLag1 = 1/(P.parameters['Tlag1']*s + 1)
    PLag2 = 1/(P.parameters['Tlag2']*s + 1)
    PLead1 = P.parameters['Tlead1']*s + 1
    PLead2 = P.parameters['Tlead2']*s + 1
    
    Ps = np.multiply(Ptheta,PGain)
    Ps = np.multiply(Ps,PLag1)
    Ps = np.multiply(Ps,PLag2)
    Ps = np.multiply(Ps,PLead1)
    Ps = np.multiply(Ps,PLead2)
    
    
    CGain = C.parameters['Kc']*np.ones_like(Ptheta)
    C1 = (1+1/(C.parameters['Ti']*s)+ (C.parameters['Td']*s/(C.parameters['alpha']*C.parameters['Td']*s+1)))
   
    
    Cs = np.multiply(CGain,C1)

    
    Ps = np.multiply(Ps,Cs)
    
    
    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Gain part
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='P(s)')
        ax_gain.semilogx(omega,20*np.log10(np.abs(np.ones_like(PGain))),label='0')
        
        gain = 20*np.log10(np.abs(Ps))
        
 
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |P| [db]')
        ax_gain.set_title('Bode plot of P')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='P(s)')
        phase = (180/np.pi)*np.unwrap(np.angle(Ps))
        
        indexc = np.where(np.round(gain,2) == 0.00)[0][0]
        wc = omega[indexc]
        indexu = np.where(np.round(phase,1) == -180.0)[0][0]
        wu = omega[indexu]
        ax_phase.vlines(wc,-180,phase[indexc])
        print('Gain margin : ',np.round(gain[indexu],5),' at ', np.round(wu,5), ' rad/s')
        print('Phase margin : ',np.round(phase[indexc],5)+180,' at ', np.round(wc,5), ' rad/s')
        
        ax_gain.vlines(wu,0,gain[indexu])
        
        
        ax_phase.axhline(y=(-180))
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle P$ [°]')
        ax_phase.legend(loc='best')
    else:
        return Ps
    
#-----------------------------------      


def GrapheMVF(): 
    nameFile = 'Open_loop_experiment_on_MV_2022-03-12-22h08.txt'

    if 'MV' in nameFile:
        ExpVariable = 'MV'
    else:    
        ExpVariable = 'DV'
    
    print(ExpVariable)    
    
    titleName = nameFile.split('.')[0]    
    data = pd.read_csv('Data/' + nameFile)
    t = data['t'].values - data['t'].values[0]
    MV = data['MV'].values
    PV = data['PV'].values
    DV = data['DV'].values
    
    if ExpVariable == 'MV':
        tstep = np.argwhere(np.diff(MV) != 0)
        tstep = tstep[0][0]
        tm = t[tstep:]
        tm = tm - tm[0]    
        MVstep = MV[tstep + 1] - MV[tstep]
        MVm = MV[tstep:]
        PVm = PV[tstep:]
        PVm = (PVm - PVm[0])/MVstep
        MVm = (MVm - MVm[0])/MVstep    
    else:    
        tstep = np.argwhere(np.diff(DV) != 0)
        tstep = tstep[0][0]
        tm = t[tstep:]
        tm = tm - tm[0]
        DVstep = DV[tstep + 1] - DV[tstep]    
        DVm = DV[tstep:]
        PVm = PV[tstep:]
        PVm = (PVm - PVm[0])/DVstep
        DVm = (DVm - DVm[0])/DVstep
        
    return tm,PVm

def GrapheDVF(): 
    nameFile = 'Open_loop_experiment_on_DV_2022-03-12-22h35.txt'

    if 'MV' in nameFile:
        ExpVariable = 'MV'
    else:    
        ExpVariable = 'DV'
    
    print(ExpVariable)    
    
    titleName = nameFile.split('.')[0]    
    data = pd.read_csv('Data/' + nameFile)
    t = data['t'].values - data['t'].values[0]
    MV = data['MV'].values
    PV = data['PV'].values
    DV = data['DV'].values
    
    if ExpVariable == 'MV':
        tstep = np.argwhere(np.diff(MV) != 0)
        tstep = tstep[0][0]
        tm = t[tstep:]
        tm = tm - tm[0]    
        MVstep = MV[tstep + 1] - MV[tstep]
        MVm = MV[tstep:]
        PVm = PV[tstep:]
        PVm = (PVm - PVm[0])/MVstep
        MVm = (MVm - MVm[0])/MVstep    
    else:    
        tstep = np.argwhere(np.diff(DV) != 0)
        tstep = tstep[0][0]
        tm = t[tstep:]
        tm = tm - tm[0]
        DVstep = DV[tstep + 1] - DV[tstep]    
        DVm = DV[tstep:]
        PVm = PV[tstep:]
        PVm = (PVm - PVm[0])/DVstep
        DVm = (DVm - DVm[0])/DVstep
        
    return tm,PVm