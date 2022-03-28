def LL_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "FO_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "FO_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+Tlead/Ts)*MV[-1]-(Tlead/Ts)*MV[-2]))
            elif method == 'EFD':
                #PV.append((1-K)*PV[-1] + K*Kp*MV[-2])
                pass
            elif method == 'TRAP':
                #PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))    
                pass
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
    else:
        PV.append(Kp*MV[-1])
        
        
        
        
def PID_RT(SP,PV,Man,MV_Man,MV_FF,Kc,Ti,Td,alpha, Ts,MVMin,MVMax,MV,MV_P,MV_I,MV_D, E ,ManFF = False, PV_Init = 0,E_init = 0, Method = 'EBD-EBD'):  
    
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
        

    if ((MV_I[-1] + MV_D[-1] + MV_P[-1] + MV_FF[-1]) > MVMax) : 
        MV.append(MVMax - MV_P[-1])
    elif ((MV_I[-1] + MV_D[-1] + MV_P[-1] + MV_FF[-1]) < MVMin ) : 
        MV.append(MVMin - MV_P[-1])
    else :
        MV.append(MV_I[-1] + MV_D[-1] + MV_P[-1] + MV_FF[-1])
    
        

def IMC_Tuning(K, Tlag1, Tlag2=0,theta=0,gamma = 0.5) : 
    Tc = gamma*Tlag1   ## Calcul ? 
    
    
    KcK = (Tlag1+Tlag2)/(theta + Tc)
    
    Kc = KcK/K
    
    Ti = (Tlag1 + Tlag2)
    
    Td = (Tlag1*Tlag2)/(Tlag1+Tlag2)
    
    return (Kc,Ti,Td)