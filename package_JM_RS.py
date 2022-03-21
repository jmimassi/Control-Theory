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
        
        
        
        
def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,alpha, Ts,MVMin,MVMax,MV,MV_P,MV_I,MV_D, E ,MV_FF = False, PVInit = 0,E_init = 0,MV_P_Init = 0,MV_I_Init = 0,MV_D_Init = 0, Method = 'EBD-EBD'):  
    
    if len(E) == 0 : 
        E.append(E_init)
    else : 
        E.append(SP[-1]-PV[-1])
    
    #PID 
    ## Proportional action
    if len(MV_P) == 0 : 
        MV_P.append(MV_P_Init)
    else : 
        MV_P.append(Kc*E[-1])
    
    ## Integral action 
    
    if len(MV_I) == 0 : 
        MV_I.append(MV_I_Init)
    else : 
        if Ti > 0 : 
            MV_I.append(MV_I[-1]+((Kc*Ts)/Ti)*(E[-1]+E[-2]))
        else : 
            MV_I = 0
    
    ## Derivative action
    Tfd = alpha*Td
    if len(MV_D) == 0 : 
        MV_D.append(MV_D_Init)
    else : 
        if Ti > 0 : 
            MV_D.append(((Tfd-Ts/2)/(Tfd+Ts/2))*MV_D[-1]+((Kc*Td)/(Tfd+Ts/2))*(E[-1]+E[-2]))
        else : 
            MV_D = 0
            
    if Man : 
        MV_I.append(MV_Man-MV_P-MV_D-MV_FF)
    
    MV.append(MV_I[-1] + MV_D[-1] + MV_P[-1] + MV_FF)
    
        
        