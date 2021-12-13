import sys
import numpy as np
import scipy.integrate

OFile = 'Outputs/scan_005.data'

# Parameters
species = 7
trep = 1E-1 # s
dt = 1E-12 # maximum timestep for solving the system
tfin = 1E+02 # Final Time
xini = np.zeros(species)
xini[0] = 1.0
xini[1] = 1.0
xini[2] = 0.1
# Model at T=298.15 K
def model(t,x):
    dxdt = np.zeros(7)
    
    #Constants
    k00 = 7.2645705800e+09
    k01 = 2.1243982518e+11
    k02 = 7.2645705800e+09
    k03 = 2.1243982518e+11
    k04 = 2.9372240517e-06
    k05 = 1.0044101389e-07
    k06 = 8.5894046618e-05
    k07 = 1.0044101389e-07
    k08 = 5.3717895536e+04
    k09 = 9.9335846365e+03
    k10 = 2.1243982518e+11
    k11 = 4.5937881793e+07
    
    #Ratelaws
    r00 = k00*x[0]*x[1]
    r01 = k01*x[3]
    r02 = k02*x[2]*x[3]
    r03 = k03*x[4]
    r04 = k04*x[0]*x[1]
    r05 = k05*x[6]
    r06 = k06*x[3]
    r07 = k07*x[6]
    r08 = k08*x[4]
    r09 = k09*x[5]
    r10 = k10*x[5]
    r11 = k11*x[2]*x[6]
    
    #MassBalances
    dxdt[0] = -r00+r01-r04+r05
    dxdt[1] = -r00+r01-r04+r05
    dxdt[2] = -r02+r03+r10-r11
    dxdt[3] = +r00-r01-r02+r03-r06+r07
    dxdt[4] = +r02-r03-r08+r09
    dxdt[5] = +r08-r09-r10+r11
    dxdt[6] = +r04-r05+r06-r07+r10-r11
    
    return dxdt
    
def jacobian(t,x):
    Jac = np.zeros(shape=(7,7))
    
    #Constants
    k00 = 7.2645705800e+09
    k01 = 2.1243982518e+11
    k02 = 7.2645705800e+09
    k03 = 2.1243982518e+11
    k04 = 2.9372240517e-06
    k05 = 1.0044101389e-07
    k06 = 8.5894046618e-05
    k07 = 1.0044101389e-07
    k08 = 5.3717895536e+04
    k09 = 9.9335846365e+03
    k10 = 2.1243982518e+11
    k11 = 4.5937881793e+07
    
    #Non-zero Elements
    Jac[0,0] = -k00*x[1]-k04*x[1]
    Jac[0,1] = -k00*x[0]-k04*x[0]
    Jac[0,3] = +k01
    Jac[0,6] = +k05
    Jac[1,0] = -k00*x[1]-k04*x[1]
    Jac[1,1] = -k00*x[0]-k04*x[0]
    Jac[1,3] = +k01
    Jac[1,6] = +k05
    Jac[2,2] = -k02*x[3]-k11*x[6]
    Jac[2,3] = -k02*x[2]
    Jac[2,4] = +k03
    Jac[2,5] = +k10
    Jac[2,6] = -k11*x[2]
    Jac[3,0] = +k00*x[1]
    Jac[3,1] = +k00*x[0]
    Jac[3,2] = -k02*x[3]
    Jac[3,3] = -k01-k02*x[2]-k06
    Jac[3,4] = +k03
    Jac[3,6] = +k07
    Jac[4,2] = +k02*x[3]
    Jac[4,3] = +k02*x[2]
    Jac[4,4] = -k03-k08
    Jac[4,5] = +k09
    Jac[5,2] = +k11*x[6]
    Jac[5,4] = +k08
    Jac[5,5] = -k09-k10
    Jac[5,6] = +k11*x[2]
    Jac[6,0] = +k04*x[1]
    Jac[6,1] = +k04*x[0]
    Jac[6,2] = -k11*x[6]
    Jac[6,3] = +k06
    Jac[6,5] = +k10
    Jac[6,6] = -k05-k07-k11*x[2]
    
    return Jac
    
t = np.arange(0,tfin,trep)
# Time indexes and Out predefinition
solution = scipy.integrate.solve_ivp(fun=model, jac=jacobian, y0=xini,
                                     t_span=(0,tfin), t_eval=t,
                                     method='LSODA', max_step=min(dt,trep),
                                     rtol=1E-6,atol=1E-12)

if not solution.success:
    print(solution.message)
else:
    print(f"""
          nfev={solution.nfev}
          njev={solution.njev}
          nlu={solution.nlu}
          status={solution.status}
          success={solution.success}
          """)
    x = np.zeros(shape=(len(solution.t),species+1))
    x[:,0] = solution.t
    x[:,1:] = solution.y[:,:].transpose()
    np.savetxt(OFile,x,delimiter='\t')
