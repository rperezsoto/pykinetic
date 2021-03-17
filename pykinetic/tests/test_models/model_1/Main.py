import sys
import numpy as np
import scipy.integrate

OFile = 'data.txt'

# alias for ode function
odeint = scipy.integrate.odeint

# Calculations Memory Constraint
MaxMem = 28*1E6

# Parameters
species = 8
dt = 1E-12 # Timestep s
trep = 1E-11 # s
tfin = 1E-9 # Final Time
xini = np.zeros(species)
xini[0] = 0.5
xini[1] = 0.5
xini[2] = 1
# Model at T=298.15 K
def model(x,t):
    dxdt = np.zeros(8)
    
    #Constants
    k00 = 2.1243982516e+11
    k01 = 3.2362038820e-39
    k02 = 2.1243982516e+11
    k03 = 3.2362038820e-39
    k04 = 2.1243982521e+11
    k05 = 2.6769180696e-24
    k06 = 2.1243982521e+11
    k07 = 6.6268592390e-26
    k08 = 2.1243982521e+11
    k09 = 6.6268592390e-26
    k10 = 2.1243982521e+11
    k11 = 2.6769180696e-24
    k12 = 2.1243982521e+11
    k13 = 5.2590657658e+09
    k14 = 2.1243982521e+11
    k15 = 5.2590657658e+09
    
    #Ratelaws
    r00 = k00*x[0]*x[2]
    r01 = k01*x[3]
    r02 = k02*x[1]*x[2]
    r03 = k03*x[4]
    r04 = k04*x[0]*x[3]
    r05 = k05*x[5]
    r06 = k06*x[1]*x[3]
    r07 = k07*x[6]
    r08 = k08*x[0]*x[4]
    r09 = k09*x[6]
    r10 = k10*x[1]*x[4]
    r11 = k11*x[7]
    r12 = k12*x[1]*x[5]
    r13 = k13*x[0]*x[6]
    r14 = k14*x[0]*x[7]
    r15 = k15*x[1]*x[6]
    
    #MassBalances
    dxdt[0] = -r00+r01-r04+r05-r08+r09+r12-r13-r14+r15
    dxdt[1] = -r02+r03-r06+r07-r10+r11-r12+r13+r14-r15
    dxdt[2] = -r00+r01-r02+r03
    dxdt[3] = +r00-r01-r04+r05-r06+r07
    dxdt[4] = +r02-r03-r08+r09-r10+r11
    dxdt[5] = +r04-r05-r12+r13
    dxdt[6] = +r06-r07+r08-r09+r12-r13+r14-r15
    dxdt[7] = +r10-r11-r14+r15
    
    return dxdt
    
def Jacobian(x,t):
    dxdt = np.zeros(shape=(8,8))
    
    #Constants
    k00 = 2.1243982516e+11
    k01 = 3.2362038820e-39
    k02 = 2.1243982516e+11
    k03 = 3.2362038820e-39
    k04 = 2.1243982521e+11
    k05 = 2.6769180696e-24
    k06 = 2.1243982521e+11
    k07 = 6.6268592390e-26
    k08 = 2.1243982521e+11
    k09 = 6.6268592390e-26
    k10 = 2.1243982521e+11
    k11 = 2.6769180696e-24
    k12 = 2.1243982521e+11
    k13 = 5.2590657658e+09
    k14 = 2.1243982521e+11
    k15 = 5.2590657658e+09
    
    #Non-zero Elements
    Jac[0,0] = -k00*x[2]-k04*x[3]-k08*x[4]-k13*x[6]-k14*x[7]
    Jac[0,1] = +k12*x[5]+k15*x[6]
    Jac[0,2] = -k00*x[0]
    Jac[0,3] = +k01-k04*x[0]
    Jac[0,4] = -k08*x[0]
    Jac[0,5] = +k05+k12*x[1]
    Jac[0,6] = +k09-k13*x[0]+k15*x[1]
    Jac[0,7] = -k14*x[0]
    Jac[1,0] = +k13*x[6]+k14*x[7]
    Jac[1,1] = -k02*x[2]-k06*x[3]-k10*x[4]-k12*x[5]-k15*x[6]
    Jac[1,2] = -k02*x[1]
    Jac[1,3] = -k06*x[1]
    Jac[1,4] = +k03-k10*x[1]
    Jac[1,5] = -k12*x[1]
    Jac[1,6] = +k07+k13*x[0]-k15*x[1]
    Jac[1,7] = +k11+k14*x[0]
    Jac[2,0] = -k00*x[2]
    Jac[2,1] = -k02*x[2]
    Jac[2,2] = -k00*x[0]-k02*x[1]
    Jac[2,3] = +k01
    Jac[2,4] = +k03
    Jac[3,0] = +k00*x[2]-k04*x[3]
    Jac[3,1] = -k06*x[3]
    Jac[3,2] = +k00*x[0]
    Jac[3,3] = -k01-k04*x[0]-k06*x[1]
    Jac[3,5] = +k05
    Jac[3,6] = +k07
    Jac[4,0] = -k08*x[4]
    Jac[4,1] = +k02*x[2]-k10*x[4]
    Jac[4,2] = +k02*x[1]
    Jac[4,4] = -k03-k08*x[0]-k10*x[1]
    Jac[4,6] = +k09
    Jac[4,7] = +k11
    Jac[5,0] = +k04*x[3]+k13*x[6]
    Jac[5,1] = -k12*x[5]
    Jac[5,3] = +k04*x[0]
    Jac[5,5] = -k05-k12*x[1]
    Jac[5,6] = +k13*x[0]
    Jac[6,0] = +k08*x[4]-k13*x[6]+k14*x[7]
    Jac[6,1] = +k06*x[3]+k12*x[5]-k15*x[6]
    Jac[6,3] = +k06*x[1]
    Jac[6,4] = +k08*x[0]
    Jac[6,5] = +k12*x[1]
    Jac[6,6] = -k07-k09-k13*x[0]-k15*x[1]
    Jac[6,7] = +k14*x[0]
    Jac[7,0] = -k14*x[7]
    Jac[7,1] = +k10*x[4]+k15*x[6]
    Jac[7,4] = +k10*x[1]
    Jac[7,6] = +k15*x[1]
    Jac[7,7] = -k11-k14*x[0]
    
    return dxdt
    
# Calculation handling a maximum matrix size
if species*tfin/dt > MaxMem:
	tfin2 = MaxMem*dt/species
	ti = 0
else:
	tfin2 = tfin
	ti = 0
t = np.arange(0,tfin2+dt,dt)
xini2 = xini
# Time indexes and Out predefinition
Out_index = []
t_old = -(trep+1.0)
for i in range(len(t)):
	if t[i] - t_old >= trep:
		Out_index.append(i)
		t_old = t[i]
Out = ['' for i in Out_index]
while ti+tfin2 <= tfin:
	print('Current interval: [{},{}] s'.format(ti,ti+tfin2))
	x = odeint(f,xini2,t,Dfun=Jacobian,rtol=1E-6,atol=1E-12)
	# Output Writing
	for j,i in enumerate(Out_index):
		Row = '\t'.join(map(str,x[i,:].tolist()))
		Out[j] = '{}\t{}'.format(ti+t[i],Row)
		t_old = t[i]
	with open('{}.txt'.format(OFile),'a') as F:
		F.write('\n'.join(Out))
		F.write('\n')
	xini2 = x[-1,:]
	ti += tfin2
#
