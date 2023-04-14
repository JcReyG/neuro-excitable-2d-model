from upDn_baseCode import *
import matplotlib 
matplotlib.use('macosx')

eCharge=1.60217733e-19 # Coulombs; 
kBoltzmann=1.38065812e-20 #mJ/K
zeroT=273.15 #deg Kelvin
TCelcius = 36
v_T = kBoltzmann * (zeroT + TCelcius)/ eCharge
C_m = 20.0; vTCm= v_T * C_m
print(r'v_T=%g, C_m=%g, v_T C_m= %g'%(v_T,C_m, vTCm))
#
voltages = {'u_w':-10.0/v_T, 'u_m': -20.0/v_T,'u_U':60.0/v_T, 'u_D':-90.0/v_T, 'u_UD':-80.0/v_T}
biases = {'b_w':0.65, 'b_U':0.5, 'b_D':0.3, 'b_UD':0.1, 'g_m':4, 'g_w':2.2, 'kappa_w':0.35}
rates = {'r_w':0.25, 'a_F': 0*100 / vTCm, 'r_U':1, 'r_D':1, 'r_UD':1e-4, 'a_U': 4, 'a_D': 6, 'a_UD':800}
numerics = {'timeMin': -0.0, 'timeMax':300.0, 'timeStep':1/40.0, 'ic': np.array([0.0001, -60.0/v_T]),\
            'uMin':-90/v_T,'uMax':40/v_T, 'wMin':0,'wMax':1, 'uStep':0.1/v_T,'wStep':0.01}
p = {'v_T': v_T, 'C_m':C_m, 'vTCm': v_T * C_m}
p= {**p, **voltages, **biases, **rates, **numerics}
#
def setBaseParameters(upDn):
    print('Resetting parameters...')
    upDn.uRange = np.arange(upDn.pars['uMin'], upDn.pars['uMax'], upDn.pars['uStep'])
    upDn.wRange = np.arange(upDn.pars['wMin'], upDn.pars['wMax'], upDn.pars['wStep'])
    upDn.pars['vMin'] = upDn.pars['uMin'] * upDn.pars['v_T']
    upDn.pars['vMax'] = upDn.pars['uMax'] * upDn.pars['v_T']
    upDn.pars['u_m']= -16.0/v_T; upDn.pars['g_m']= 4;  
    upDn.pars['u_w']= -10.0/v_T; upDn.pars['g_w']= 2.2; upDn.pars['b_w'] = 0.65; upDn.pars['kappa_w']=0.35; upDn.pars['r_w'] = 0.25
    upDn.pars['b_D']= 0.3; #upDn.pars['b_UD']= 0.8
    upDn.pars['a_F'] = 0 / vTCm; # 70 pA en rheobase 
    upDn.pars['a_U'] = 4.; upDn.pars['a_D'] = 6; upDn.pars['a_UD'] = 800
    upDn.pars['timeMax'] = 500.0; upDn.pars['timeStep']=1/60.0; 
    upDn.pars['ic'] = np.array([-80.0/upDn.pars['v_T'], 0.001])
    return upDn

upDn = UD(params= p, variables=['u','w'] )
upDn = setBaseParameters(upDn)

# --------------------------
def fullIClamp_range(ampMin=-20, ampStep=20, timeMax=1000, dvdtThresh = 80):
    upDn.pars['a_F'] = 0; upDn.pars['timeMax'] = 200.0; 
    upDn.updateFunctions()
    uOrbit, wOrbit = upDn.trayectory_Autonomous(upDn.duw_)
    vOrbit = uOrbit * upDn.pars['v_T']
    upDn.pars['ic'] = np.array([vOrbit[-1],wOrbit[-1]])
    print('Found steady state near (%g,%g)'%(vOrbit[-1],wOrbit[-1]))
    #
    upDn.pars['a_F'] = ampMin/upDn.pars['vTCm']
    upDn.pars['timeMax'] = timeMax; 
    iBase = 1e6; iTop = 0; nSI = 2
    nSI = list()
    amps = list()
    while (ampStep>0):
        amps.append(upDn.pars['a_F']*upDn.pars['vTCm'] + ampStep)
        upDn.updateFunctions()
        upDn.uOrbit, upDn.wOrbit = upDn.trayectory_Autonomous(upDn.duw_)
        upDn.dvdt = np.zeros(len(upDn.uOrbit))
        upDn.dvdt[1:] = upDn.pars['v_T'] *(upDn.uOrbit[1:]-upDn.uOrbit[:-1])/upDn.pars['timeStep'] 
        print('Max dv/dt = %g V/s'% upDn.dvdt.max())
        thr = (upDn.dvdt.max() > dvdtThresh)
        if  (thr>0):
            iBase = np.minimum([iBase, upDn.pars['a_F']*upDn.pars['vTCm']])
            si = spikeInds(upDn.dvdt,dvdtThresh) 
            nSI.append(len(si))
            iTop = np.maximum([iBase, upDn.pars['a_F']*upDn.pars['vTCm']]) 
            print('(iBase,iTop) = (%g,%g) '% (iBase,iTop))
        else : 
            nSI.append(0)
        upDn.pars['a_F'] = upDn.pars['a_F'] + ampStep/upDn.pars['vTCm']
        
        print(len(nSI), nSI)
        print(amps)

        if (len(nSI)>2) : 
            ampStep =  1 - np.int32(nSI[-1] < nSI[-2])
        
    return iBase, iTop 

# -------------------------

upDn.pars['a_D'] = 8
iBase, iTop = fullIClamp_range(ampMin=-20, ampStep=20, timeMax=1000, dvdtThresh = 80)

# -------------
# Fig 1. Two kinds of partitions of the phase plane 
# Needs fixed points by type and different slopes
# -------------
#upDn.resetParameters()
def figure1_partitionsPP(aDs = np.array([5, 15]),wMin=0,wMax=1,nWPts=100, vMin=-90,vMax=30,nVPts=200,saveFig=0):
    f = pl.figure(figsize=(11,4)); pl.ioff(); 
    ax = list(); rows = 1; cols= len(aDs);
    for n in range(rows*cols):
        upDn.pars['a_D']= aDs[n]
        ax.append(f.add_subplot(rows,cols,n+1))  
        upDn.phasePlane(ax=ax[n], W = np.linspace(wMin,wMax,nWPts), U = np.linspace(vMin,vMax,nVPts)/upDn.pars['v_T'], \
            wNullLabel=r'$\partial_t w = 0$', vNullLabel=r'$\partial_t v = 0$', plotNullClines=1)
        ax[n].set_ylim(upDn.pars['vMin'],upDn.pars['vMax'])
        ax[n].set_xlim(-0.01,0.55)
        ax[n].set_xlabel(r'$w$')
        ax[n].set_title(r'$N_{KD}$=%d'%(upDn.N_D.subs(upDn.pars)))
    ax[0].set_ylabel(r'$v$ (mV)')
    pl.subplots_adjust(left=0.1, right=0.95, bottom=0.13,  top=0.9, wspace=0.1, hspace=0.25)
    pl.ion(); pl.draw(); pl.show()
    if saveFig>0:
        figName1 = './twoTypesPP_NKD_2131_7992_labeled.jpg'
        figName2 = './twoTypesPP_NKD_2131_7992_labeled.png'
        f.savefig(figName1)
        f.savefig(figName2)
    return f


figure1= figure1_partitionsPP(aDs = np.array([5, 15]),wMin=0,wMax=1,nWPts=100, vMin=-90,vMax=30,nVPts=200,saveFig=1)

# -------------
# Fig 2. Excitability dynamics from different initial conditions, several phenotype2
# Pass
# -------------
# 
def plotProfile_v_twdvdt(ax3p, orbits,markers = ['-','--',':','-'],colors = ['blue', 'black', 'red', 'cyan']):
    print('Generating plots of the time course, phase plane, and $v$ vs $dv/dt$ with trajectories from different initial conditions')
    nIcs = len(orbits)
    for n in range(nIcs):
        alp = 0.4 + n/(4+nIcs)
        lwi = (n+4)/2
        ax3p[0].plot(orbits[n]['timeSamples'], orbits[n]['vOrbit'], markers[n], lw=lwi, color=colors[n], alpha=alp, label =r'$(t, v)$')
        upDn.phasePlane(ax=ax3p[1], W = np.linspace(0,1,50), U = np.linspace(upDn.pars['vMin'],upDn.pars['vMax'],200)/upDn.pars['v_T'], \
            wNullLabel=r'', vNullLabel=r'', plotNullClines=1)
        ax3p[1].plot(orbits[n]['wOrbit'], orbits[n]['vOrbit'], markers[n], lw=lwi, color=colors[n], alpha=alp, label =r'$(w,v)$')
        ax3p[1].set_xlim(-0.01,0.55)
        ax3p[2].plot(orbits[n]['dvdt'], orbits[n]['vOrbit'], markers[n], lw=lwi, color=colors[n], alpha=alp, label =r'$( \partial_t v, v)$')
        ax3p[0].set_xlim(-1,20)

    ax3p[0].set_ylabel(r'$v$ (mV)'); 
    ax3p[0].set_xlabel(r'$t$ (ms)')
    ax3p[1].set_xlabel(r'$w$'); 
    ax3p[2].set_xlabel(r'$\partial_t v$ (V/s)')
    for m in range(3): 
        ax3p[m].set_ylim(upDn.pars['vMin'],upDn.pars['vMax'])
        #ax3p[m].legend(loc='lower right')
    return ax3p

#
def figure2_profile_v_twdvdt(vics,aDs = np.arange(4,20,1),saveFigs = 0,markers = ['-','--',':','-'],colors = ['blue', 'black', 'red', 'cyan']):
    nIcs0 = len(vics)
    ics = np.zeros((nIcs0,2),'float64')
    ics[:,0] = vics
    ics[:,1] = upDn.pars['ic'][1]
    #
    nIcs = len(vics)
    profiles = list()
    figures = list()
    for nn in range(len(aDs)):
        upDn.pars['a_D'] = aDs[nn]
        nkd = np.int32(aDs[nn] * upDn.pars['vTCm'])
        orbits = upDn.orbitsFromICs(ics)
        f = pl.figure(figsize=(11,4)); pl.ioff(); 
        f.suptitle(r'$N_{KD}$= %d'% nkd)
        print(r'$N_{KD}$= %d'% nkd)
        ax = list(); rows = 1; cols= 3;
        for m in range(rows*cols): ax.append(f.add_subplot(rows,cols,m+1))  
        ax = plotProfile_v_twdvdt(ax, orbits,markers,colors)
        pl.subplots_adjust(left=0.075, right=0.98, bottom=0.12,  top=0.85, wspace=0.15, hspace=0.25)
        pl.ion(); pl.draw(); pl.show()
        figName = r'./ePhysPhenotypes_NKD%d.png'%nkd
        if len(figName)>0: f.savefig(figName, transparent=False)
        profiles.append(orbits); figures.append(f)
    return profiles,figures

vics = np.array([-40, -50, -60,-90])/upDn.pars['v_T']
profiles, figures = figure2_profile_v_twdvdt(vics,aDs = np.arange(4,20,1),saveFigs = 0,markers = ['-','--',':','-'],colors = ['blue', 'black', 'red', 'cyan'])

# -------------
# Fig. 3 Bifurcation diagrams for I_F, all different values for a a_K in same graph
# -------------
#
aDs = np.arange(4,30,1)
nkds = aDs[n] * upDn.pars['vTCm']
print(aDs,'\n', nkds)

naDs = len(aDs); nkds=list()
startTime = time.time()
upDn.updateFunctions()
startTime = time.time()
uStars = np.hstack([np.linspace(-90,-20,50),np.linspace(-20,0,50)])/upDn.pars['v_T']
uStars = np.linspace(-90,0,100)/upDn.pars['v_T']
wStars = upDn.w_inf_(uStars)
fps = np.vstack([uStars,wStars]).transpose()
aDs = np.arange(4,20,1)
nPheno = len(aDs)
fpsTL = upDn.cod1SecondParameterVariation(secParName='a_D', secParVals = aDs, fps= fps)
print('Took %d seconds to calculate the fixed points and their types'%(time.time()-startTime))
#
f = pl.figure(figsize=(9,5)); pl.ioff()
axBif = f.add_subplot(111)
for n  in range(naDs):
    upDn.pars['a_D'] = aDs[n]
    tStr = r'$N_{KD}$ = %d'%(nkds[n])
    axBif = upDn.bifurcationDiagram_Cod1(axBif, fpsTL[n], coordinate=0, fpScaleFactor=upDn.pars['v_T'],  \
                parScaleFactor= upDn.pars['vTCm'], xLabel=r'$I_F$ (pA)', yLabel=r'$v_{*}$ (mV)')
        
axBif.set_ylim(uStars.min()*upDn.pars['v_T'],uStars.max()*upDn.pars['v_T'])
f.suptitle(tStr)
pl.subplots_adjust(left=0.1, right=0.95, bottom=0.1,  top=0.9, wspace=0.1, hspace=0.25)
pl.ion(); pl.draw();pl.show()
figName = r'./ePhysPhenotypes_Bifurcation_NKDs%d-%d.png'%(nkds[0],nkds[-1])
if len(figName)>0: f.savefig(figName, transparent=True)
print('Took %d seconds to plot these bifurcation diagrams'%(time.time()-startTime))


# -------------
# Fig 4. Bifurcation diagrams for different phenotypes
# -------------
aDs = np.arange(4,20,5)
nPheno = len(aDs)
startTime = time.time()
upDn.updateFunctions()
startTime = time.time()
uStars = np.hstack([np.linspace(-90,-20,50),np.linspace(-20,0,50)])/upDn.pars['v_T']
uStars = np.linspace(-90,0,100)/upDn.pars['v_T']
wStars = upDn.w_inf_(uStars)
fps = np.vstack([uStars,wStars]).transpose()
aDs = np.arange(4,20,1)
nPheno = len(aDs)
fpsTL = upDn.cod1SecondParameterVariation(secParName='a_D', secParVals = aDs, fps= fps)
#
naDs = len(aDs)
for n  in range(naDs):
    upDn.pars['a_D'] = aDs[n]
    nkd = np.int32 ( aDs[n] * upDn.pars['vTCm'])
    tStr = r'$N_{KD}$ = %d'%(nkd)
    f = pl.figure(figsize=(9,5)); pl.ioff()
    axBif = f.add_subplot(111)
    axBif = upDn.bifurcationDiagram_Cod1(axBif, fpsTL[n], coordinate=0, fpScaleFactor=upDn.pars['v_T'],  \
                parScaleFactor= upDn.pars['vTCm'], xLabel=r'$I_F$ (pA)', yLabel=r'$v_{*}$ (mV)')
        
    axBif.set_ylim(uStars.min()*upDn.pars['v_T'],uStars.max()*upDn.pars['v_T'])
    f.suptitle(tStr)
    pl.subplots_adjust(left=0.1, right=0.95, bottom=0.1,  top=0.9, wspace=0.1, hspace=0.25)
    pl.ion(); pl.draw();pl.show()
    figName = r'./ePhysPhenotypes_Bifurcation_NKD%d.png'%nkd
    if len(figName)>0: f.savefig(figName, transparent=False)
print('Took %d seconds to calculate these bifurcation diagrams'%(time.time()-startTime))


# -------------
# Fig 5. I-clamp for two different profiles with bifurcation diagrams for I_F
# -------------
aDs = np.array([4,8])
nDs = np.int32 ( aDs * upDn.pars['vTCm'])
nPheno = len(aDs)
# 
upDn.pars['timeMin'] = -10
upDn.pars['timeMax'] = 200
timeStimStart=0; timeStimStop=200
#
startTime = time.time()
iClampVolts = list()
iStimList = np.array([[32,33,435,436],[80,81,2210,2215]])/upDn.pars['vTCm']
for nn in range(nPheno):
    #iStim = np.arange(ampMin,ampMax,ampStep)/upDn.pars['vTCm']
    upDn.pars['a_D'] = aDs[nn]
    print('n_D = %g'% nDs[nn])
    upDn.updateFunctions() 
    uOrbit,wOrbit = upDn.trayectory_Autonomous(upDn.duw_)
    upDn.pars['ic'] = np.array([uOrbit[-1],wOrbit[-1]])
    print('Found steady state near (%g,%g)'%(upDn.pars['ic'][0]*upDn.pars['v_T'],upDn.pars['ic'][1]))
    upDn.pars['timeStart'] = -10
    iStim = iStimList[nn]
    iList = upDn.iClampSquareStims(iLevels=iStim, timeStimStart=timeStimStart, timeStimStop=timeStimStop)
    #
    upDn.updateFunctions()
    vOrbits = upDn.iClamp(iList) 
    iClampVolts.append(vOrbits)
print('Calculated responses to I-Clamp for N_KD en ', nDs)

hShift = 0; vShift = 0; tMin = -50; tMax = upDn.pars['timeMax']
widths = [4,2,3,1]
colors = ['black','blue', 'orange', 'black']
transp = [0.4, 0.7, 1, 0.6]
for nn in range(nPheno):
    f = pl.figure(figsize=(9,5)); pl.ioff()
    axIC = f.add_subplot(111)   
    vOrbits = iClampVolts[nn] 
    nCommands = len(iStimList[nn])
    for n in range(nCommands):
        axIC.plot(upDn.timeSamples - n*hShift, vOrbits[n]*upDn.pars['v_T'] + n*vShift, color = colors[n], lw = widths[n], alpha = transp[n], label = r'$I_F=%g$ pA'%(iStimList[nn][n]*upDn.pars['vTCm']))
        axIC.legend(loc='upper right')
    axIC.set_ylabel(r'$v$ (mV)')
    tStr = r'$N_{KD}$ = %d'%(nDs[nn])
    f.suptitle(tStr)
    pl.subplots_adjust(left=0.1, right=0.95, bottom=0.1,  top=0.9, wspace=0.1, hspace=0.25)
    pl.ion(); pl.draw();pl.show()
    figName = r'./ePhysPhenotypes_IClamp_NKD%d.png'%nkd
    if len(figName)>0: f.savefig(figName, transparent=False)
    print('Took %d seconds to perform this I-clamp'%(time.time()-startTime))


