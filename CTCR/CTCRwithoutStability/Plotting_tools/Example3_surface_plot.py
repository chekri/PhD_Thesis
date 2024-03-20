import sys
import numpy as np
from numpy import array
import string
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.tri as mtri
import matplotlib.tri as mtri
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy.optimize import newton
from numpy import interp

fiLDame=sys.argv[1]
f=open(fiLDame,"r")

PAR=[]

alf1=0.0
alf0=0.0
alf10=0.0
b0=0.0
b10=0.0
alf0old=0.0
alf1old=0.0
LD=[]
le=0.2
No=48/7 + 1
Nl=No
oldLAB=0
LAB=0
lam=[]
bif=[]
BIF=[]
THT=[]
while(1):
    bet1=[]
    par=[]
    A=f.readline()
    if (A==''):
        break
    K=A.split()
    LAB=int(K[3])
    K1= int(A[41:44]) 
    K2 =int(A[54:58])  
    for i in range(0,int(K[8])-2):
        B=f.readline()
        Bb=B.split()
        if (i%Nl==5 and i< (int(K[6])*Nl)):    
            q4=(float(Bb[0]))
            if i==5:
                mu01=(float(Bb[3]))
                mu02=(float(Bb[4]))
                mu03=(float(Bb[5]))
                mu04=(float(Bb[6]))    
    if LAB < oldLAB:
        BIF.append(bif)
        THT.append(lam)
        bif=[]
        tht=[]
        lam=[]
        LD.append(ld)


    
    oldLAB=LAB
    B=f.readline()
    Bb=B.split()
    ld=float(Bb[6])
    B=f.readline()
    Bb=B.split()
    
    th=float(Bb[1])
    lam.append(th)
    bif.append( 0.5*mu03*np.cos(th/2.0) - 0.5*mu04*np.sin(th/2.0))
    


l=len(THT)

############################################################################################################################################################
i=1
while(i<l):
    if abs(THT[i][0])>1.0:
        del THT[i] 
        del BIF[i]
        i=i-1
    l=len(THT)
    i=i+1
############################################################################################################################################################    

i=2
while(i>1):
    if THT[0][0]==0.0 and THT[0][1]==0.0:
        del THT[0][0] 
        del BIF[0][0]
    else:
        break 

#Remove the first array 
############################################################################################################################################################        
for j in range(1,len(THT)):
    if np.amin(THT[j]) < -np.pi :
        for o in range(0,len(THT[j])):
            THT[j][o]=2*np.pi + THT[j][o]

l=len(THT)

############################################################################################################################################################    
B=0
L=0
LE=len(LD)
while (L<LE):
    A=LD[L]
    if A==B:
        del LD[L]
        L=L-1
    LE=len(LD)
    B=A    
    L=L+1

XX=[]
YY=[]
ZZ=[]    
maxlen=0
maxindex=0


for i in range(0,len(LD)):
    if (len(THT[i])> maxlen):
        maxlen=len(THT[i])
        maxindex=i
    for j in range(len(THT[i])):
        XX.append(LD[i])
        YY.append(THT[i][j])
        ZZ.append(BIF[i][j])
        
dist=0.0
length=0.0
Arclength=[]
Arclength.append(length)

############################################################################################################################################################    
for i in range(0,maxlen-1):
    dist=np.sqrt( (THT[maxindex][i+1]-THT[maxindex][i])**2  + (BIF[maxindex][i+1]-BIF[maxindex][i])**2 )
    length=length + dist
    Arclength.append(length)

Arclength=np.array(Arclength)/length
ARCLENGTH=[]
length=0.0

for i in range(0, len(LD)):
    ARK=[]
    length=0.0
    dist=0.0
    ARK.append(length)
    for j in range(0,len(THT[i])-1):
        dist=np.sqrt( (THT[i][j+1]-THT[i][j])**2  + (BIF[i][j+1]-BIF[i][j])**2 )
        length=length + dist
        ARK.append(length)
    Arklength=np.array(ARK)/length    
    ARCLENGTH.append(Arklength)

    
THTnew=THT[maxindex]

THTN=[]
BIFN=[]
API=[]

for i in range(0,len(LD)):
    Ff0= interp1d(ARCLENGTH[i],THT[i])
    Ff1= interp1d(ARCLENGTH[i],BIF[i])    

    THTnew=Ff0(Arclength)
    UNFnew=Ff1(Arclength)
    THTN.append(THTnew)
    BIFN.append(UNFnew)
    idx = np.argwhere(np.diff(np.sign(THTnew - np.pi))).flatten()
    AP=[]
    for id in idx:
        P=interp(np.pi,[THTnew[id], THTnew[id+1]],[Arclength[id],Arclength[id+1]])
        AP.append(Ff1(P))
        
    API.append(AP)
    
THTPLT=[]
BIFPLT=[]
BIF1PLT=[]
XX=[]

idx = np.argwhere(np.diff(np.sign(THTN[0] - np.pi))).flatten()

#############################################################################################################
for i in range(0,maxlen):
    Y=[]
    Z=[]    
    for j in range(0,len(LD)):
        Y.append(THTN[j][i])
        Z.append(BIFN[j][i])
    THTPLT.append(Y)
    BIFPLT.append(Z)

for i in range(0,len(LD)):
    Z=[]
    for j in range(0,maxlen):
        Z.append(BIFN[i][j])
    BIF1PLT.append(Z)


for j in range(0,maxlen):
    XX.append(LD)

XX=np.array(XX)

THTPLT=np.array(THTPLT)
BIFPLT=np.array(BIFPLT)
BIF1PLT=np.array(BIF1PLT)

########################################################################################################################################################################
fig = plt.figure()
ax=plt.axes(projection='3d')

yy, zz = np.meshgrid(np.linspace(-0.6,6.582,10),np.linspace(-0.15,0.15,40))
ax.plot_surface(XX, THTPLT, BIFPLT,cmap=cm.coolwarm,alpha=0.65)

F1=0.02#vary this value if plane at different F2 is required
xx = zz*0 + F1
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

F2=0.08#vary this value if plane at different F2 is required
xx = zz*0 + F2
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

F3=0.12#vary this value if plane at different F2 is required
xx = zz*0 + F3
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

ax.set_xlabel('$F_{2}$')
ax.set_zlabel(' $m_{3}(0)$ ')
ax.set_ylabel('$\Theta$')

ax.contour3D(XX, THTPLT, BIFPLT , zdir='x', levels=[F1],offset=F1, colors=['black'],linewidths=[2])
ax.contour3D(XX, THTPLT, BIFPLT , zdir='x', levels=[F2],offset=F2, colors=['black'],linewidths=[2])
ax.contour3D(XX, THTPLT, BIFPLT , zdir='x', levels=[F3],offset=F3, colors=['black'],linewidths=[2])
zdirs = (None, 'x', 'y', 'z', (1, 1, 0), (1, 1, 1))
ax.text(0.01, 4.7,0.15, "$F_{2}$=0.02", 'y')
ax.text(0.08, 4.7,0.15, "$F_{2}$=0.08", 'y')
ax.text(0.12, 4.7,0.15, "$F_{2}$=0.12", 'y')
ax.text(0.2, 3.1415, 6.6, "$\pi$", 'x')
plt.title('Bifurcation Surface')
plt.show()

plt.plot(THT[1],BIF[1],label='$F_{2}=0.02$')
plt.plot(THT[11],BIF[11],label='$F_{2}=0.08$')
plt.plot(THT[26],BIF[26],label='$F_{2}=0.12$')
plt.ylabel(' $m_{3}(0)$ ')
plt.xlabel('$\Theta$')
plt.title('Distinguished Bifurcation diagrams for different values of $F_{2}$')
plt.legend()
plt.grid()
plt.show()
##########################################################################
