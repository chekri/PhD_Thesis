#This python code generates Bifurcation surfaces and works only with the solutions generated by "Example1_surface.auto" and "Example1_surface_load.auto.
# It generates errors or noisy images with other solution files
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import interp1d
from numpy import interp

filename=sys.argv[1]
f=open(filename,"r")

PAR=[]
ALF1=[]
BET1=[]
ALF0=[]
Alf0=[]
Alf1=[]
Bet1=[]
Alf10=[]
ALF10=[]
alf1=0.0
alf0=0.0
alf10=0.0
b0=0.0
b10=0.0
alf0old=0.0
alf1old=0.0
LEN=[]
le=0.2
No=int(30/7) + 1
No=int(48/7) + 1
#No=9
oldLAB=0
LAB=0
BET1=[]
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
        if (i== ((int(K[6])-1)*No)+1+1):
            alf1=float(Bb[1])
            Alf1.append(alf1)

        if (i== ((int(K[6])-1)*0)+1+1):
            bet1=float(Bb[2])
            Bet1.append(bet1)
    
    if (LAB < oldLAB):
        alf1=Alf1[-1];
        del Alf1[-1];
        ALF1.append(Alf1)
        ALF0.append(Alf0)
        LEN.append(le)
        Alf1=[alf1]
        Alf0=[]
        bet1=Bet1[-1];
        del Bet1[-1];
        BET1.append(Bet1);
        Bet1=[bet1];
        alf1old=alf1;
        bet1old=bet1;

    oldLAB=LAB
    B=f.readline()

    Bb=B.split()
    Alf0.append(float(Bb[1]))
    le=float(float(Bb[0]))
    B=f.readline()

l=len(ALF0)

############################################################################################################################################################
i=1
while(i<l):
    if abs(ALF0[i][0])>1.0:
        del ALF0[i]
        del ALF1[i]
        del BET1[i]
        i=i-1
    l=len(ALF0)
    i=i+1
############################################################################################################################################################

i=2
while(i>1):
    if ALF0[0][0]==0.0 and ALF0[0][1]==0.0:
        del ALF0[0][0]
        del ALF1[0][0]
        del BET1[0][0]
    else:
        break

#Remove the first array
############################################################################################################################################################
for j in range(1,len(ALF0)):
    if np.amin(ALF0[j]) < -np.pi :
        for o in range(0,len(ALF0[j])):
            ALF0[j][o]=2*np.pi + ALF0[j][o]

l=len(ALF0)

############################################################################################################################################################
B=0
L=0
LE=len(LEN)
while (L<LE):
    A=LEN[L]
    if A==B:
        del LEN[L]
        L=L-1
    LE=len(LEN)
    B=A
    L=L+1

XX=[]
YY=[]
ZZ=[]
maxlen=0
maxindex=0
for i in range(0,len(LEN)):
    if (len(ALF0[i])> maxlen):
        maxlen=len(ALF0[i])
        maxindex=i
    for j in range(len(ALF0[i])):
        XX.append(LEN[i])
        YY.append(ALF0[i][j])
        ZZ.append(BET1[i][j])

dist=0.0
length=0.0
Arclength=[]
Arclength.append(length)

############################################################################################################################################################

for i in range(0,maxlen-1):
    dist=np.sqrt( (ALF0[maxindex][i+1]-ALF0[maxindex][i])**2  + (BET1[maxindex][i+1]-BET1[maxindex][i])**2 )
    length=length + dist
    Arclength.append(length)

Arclength=np.array(Arclength)/length
ARCLENGTH=[]
length=0.0

for i in range(0, len(LEN)):
    ARK=[]
    length=0.0
    dist=0.0
    ARK.append(length)
    for j in range(0,len(ALF0[i])-1):
        dist=np.sqrt( (ALF0[i][j+1]-ALF0[i][j])**2  + (BET1[i][j+1]-BET1[i][j])**2 )
        length=length + dist
        ARK.append(length)
    Arklength=np.array(ARK)/length
    ARCLENGTH.append(Arklength)

ALF0new=ALF0[maxindex]
ALFINDEX = range(0, maxlen)
ALF0N=[]
BET1N=[]
API=[]


for i in range(0,len(LEN)):
    Ff0= interp1d(ARCLENGTH[i],ALF0[i])
    Ff1= interp1d(ARCLENGTH[i],BET1[i])
    #Ff3= interp1d(ARCLENGTH[i],ALF0[i])
    #Use Newton Raphson method to get FSOLVE(Ff0(Arclength)-np.pi==0)
    #Get Arclength and set of alpha 1
    #The plot plot all the branches
    #interp function

    ALF0new=Ff0(Arclength)
    BET1new=Ff1(Arclength)
    ALF0N.append(ALF0new)
    BET1N.append(BET1new)
    idx = np.argwhere(np.diff(np.sign(ALF0new - np.pi))).flatten()
    AP=[]
    for id in idx:
        P=interp(np.pi,[ALF0new[id], ALF0new[id+1]],[Arclength[id],Arclength[id+1]])
        AP.append(Ff1(P))
    API.append(AP)

ALF0PLT=[]
BET1PLT=[]
BET11PLT=[]
XX=[]
idx = np.argwhere(np.diff(np.sign(ALF0N[0] - np.pi))).flatten()


#print len(API)
Br1=[]
Br2=[]
Br3=[]
for k in range(0,len(API)):
    if len(API[k])==1:
        Br1.append(API[k][0])
        Br2.append(API[k][0])
        Br3.append(API[k][0])
    else:
        Br1.append(API[k][0])
        Br2.append(API[k][1])
        Br3.append(API[k][2])




for i in range(0,maxlen):
    Y=[]
    Z=[]
    for j in range(0,len(LEN)):
        Y.append(ALF0N[j][i])
        Z.append(BET1N[j][i])
    ALF0PLT.append(Y)
    BET1PLT.append(Z)

for i in range(0,len(LEN)):
    Z=[]
    for j in range(0,maxlen):

        Z.append(BET1N[i][j])

    BET11PLT.append(Z)


for j in range(0,maxlen):
    XX.append(LEN)

XX=np.array(XX)
ALF0PLT=np.array(ALF0PLT)
BET1PLT=np.array(BET1PLT)
BET11PLT=np.array(BET11PLT)


# Plot the surface
yy, zz = np.meshgrid(np.linspace(0,6.28,10),np.linspace(-0.11,0.11,40))
L2_proj1=0.5;#change its value if plane at different L2 is required
xx = zz*0 + L2_proj1; 

fig = plt.figure()
ax=plt.axes(projection='3d')
ax.plot_surface(XX, ALF0PLT, BET1PLT,cmap=cm.coolwarm,alpha=0.65)
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)
L2_proj2=0.7;#change its value if plane at different L2 is required
xx = zz*0 + L2_proj2 
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

L2_proj3=0.9;#change its value if plane at different L2 is required
xx = zz*0 +0.9 
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

ax.set_xlabel(' $L_{2}$')
ax.set_zlabel(' $\\beta^{[2]}(0)$ ')
ax.set_ylabel('$\\alpha^{[2]}(0)$ ')

xx, zz = np.meshgrid(np.linspace(0,2.0,20),np.linspace(-0.11,0.11,40))
alf_plane=np.pi;#change its value if plane at different alpha2_{o} is required
yy= zz*0 + alf_plane;
ax.plot_surface(xx, yy, zz,color='m',alpha=0.25)

ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[L2_proj1],offset=L2_proj1, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[L2_proj2],offset=L2_proj2, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[L2_proj3],offset=L2_proj3, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='y', levels=[alf_plane],offset=alf_plane, colors=['red'],linewidths=[4])
zdirs = (None, 'x', 'y', 'z', (1, 1, 0), (1, 1, 1))
ax.text(0.5, 4.7,0.1, "$L_{2}$=0.5", 'y')
ax.text(0.7, 4.7,0.1, "$L_{2}$=0.7", 'y')
ax.text(0.9, 4.7,0.1, "$L_{2}$=0.9", 'y')
ax.text(0.2, 3.1415, 0.1, "$\pi$", 'x')
plt.title('Bifurcation surface when $ \\alpha^{[2]}_{o}$ of a two-tube CTCR \n in terms of $L_{2}$ at a given tip load')
plt.show()

##########################################################################

