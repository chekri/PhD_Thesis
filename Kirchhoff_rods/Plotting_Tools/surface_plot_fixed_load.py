import numpy as np
import sys
import numpy as np
from numpy import array
import string
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.tri as mtri
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy.optimize import newton
from numpy import interp

filename=sys.argv[1]
f=open(filename,"r")


LEN=[]
le=0.2
No=14/7 + 1
Nl=No
oldLAB=0
LAB=0
lam=[]
unf=[]
UNF=[]
THT=[]
while(1):
	A=f.readline()
	if (A==''):
		break
	K=A.split()
	LAB=int(K[3])
	for i in range(0,int(K[8])-1):
		B=f.readline()
		Bb=B.split()
		if (i%Nl==1 and i< (int(K[6])*Nl)):	
			q4=(float(Bb[0]))
			if i==1:
				mu01=(float(Bb[1]))
				mu02=(float(Bb[2]))
				mu03=(float(Bb[3]))
				mu04=(float(Bb[4]))	
						
	if LAB < oldLAB:
		UNF.append(unf)
		THT.append(lam)
		unf=[]
		tht=[]
		lam=[]
		LEN.append(le)
		
	oldLAB=LAB			
	B=f.readline();
	Bb=B.split()
	th=float(Bb[4])
	lam.append(th)
	unf.append( 0.5*mu03*np.cos(th/2.0) - 0.5*mu04*np.sin(th/2.0))
	le=float(Bb[2])

l=len(THT)


#Data for plots
Y1=UNF[10]
X1=THT[10]



############################################################################################################################################################
i=1
while(i<l):

	if abs(THT[i][0])>1.0:
		del THT[i] 
		del UNF[i]
		i=i-1
	i=i+1
############################################################################################################################################################	
i=2
while(i>1):
    if THT[0][0]==0.0 and THT[0][1]==0.0:
        del THT[0][0] 
        del UNF[0][0]
    else:
        break 

############################################################################################################################################################	

for j in range(1,len(THT)):
	if np.amin(THT[j]) < -np.pi :
		for o in range(0,len(THT[j])):
			THT[j][o]=2*np.pi + THT[j][o]

l=len(THT)


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
    if (len(THT[i])> maxlen):
        maxlen=len(THT[i])
        maxindex=i
    for j in range(len(THT[i])):
        XX.append(LEN[i])
        YY.append(THT[i][j])
        ZZ.append(UNF[i][j])
		
dist=0.0
length=0.0
Arclength=[]
Arclength.append(length)

############################################################################################################################################################	
for i in range(0,maxlen-1):
	dist=np.sqrt( (THT[maxindex][i+1]-THT[maxindex][i])**2  + (UNF[maxindex][i+1]-UNF[maxindex][i])**2 )
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
	for j in range(0,len(THT[i])-1):
		dist=np.sqrt( (THT[i][j+1]-THT[i][j])**2  + (UNF[i][j+1]-UNF[i][j])**2 )
		length=length + dist
		ARK.append(length)
	Arklength=np.array(ARK)/length	
	ARCLENGTH.append(Arklength)
	
THTnew=THT[maxindex]
ALFINDEX = range(0, maxlen)
THTN=[]
UNFN=[]
API=[]

for i in range(0,len(LEN)):
	Ff0= interp1d(ARCLENGTH[i],THT[i])
	Ff1= interp1d(ARCLENGTH[i],UNF[i])	
	#Ff3= interp1d(ARCLENGTH[i],THT[i])
	#Use Newton Raphson method to get FSOLVE(Ff0(Arclength)-np.pi==0)
	#Get Arclength and set of alpha 1
	#The plot plot all the branches
	#interp function
	
	THTnew=Ff0(Arclength)
	UNFnew=Ff1(Arclength)
	THTN.append(THTnew)
	UNFN.append(UNFnew)
	idx = np.argwhere(np.diff(np.sign(THTnew - np.pi))).flatten()
	AP=[]
	for id in idx:
		P=interp(np.pi,[THTnew[id], THTnew[id+1]],[Arclength[id],Arclength[id+1]])
		AP.append(Ff1(P))
		
	API.append(AP)
    
THTPLT=[]
UNFPLT=[]
UNF1PLT=[]
XX=[]


idx = np.argwhere(np.diff(np.sign(THTN[0] - np.pi))).flatten()


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
		Y.append(THTN[j][i])
		Z.append(UNFN[j][i])
	THTPLT.append(Y)
	UNFPLT.append(Z)

for i in range(0,len(LEN)):
	Z=[]	
	for j in range(0,maxlen):
		Z.append(UNFN[i][j])

	UNF1PLT.append(Z)


for j in range(0,maxlen):
	XX.append(LEN)

XX=np.array(XX)
THTPLT=np.array(THTPLT)
UNFPLT=np.array(UNFPLT)
UNF1PLT=np.array(UNF1PLT)

# Plot the surface
fig = plt.figure()
ax=plt.axes(projection='3d')
yy, zz = np.meshgrid(np.linspace(0,6.282,10),np.linspace(-2.0,2.0,40))
xx = zz*0 + 0.8
ax.plot_surface(XX, THTPLT, UNFPLT,color='m',alpha=0.5)
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)
xx = zz*0 + 1.6
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)
#ax.plot_surface(XX, THTPLT, UNFPLT,color='m',alpha=0.5)

xx = zz*0 + 2.5
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)
#ax.plot_surface(XX, THTPLT, UNFPLT,color='m',alpha=0.5)


ax.set_zlim(-1.8, 1.8)

ax.set_xlabel('Length $l$')
ax.set_zlabel('$m_{3}(0)$')
ax.set_ylabel('Rotation of the clamped end $\Theta$')
xx, zz = np.meshgrid(np.linspace(0,4.0,20),np.linspace(-2.0,2.0,40))
yy= zz*0 + 3.1414159

ax.plot_surface(xx, yy, zz,color='b',alpha=0.25)
ax.contour3D(XX, THTPLT,UNFPLT , zdir='x', levels=[0.8],offset=0.8, colors=['black'],linewidths=[2])
ax.contour3D(XX, THTPLT,UNFPLT , zdir='x', levels=[1.6],offset=1.6, colors=['black'],linewidths=[2])
ax.contour3D(XX, THTPLT,UNFPLT , zdir='x', levels=[2.5],offset=2.5, colors=['black'],linewidths=[2])
ax.contour3D(XX, THTPLT,UNFPLT , zdir='y', levels=[3.14159],offset=3.14159, colors=['red'],linewidths=[2])
zdirs = (None, 'x', 'y', 'z', (1, 1, 0), (1, 1, 1))
ax.text(0.8, 4.7,1.5, "$L_{2}$=0.8", 'y')
ax.text(1.6, 4.7, 1.5, "$L_{2}$=1.6", 'y')
ax.text(2.5, 4.7, 1.5, "$L_{2}$=2.5", 'y')
ax.text(0.2, 3.1415, 0.75, "$\Theta=\pi$", 'x')
plt.show()
# Plot the surface
fig = plt.figure()
ax=plt.axes(projection='3d')
ax.contour3D(XX, THTPLT,UNFPLT , zdir='y', levels=[3.14159],offset=3.14159, colors=['red'],linewidths=[4])
ax.plot_surface(XX, THTPLT, UNFPLT,cmap=cm.coolwarm,alpha=0.850)
ax.plot_wireframe(XX, THTPLT, UNFPLT,rstride=2, cstride=2, color='grey',alpha=0.25,linewidths=[1] )


xx, zz = np.meshgrid(np.linspace(0,4.0,20),np.linspace(-2.0,2.0,40))
yy= zz*0 + 3.1414159
ax.plot_surface(xx, yy, zz,color='m',alpha=0.15)

ax.contour3D(XX, THTPLT,UNFPLT , zdir='x', levels=[0.8],offset=0.8, colors=['blue'],linewidths=[3])
ax.contour3D(XX, THTPLT,UNFPLT , zdir='x', levels=[1.6],offset=1.6, colors=['blue'],linewidths=[3])
ax.contour3D(XX, THTPLT,UNFPLT , zdir='x', levels=[2.5],offset=2.5, colors=['blue'],linewidths=[3])
ax.contour3D(XX, THTPLT,UNFPLT , zdir='y', levels=[3.14159],offset=3.14159, colors=['red'],linewidths=[2])

yy, zz = np.meshgrid(np.linspace(0,6.282,10),np.linspace(-2.0,2.0,40))
xx = zz*0 + 0.8

xx = zz*0 + 0.8
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)
xx = zz*0 + 1.6
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)
xx = zz*0 + 2.5
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

ax.set_xlabel('Length $l$')
ax.set_zlabel('$m_{3}(0)$')
ax.set_ylabel('Rotation of the clamped end $\Theta$')
ax.text(0.8, 4.7,1.5, "$L_{2}$=0.8", 'y')
ax.text(1.6, 4.7, 1.5, "$L_{2}$=1.6", 'y')
ax.text(2.5, 4.7, 1.5, "$L_{2}$=2.5", 'y')
ax.text(0.2, 3.1415, 0.75, "$\Theta=\pi$", 'x')
plt.show()

##########################################################################
