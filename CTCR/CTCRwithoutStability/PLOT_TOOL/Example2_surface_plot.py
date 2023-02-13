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

filename=sys.argv[1]
f=open(filename,"r")

PAR=[]
ALF1=[]
BET1=[]
ALF0=[]
Alf0=[]
Alf1=[]
Bet1=[]
Bet2=[]
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
No=48/7 + 1
oldLAB=0
LAB=0
BET1=[]
BET2=[]
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

		if (i== ((int(K[6])-1)*0)+6):
			bet1=float(Bb[5])
			Bet1.append(bet1)
			bet2=float(Bb[6])
			Bet2.append(bet2)

	if LAB < oldLAB:
		alf1=Alf1[-1];
		del Alf1[-1];
		ALF1.append(Alf1)
		ALF0.append(Alf0)
		LEN.append(le)
		Alf1=[alf1]
		Alf0=[]

		bet1=Bet1[-1];
		del Bet1[-1];
		BET1.append(Bet1)
		BET2.append(Bet2)

		Bet1=[bet1]
		Bet2=[bet2]
		

		alf1old=alf1
		bet1old=bet1
		bet2old=bet2


	
	oldLAB=LAB
	B=f.readline()

	Bb=B.split()

	Alf0.append(float(Bb[1]))
        B=f.readline()

	Bb=B.split()
	le=float(float(Bb[4]))


l=len(ALF0)
print 'Collected arrays'
print 'Check lines'
#print ALF1
############################################################################################################################################################
i=1
while(i<l):
	#print ALF0[i][0]
	if abs(ALF0[i][0])>1.0:
		del ALF0[i] 
		del ALF1[i]
		del BET1[i]
		del BET2[i]
		i=i-1
	l=len(ALF0)
	i=i+1
############################################################################################################################################################	
print 'Input Angles'
print 'Angle at the end'



print 'Finished step 1'
i=2
while(i>1):
    if ALF0[0][0]==0.0 and ALF0[0][1]==0.0:
        del ALF0[0][0] 
        del ALF1[0][0]
	del BET1[0][0]
	del BET2[0][0]
        #print '2',i
    else:
        #print 'c'
        break 

print 'Finished step 2'


#Remove the first array 
############################################################################################################################################################	
print len(LEN)	
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
	#print LE
	B=A	
	L=L+1
print 'lengths'

XX=[]
YY=[]
ZZ=[]	
ZZ2=[]	
maxlen=0
maxindex=0
print len(LEN)
for i in range(0,len(LEN)):
        if (len(ALF0[i])> maxlen):
		maxlen=len(ALF0[i])
		maxindex=i
	for j in range(len(ALF0[i])):
		XX.append(LEN[i])
		YY.append(ALF0[i][j])
		ZZ.append(BET1[i][j])
		ZZ2.append(BET2[i][j])

		
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

print len(LEN)
for i in range(0,len(LEN)):

	Ff0= interp1d(ARCLENGTH[i],ALF0[i])
	Ff1= interp1d(ARCLENGTH[i],BET1[i])	

	ALF0new=Ff0(Arclength)
	BET1new=Ff1(Arclength)
	ALF0N.append(ALF0new)
	BET1N.append(BET1new)
	idx = np.argwhere(np.diff(np.sign(ALF0new - np.pi))).flatten()
	AP=[]
	for id in idx:
		P=interp(np.pi,[ALF0new[id], ALF0new[id+1]],[Arclength[id],Arclength[id+1]])
		#print P
		AP.append(Ff1(P))
		

	API.append(AP)
    
ALF0PLT=[]
BET1PLT=[]
BET11PLT=[]
XX=[]

idx = np.argwhere(np.diff(np.sign(ALF0N[0] - np.pi))).flatten()
	

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
########################################################################################################################################################################
# Plot the surface
fig = plt.figure()
ax=plt.axes(projection='3d')

yy, zz = np.meshgrid(np.linspace(0,6.28,10),np.linspace(-0.1,0.1,40))
xx = zz*0 + 3.1415/3

ax.plot_surface(XX, ALF0PLT, BET1PLT,cmap=cm.coolwarm,alpha=0.85)
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

xx = zz*0 + 2*3.1415/3
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

xx = zz*0 + 3.1415
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

xx = zz*0 + 4*3.1415/3
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

xx = zz*0 + 5*3.1415/3
ax.plot_surface(xx, yy, zz,color='honeydew',alpha=0.25)

ax.set_xlabel("$\\alpha^{[3]}(0)$")
ax.set_zlabel("$\\beta^{[2]}(0)$ ")
ax.set_ylabel("$\\alpha^{[2]}(0)$")
xx, zz = np.meshgrid(np.linspace(0,6.28,20),np.linspace(-0.1,0.1,40))
yy= zz*0 + 3.1414159

ax.plot_surface(xx, yy, zz,color='m',alpha=0.15)
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[3.1415/3],offset=3.1415/3, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[2*3.1415/3],offset=2*3.1415/3, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[3*3.1415/3],offset=3*3.1415/3, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[4*3.1415/3],offset=4*3.1415/3, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='x', levels=[5*3.1415/3],offset=5*3.1415/3, colors=['black'],linewidths=[4])
ax.contour3D(XX, ALF0PLT,BET1PLT , zdir='y', levels=[3.0],offset=3.0, colors=['red'],linewidths=[4])

zdirs = (None, 'x', 'y', 'z', (1, 1, 0), (1, 1, 1))
ax.text(3.14/3, 4.7,0.09, "$\pi/3$", 'y')
ax.text(2*3.14/3, 4.7, 0.09, "$2\pi/3$", 'y')
ax.text(3.14, 4.7, 0.09, "$\pi$", 'y')
ax.text(4*3.14/3, 4.7,0.09, "$4\pi/3$", 'y')
ax.text(5*3.14/3, 4.7, 0.09, "$5\pi/3$", 'y')
ax.text(0.2, 3.1415, 0.09, "$\pi$", 'x')
plt.title('Bifurcation surface of a Three-tube CTCR, \n when its tubes are rotated relative to each other')
plt.show()
##########################################################################
