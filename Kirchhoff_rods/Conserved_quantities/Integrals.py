import numpy as np
import matplotlib.pyplot as plt
import sys


def Dir1(q1,q2,q3,q4):
	return([q1**2 - q2**2 - q3**2 + q4**2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4)])

def Dir2(q1,q2,q3,q4):
	return([2*(q1*q2 - q3*q4),-q1**2 + q2**2 - q3**2 + q4**2, 2*(q2*q3 + q1*q4)])

def Dir3(q1,q2,q3,q4):
	return([2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), -q1**2 - q2**2 + q3**2 + q4**2])
	
def nor(q1,q2,q3,q4):
	return(q1**2 + q2**2 + q3**2 + q4**2)
	
def Moment(Mu,Q):
	m1= ( Mu[0]*Q[3] + Mu[1]*Q[2] - Mu[2]*Q[1] - Mu[3]*Q[0])/2  # Muq + Mu q + Mu + Mu 
	m2= (-Mu[0]*Q[2] + Mu[1]*Q[3] + Mu[2]*Q[0] - Mu[3]*Q[1])/2
	m3= ( Mu[0]*Q[1] - Mu[1]*Q[0] + Mu[2]*Q[3] - Mu[3]*Q[2])/2
	return([m1,m2,m3])


il=1
int1=[]
Ham=[]
norm=[]
filename=sys.argv[1]
f=open(filename,"r")
INT1=[]
HAM=[]
NOR=[]
S=[]
LP=[]
ii=0
th=2*np.pi
Nl=14/7+1
K1=1.0
K3=1.0/1.3
uhat1=1.5
uhat2=0.0
uhat3=0.0
while(1):
	s=[]
	x1=[]
	y1=[]
	z1=[]
	int1=[]
	A=f.readline()
	if (A==''):
		break
	K=A.split()
	ii=ii+1
	if int(K[2])==5:
		LP.append(ii)
	for i in range(0,int(K[8])-2):
		B=f.readline()
		Bb=B.split()
		if (i%Nl==0 and i< (int(K[6])*Nl)):	
			s.append(float(Bb[0]))
			x1=(float(Bb[1]))
			y1=(float(Bb[2]))
			z1=(float(Bb[3]))
			q1=float(Bb[4]);
			q2=float(Bb[5]);
			q3=float(Bb[6]);
		if (i%Nl==1 and i< (int(K[6])*Nl)):	
			q4=float(Bb[0]);
			mu1=float(Bb[1]);
			mu2=float(Bb[2]);
			mu3=float(Bb[3]);
			mu4=float(Bb[4]);
			n1=float(Bb[5]);
			n2=float(Bb[6]);

		if (i%Nl==2 and i< (int(K[6])*Nl)):	
			n3=float(Bb[0]);
			d3=Dir3(q1,q2,q3,q4)
			M=Moment([mu1,mu2,mu3,mu4],[q1,q2,q3,q4]);
			Hm =  M[0]*M[0]/(2*K1) + M[1]*M[1]/(2*K1) + M[2]*M[2]/(2*K3)  + M[0]*uhat1 + M[1]*uhat2 + M[2]*uhat3  + n1*d3[0] + n2*d3[1] + n3*d3[2]
			int1.append(mu1*q1 + mu2*q2 + mu3*q3 + mu4*q4 + 2.0*(x1*n1 + y1*n2 + z1*n3))
			norm.append(q1*q1 + q2*q2 + q3*q3 + q4*q4 -1)
			Ham.append(Hm)

	B=f.readline();
	Bb=B.split()
	B=f.readline();
	INT1.append(int1)
	HAM.append(Ham)
	NOR.append(norm)
	S.append(s)
	int1=[]
	s=[]
	Ham=[]
	norm=[]



#Give a label of the equilibrum from the s.filename to plot the corresponding conserved quantities
k=23
plt.plot(S[k],INT1[k])
plt.xlabel("Arclength of the rod s normalized in [0,1]")
plt.ylabel("$\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n}$")
plt.xlim([0,1])
plt.grid()
plt.show()

plt.plot(S[k],HAM[k])
plt.xlabel("Arclength of the rod s normalized in [0,1]")
plt.ylabel("Hamiltonian H")
plt.xlim([0,1])
plt.grid()
plt.show()

plt.plot(S[k],NOR[k])
plt.xlabel("Arclength of the rod s normalized in [0,1]")
plt.ylabel("$\mathbf{q} \cdot \mathbf{q} -1$")
plt.xlim([0,1])
plt.grid()
plt.show()

