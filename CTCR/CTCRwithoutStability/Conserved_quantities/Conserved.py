import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import array
import string
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.ticker as mticker  
import sys

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


filename=sys.argv[1]
f=open(filename,"r")
rib=0.1

#f=open("s.full","r")


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

def LMoment(Mu,Q):
    m1= ( Mu[0]*Q[3] - Mu[1]*Q[2] + Mu[2]*Q[1] - Mu[3]*Q[0])/2  # Muq + Mu q + Mu + Mu 
    m2= ( Mu[0]*Q[2] + Mu[1]*Q[3] - Mu[2]*Q[0] - Mu[3]*Q[1])/2
    m3= (-Mu[0]*Q[1] + Mu[1]*Q[0] + Mu[2]*Q[3] - Mu[3]*Q[2])/2
    return([m1,m2,m3])

def Rodrigues(K,V,th):
    A=np.multiply(np.cos(th),V) 
    B=np.multiply(np.cross(K,V),np.sin(th))
    C=np.multiply(np.dot(K,V)*(1- np.cos(th)),K)
    D=[A[0]+B[0]+C[0], A[1]+B[1]+C[1], A[2]+B[2]+C[2] ]
    return(D)
 
A1=0.03487
A2=0.09422
A3=3.4
C1=A1/1.3
C2=A2/1.3
C3=A3/1.3


uhat11=2.8571
uhat12=0.0
uhat13=0.0

uhat21=1.4286
uhat22=0.0
uhat23=0.0

uhat31=0.625
uhat32=0.0
uhat33=0.0

HAMIL1=[]
HAMIL2=[]
HAMIL3=[]
INT=[]
INT1=[]
INT2=[]
INT3=[]
QN1=[]
QN2=[]
QN3=[]
INT3_1=[]
INT3_2=[]
MOM=[]
ii=0
th=2*np.pi
Nl=int(48/7)+1
TX=[]
TY=[]
TZ=[]
S=[]
LP=[]
Or=[]
QQ=[]
Lengths=[];
while(1):
    ######
    s=[]
    int1_2=[]

    int2_2=[]
    int3_2=[]
    int3_1=[]
    HAM1=[]
    HAM2=[]
    HAM3=[]
    qnorm1=[];
    qnorm2=[];
    qnorm3=[];

    A=f.readline()
    if (A==''):
        break
    #print A
    K=A.split()
    ii=ii+1
    if int(K[2])==5:
        LP.append(ii)
    for i in range(0,int(K[8])-2):
        B=f.readline()
        Bb=B.split()
        if (i%Nl==0 and i< (int(K[6])*Nl)):    
            s.append(float(Bb[0]))
            xx2=float(Bb[1])
            yy2=float(Bb[2])
            zz2=float(Bb[3])
            qq21=float(Bb[4])
            qq22=float(Bb[5])
            qq23=float(Bb[6])
        elif (i%Nl==1 and i< (int(K[6])*Nl)):        
            qq24=float(Bb[0])
            mu21=float(Bb[1])
            mu22=float(Bb[2])
            mu23=float(Bb[3])
            mu24=float(Bb[4])
            N21=float(Bb[5])
            N22=float(Bb[6])
        elif (i%Nl==2 and i< (int(K[6])*Nl)):    
            N23=float(Bb[0])
            alf=float(Bb[1])
            bet=float(Bb[2])
            xx1=float(float(Bb[3]));
            yy1=float(float(Bb[4]));
            zz1=float(float(Bb[5]));
            qq11=float(Bb[6])
            
        elif (i%Nl==3 and i< (int(K[6])*Nl)):    
            qq12=float(Bb[0])
            qq13=float(Bb[1])
            qq14=float(Bb[2])
            mu11=float(Bb[3])
            mu12=float(Bb[4])
            mu13=float(Bb[5])
            mu14=float(Bb[6])
    
        elif (i%Nl==4 and i< (int(K[6])*Nl)):    
            N11=float(Bb[0])
            N12=float(Bb[1])
            N13=float(Bb[2])
            xx3=float(Bb[3])
            yy3=float(Bb[4])
            zz3=float(Bb[5])
            qq31=float(Bb[6])
        
        
        elif (i%Nl==5 and i< (int(K[6])*Nl)):    
            qq32=float(Bb[0])
            qq33=float(Bb[1])
            qq34=float(Bb[2])
            mu31=float(Bb[3])
            mu32=float(Bb[4])
            mu33=float(Bb[5])
            mu34=float(Bb[6])
               
        elif (i%Nl==6 and i< (int(K[6])*Nl)):    
            N31=float(Bb[0])
            N32=float(Bb[1])
            N33=float(Bb[2])
            Alf2=float(Bb[3])
            Alf3=float(Bb[4])
            Bet2=float(Bb[5])
            Bet3=float(Bb[6])

        
            Tang= Dir3(qq21,qq22,qq23,qq24)
            M2=Moment([mu21,mu22,mu23,mu24],[qq21,qq22,qq23,qq24]);
            phi= (C1+C2)*bet/(C1*C2)  - (M2[2]/C1) + (uhat23-uhat13); 
            util1=(A1*uhat11 + A2*(uhat21*np.cos(alf) - uhat22*np.sin(alf) ) )/(A1+A2);
            util2=(A1*uhat12 + A2*(uhat21*np.sin(alf) + uhat22*np.cos(alf) ) )/(A1+A2);
            util3= (C1*uhat13 + C2*(uhat23 - phi))/(C1+C2);
            um = M2[0]*util1 + M2[1]*util2 + M2[2]*util3;
            mkm= M2[0]*M2[0]/(A1+A2) + M2[1]*M2[1]/(A1+A2) + M2[2]*M2[2]/(C1+C2);
            nd3= N21*Tang[0] + N22*Tang[1] + N23*Tang[2];      
            phibe= bet*phi;
            Cab= (uhat11*A1*uhat11 + uhat12*A1*uhat12 + uhat13*C1*uhat13 + uhat21*A2*uhat21 + uhat22*A2*uhat22 + (uhat23-phi)*C2*(uhat23-phi)) - (util1*(A1+A2)*util1  + util2*(A1+A2)*util2 + util3*(C1+C2)*util3);
            HAM2.append( nd3 + um + phibe + mkm/2 - Cab/2 );
    
            M1=Moment([mu11,mu12,mu13,mu14],[qq11,qq12,qq13,qq14]);
            nd3=np.dot([N11,N12,N13],Dir3(qq11,qq12,qq13,qq14))
            Hm1= (M1[0]*M1[0]/A1 + M1[1]*M1[1]/A1 + M1[2]*M1[2]/C1)/2 + M1[0]*uhat11 + M1[1]*uhat12 +  M1[2]*uhat13 + nd3;
            HAM1.append(Hm1);

            int2_2.append(mu21*qq21+ mu22*qq22 + mu23*qq23 + mu24*qq24 + 2*(xx2*N21 + yy2*N22 + zz2*N23));
            int1_2.append(mu11*qq11+ mu12*qq12 + mu13*qq13 + mu14*qq14 + 2*(xx1*N11 + yy1*N12 + zz1*N13));
            int3_2.append(mu31*qq31+ mu32*qq32 + mu33*qq33 + mu34*qq34 + 2*(xx3*N31 + yy3*N32 + zz3*N33));

            qnorm1.append(qq11*qq11+ qq12*qq12 + qq13*qq13 + qq14*qq14 -1.0)
            qnorm2.append(qq21*qq21+ qq22*qq22 + qq23*qq23 + qq24*qq24-1.0)
            qnorm3.append(qq31*qq31+ qq32*qq32 + qq33*qq33 + qq34*qq34 -1.0)
            

            #3-tube overlap:
            M3=Moment([mu31,mu32,mu33,mu34],[qq31,qq32,qq33,qq34]);
            m1=M3[0];
            m2=M3[1];
            m3=M3[2];
            Tang3= Dir3(qq31,qq32,qq33,qq34);
            NNND3= N31*Tang3[0] + N32*Tang3[1] + N33*Tang3[2]; 

            phi2= (Bet2 + Bet3)/C1 + Bet2/C2 - M3[2]/C1 + uhat23 - uhat13;
            phi3= (Bet2 + Bet3)/C1 + Bet3/C3 - M3[2]/C1 + uhat33 - uhat13;

            util31=(A1*uhat11 + A2*(uhat21*np.cos(Alf2) - uhat22*np.sin(Alf2) ) + A3*(uhat31*np.cos(Alf3) - uhat32*np.sin(Alf3) ) )/(A1 + A2 + A3);
            util32=(A1*uhat12 + A2*(uhat21*np.sin(Alf2) + uhat22*np.cos(Alf2) ) + A3*(uhat31*np.sin(Alf3) + uhat32*np.cos(Alf3) ) )/(A1 + A2 + A3);
            util33= (C1*uhat13 + C2*(uhat23 - phi2) + C3*(uhat33 - phi3))/(C1+C2+C3);
            
            UM3 = M3[0]*util31 + M3[1]*util32 + M3[2]*util33;
            MKM3= M3[0]*M3[0]/(A1+A2+A3) + M3[1]*M3[1]/(A1+A2+A3) + M3[2]*M3[2]/(C1+C2+C3);
            PHBE3= Bet3*phi3 + Bet2*phi2;

            CAB3= (A1*uhat11**2 + A1*uhat12**2 + C1*uhat13**2 + A2*uhat21**2 + A2*uhat22**2 + C2*(uhat23-phi2)**2 + A3*uhat31**2 + A3*uhat32**2 + C3*(uhat33-phi3)*2 ) - ( (A1+A2+A3)*util31**2  + (A1+A2+A3)*util32**2 + (C1+C2+C3)*util33**2 );

            Ham3=(1/(2*(A1 + A2 + A3)*C1*C2*C3))*(A1*(Bet3**2)*C1*C2 + A2*(Bet3**2)*C1*C2 + A3*(Bet3**2)*C1*C2 + A1*(Bet2**2)*C1*C3 + A2*(Bet2**2)*C1*C3 + A3*(Bet2**2)*C1*C3 + A1*(Bet2**2)*C2*C3 + A2*(Bet2**2)*C2*C3 + A3*(Bet2**2)*C2*C3 + 2*A1*Bet2*Bet3*C2*C3 +2*A2*Bet2*Bet3*C2*C3 + 2*A3*Bet2*Bet3*C2*C3 + A1*(Bet3**2)*C2*C3 + A2*(Bet3**2)*C2*C3 + A3*(Bet3**2)*C2*C3 + C1*C2*C3*m1**2 + C1*C2*C3*m2**2 - 2*A1*Bet2*C2*C3*m3 - 2*A2*Bet2*C2*C3*m3 -2*A3*Bet2*C2*C3*m3 - 2*A1*Bet3*C2*C3*m3 - 2*A2*Bet3*C2*C3*m3 - 2*A3*Bet3*C2*C3*m3 + A1*C2*C3*(m3**2) + A2*C2*C3*(m3**2) + A3*C2*C3*(m3**2) + 2*A1*C1*C2*C3*m1*uhat11 - A1*A2*C1*C2*C3*(uhat11**2) - A1*A3*C1*C2*C3*(uhat11**2) -    2*A1*Bet2*C1*C2*C3*uhat13 - 2*A2*Bet2*C1*C2*C3*uhat13 - 2*A3*Bet2*C1*C2*C3*uhat13 - 2*A1*Bet3*C1*C2*C3*uhat13 - 2*A2*Bet3*C1*C2*C3*uhat13 -2*A3*Bet3*C1*C2*C3*uhat13 +    2*A1*C1*C2*C3*m3*uhat13 + 2*A2*C1*C2*C3*m3*uhat13 + 2*A3*C1*C2*C3*m3*uhat13 -A1*A2*C1*C2*C3*(uhat21**2) - A2*A3*C1*C2*C3*(uhat21**2) + 2*A1*Bet2*C1*C2*C3*uhat23 +    2*A2*Bet2*C1*C2*C3*uhat23 + 2*A3*Bet2*C1*C2*C3*uhat23 - A1*A3*C1*C2*C3*(uhat31**2) -A2*A3*C1*C2*C3*(uhat31**2) + 2*A1*Bet3*C1*C2*C3*uhat33 + 2*A2*Bet3*C1*C2*C3*uhat33 +   2*A3*Bet3*C1*C2*C3*uhat33 + 2*A2*C1*C2*C3*(m1 + A1*uhat11)*uhat21*np.cos(Alf2) + 2*A2*A3*C1*C2*C3*uhat21*uhat31*np.cos(Alf2 - Alf3) + 2*A3*C1*C2*C3*m1*uhat31*np.cos(Alf3) + 2*A1*A3*C1*C2*C3*uhat11*uhat31*np.cos(Alf3) + 2*A2*C1*C2*C3*m2*uhat21*np.sin(Alf2) + 2*A3*C1*C2*C3*m2*uhat31*np.sin(Alf3))
            Ham3= Ham3 + np.dot(Tang3,[N31,N32,N33]);
            HAM3.append(Ham3)        
    B=f.readline();
    Bb=B.split()
    Or.append(float(Bb[5]));
    L2=float(Bb[0])
    L1=float(Bb[3])

    

    B=f.readline()
    Bb=B.split()
    L3=float(Bb[2])

    Lengths.append([L1,L2,L3])
    S.append(s);
    HAMIL1.append(HAM1)
    HAMIL2.append(HAM2);
    HAMIL3.append(HAM3)
    
    INT1.append(int1_2);
    INT2.append(int2_2);
    INT3.append(int3_2);

    QN1.append(qnorm1);
    QN2.append(qnorm2);
    QN3.append(qnorm3);

print(LP)

k=40
print("Successfully Loaded solutions")
grid = plt.GridSpec(5, 1, wspace=0.4, hspace=1.0)

plt.subplot(grid[0,0])

plt.title("Hamiltonian along the CTCR")
plt.plot(Lengths[k][0] + Lengths[k][1] + np.multiply(Lengths[k][2],S[k]),HAMIL1[k])
plt.xlabel("Arclength s")
plt.ylabel('Hamiltonian  \n 1-tube section')
plt.grid()

plt.subplot(grid[2,0])
plt.plot(Lengths[k][0] + np.multiply(Lengths[k][1],S[k]),HAMIL2[k])
plt.xlabel("Arclength s")
plt.ylabel('Hamiltonian  \n 2-tube section')
plt.grid()
plt.title(" ")

plt.subplot(grid[4,0])
plt.plot( np.multiply(Lengths[k][0],S[k]),HAMIL3[k])
plt.xlabel("Arclength s")
plt.ylabel('Hamiltonian  \n 3-tube section')
plt.grid()
plt.title(" ")
plt.show()


grid = plt.GridSpec(5, 1, wspace=0.4, hspace=1.0)
plt.subplot(grid[0,0])

plt.title("$\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n}$ along the CTCR")
plt.plot(Lengths[k][0] + Lengths[k][1] + np.multiply(Lengths[k][2],S[k]),INT1[k])
plt.xlabel("Arclength s")
plt.ylabel('$\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n}$  \n 3 tube section')
plt.grid("True")

plt.subplot(grid[2,0])
plt.plot(Lengths[k][0] + np.multiply(Lengths[k][1],S[k]),INT2[k])
plt.xlabel("Arclength s")
plt.ylabel('$\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n}$  \n 2 tube section')
plt.grid("True")

plt.subplot(grid[4,0])
plt.plot(np.multiply(Lengths[k][0],S[k]),INT3[k])
plt.xlabel("Arclength s")
plt.ylabel('$\mu \cdot \mathbf{q} + 2 \mathbf{r} \cdot \mathbf{n}$  \n 1 tube section')
plt.grid("True")
plt.show()


grid = plt.GridSpec(5, 1, wspace=0.4, hspace=1.0)
plt.subplot(grid[0,0])
plt.title("$\mathbf{q} \cdot \mathbf{q} -1 $ along the CTCR")
plt.plot(Lengths[k][0] + Lengths[k][1] + np.multiply(Lengths[k][2],S[k]),QN1[k])
plt.xlabel("Arclength")
plt.ylabel('$\mathbf{q} \cdot \mathbf{q} -1$ \n 3 tube section')
plt.grid("True")

plt.subplot(grid[2,0])
plt.plot(Lengths[k][0] + np.multiply(Lengths[k][1],S[k]),QN2[k])
plt.xlabel("Arclength")
plt.ylabel('$\mathbf{q} \cdot \mathbf{q} -1$ \n 2 tube section')
plt.grid("True")

plt.subplot(grid[4,0])
plt.plot(np.multiply(Lengths[k][0],S[k]),QN3[k])
plt.xlabel("Arclength")
plt.ylabel('$\mathbf{q} \cdot \mathbf{q} -1$ \n 1 tube section')
plt.grid("True")
plt.show()

