import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np
from numpy import array
import string
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


filename=sys.argv[1]
f=open(filename,"r")
PQ1=[]
PQ2=[]
PQ3=[]
SOLL=[]
SOL=[]
DDX=[]
QQ=[]
Lengths=[]
    
ii=0
BR=1
dX1=[]
dY1=[]
dZ1=[]
Pq1=[]
Pq2=[]
Pq3=[]
BP=[]
LP=[]
DX=[]
DX1=[]
DY1=[]
DY2=[]
DY3=[]
Q=[]
DDY1=[]
DDY2=[]
DDY3=[]
DDDY1=[]
DDDY2=[]
DDDY3=[]

SOLL=[]

No=198/7 + 1
while(1):
    dalf=[]
    Dt=[]
    A=f.readline()
    #print A
    if A=='': #and f.tell() == endLocation:
        break
    K=A.split()
    #print K;
    if K==[]:
        break
    #print int(K[8]);
    Sol=[]
    #print K[8]
    for i in range(0,int(K[8])-2):
        #print i
        B=f.readline()
        Bb=B.split()
        #print B
        Bb=np.array(Bb)
        Bb=Bb.astype(np.float)
        dy1=[]
        dy2=[]
        dy3=[]
        if ((i%No ==0) and (i != 0) and  i <= (int(K[6])*No)):
            s=Sol[0]
            q1=Sol[4]
            q2=Sol[5]
            q3=Sol[6]
            q4=Sol[7]

            qq1=Sol[20]
            qq2=Sol[21]
            qq3=Sol[22]
            qq4=Sol[23]
            
            q31=Sol[34]
            q32=Sol[35]
            q33=Sol[36]
            q34=Sol[37]
            nor=q1*q1 + q2*q2 + q3*q3 + q4*q4
            nor1=qq1*qq1 + qq2*qq2 + qq3*qq3 + qq4*qq4
            nor3=q31*q31 + q32*q32 + q33*q33 + q34*q34
            q=[s,q1,q2,q3,q4,qq1,qq2,qq3,qq4,q31,q32,q33,q34]
            #print s
            for m in range(0,5):
                for i2 in range(1,5):
                    dy2.append(Sol[48 + 30*m +i2])
                    dy1.append(Sol[48 + 30*m + 10 + i2])
                    dy3.append(Sol[48 + 30*m + 18 + i2])
                dy2.append(Sol[48 + 30*m + 9])
                dy3.append(Sol[48 + 30*m + 27])
                dy3.append(Sol[48 + 30*m + 28])
                DY2.append(dy2)
                DY1.append(dy1)
                DY3.append(dy3)
                dy1=[]
                dy2=[] 
                dy3=[]
            DDY3.append(DY3)
            del DY2[-1]
            DDY2.append(DY2)
            del DY1[-1]
            del DY1[-1]
            DDY1.append(DY1)
            Q.append(q)    
            SOL.append(Sol)
            Sol=[]
            DX=[]
            DY1=[]
            DY2=[]
            DY3=[]
        Sol=np.concatenate([Sol,Bb])
    B=f.readline()
    #print "B"

    B=B.split()
  
    L2=float(B[0])
    L1=float(B[3])
    
    B=f.readline()
    B=B.split()
    L3=float(B[2])
    Lengths.append([L1,L2,L3])
    DX1.append(L3) 
    DX1.append(L2)  #Double overlap
    DX1.append(L1)  #Single overlap
    DDDY2.append(DDY2)
    DDDY1.append(DDY1)
    DDDY3.append(DDY3)
    
    QQ.append(Q)                               
    SOLL.append(SOL)

    DDY1=[]
    DDY2=[]
    DDY3=[]
    Q=[]                                               
print "Successfully Loaded solutions";
#print Pq1

DT3=[]
DT2=[]
DT1=[]
DET2=[]
DET1=[]
DET3=[]
St=[]
St1=[]
QQ[0][1];
print len(QQ)
print len(QQ[0])
for i in range(0,len(DDDY1)):
    for j in range(0,len(DDDY1[i])):
###############################################################        
        q1=QQ[i][j][1]
        q2=QQ[i][j][2]
        q3=QQ[i][j][3]
        q4=QQ[i][j][4]
        Pr=[[q4,q3,-q2,-q1,0],[-q3,q4,q1,-q2,0],[q2,-q1,q4,-q3,0],[0,0,0,0,1.0]] 
        PM=np.dot(DDDY2[i][j],np.transpose(Pr))
        #print PM
        Det2=np.linalg.det(PM)
###############################################################
        qq1=QQ[i][j][5]
        qq2=QQ[i][j][6]
        qq3=QQ[i][j][7]
        qq4=QQ[i][j][8]
        nor1=qq1*qq1 + qq2*qq2 + qq3*qq3 + qq4*qq4
        
        Pr1=[[qq4,qq3,-qq2,-qq1],[-qq3,qq4,qq1,-qq2],[qq2,-qq1,qq4,-qq3]] 
        PM1=np.dot(DDDY1[i][j],np.transpose(Pr1))
        #print PM1
        Det1=np.linalg.det(PM1)
##############################################################
        q31=QQ[i][j][9]
        q32=QQ[i][j][10]
        q33=QQ[i][j][11]
        q34=QQ[i][j][12]
        nor3=q31*q31 + q32*q32 + q33*q33 + q34*q34        
        Pr3=[[q34,q33,-q32,-q31,0,0],[-q33,q34,q31,-q32,0,0],[q32,-q31,q34,-q33,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]] 
        PM3=np.dot(DDDY3[i][j],np.transpose(Pr3))
        Det3=np.linalg.det(PM3)
        DET3.append([Lengths[i][2]*QQ[i][j][0],Det3])
        DET2.append([Lengths[i][2] + Lengths[i][1]*QQ[i][j][0],Det2])
        DET1.append([Lengths[i][2] + Lengths[i][1] + (Lengths[i][0]*QQ[i][j][0]),Det1])
    DT2.append(DET2)
    DT1.append(DET1)
    DT3.append(DET3)
    DET2=[]
    DET1=[]
    DET3=[]
print "Successful Computation!!!"


Datax=[]
Datay=[]
Dataxx=[]
Datayy=[]
Dataxxx=[]
Datayyy=[]

for k in range(0,len(DT1),1):
    plotx=[]
    ploty=[]
    plotxx=[]
    plotyy=[]
    plotxxx=[]
    plotyyy=[]
    for i in range(0,len(DT1[k])):
        plotx.append(DT2[k][i][0])
        ploty.append(DT2[k][i][1])
        plotxx.append(DT1[k][i][0])
        plotyy.append(DT1[k][i][1])
        plotxxx.append(DT3[k][i][0])
        plotyyy.append(DT3[k][i][1])
    Datax.append(plotx)
    Datay.append(ploty)
    Dataxx.append(plotxx)
    Datayy.append(plotyy)
    Dataxxx.append(plotxxx)
    Datayyy.append(plotyyy)
    plt.plot(plotx,ploty,'r')
    plt.plot(plotxx,plotyy,'g')
    plt.plot(plotxxx,plotyyy,'b')
plt.grid()
plt.show()


#Choose a solution
k=24
plotx=[]
ploty=[]
plotxx=[]
plotyy=[]
plotxxx=[]
plotyyy=[]
for i in range(0,len(DT2[k])):
    plotx.append(DT2[k][i][0])
    ploty.append(DT2[k][i][1])
    plotxx.append(DT1[k][i][0])
    plotyy.append(DT1[k][i][1])
    plotxxx.append(DT3[k][i][0])
    plotyyy.append(DT3[k][i][1])
plt.xlabel('Length of the robot (X100 mm)')
plt.ylabel('Determinant of the stability matrix')
plt.plot(plotx,ploty,'r')
plt.plot(plotxx,plotyy,'g')
plt.plot(plotxxx,plotyyy,'b')
plt.grid()
plt.show()




