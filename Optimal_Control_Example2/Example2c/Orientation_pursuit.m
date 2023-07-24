function [Rt,Jac]=Orientation_pursuit(Ori,U,Tr,P)
B0=[0,0];
nu=0.0;
mu=1.0;


qend1=Tr.y(20,end);
qend2=Tr.y(21,end);
qend3=Tr.y(22,end);
qend4=Tr.y(23,end);

d31=2*(qend1*qend3 + qend2*qend4);
d32=2*(qend2*qend3 - qend1*qend4);
d33=-qend1*qend1 - qend2*qend2 + qend3*qend3 + qend4*qend4;

Rt= mu*dot([d31,d32,d33],Ori);
IC=Tr.y(31:48,1);
e=1e-04;
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U+diff(i,:);
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',UU);
    
    qend1=p1(end,4);
    qend2=p1(end,5);
    qend3=p1(end,6);
    qend4=p1(end,7);

    d31=2*(qend1*qend3 + qend2*qend4);
    d32=2*(qend2*qend3 - qend1*qend4);
    d33=-qend1*qend1 - qend2*qend2 + qend3*qend3 + qend4*qend4;

    tar_diff= mu*dot([d31,d32,d33],Ori);
    Jac(i)=(tar_diff-Rt)/e;
    Dbdc(i,:)=[p2(end,16)/e,p3(end,18)/e];   
end

[y13,y12,y11,t30,t20,t10]=IVP_trajectory(IC'+[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e,0,0],U);%MU1
Dbdy0(1,:)=([y12(end,16),y13(end,18)]-B0 )/e;
qend1=y11(end,4);
qend2=y11(end,5);
qend3=y11(end,6);
qend4=y11(end,7);
d31=2*(qend1*qend3 + qend2*qend4);
d32=2*(qend2*qend3 - qend1*qend4);
d33=-qend1*qend1 - qend2*qend2 + qend3*qend3 + qend4*qend4;

tar_diff= mu*dot([d31,d32,d33],Ori) ;
Jac_y0(1)=(tar_diff-Rt)/e;

[y23,y22,y21,t30,t20,t10]=IVP_trajectory(IC'+[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e],U);%Mu2
Dbdy0(2,:)=([y22(end,16),y23(end,18)]- B0)/e;

qend1=y21(end,4);
qend2=y21(end,5);
qend3=y21(end,6);
qend4=y21(end,7);
d31=2*(qend1*qend3 + qend2*qend4);
d32=2*(qend2*qend3 - qend1*qend4);
d33=-qend1*qend1 - qend2*qend2 + qend3*qend3 + qend4*qend4;

tar_diff= mu*dot([d31,d32,d33],Ori) ;
Jac_y0(2)=(tar_diff-Rt)/e;

Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));

