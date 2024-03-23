function [NOR] = Integrals()
%This function computes the integrals of the CTCR system
%   Detailed explanation goes here
u=[rand(1),rand(1),rand(1),rand(1),rand(1),rand(1)]
%Load=[-0.2,0.4,-0.4];
Load=[-1+2*rand(1),-1+2*rand(1),-1+2*rand(1)]
     A1=1.0;
     A2=1.2;
     A3=1.4;
     C1=A1/1.3;
     C2=A2/1.3;
     C3=A3/1.3;
Uh1=[0.5,0,0];
Uh2=[0.8,0,0];
Uh3=[1.0,0,0];

uhat21=Uh2(1);
uhat22=Uh2(2);
uhat23=Uh2(3);
uhat11=Uh1(1);
uhat12=Uh1(2);
uhat13=Uh1(3);
uhat31=Uh3(1);
uhat32=Uh3(2);
uhat33=Uh3(3);

Trj=Trajectory(u,Load,[]);
Q2=[Trj.y(4,:);Trj.y(5,:);Trj.y(6,:);Trj.y(7,:)];
MU2=[Trj.y(8,:);Trj.y(9,:);Trj.y(10,:);Trj.y(11,:)];
R2=[Trj.y(1,:);Trj.y(2,:);Trj.y(3,:)];
N2=[Trj.y(12,:);Trj.y(13,:);Trj.y(14,:)];

S=Trj.x;
NOR2=[];
INT2=[];
Ham2=[];
for k=1:size(Q2,2)
    NOR2(k)=dot(Q2(:,k),Q2(:,k));
    INT2(k)=dot(MU2(:,k),Q2(:,k)) + 2*dot(R2(:,k),N2(:,k)); 
    [m1,m2,m3]=Local_Moment(MU2(:,k),Q2(:,k));
    [D31,D32,D33]=director_D3(Q2(:,k));
    bet=Trj.y(16,k);
    alf=Trj.y(15,k);
     phi= (C1+C2)*bet/(C1*C2)  - (m3/C1) + (uhat23-uhat13) ;
     util1=(A1*uhat11 + A2*(uhat21*cos(alf) - uhat22*sin(alf) ) )/(A1+A2);
     util2=(A1*uhat12 + A2*(uhat21*sin(alf) + uhat22*cos(alf) ) )/(A1+A2);
     util3= (C1*uhat13 + C2*(uhat23 - phi))/(C1+C2);
         
     um = m1*util1 + m2*util2 + m3*util3;
     mkm = (m1*m1 + m2*m2)/(A1+A2) + m3*m3/(C1+C2);
     nd3= dot(N2(:,k),[D31,D32,D33]); 
     phibe= bet*phi;
     Cab= (uhat11*A1*uhat11 + uhat12*A1*uhat12 + uhat13*C1*uhat13 + uhat21*A2*uhat21 + uhat22*A2*uhat22 + (uhat23-phi)*C2*(uhat23-phi)) - (util1*(A1+A2)*util1  + util2*(A1+A2)*util2 + util3*(C1+C2)*util3);
            
            
     Ham2(k)= nd3 + um + phibe + mkm/2 - Cab/2 ;
    
    
end


I=16;
Q1=[Trj.y(I+4,:);Trj.y(I+5,:);Trj.y(I+6,:);Trj.y(I+7,:)];
MU1=[Trj.y(I+8,:);Trj.y(I+9,:);Trj.y(I+10,:);Trj.y(I+11,:)];
R1=[Trj.y(I+1,:);Trj.y(I+2,:);Trj.y(I+3,:)];
N1=[Trj.y(I+12,:);Trj.y(I+13,:);Trj.y(I+14,:)];
NOR1=[];
INT1=[];
Ham1=[];
for k=1:size(Q1,2)
    NOR1(k)=norm(Q1(:,k));
    INT1(k)=dot(MU1(:,k),Q1(:,k)) + 2*dot(R1(:,k),N1(:,k)); 
     [m1,m2,m3]=Local_Moment(MU1(:,k),Q1(:,k));
     [D31,D32,D33]=director_D3(Q1(:,k));     
    
     Ham1(k)=((m1*m1 + m2*m2)/A1 + (m3*m3)/C1)/2 + dot([m1,m2,m3],Uh1) + dot(N1(:,k),[D31,D32,D33]);
    
end

I=30;
Q3=[Trj.y(I+4,:);Trj.y(I+5,:);Trj.y(I+6,:);Trj.y(I+7,:)];
MU3=[Trj.y(I+8,:);Trj.y(I+9,:);Trj.y(I+10,:);Trj.y(I+11,:)];
R3=[Trj.y(I+1,:);Trj.y(I+2,:);Trj.y(I+3,:)];
N3=[Trj.y(I+12,:);Trj.y(I+13,:);Trj.y(I+14,:)];
NOR3=[];
INT3=[];
Ham3=[];
INT_MOM1=[];
INT_MOM2=[];
INT_MOM3=[];
for k=1:size(Q3,2)
    NOR3(k)=norm(Q3(:,k));
    INT3(k)=dot(MU3(:,k),Q3(:,k)) + 2*dot(R3(:,k),N3(:,k)); 
    
    [m1,m2,m3]=Local_Moment(MU3(:,k),Q3(:,k));
    [M1,M2,M3]=Global_Moment(MU3(:,k),Q3(:,k));
    INT_MOM=[M1,M2,M3]+ cross(R3(:,k)',N3(:,k)');
    INT_MOM1(k)=INT_MOM(1);
    INT_MOM2(k)=INT_MOM(2);
    INT_MOM3(k)=INT_MOM(3);
    [D31,D32,D33]=director_D3(Q3(:,k));

    alf2=Trj.y(I+15,k);  %alf2
    alf3=Trj.y(I+17,k);  %alf2
    bet2=Trj.y(I+16,k);  %alf2
    bet3=Trj.y(I+18,k);  %alf2
    
     ab1=-C2*((bet2 + bet3)/C1 - m3/C1  + bet2/C2 + uhat23- uhat13 - uhat23); 
     ab2=-C3*((bet2 + bet3)/C1 - m3/C1  + bet3/C3 + uhat33- uhat13 - uhat33);
    

    %alf3  alf3
    %bet2  bet2
    %bet2  bet3
    phi2=(bet2 + bet3)/C1 - m3/C1  + bet2/C2 + uhat23- uhat13;
    phi3=(bet2 + bet3)/C1 - m3/C1  + bet3/C3 + uhat33- uhat13;
    
    util1=(A1*uhat11 + A2*(uhat21*cos(alf2)- uhat22*sin(alf2) ) + A3*(uhat31*cos(alf3)- uhat32*sin(alf3)))/(A1+A2+A3);
	util2=(A1*uhat12 + A2*(uhat21*sin(alf2)+ uhat22*cos(alf2) ) + A3*(uhat31*sin(alf3)+ uhat32*cos(alf3)))/(A1+A2+A3);
	util3= (ab1 + ab2 + C1*uhat13)/(C1+C2+C3);

     um = m1*util1 + m2*util2 + m3*util3;
     mkm = (m1*m1 + m2*m2)/(A1+A2+A3) + m3*m3/(C1+C2+C3);
     nd3= dot(N3(:,k),[D31,D32,D33]); 
     phibe= bet2*phi2  + bet3*phi3;
     Cab= uhat11*A1*uhat11 + uhat12*A1*uhat12 + uhat13*C1*uhat13 + uhat21*A2*uhat21 + uhat22*A2*uhat22 + (uhat23-phi2)*C2*(uhat23-phi2) + uhat31*A3*uhat31 + uhat32*A3*uhat32 + (uhat33-phi3)*C3*(uhat33-phi3)  - ( (util1*(A1+A2+A3)*util1)  + (util2*(A1+A2+A3)*util2) + (util3*(C1+C2+C3)*util3));
    % Hm=(1/(2*(A1 + A2 + A3)*C1*C2*C3))*(A1*(bet3^2)*C1*C2 + A2*(bet3^2)*C1*C2 + A3*(bet3^2)*C1*C2 + A1*(bet2^2)*C1*C3 + A2*(bet2^2)*C1*C3 + A3*(bet2^2)*C1*C3 + A1*(bet2^2)*C2*C3 +  A2*(bet2^2)*C2*C3 + A3*(bet2^2)*C2*C3 + 2*A1*bet2*bet3*C2*C3 + 2*A2*bet2*bet3*C2*C3 + 2*A3*bet2*bet3*C2*C3 + A1*(bet3^2)*C2*C3 + A2*(bet3^2)*C2*C3 + A3*(bet3^2)*C2*C3 -  A1*A2*C1*C2*C3*(uhat11^2) - A1*A3*C1*C2*C3*(uhat11^2) - 2*A1*bet2*C1*C2*C3*uhat13 - 2*A2*bet2*C1*C2*C3*uhat13 - 2*A3*bet2*C1*C2*C3*uhat13 - 2*A1*bet3*C1*C2*C3*uhat13 -  2*A2*bet3*C1*C2*C3*uhat13 - 2*A3*bet3*C1*C2*C3*uhat13 - A1*A2*C1*C2*C3*(uhat21^2) - A2*A3*C1*C2*C3*(uhat21^2) + 2*A1*bet2*C1*C2*C3*uhat23 + 2*A2*bet2*C1*C2*C3*uhat23 + 2*A3*bet2*C1*C2*C3*uhat23 - A1*A3*C1*C2*C3*(uhat31^2) - A2*A3*C1*C2*C3*(uhat31^2) + 2*A1*bet3*C1*C2*C3*uhat33 + 2*A2*bet3*C1*C2*C3*uhat33 + 2*A3*bet3*C1*C2*C3*uhat33 +  2*A1*A2*C1*C2*C3*uhat11*uhat21*cos(alf2) + 2*A2*A3*C1*C2*C3*uhat21*uhat31*cos(alf2 - alf3) + 2*A1*A3*C1*C2*C3*uhat11*uhat31*cos(alf3));
     %Ham3(k)= nd3 + um + phibe + mkm/2 - Cab/2 ;
     
	% Hm= A1*(bet3^2)*C1*C2 + A2*(bet3^2)*C1*C2 + A3*(bet3^2)*C1*C2 + A1*(bet2^2)*C1*C3 + A2*(bet2^2)*C1*C3 + A3*(bet2^2)*C1*C3 + A1*(bet2^2)*C2*C3 + A2*(bet2^2)*C2*C3 + A3*(bet2^2)*C2*C3 + 2*A1*bet2*bet3*C2*C3 + 2*A2*bet2*bet3*C2*C3 + 2*A3*bet2*bet3*C2*C3 + A1*(bet3^2)*C2*C3 + A2*(bet3^2)*C2*C3 + A3*(bet3^2)*C2*C3 + C1*C2*C3*(m1^2) +  C1*C2*C3*(m2^2) -  2*A1*bet2*C2*C3*m3 - 2*A2*bet2*C2*C3*m3 - 2*A3*bet2*C2*C3*m3 -  2*A1*bet3*C2*C3*m3 - 2*A2*bet3*C2*C3*m3 - 2*A3*bet3*C2*C3*m3 +  A1*C2*C3*(m3^2) + A2*C2*C3*(m3^2) +   A3*C2*C3*(m3^2) + 2*A1*C1*C2*C3*m1*uhat11 - A1*A2*C1*C2*C3*(uhat11^2) -  A1*A3*C1*C2*C3*(uhat11^2) - 2*A1*bet2*C1*C2*C3*uhat13 -  2*A2*bet2*C1*C2*C3*uhat13 -  2*A3*bet2*C1*C2*C3*uhat13 - 2*A1*bet3*C1*C2*C3*uhat13 -  2*A2*bet3*C1*C2*C3*uhat13 - 2*A3*bet3*C1*C2*C3*uhat13 +  2*A1*C1*C2*C3*m3*uhat13 + 2*A2*C1*C2*C3*m3*uhat13 +  2*A3*C1*C2*C3*m3*uhat13 - A1*A2*C1*C2*C3*(uhat21^2) -  A2*A3*C1*C2*C3*(uhat21^2) + 2*A1*bet2*C1*C2*C3*uhat23 +  2*A2*bet2*C1*C2*C3*uhat23 + 2*A3*bet2*C1*C2*C3*uhat23 - A1*A3*C1*C2*C3*(uhat31^2) - A2*A3*C1*C2*C3*(uhat31^2) +  2*A1*bet3*C1*C2*C3*uhat33 + 2*A2*bet3*C1*C2*C3*uhat33 +  2*A3*bet3*C1*C2*C3*uhat33 +  2*A2*C1*C2*C3*(m1 + A1*uhat11)*uhat21*cos(alf2) +  2*A2*A3*C1*C2*C3*uhat21*uhat31*cos(alf2 - alf3) + 2*A3*C1*C2*C3*m1*uhat31*cos(alf3) +  2*A1*A3*C1*C2*C3*uhat11*uhat31*cos(alf3) + 2*A2*C1*C2*C3*m2*uhat21*sin(alf2) + 2*A3*C1*C2*C3*m2*uhat31*sin(alf3);
    
     Hm=(1/(2*(A1 + A2 + A3)*C1*C2*C3))*( A1*(bet3^2)*C1*C2 + A2*(bet3^2)*C1*C2 + A3*(bet3^2)*C1*C2 + A1*(bet2^2)*C1*C3 + A2*(bet2^2)*C1*C3 + A3*(bet2^2)*C1*C3 + A1*(bet2^2)*C2*C3 + A2*(bet2^2)*C2*C3 + A3*(bet2^2)*C2*C3 + 2*A1*bet2*bet3*C2*C3 +2*A2*bet2*bet3*C2*C3 + 2*A3*bet2*bet3*C2*C3 + A1*(bet3^2)*C2*C3 + A2*(bet3^2)*C2*C3 + A3*(bet3^2)*C2*C3 + C1*C2*C3*m1^2 + C1*C2*C3*m2^2 - 2*A1*bet2*C2*C3*m3 - 2*A2*bet2*C2*C3*m3 -2*A3*bet2*C2*C3*m3 - 2*A1*bet3*C2*C3*m3 - 2*A2*bet3*C2*C3*m3 - 2*A3*bet3*C2*C3*m3 + A1*C2*C3*(m3^2) + A2*C2*C3*(m3^2) + A3*C2*C3*(m3^2) + 2*A1*C1*C2*C3*m1*uhat11 - A1*A2*C1*C2*C3*(uhat11^2) - A1*A3*C1*C2*C3*(uhat11^2) -    2*A1*bet2*C1*C2*C3*uhat13 - 2*A2*bet2*C1*C2*C3*uhat13 - 2*A3*bet2*C1*C2*C3*uhat13 - 2*A1*bet3*C1*C2*C3*uhat13 - 2*A2*bet3*C1*C2*C3*uhat13 -2*A3*bet3*C1*C2*C3*uhat13 +    2*A1*C1*C2*C3*m3*uhat13 + 2*A2*C1*C2*C3*m3*uhat13 + 2*A3*C1*C2*C3*m3*uhat13 -A1*A2*C1*C2*C3*(uhat21^2) - A2*A3*C1*C2*C3*(uhat21^2) + 2*A1*bet2*C1*C2*C3*uhat23 +    2*A2*bet2*C1*C2*C3*uhat23 + 2*A3*bet2*C1*C2*C3*uhat23 - A1*A3*C1*C2*C3*(uhat31^2) -A2*A3*C1*C2*C3*(uhat31^2) + 2*A1*bet3*C1*C2*C3*uhat33 + 2*A2*bet3*C1*C2*C3*uhat33 +   2*A3*bet3*C1*C2*C3*uhat33 + 2*A2*C1*C2*C3*(m1 + A1*uhat11)*uhat21*cos(alf2) + 2*A2*A3*C1*C2*C3*uhat21*uhat31*cos(alf2 - alf3) + 2*A3*C1*C2*C3*m1*uhat31*cos(alf3) + 2*A1*A3*C1*C2*C3*uhat11*uhat31*cos(alf3) + 2*A2*C1*C2*C3*m2*uhat21*sin(alf2) + 2*A3*C1*C2*C3*m2*uhat31*sin(alf3));
     Ham3(k)=Hm + nd3;
end

figure(1)
plot(S,NOR1-1.0,'r-o',S,NOR2-1.0,'b-o',S,NOR3-1.0,'k-o')
grid on

figure(2)
plot(S,INT1,'r-o',S,INT2,'b-o',S,INT3,'k-o')
grid on

figure(3)
subplot(3,1,1);
plot(S,Ham1- min(Ham1),'r-o')
grid on

subplot(3,1,2); 
plot(S,Ham2- min(Ham2),'b-o')
grid on

subplot(3,1,3); 
plot(S,Ham3- min(Ham3),'k-o')
grid on

figure(4)
subplot(3,1,1);
plot(S,INT_MOM1)
grid on

subplot(3,1,2); 
plot(S,INT_MOM2)
grid on

subplot(3,1,3); 
plot(S,INT_MOM3)
grid on


function [m1,m2,m3]=Local_Moment(Mu,Q)
    B1=[0 ,0, 0, 1; 0, 0, 1, 0; 0,-1, 0, 0;-1, 0, 0, 0];
    B2=[0 ,0,-1, 0; 0, 0, 0, 1; 1, 0, 0, 0; 0,-1, 0, 0];
    B3=[0 ,1, 0, 0;-1, 0, 0, 0; 0, 0, 0, 1; 0, 0, -1, 0];
    m1=dot(Mu,(B1*Q))/2;
    m2=dot(Mu,(B2*Q))/2;
    m3=dot(Mu,(B3*Q))/2;

 function [m1,m2,m3]=Global_Moment(Mu,Q)
    F1=[0 ,0, 0, 1; 0, 0, -1, 0; 0, 1, 0, 0;-1, 0, 0, 0];
    F2=[0 ,0, 1, 0; 0, 0, 0, 1;-1, 0, 0, 0; 0,-1, 0, 0];
    F3=[0 ,-1, 0, 0;1, 0, 0, 0; 0, 0, 0, 1; 0, 0, -1, 0];
    m1=dot(Mu,(F1*Q))/2;
    m2=dot(Mu,(F2*Q))/2;
    m3=dot(Mu,(F3*Q))/2;
   
    
    
function [d31,d32,d33]=director_D3(Q)
d31=2*(Q(1)*Q(3) + Q(2)*Q(4));
d32=2*(-Q(1)*Q(4) + Q(2)*Q(3));
d33=-Q(1)*Q(1) - Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4);
