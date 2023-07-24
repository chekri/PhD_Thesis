function [sum,jac] = Objectivefunction(Y)
R_init=Y(1:6)
y=Y;
y(1:6)=[];
Rtar=;
R_init=;
UFL=
P=[0,0,0];
TRF=
TRF2=
Objective(y,Rtar,R_init,Ori,UFL,USt,P,TRF,TRF2)


function [sum,jac] = Objective(y,Rtar,R_init,Ori,UFL,USt,P,TRF,TRF2) %{RF3,RF2,RF1,t3,t2,t1}
y=[USt,y]; 
%R_init=[-0.378,0.0015,1.44];
%R_init=[1.436,0.947,-1.263];

sum=0;
En=0;
total_path=0;
sweep_area=0;
reach_along=0;
N=floor(size(y,2)/12);
h=1/N;
p1=5.0; % 0.5 is good
%p2=15.0;
p2=0.0;%40;%40;%10.0;
p3=50.0;
p4=400;
rp1=1;
Tr_prev=TRF;
U_prev=USt;
T=1;
Rpath=R_init*((N-T)/N) + Rtar*(T/N);
%Rpath=R_init*((N-T)/N)^2 + Rtar*(T/N)^2;


En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*(y(10)^2 + y(11)^2 + y(12)^2) )*h/2;
prevsol=[];

U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
[reach,jac_reach]=Reach_targetandorientation(Rpath,Ori,U_current,Tr_current);
[sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
sweep_area= sweep*h; 
reach_along=  reach*h; 

jac(1:6,1) = p1*[y(7)*h;y(8)*h;y(9)*h;rp1*y(10)*h;rp1*y(11)*h;rp1*y(12)*h];
jac(7:12,1) = p2*h*jac_sweep'+ p3*h*jac_reach';

Tr_prev=Tr_current;
U_prev=U_current;
prevsol=Tr_current;
for T=(2:N)
    En = En + (y(12*(T-1) + 7)^2 + y(12*(T-1) + 8)^2 + y(12*(T-1)+ 9)^2 + rp1*(y(12*(T-1) + 10)^2 + y(12*(T-1) + 11)^2 + y(12*(T-1) + 12)^2) )*h/2;
    %Rpath=R_init*((N-T)/N)^2 + Rtar*(T/N)^2;
    Rpath=R_init*((N-T)/N) + Rtar*(T/N);
    U_current=[y(12*(T-1)+13),y(12*(T-1)+14),y(12*(T-1)+15),y(12*(T-1)+16),y(12*(T-1)+17),y(12*(T-1)+18)];
    Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
    [reach,jac_reach]=Reach_targetandorientation(Rpath,Ori,U_current,Tr_current);
    [sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
    sweep_area=sweep_area + sweep*h; 
    reach_along= reach_along + reach*h;
    jac(12*(T-1)+1:12*(T-1)+6,1) = p1*[y(12*(T-1) + 7)*h;y(12*(T-1) + 8)*h;y(12*(T-1) + 9)*h;rp1*y(12*(T-1) + 10)*h;rp1*y(12*(T-1) + 11)*h;rp1*y(12*(T-1) + 12)*h];
    jac(12*(T-1)+7:12*(T-1)+12,1) = p2*h*jac_sweep'+ p3*h*jac_reach';
    %jac(12*(T-1)+7:12*(T-1)+12,1) = p3*h*jac_reach';
    Tr_prev=Tr_current;
    U_prev=U_current;
    prevsol=Tr_current;
end

%[reach,jac_reach]=Reach_targetandorientation(Rtar,Ori,[y(12*(N-1)+13),y(12*(N-1)+14),y(12*(N-1)+15),y(12*(N-1)+16),y(12*(N-1)+17),y(12*(N-1)+18)],Tr_current);
reach_along= reach_along - reach*h/2; 
sweep_area=sweep_area - sweep*h/2;
reach;
sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) )/2 + p4*jac_reach';