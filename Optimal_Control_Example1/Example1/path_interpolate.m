function  path_interpolate(Y1,Y2)
Rtar=[0.4,-0.0,1.0];
%%%%%%%%%%%%%%%%%%%%
%Specify Tip Load
F=[0.1,0.0,0.11];
%F=[0,0,0]
%Weights of objective functions
pp=1.0;
p1=5*pp; % 0.5 is good
p2=100.0*pp;%0,50,100, 200.0,500.0;
p3=0.0;
p4=800*pp;

%Specify Follow the Leader controls of the working environment
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial configuration
%Slight perturbation from Follow-the Leader
USt=UFL + [0,0.1,0.2,0.0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference "Follow the leader" Curve
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.99+0.4,0.6,0.5];
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1}

UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %TH1 TH2 L1 L2
% The state at t=0
Tr_initial=Trajectory(USt,F,[]);
S=Tr_initial.y;
%The Robot tip at t=0
R_init=[S(17,end),S(18,end),S(19,end)];


for k=0:0.1:1
    Y=k*(Y1) + (1-k)*Y2;
    Y(1:6)=[];
    Objective(Y,Rtar,R_init,UFL,USt,F,Tr_initial,TRF2,p1,p2,p3,p4)
end



function [sum,jac] = Objective(y,Rtar,R_init,UFL,USt,F,TRF,TRF2,p1,p2,p3,p4) 
y=[USt,y]; 
En=0;
%Number of steps
N=floor(size(y,2)/12);
h=1/N;
rp1=1;
%Previous solution provided as an initial guess to the bvp to get its quicker convergence.
Tr_prev=TRF;
U_prev=USt;
T=1;
%Rpath for follwing a path 
%Straight line between initial point and target
Rpath=R_init*((N-T)/N) + Rtar*(T/N);
jac(7:12,1)=[0,0,0,0,0,0];
En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*(y(10)^2 + y(11)^2 + y(12)^2) )*h/2;
prevsol=Tr_prev;
U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,F,prevsol); %TH1 TH2 L1 L2
[reach,jac_reach]=Reach_target(Rpath,U_current,Tr_current);
[sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
sweep_area= sweep*h; 
reach_along=  reach*h; 
jac(1:6,1) = p1*[y(7)*h;y(8)*h;y(9)*h;rp1*y(10)*h;rp1*y(11)*h;rp1*y(12)*h];
jac(7:12,1) = p2*h*jac_sweep'+ p3*h*jac_reach';
Tr_prev=Tr_current;
U_prev=U_current;
prevsol=Tr_current;
for T=(1:N-1)
    En = En + (y(12*T + 7)^2 + y(12*T + 8)^2 + y(12*T + 9)^2 + rp1*(y(12*T + 10)^2 + y(12*T + 11)^2 + y(12*T + 12)^2) )*h/2;
    Rpath=R_init*((N-(T+1))/N) + Rtar*((T+1)/N);
    U_current=[y(12*T+12+1),y(12*T+12+2),y(12*T+12+3),y(12*T+12+4),y(12*T+12+5),y(12*T+12+6)];
    Tr_current=Trajectory(U_current,F,prevsol); %TH1 TH2 L1 L2
    [reach,jac_reach]=Reach_target(Rpath,U_current,Tr_current);
    [sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
    sweep_area=sweep_area + sweep*h; 
    reach_along= reach_along + reach*h;
    jac(12*T+1:12*T+6,1) = p1*[y(12*T + 7)*h;y(12*T + 8)*h;y(12*T + 9)*h;rp1*y(12*T + 10)*h;rp1*y(12*T + 11)*h;rp1*y(12*T + 12)*h];
    jac(12*T+7:12*T+12,1) = p2*h*jac_sweep'+ p3*h*jac_reach';  
    Tr_prev=Tr_current;
    U_prev=U_current;
    prevsol=Tr_current;
end
[reach,jac_reach]=Reach_target2(Rtar,[y(12*(N-1)+13),y(12*(N-1)+14),y(12*(N-1)+15),y(12*(N-1)+16),y(12*(N-1)+17),y(12*(N-1)+18)],Tr_current,F);
%[reach,jac_reach]=Reach_target(Rtar,[y(12*(N-1)+13),y(12*(N-1)+14),y(12*(N-1)+15),y(12*(N-1)+16),y(12*(N-1)+17),y(12*(N-1)+18)],Tr_current)

reach_along= reach_along - reach*h/2; 
sweep_area=sweep_area - sweep*h/2;

sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) )/2 + p4*jac_reach';