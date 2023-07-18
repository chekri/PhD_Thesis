function [a]=Optimize_path_alongpath(y_init)

%%% Target point
Rtar=[0.4,-0.0,1.0]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Penalization weights
pp=1.0
p1=5*pp; % 0.5 is good
p2=400.0*pp;%20.0;
p3=0*pp;
p4=800*pp;
%% Tip Load
Load=[0.1,0.15,0.12];

%% Reference 'Follow the Leader' controls
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
Tr_current=Trajectory(UFL,Load,[]); %TH1 TH2 L1 L2
%% Initial Configuration
USt=UFL + [0,0.1,0.2,0.0,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3d shape of the 'Follow the leader' curve
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4+1,0.6,0.5];
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=Trajectory(USt,Load,[]).y;
R_init=[S(17,end),S(18,end),S(19,end)];
T_steps=10;
lam=1;
if size(y_init,2)==0
    Yo=[];
    for T=(1:T_steps)
        Yo=[Yo,[0,0,0,0,0,0,USt(1)+0.0*rand(1)*pi,USt(2)+0.0*(rand(1))*pi,USt(3)+0.00*pi*rand(1),USt(4)+0.00*(rand(1)-1)*0.3,USt(5)+0.00*(rand(1)-1)*0.3,USt(6)+0.00*(rand(1)-1)*0.3]];
    end
else 
    Yo=y_init;
    Yo(1:6)=[];
    for T=(1:size(Yo,2))
        Yo(T)=Yo(T)+ 0.00*rand(1)
    end
end


Aa=[]
b=[];
Aeq=[];
beq=[];

lb=-inf*ones(1,size(Yo,2));
ub= inf*ones(1,size(Yo,2));

Aeq=zeros(12*(T_steps));
h=1/T_steps;
for T=[1:T_steps]
 Aeq(12*(T-1)+1,12*(T-1)+1)=1;
 Aeq(12*(T-1)+1,12*(T-1)+7)=1*h;
 Aeq(12*(T-1)+1,12*T+1)=-1;
% Aeq(12*(T-1)+1,12*T+7)=1*h/2;
 
 Aeq(12*(T-1)+2,12*(T-1)+2)=1;
 Aeq(12*(T-1)+2,12*(T-1)+8)=1*h;
 Aeq(12*(T-1)+2,12*T+2)=-1;
% Aeq(12*(T-1)+2,12*T+8)=1*h/2;
 
 Aeq(12*(T-1)+3,12*(T-1)+3)=1;
 Aeq(12*(T-1)+3,12*(T-1)+9)=1*h;
 Aeq(12*(T-1)+3,12*T+3)=-1;
% Aeq(12*(T-1)+3,12*T+9)=1*h/2;
 
 Aeq(12*(T-1)+4,12*(T-1)+4)=1;
 Aeq(12*(T-1)+4,12*(T-1)+10)=1*h;
 Aeq(12*(T-1)+4,12*T+4)=-1;
% Aeq(12*(T-1)+4,12*T+10)=1*h/2;
 
 Aeq(12*(T-1)+5,12*(T-1)+5)=1;
 Aeq(12*(T-1)+5,12*(T-1)+11)=1*h;
 Aeq(12*(T-1)+5,12*T+5)=-1;
% Aeq(12*(T-1)+5,12*T+11)=1*h/2;
 
 Aeq(12*(T-1)+6,12*(T-1)+6)=1;
 Aeq(12*(T-1)+6,12*(T-1)+12)=1*h;
 Aeq(12*(T-1)+6,12*T+6)=-1;
 %Aeq(12*(T-1)+6,12*T+12)=1*h/2;
 
 lb(12*(T-1)+10)=0;
 lb(12*(T-1)+11)=0;
 lb(12*(T-1)+12)=0;

 ub(12*(T-1)+10)=0.8;
 ub(12*(T-1)+11)=0.8;
 ub(12*(T-1)+12)=0.8;
end

beq=zeros(12*(T_steps),1);
beq(1)=-USt(1);
beq(2)=-USt(2);
beq(3)=-USt(3);
beq(4)=-USt(4);
beq(5)=-USt(5);
beq(6)=-USt(6);
Aeq(:,1)=[];
Aeq(:,1)=[];
Aeq(:,1)=[];
Aeq(:,1)=[];
Aeq(:,1)=[];
Aeq(:,1)=[];



Yo
options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'Diagnostics','on','MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',true,'StepTolerance',0.2*1e-09,'OptimalityTolerance',0.2*1e-09);
nonlcon=[];
Objective(Yo,Rtar,R_init,UFL,USt,Load,Tr_current,TRF2,p1,p2,p3,p4)
fun3=@(U)Objective(U,Rtar,R_init,UFL,USt,Load,Tr_current,TRF2,p1,p2,p3,p4)

tic
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,Aa,b,Aeq,beq,lb,ub,nonlcon,options);
toc
disp('Time for 10 time steps ')
y10=[USt,x];
a=y10;
Objective(x,Rtar,R_init,UFL,USt,Load,Tr_current,TRF2,p1,p2,p3,p4)

function [sum,jac] = Objective(y,Rtar,R_init,UFL,USt,P,TRF,TRF2,p1,p2,p3,p4) 
y=[USt,y]; 
sum=0;
En=0;
total_path=0;
sweep_area=0;
reach_along=0;
Load=P;
N=floor(size(y,2)/12);
h=1/N;
rp1=1.0;
Tr_prev=TRF;
U_prev=USt;
T=1;
Rpath=R_init*((N-T)/N) + Rtar*(T/N);
%Rpath=R_init*((N-T)/N)^2 + Rtar*(T/N)^2;
jac(7:12,1)=[0,0,0,0,0,0];

En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*(y(10)^2 + y(11)^2 + y(12)^2) )*h/2;
prevsol=[];

U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
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
    Tr_current=Trajectory(U_current,Load,prevsol); %TH1 TH2 L1 L2
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
[reach,jac_reach]=Reach_target2(Rtar,[y(12*(N-1)+13),y(12*(N-1)+14),y(12*(N-1)+15),y(12*(N-1)+16),y(12*(N-1)+17),y(12*(N-1)+18)],Tr_current,Load);

reach_along= reach_along - reach*h/2; 
sweep_area=sweep_area - sweep*h/2;
sweep_area
reach
sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) )/2 + p4*jac_reach';