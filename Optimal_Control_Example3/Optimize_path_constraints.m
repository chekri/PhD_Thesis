function [a]=Optimize_path_constraints(y_init)
initialize() 
%R_init=[1.436,0.947,-1.263];
Rtar=[-0.174953,0.32049,1.44671]
Rtar=[-0.2,0.44,1.1]
Rtar=[-0.18,0.32,1.446]
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4+1,0.6,0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tr_current=Trajectory(UFL,[0,0,0],[]); %TH1 TH2 L1 L2
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4+1,0.6,0.5];
USt=UFL + [0,0.2,1.3,0.0,0,0];
%USt=[0.5,2.0,3.4,0.5,0.6,0.5]
USt=[0.0,0.2,1.0,0.75,0.6,0.5];
USt=[0.5,0.5+2.1415,0.75+3.1415,0.4,0.6,0.5];
P=[0,0,0];
S=Trajectory(USt,[0,0,0]).y;
R_init=[S(17,end),S(18,end),S(19,end)];
Rtar=R_init;
T_steps=5;
lam=1;
if size(y_init,2)==0
    Yo=[0,0,0,0,0,0];
    for T=(1:T_steps)
        %Yo=[Yo,[USt(1)+0.1*rand(1)*pi,USt(2)+0.1*(rand(1))*pi,USt(3)+0.10*pi*rand(1),USt(4)+0.1*(rand(1)-1)*0.3,USt(5)+0.10*(rand(1)-1)*0.3,USt(6)+0.10*(rand(1)-1)*0.3,0,0,0,0,0,0]];
        Yo=[Yo,[USt(1),USt(2),USt(3),USt(4),USt(5),USt(6),0,0,0,0,0,0]];

    end
else 
    Yo=y_init;
    Yo(1:6)=[];
    for T=(1:size(Yo,2))
        Yo(T)=Yo(T)+ 0.1*rand(1)
    end
end


Aa=[]
b=[];
Aeq=[];
beq=[];

lb=-inf*ones(1,size(Yo,2));
ub= inf*ones(1,size(Yo,2));

Aeq=zeros(12*(T_steps + 1));
h=1/T_steps;
for T=[1:T_steps]
 Aeq(12*(T-1)+1,12*(T-1)+1)=1;
 Aeq(12*(T-1)+1,12*(T-1)+7)=1*h/2;
 Aeq(12*(T-1)+1,12*T+1)=-1;
 Aeq(12*(T-1)+1,12*T+7)=1*h/2;
 
 Aeq(12*(T-1)+2,12*(T-1)+2)=1;
 Aeq(12*(T-1)+2,12*(T-1)+8)=1*h/2;
 Aeq(12*(T-1)+2,12*T+2)=-1;
 Aeq(12*(T-1)+2,12*T+8)=1*h/2;
 
 Aeq(12*(T-1)+3,12*(T-1)+3)=1;
 Aeq(12*(T-1)+3,12*(T-1)+9)=1*h/2;
 Aeq(12*(T-1)+3,12*T+3)=-1;
 Aeq(12*(T-1)+3,12*T+9)=1*h/2;
 
 Aeq(12*(T-1)+4,12*(T-1)+4)=1;
 Aeq(12*(T-1)+4,12*(T-1)+10)=1*h/2;
 Aeq(12*(T-1)+4,12*T+4)=-1;
 Aeq(12*(T-1)+4,12*T+10)=1*h/2;
 
 Aeq(12*(T-1)+5,12*(T-1)+5)=1;
 Aeq(12*(T-1)+5,12*(T-1)+11)=1*h/2;
 Aeq(12*(T-1)+5,12*T+5)=-1;
 Aeq(12*(T-1)+5,12*T+11)=1*h/2;
 
 Aeq(12*(T-1)+6,12*(T-1)+6)=1;
 Aeq(12*(T-1)+6,12*(T-1)+12)=1*h/2;
 Aeq(12*(T-1)+6,12*T+6)=-1;
 Aeq(12*(T-1)+6,12*T+12)=1*h/2;
 
 lb(12*(T-1)+10)=0.001;
 lb(12*(T-1)+11)=0.001;
 lb(12*(T-1)+12)=0.001;

 ub(12*(T-1)+10)=1.5;
 ub(12*(T-1)+11)=1.5;
 ub(12*(T-1)+12)=1.5;
end

beq=zeros(12*(T_steps)+6,1);
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

Aeq(end,:)=[];
Aeq(end,:)=[];
Aeq(end,:)=[];
Aeq(end,:)=[];
Aeq(end,:)=[];
Aeq(end,:)=[];
options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'Diagnostics','on','MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',false);
nonlcon=[];
nonlcon=@(u)stay_atpoint(u,Rtar)
Objective(Yo,Rtar,R_init,UFL,USt,P,Tr_current,TRF2)
fun3=@(U)Objective(U,Rtar,R_init,UFL,USt,P,Tr_current,TRF2)

tic
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,Aa,b,Aeq,beq,lb,ub,nonlcon,options);
toc
disp('Time for 10 time steps ')
y10=[USt,x];
end_c=y10(12*T_steps + 1:12*T_steps + 6)
Tr=Trajectory(end_c,[0,0,0],[])
a=y10;

function [sum,jac] = Objective(y,Rtar,R_init,UFL,USt,P,TRF,TRF2) %{RF3,RF2,RF1,t3,t2,t1}
y=[USt,y]; 
%R_init=[-0.378,0.0015,1.44];
%R_init=[1.436,0.947,-1.263];
sum=0;
En=0;
total_path=0;
sweep_area=0;
reach_along=0;
N=floor(size(y,2)/12)-1;
tar_ori=[0,0.707,0.707];
tar_ori=[0.0,0.0,1.0];
h=1/N;
p1=1; % 0.5 is good
p2=0.0;
p3=0.0;%400;
p4=4;
rp1=1;
Tr_prev=TRF;
U_prev=USt;
T=1;
Rpath=R_init*((N-T)/N) + Rtar*(T/N);
%Rpath=R_init*((N-T)/N)^2 + Rtar*(T/N)^2;
jac(7:12,1)=[0,0,0,0,0,0];
En = En + (y(12 + 7)^2 + y(12 + 8)^2 + y(12 + 9)^2 + rp1*(y(12 + 10)^2 + y(12 + 11)^2 + y(12+ 12)^2) )*h/2;
prevsol=[];
U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
[reach,jac_reach]=Reach_target(Rpath,[y(12*T+1),y(12*T+2),y(12*T+3),y(12*T+4),y(12*T+5),y(12*T+6)],Tr_current);
[sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);

sweep_area= sweep*h; 
reach_along=  reach*h; 

jac(12+1:12+6,1) = p1*[y(12 + 7)*h;y(12 + 8)*h;y(12 + 9)*h;rp1*y(12 + 10)*h;rp1*y(12 + 11)*h;rp1*y(12 + 12)*h];
jac(12*(T-1)+7:12*(T-1)+12,1) = p2*h*jac_sweep'+ p3*h*jac_reach';
%jac(12*(T-1)+7:12*(T-1)+12,1) = p3*h*jac_reach';
Tr_prev=Tr_current;
U_prev=U_current;
prevsol=Tr_current;
for T=(2:N)
    En = En + (y(12*T + 7)^2 + y(12*T + 8)^2 + y(12*T + 9)^2 + rp1*(y(12*T + 10)^2 + y(12*T + 11)^2 + y(12*T + 12)^2) )*h/2;
    %Rpath=R_init*((N-T)/N)^2 + Rtar*(T/N)^2;
    %Rpath=R_init*((N-T)/N) + Rtar*(T/N);
    U_current=[y(12*T+1),y(12*T+2),y(12*T+3),y(12*T+4),y(12*T+5),y(12*T+6)];
    Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
    %[reach,jac_reach]=Reach_target(Rpath,U_current,Tr_current);
    [sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
    sweep_area=sweep_area + sweep*h; 
    %reach_along= reach_along + reach*h; 
    jac(12*T+1:12*T+6,1) = p1*[y(12*T + 7)*h;y(12*T + 8)*h;y(12*T + 9)*h;rp1*y(12*T + 10)*h;rp1*y(12*T + 11)*h;rp1*y(12*T + 12)*h];
    jac(12*(T-1)+7:12*(T-1)+12,1) = p2*h*jac_sweep';%+ p3*h*jac_reach';
    %jac(12*(T-1)+7:12*(T-1)+12,1) =  p3*h*jac_reach';
    Tr_prev=Tr_current;
    U_prev=U_current;
    prevsol=Tr_current;
end

En = En - (y(12*N + 7)^2 + y(12*N + 8)^2 + y(12*N + 9)^2 + rp1*(y(12*N + 10)^2 +  y(12*N + 11)^2 + y(12*N + 12)^2))*h/4 ;
En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*( y(10)^2  + y(11)^2 + y(12)^2) )*h/4;
%reach_along= reach_along - reach*h/2;
[reach_tar_ori,jac_reach]=Orient(Rtar,tar_ori,U_current,Tr_current);
 
%sweep_area=sweep_area - sweep*h/2;
sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach_tar_ori;


jac(1:6,1) = p1*[y(7)*h/2;y(8)*h/2;y(9)*h/2;rp1*y(10)*h/2;rp1*y(11)*h/2;rp1*y(12)*h/2];
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) ) - p3/2*h*jac_reach' + p4*jac_reach';
jac(12*N+1:12*N+6,1)=jac(12*N+1:12*N+6,1)/2;


function [gneq,geq,gradneq,gradceq]=stay_atpoint(y,Rtar)
    tol=1e-01;
    N=floor(size(y,2)/12)-1;
    Jac_neq=zeros(size(y,2),N+1);
    gneq=zeros(1,N+1);
    prevsol=[];
    for T=(1:N)
        %global A;
        U_cur=[y(12*T+1-6),y(12*T+2-6),y(12*T+3-6),y(12*T+4-6),y(12*T+5-6),y(12*T+6-6)];
        Tr_current=Trajectory(U_cur,[0,0,0],prevsol); %TH1 TH2 L1 L2
        [reach,jac_reach]=Stay_at_target(Rtar,U_cur,Tr_current);
        gneq(T)=reach-tol;
        Jac_neq(12*T-5:12*T,T)=jac_reach;
        prevsol=Tr_current;
    end

    U_cur=[y(12*N+7),y(12*N+8),y(12*N+9),y(12*N+10),y(12*N+11),y(12*N+12)];
    Tr_current=Trajectory(U_cur,[0,0,0],prevsol); %TH1 TH2 L1 L2
    [reach,jac_reach]=Stay_at_target(Rtar,U_cur,Tr_current);
    gneq(T)=reach-tol;
    Jac_neq(12*N+7:12*N+12,N)=jac_reach;
    gradneq=Jac_neq;
    geq=[];
    gradceq=[];

