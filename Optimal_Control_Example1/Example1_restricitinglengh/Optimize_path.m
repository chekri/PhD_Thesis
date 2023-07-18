function [a]=Optimize_path(y_init)
initialize() 
Rtar=[-0.174953,0.32049,1.44671]
Rtar=[-0.2,0.44,1.1]
R_init=[-0.348,0.015,1.4]
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4+1,0.6,0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tr_current=Trajectory(UFL,[0,0,0],[]); %TH1 TH2 L1 L2
UFL=[0.,0.+3.1415,0.+3.1415,0.4,0.6,0.5];
USt=UFL + [0,0.2,0,0.0,0,0];
USt=[0.5,0.5,1.0,0.75,0.86,0.75]
P=[0,0,0]
T_steps=20;
lam=1;
if size(y_init,2)==0
    Yo=[0,0,0,0,0,0];
    for T=(1:T_steps)
        Yo=[Yo,[USt(1)+0.0*rand(1)*pi,USt(2)+0.0*(rand(1))*pi,USt(3)+0.00*pi*rand(1),USt(4)+0.00*(rand(1)-1)*0.3,USt(5)+0.00*(rand(1)-1)*0.3,USt(6)+0.00*(rand(1)-1)*0.3,0,0,0,0,0,0]];
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
 
 lb(12*(T-1)+10)=0;
 lb(12*(T-1)+11)=0;
 lb(12*(T-1)+12)=0;

 ub(12*(T-1)+10)=1.2;
 ub(12*(T-1)+11)=1.2;
 ub(12*(T-1)+12)=1.2;
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

options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'Diagnostics','on','MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',true);
nonlcon=[];
Objective(Yo,Rtar,UFL,USt,P,Tr_current,TRF2)
fun3=@(U)Objective(U,Rtar,UFL,USt,P,Tr_current,TRF2)

tic
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,Aa,b,Aeq,beq,lb,ub,nonlcon,options);
toc
disp('Time for 10 time steps ')
y10=[USt,x];
end_c=y10(12*T_steps + 1:12*T_steps + 6)
Tr=Trajectory(end_c,[0,0,0],[])
a=y10;

function [sum,jac] = Objective(y,Rtar,UFL,USt,P,TRF,TRF2) %{RF3,RF2,RF1,t3,t2,t1}
y=[USt,y]; 
sum=0;
En=0;
total_path=0;
sweep_area=0;
N=floor(size(y,2)/12)-1;
h=1/N;
p1=0.5; % 0.5 is good
p2=0.0;
p3=500;
R_init=[-0.348,0.015,1.4];
rp1=1;
Tr_prev=TRF;
U_prev=USt;
jac(7:12,1)=[0,0,0,0,0,0];
En = En + (y(12 + 7)^2 + y(12 + 8)^2 + y(12 + 9)^2 + rp1*(y(12 + 10)^2 + y(12 + 11)^2 + y(12+ 12)^2) )*h/2;
global A;
prevsol=A;
U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
[sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
sweep_area=sweep_area + sweep*h; 

jac(12+1:12+6,1) = p1*[y(12 + 7)*h;y(12 + 8)*h;y(12 + 9)*h;rp1*y(12 + 10)*h;rp1*y(12 + 11)*h;rp1*y(12 + 12)*h];
jac(12*(1-1)+7:12*(1-1)+12,1) =  p2*h*jac_sweep';

Tr_prev=Tr_current;
U_prev=U_current;

for T=(2:N)
    En = En + (y(12*T + 7)^2 + y(12*T + 8)^2 + y(12*T + 9)^2 + rp1*(y(12*T + 10)^2 + y(12*T + 11)^2 + y(12*T + 12)^2) )*h/2;
    global A;
    prevsol=A;
    U_current=[y(12*T+1),y(12*T+2),y(12*T+3),y(12*T+4),y(12*T+5),y(12*T+6)];
    Tr_current=Trajectory(U_current,[0,0,0],prevsol); %TH1 TH2 L1 L2
    [sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
    sweep_area=sweep_area + sweep*h; 
    jac(12*T+1:12*T+6,1) = p1*[y(12*T + 7)*h;y(12*T + 8)*h;y(12*T + 9)*h;rp1*y(12*T + 10)*h;rp1*y(12*T + 11)*h;rp1*y(12*T + 12)*h];
    jac(12*(T-1)+7:12*(T-1)+12,1) = p2*h*jac_sweep';
    Tr_prev=Tr_current;
    U_prev=U_current;
end

En = En - (y(12*N + 7)^2 + y(12*N + 8)^2 + y(12*N + 9)^2 + rp1*(y(12*N + 10)^2 +  y(12*N + 11)^2 + y(12*N + 12)^2))*h/4 ;
En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*( y(10)^2  + y(11)^2 + y(12)^2) )*h/4;
[reach,jac_reach]=Reach_target(Rtar,[y(12*N+1),y(12*N+2),y(12*N+3),y(12*N+4),y(12*N+5),y(12*N+6)],Tr_current);
reach;
En;
sweep_area;
total_path;
reach;
sum = p1*En + p2*sweep_area + p3*reach;
jac(1:6,1) = p1*[y(7)*h/2;y(8)*h/2;y(9)*h/2;rp1*y(10)*h/2;rp1*y(11)*h/2;rp1*y(12)*h/2];
jac(12*(N-1)+7:12*(N-1)+12,1)=jac(12*(N-1)+7:12*(N-1)+12,1)/2 + p3*jac_reach';
jac(12*N+1:12*N+6,1)=jac(12*N+1:12*N+6,1)/2;