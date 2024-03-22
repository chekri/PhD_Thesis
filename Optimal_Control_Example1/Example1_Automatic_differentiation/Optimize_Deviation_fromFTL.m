function [Output]=Optimize_Deviation_fromFTL(y_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output : Optimization solution from an initial soluton for the given set 
% of parameters, along with Target point, Tip Load and "Follow the Leader" 
% control values, useful for post-processing.Run Optimize_Deviation_fromFTL([]) 
% in command window to start optimization. This optimization framework
% utilizes "Deviation from FTL" as a measure of covered volume (Chapter 9 of thesis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Speicfy Target point
Rtar=[0.4,-0.0,1.0];
%%%%%%%%%%%%%%%%%%%%
%Specify Tip Load
F=[0.1,0.15,0.12];
%F=[0,0,0]
%Time steps
T_steps=10;
%Weights of objective functions
pp=1.0;
p1=5*pp; % 0.5 is good
p2=100.0*pp;%0,50,100, 200.0,500.0;
p3=0;
p4=800*pp;
%CTCR properties are specified in the Trajectory.m file.
%Specify Follow the Leader controls of the working environment
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial configuration
%Slight perturbation from Follow-the Leader
%Can be specified any value in principle
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
TRF2={RF3,RF2,RF1,t3,t2,t1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %TH1 TH2 L1 L2
% The state at t=0
Tr_initial=Trajectory(USt,F,[]);
S=Tr_initial.y;
%The Robot tip at t=0
R_init=[S(17,end),S(18,end),S(19,end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If initial solution is not provided, take the configuration at initial
%time as configrations at all time steps.
%If it is provided, take it as inital guess and interpolate over the
% required time steps.

if size(y_init,2)==0
    Yo=[];
    for T=(1:T_steps)
        Yo=[Yo,[0,0,0,0,0,0,USt(1)+0.0*rand(1)*pi,USt(2)+0.0*(rand(1))*pi,USt(3)+0.00*pi*rand(1),USt(4)+0.00*(rand(1)-1)*0.3,USt(5)+0.00*(rand(1)-1)*0.3,USt(6)+0.00*(rand(1)-1)*0.3]];
           end
else 
    Yo=y_init;
    N=T_steps;
    No=floor(size(Yo,2)/12);
    y5=Yo;
    a1=Yo(12*(0:No)+1);
    a1=interp1(linspace(0,1,No+1),a1,linspace(0,1,N+1));
    
    a2=y5(12*(0:No)+2);
    a2=interp1(linspace(0,1,No+1),a2,linspace(0,1,N+1));
    a3=y5(12*(0:No)+3);
    a3=interp1(linspace(0,1,No+1),a3,linspace(0,1,N+1));
    a4=y5(12*(0:No)+4);
    a4=interp1(linspace(0,1,No+1),a4,linspace(0,1,N+1));
    a5=y5(12*(0:No)+5);
    a5=interp1(linspace(0,1,No+1),a5,linspace(0,1,N+1));
    a6=y5(12*(0:No)+6);
    a6=interp1(linspace(0,1,No+1),a6,linspace(0,1,N+1));
    
    a7=y5(12*(0:No-1)+7);
    a7=interp1(linspace(0,1,No),a7,linspace(0,1,N));
    a8=y5(12*(0:No-1)+8);
    a8=interp1(linspace(0,1,No),a8,linspace(0,1,N));
    a9=y5(12*(0:No-1)+9);
    a9=interp1(linspace(0,1,No),a9,linspace(0,1,N)); 
    a10=y5(12*(0:No-1)+10);
    a10=interp1(linspace(0,1,No),a10,linspace(0,1,N));
    a11=y5(12*(0:No-1)+11);
    a11=interp1(linspace(0,1,No),a11,linspace(0,1,N));
    a12=y5(12*(0:No-1)+12);
    a12=interp1(linspace(0,1,No),a12,linspace(0,1,N));
        
    Yo=[];
    for T=(1:N)
        Yo=[Yo,[a7(T),a8(T),a9(T),a10(T),a11(T),a12(T),a1(T+1),a2(T+1),a3(T+1),a4(T+1),a5(T+1),a6(T+1)]];
    end
end

Aeq=zeros(12*(T_steps));
lb=-inf*ones(1,size(Yo,2));
ub= inf*ones(1,size(Yo,2));
h=1/T_steps;
%Linear Constraints
for T=[1:T_steps]
 Aeq(12*(T-1)+1,12*(T-1)+1)=1;
 Aeq(12*(T-1)+1,12*(T-1)+7)=1*h;
 Aeq(12*(T-1)+1,12*T+1)=-1;

 
 Aeq(12*(T-1)+2,12*(T-1)+2)=1;
 Aeq(12*(T-1)+2,12*(T-1)+8)=1*h;
 Aeq(12*(T-1)+2,12*T+2)=-1;

 Aeq(12*(T-1)+3,12*(T-1)+3)=1;
 Aeq(12*(T-1)+3,12*(T-1)+9)=1*h;
 Aeq(12*(T-1)+3,12*T+3)=-1;

 
 Aeq(12*(T-1)+4,12*(T-1)+4)=1;
 Aeq(12*(T-1)+4,12*(T-1)+10)=1*h;
 Aeq(12*(T-1)+4,12*T+4)=-1;

 Aeq(12*(T-1)+5,12*(T-1)+5)=1;
 Aeq(12*(T-1)+5,12*(T-1)+11)=1*h;
 Aeq(12*(T-1)+5,12*T+5)=-1;

 
 Aeq(12*(T-1)+6,12*(T-1)+6)=1;
 Aeq(12*(T-1)+6,12*(T-1)+12)=1*h;
 Aeq(12*(T-1)+6,12*T+6)=-1;

 %Bounds on the states
 lb(12*(T-1)+10)=0;
 lb(12*(T-1)+11)=0;
 lb(12*(T-1)+12)=0;

 ub(12*(T-1)+10)=0.7;
 ub(12*(T-1)+11)=0.7;
 ub(12*(T-1)+12)=0.7;
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

%Optimization controls
options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',true,'StepTolerance',1e-06,'OptimalityTolerance',1e-06,'ConstraintTolerance',1e-05);
%options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true);
Aeq*Yo' - beq ;

nonlcon=[];
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
fun3=@(U)Objective(U,Rtar,R_init,UFL,USt,F,Tr_initial,TRF2,p1,p2,p3,p4)
tic

Objective(Yo,Rtar,R_init,UFL,USt,F,Tr_initial,TRF2,p1,p2,p3,p4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve the optimization problem.
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,[],[],Aeq,beq,lb,ub,nonlcon,options);
disp('Time for the optimization: ')
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specifies the value of optimized objective
Objective(x,Rtar,R_init,UFL,USt,F,Tr_initial,TRF2,p1,p2,p3,p4)
Output={[USt,x],Rtar,F,UFL};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function provides objctive function and its jacobian
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
%Rpath for following a path 
%Straight line between initial point and target
jac(7:12,1)=[0,0,0,0,0,0];
En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*(y(10)^2 + y(11)^2 + y(12)^2) )*h/2;
prevsol=Tr_prev;
U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,F,prevsol); %TH1 TH2 L1 L2
[J,Fc]=Collocation_Equation(Tr_current,U_current,F);

[sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2,J,Fc);
sweep_area= sweep*h; 
Tr_prev=Tr_current;
U_prev=U_current;

jac(1:6,1) = p1*[y(7)*h;y(8)*h;y(9)*h;rp1*y(10)*h;rp1*y(11)*h;rp1*y(12)*h];
jac(7:12,1) = p2*h*jac_sweep;
prevsol=Tr_current;
for T=(1:N-1)
    En = En + (y(12*T + 7)^2 + y(12*T + 8)^2 + y(12*T + 9)^2 + rp1*(y(12*T + 10)^2 + y(12*T + 11)^2 + y(12*T + 12)^2) )*h/2;
    U_current=[y(12*T+12+1),y(12*T+12+2),y(12*T+12+3),y(12*T+12+4),y(12*T+12+5),y(12*T+12+6)];
    Tr_current=Trajectory(U_current,F,prevsol); %TH1 TH2 L1 L2
    [J,Fc]=Collocation_Equation(Tr_current,U_current,F); 
    [sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2,J,Fc);
    sweep_area=sweep_area + sweep*h;
    jac(12*T+1:12*T+6,1) = p1*[y(12*T + 7)*h;y(12*T + 8)*h;y(12*T + 9)*h;rp1*y(12*T + 10)*h;rp1*y(12*T + 11)*h;rp1*y(12*T + 12)*h];
    jac(12*T+7:12*T+12,1) = p2*h*jac_sweep';  
    prevsol=Tr_current;
 
end
[reach,jac_reach]=Reach_target(Rtar,Tr_current,J,Fc);
sweep_area=sweep_area - sweep*h/2;
sum = p1*En + p2*sweep_area  + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) )/2 + p4*jac_reach';