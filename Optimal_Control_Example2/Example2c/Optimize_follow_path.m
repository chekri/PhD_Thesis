function [Output]=Optimize_follow_path(y_init)
%Run Optimize_path_alongpath([]) in command window to start the
%optimization
%(Example 2 in the Chapter 10). 
%Navigate the CTCR such that it follows a path and maintains a Fixed
%Orientation. Orientation is maintained through constraints.
%Can be run from the provided solution y_init if available, otherwise start
%from []
%Output: Array of collection of control parameters, control parameter vectors at the mesh points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify Target point
Rtar= [-0.565, 0.2324, 1.45];

%Specify Follow the Leader controls of the working environment
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];

%Specify Tip Load
P=[0.1,0.01,0.12];
%Specify the mesh size: Number of time steps
T_steps=20;

%Specify the tolerance angle here (in degrees)
d_angle_tol=10;
%CTCR parameters such as the stiffnesses, intrinsic curvatures of tubes are
%specified in Trajectory.m file
%Penalizing Weights of the objective functions
p1=5.0; % 0.5 is good
p2=0.0;%40;%40;%10.0;
p3=40.0;
p4=400;

% Penaliing Weigths for prioritizing between Position and Orientation pursuit
nu=1.0; % 0.5 % For position
mu=0.0;%-2.5 % For orientation. maximizes dot product and that's why it is negative 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the shape of Follow The Leader curve with the specified UFL
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tr_current=Trajectory(UFL,P,[]); %TH1 TH2 L1 L2
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
USt=UFL + [0,0.1,0.2,0.0,0,0];
USt=[0.8267,3.4301,3.5024,0.1967,0.2684,0.3262];

S=Trajectory(USt,P,[]).y;
R_init=[S(17,end),S(18,end),S(19,end)];

%Orientation of the CTCR tip: quaternions
qend1=S(20,end);
qend2=S(21,end);
qend3=S(22,end);
qend4=S(23,end);

%Orientation of the CTCR tip tangent: d3-axis
d31=2*(qend1*qend3 + qend2*qend4);
d32=2*(qend2*qend3 - qend1*qend4);
d33=-qend1*qend1 - qend2*qend2 + qend3*qend3 + qend4*qend4;
Ori=[d31,d32,d33];

%Inflated UFL for computing the Deviation from FTL objective.(1 unit is added to the one tube section's length)
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4+1,0.6,0.5];

S=Trajectory(USt,P,[]).y;
R_init=[S(17,end),S(18,end),S(19,end)];

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
        Yo(T)=Yo(T)+ 0.1*rand(1)
    end
end

lb=-inf*ones(1,size(Yo,2));
ub= inf*ones(1,size(Yo,2));
Aeq=zeros(12*(T_steps));
h=1/T_steps;
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

 %Bounds on the length control paramters
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

options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'Diagnostics','on','MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',true,'StepTolerance',1e-12,'OptimalityTolerance',1e-05);
nonlcon=@(U)Orientation_Constraint(U,Ori,P,d_angle_tol);
fun3=@(U)Objective(U,Rtar,R_init,Ori,UFL,USt,P,Tr_current,TRF2,p1,p2,p3,p4,mu,nu)

tic
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,[],[],Aeq,beq,lb,ub,nonlcon,options);
disp('Time taken for the optimization:')
toc
Output={[USt,x],Rtar,P,UFL};
Objective(x,Rtar,R_init,Ori,UFL,USt,P,Tr_current,TRF2,p1,p2,p3,p4,mu,nu)

function [sum,jac] = Objective(y,Rtar,R_init,Ori,UFL,USt,P,TRF,TRF2,p1,p2,p3,p4,mu,nu) %{RF3,RF2,RF1,t3,t2,t1}
y=[USt,y]; 
sum=0;
En=0;
total_path=0;
sweep_area=0;
reach_along=0;
N=floor(size(y,2)/12);
h=1/N;

rp1=1;
Tr_prev=TRF;
U_prev=USt;
T=1;
Rpath=R_init*((N-T)/N) + Rtar*(T/N);
En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*(y(10)^2 + y(11)^2 + y(12)^2) )*h/2;
prevsol=[];

U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,P,prevsol); %TH1 TH2 L1 L2
[reach,jac_reach]=Reach_targetandorientation(Rpath,Ori,U_current,Tr_current,P,mu,nu);
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
    Rpath=R_init*((N-T)/N) + Rtar*(T/N);
    U_current=[y(12*(T-1)+13),y(12*(T-1)+14),y(12*(T-1)+15),y(12*(T-1)+16),y(12*(T-1)+17),y(12*(T-1)+18)];
    Tr_current=Trajectory(U_current,P,prevsol); 
    [reach,jac_reach]=Reach_targetandorientation(Rpath,Ori,U_current,Tr_current,P,mu,nu);
    [sweep,jac_sweep]=Deviation_from_FTL(U_current,Tr_current,UFL,TRF2);
    sweep_area=sweep_area + sweep*h; 
    reach_along= reach_along + reach*h;
    jac(12*(T-1)+1:12*(T-1)+6,1) = p1*[y(12*(T-1) + 7)*h;y(12*(T-1) + 8)*h;y(12*(T-1) + 9)*h;rp1*y(12*(T-1) + 10)*h;rp1*y(12*(T-1) + 11)*h;rp1*y(12*(T-1) + 12)*h];
    jac(12*(T-1)+7:12*(T-1)+12,1) = p2*h*jac_sweep'+ p3*h*jac_reach';
    Tr_prev=Tr_current;
    U_prev=U_current;
    prevsol=Tr_current;
end
reach_along= reach_along - reach*h/2; 
sweep_area=sweep_area - sweep*h/2;
reach;
sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) )/2 + p4*jac_reach';

function [gneq,geq,gradneq,gradceq]=Orientation_Constraint(y,Ori_tar,P,d_angle_tol)
% This function generates the nonlinear constraints 
% and their gradients with respect to control paramters.
% Ori_tar     - Target Orientation
% d_angle_tol - Tolerance angle
% P           - Tip load
% y           - Array of control paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tol1=cosd(d_angle_tol);   %Tolerance angle
    N=floor(size(y,2)/12)-1;
    Jac_neq=zeros(size(y,2),N+1);
    gneq=zeros(1,N+1);
    prevsol=[];
    for T=(1:N+1)
        U_cur=[y(12*T+1-6),y(12*T+2-6),y(12*T+3-6),y(12*T+4-6),y(12*T+5-6),y(12*T+6-6)];
        Tr_current=Trajectory(U_cur,P,prevsol); %TH1 TH2 L1 L2
        [align,jac_align]=Orientation_pursuit(Ori_tar,U_cur,Tr_current,P);
        gneq(T)=tol1 - align;
        Jac_neq(12*T-5:12*T,T)=-jac_align;
        prevsol=Tr_current;
    end
    gradneq=Jac_neq;
    geq=[];
    gradceq=[];