function [Output]=Optimize_stay_at_point(y_init)
%Run Optimize_stay_at_point([]) in command window to start the
%optimization.
%Can be run from the provided solution y_init if available, otherwise start
%from [].(Example 3 "Fixed Tip Position with Adjustable Orientation" in Chapter 10)
%Output: Array of collection of control parameters, control parameter vectors at the mesh points
% Run "Animate(output_file)" to get the relevant animations and plots,
% where output_file is the name of the optimization solution file (output of this program)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numberof mesh-intervals in time
T_steps=10;
%The system paramters of CTCR such as the intrinisic curvatures and stiffnesses of the
%tubes can be varied in Trajectory.m file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Specify penalizing weights
p1=5; % 0.5 is good % Penalizing weight of Regularization 
p2=0.0;             % Penalizing weight for Covered Volume 
p3=100.0;%200;      % Penalizing weight for deviation from prescribed path
p4=800; %200;       % Penalizing weight for deviation from Target at final time
% For prioritizing between the tip reaching the target position or reaching the
% target orientation.
mu=1.0;     %Weightage for path pursuit (position) along the path
nu=0.0;     %Weightage for path pursuit (orientation)along the path
mu_1=0.55;  %Weightage for position pursuit at the final time
nu_1=-1.0;  %Weightage for orientation pursuit at the final time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Specify tip load
P=[0.1,0.1,0.02];

%%Target orientation
%Ori=[-1.0,0.0,1.0];
Ori=[-0.1,0.0,1.0];
Ori=Ori/norm(Ori);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Specify initial configuration. (Configuration at t=0) 
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
USt=[0.5,0.5+2.1415,0.75+3.1415,0.4,0.6,0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Follow the Leader curve for a CTCR
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1};
Tr_current=Trajectory(UFL,[0,0,0],[]); %TH1 TH2 L1 L2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=Trajectory(USt,P,[]).y;
R_init=[S(17,end),S(18,end),S(19,end)];
Rtar=R_init;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if size(y_init,2)==0
    Yo=[];
    for T=(1:T_steps)
        Yo=[Yo,[0,0,0,0,0,0,USt(1)+0.0*rand(1)*pi,USt(2)+0.0*(rand(1))*pi,USt(3)+0.00*pi*rand(1),USt(4)+0.00*(rand(1)-1)*0.3,USt(5)+0.001*(rand(1)-1)*0.3,USt(6)+0.00*(rand(1)-1)*0.3]];
    end
else 
    Yo=y_init;
    Yo(1:6)=[];
    for T=(1:size(Yo,2))
        Yo(T)=Yo(T)+ 0.1*rand(1)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Linear constraints
Aa=[];
b=[];

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

 %% Lower and upper bounds on the lengths of the sections.
 lb(12*(T-1)+10)=0.001;
 lb(12*(T-1)+11)=0.001;
 lb(12*(T-1)+12)=0.001;

 ub(12*(T-1)+10)=0.9;
 ub(12*(T-1)+11)=0.9;
 ub(12*(T-1)+12)=0.9;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimization problem
options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'Diagnostics','on','MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',true,'StepTolerance',1e-12,'OptimalityTolerance',1e-06);
nonlcon=[];
Objective(Yo,Rtar,R_init,Ori,UFL,USt,P,Tr_current,TRF2,p1,p2,p3,p4,mu,nu,mu_1,nu_1);
fun3=@(U)Objective(U,Rtar,Rtar,Ori,UFL,USt,P,Tr_current,TRF2,p1,p2,p3,p4,mu,nu,mu_1,nu_1)

tic
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,Aa,b,Aeq,beq,lb,ub,nonlcon,options);
%Minimize the objective function
disp('Time taken for the optimization: ')
toc
Output={[USt,x],Rtar,P,UFL,Ori};
Objective(x,Rtar,R_init,Ori,UFL,USt,P,Tr_current,TRF2,p1,p2,p3,p4,mu,nu,mu_1,nu_1);

%Optimization framework. Overall Objective function
function [sum,jac] = Objective(y,Rtar,R_init,Ori,UFL,USt,P,TRF,TRF2,p1,p2,p3,p4,mu,nu,mu_1,nu_1) %{RF3,RF2,RF1,t3,t2,t1}
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
reach_along=  reach*h; 

jac(1:6,1) = p1*[y(7)*h;y(8)*h;y(9)*h;rp1*y(10)*h;rp1*y(11)*h;rp1*y(12)*h];
jac(7:12,1) = p3*h*jac_reach';

Tr_prev=Tr_current;
U_prev=U_current;
prevsol=Tr_current;
for T=(2:N)
    En = En + (y(12*(T-1) + 7)^2 + y(12*(T-1) + 8)^2 + y(12*(T-1)+ 9)^2 + rp1*(y(12*(T-1) + 10)^2 + y(12*(T-1) + 11)^2 + y(12*(T-1) + 12)^2) )*h/2;
    Rpath=R_init*((N-T)/N) + Rtar*(T/N);
    U_current=[y(12*(T-1)+13),y(12*(T-1)+14),y(12*(T-1)+15),y(12*(T-1)+16),y(12*(T-1)+17),y(12*(T-1)+18)];
    Tr_current=Trajectory(U_current,P,prevsol); %TH1 TH2 L1 L2
    [reach,jac_reach]=Reach_targetandorientation(Rpath,Ori,U_current,Tr_current,P,mu,nu);
    reach_along= reach_along + reach*h;
    jac(12*(T-1)+1:12*(T-1)+6,1) = p1*[y(12*(T-1) + 7)*h;y(12*(T-1) + 8)*h;y(12*(T-1) + 9)*h;rp1*y(12*(T-1) + 10)*h;rp1*y(12*(T-1) + 11)*h;rp1*y(12*(T-1) + 12)*h];
    jac(12*(T-1)+7:12*(T-1)+12,1) = p3*h*jac_reach';
    Tr_prev=Tr_current;
    U_prev=U_current;
    prevsol=Tr_current;
end
reach_along= reach_along - reach*h/2; 
[reach,jac_reach]=Reach_targetandorientation(Rpath,Ori,U_current,Tr_current,P,mu_1,nu_1);

sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) )/2 + p4*jac_reach';