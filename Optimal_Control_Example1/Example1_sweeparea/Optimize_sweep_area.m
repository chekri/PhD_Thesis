function [Output]=Optimize_sweep_area(y_init)
%Run Optimize_sweep_area([]) in command window to start the
%optimization. This function utilizes swept area by CTCR as a Covered
%Volume measure.(Chapter 9 of the thesis)
%Output: Array of collection of control parameters, control paramter vectors at the mesh points
%%% Target point
%Specify Target point
Rtar=[0.4,-0.0,1.0];
%Specify Tip Load
P=[0.1,0.15,0.12];
%P=[0.,0.,0.];
%Time steps
T_steps=10;
%Weights of objective functions
pp=1;
p1=5*pp; % 0.5 is good
p2=50.0*pp;%20.0; 750
p3=0;
p4=1000*pp;
rp1=1;
%Specify Follow the Leader controls to reach initial point
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
%Initial controls
%Slight perturbation from Follow-the Leader
USt=UFL + [0,0.1,0.2,0.0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference Follow the leader Curve
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TRF2={RF3,RF2,RF1,t3,t2,t1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %TH1 TH2 L1 L2
% The state at t=0
Tr_initial=Trajectory(USt,P,[]);
S=Tr_initial.y;
%The Robot tip at t=0
R_init=[S(17,end),S(18,end),S(19,end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If initial solution is not provided, take the configuration at initial
%time as configrations at all time steps.
%If it is provided, take it as inital guess.
if size(y_init,2)==0
    Yo=[];
    for T=(1:T_steps)
        Yo=[Yo,[0,0,0,0,0,0,USt(1)+0.0*rand(1)*pi,USt(2)+0.0*(rand(1))*pi,USt(3)+0.00*pi*rand(1),USt(4)+0.00*(rand(1)-1)*0.3,USt(5)+0.00*(rand(1)-1)*0.3,USt(6)+0.00*(rand(1)-1)*0.3]];
           end
else 
    Yo=y_init;
    Yo(1:6)=[];
    for T=(1:size(Yo,2))
        Yo(T)=Yo(T)
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

%Optimization controls
options  = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'Diagnostics','on','MaxFunctionEvaluations',5000,'SpecifyConstraintGradient',true,'StepTolerance',1e-10,'OptimalityTolerance',1e-04);
nonlcon=[];
Objective(Yo,Rtar,R_init,UFL,USt,P,Tr_initial,TRF2,p1,p2,p3,p4)
fun3=@(U)Objective(U,Rtar,R_init,UFL,USt,P,Tr_initial,TRF2,p1,p2,p3,p4)
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve the optimization problem.
[x,fval,exitflag,output,lambda,grad,Hess] = fmincon(fun3,Yo,[],[],Aeq,beq,lb,ub,nonlcon,options);
disp('Time taken for the optimization')
toc
Output={[USt,x],Rtar,P,UFL};
Objective(x,Rtar,R_init,UFL,USt,P,Tr_initial,TRF2,p1,p2,p3,p4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sum,jac] = Objective(y,Rtar,R_init,UFL,USt,Load,TRF,TRF2,p1,p2,p3,p4,rp1) %{RF3,RF2,RF1,t3,t2,t1}
y=[USt,y]; 
En=0;
rp1=1;
N=floor(size(y,2)/12);
h=1/N;
Tr_prev=TRF;
U_prev=USt;
T=1;
Rpath=R_init*((N-T)/N) + Rtar*(T/N);
jac(7:12,1)=[0,0,0,0,0,0];
En = En + (y(7)^2 + y(8)^2 + y(9)^2 + rp1*(y(10)^2 + y(11)^2 + y(12)^2) )*h/2;
prevsol=Tr_prev;
U_current=[y(12+1),y(12+2),y(12+3),y(12+4),y(12+5),y(12+6)];
Tr_current=Trajectory(U_current,Load,prevsol); %TH1 TH2 L1 L2
[reach,jac_reach]=Reach_target(Rpath,U_current,Tr_current);
[sweep,jac_sweep1,jac_sweep2]=Sweep_area(Tr_current,Tr_prev,U_current,U_prev,Load);
sweep_area= sweep; 
reach_along=  reach*h; 
jac(1:6,1) = p1*[y(7)*h;y( 8)*h;y(9)*h;rp1*y(10)*h;rp1*y(11)*h;rp1*y( 12)*h];
jac(7:12,1) = p2*jac_sweep1'+ p3*h*jac_reach';
Tr_prev=Tr_current;
U_prev=U_current;
prevsol=Tr_current;
for T=(2:N)
    En = En + (y(12*(T-1) + 7)^2 + y(12*(T-1) + 8)^2 + y(12*(T-1) + 9)^2 + rp1*(y(12*(T-1) + 10)^2 + y(12*(T-1) + 11)^2 + y(12*(T-1) + 12)^2) )*h/2;
    Rpath=R_init*((N-T)/N) + Rtar*(T/N);
    U_current=[y(12*T+1),y(12*T+2),y(12*T+3),y(12*T+4),y(12*T+5),y(12*T+6)];
    Tr_current=Trajectory(U_current,Load,prevsol); %TH1 TH2 L1 L2
    [reach,jac_reach] = Reach_target(Rpath,U_current,Tr_current);
    [sweep,jac_sweep1,jac_sweep2] = Sweep_area(Tr_current,Tr_prev,U_current,U_prev,Load);
    sweep_area=sweep_area + sweep; 
    reach_along= reach_along + reach*h;
    
    jac(12*(T-1)+1:12*(T-1)+6,1) = p1*[y(12*(T-1) + 7)*h;y(12*(T-1) + 8)*h;y(12*(T-1) + 9)*h;rp1*y(12*(T-1) + 10)*h;rp1*y(12*(T-1) + 11)*h;rp1*y(12*(T-1) + 12)*h];
    jac(12*(T-1)+7:12*(T-1)+12,1) = p2*jac_sweep1'+ p3*h*jac_reach';
    jac(12*(T-2)+7:12*(T-2)+12,1) =  jac(12*(T-1)+7:12*(T-1)+12,1) + p2*jac_sweep2';
    
    Tr_prev=Tr_current;
    U_prev=U_current;
    prevsol=Tr_current;
end
[reach,jac_reach]=Reach_target2(Rtar,U_current,Tr_current,Load);
reach_along= reach_along - reach*h/2; 
sum = p1*En + p2*sweep_area + p3*reach_along + p4*reach;
jac(12*(N-1)+7:12*(N-1)+12,1)=(jac(12*(N-1)+7:12*(N-1)+12,1) ) + p4*jac_reach';
