function [earea,Jac,Jac1]=Sweep_area_non(Tr_Current,Tr_Prev,U_Current,U_Prev,Load)
% Output:   This function evaluates the approximate area (sum of the Trapezoidal areas),
%           traced between the Current configuration and the previous
%           Configuration and its its gradient with respect to control 
%           parameters using Finite difference approach ("Using IVPs" section 
%           in chapter 9 of the thesis)
%           The gradient components corresponding to Length parameters are
%           set to zero in this particular case to restrict the function's
%           dependence on them.
%Arguments: U_Current  - Control Paramters at current time step
%           Tr_Current - Array of position vectors of CTCR at current time step.
%           U_Prev     - Control Paramters at previous time step
%           Tr_Prev    - Array of position vectors of CTCR at previous time step
%           Load       - Tip Load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC_Prev=Tr_Prev.y(31:48,1);        % State variables of Previous CTCR configuration at the end s=0
IC_Current=Tr_Current.y(31:48,1);  % State variables of current CTCR configuration at the end s=0

Steps=50; %Spatial discretization
e=1e-04;  % perturbation parameter for evaluating internal derivatives through finite differences

%Array of arclengths of CTCR
Arc_c=[U_Current(6)*(Tr_Current.x), U_Current(6)+ U_Current(5)*(Tr_Current.x), U_Current(6)+ U_Current(5)+ U_Current(4)*(Tr_Current.x)];
Arc_p=[U_Prev(6)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)*(Tr_Prev.x), U_Prev(6)+ U_Prev(5)+ U_Prev(4)*(Tr_Prev.x)];

%Array of Position vectors of curent CTCR. Appended along the sections
Cur_x=[Tr_Current.y(31,:),Tr_Current.y(1,:),Tr_Current.y(17,:)];
Cur_y=[Tr_Current.y(32,:),Tr_Current.y(2,:),Tr_Current.y(18,:)];
Cur_z=[Tr_Current.y(33,:),Tr_Current.y(3,:),Tr_Current.y(19,:)];

%Remove the repeated term. The term at the boundary is repeated twice when
%coupled. i.e., [.......term(s=l),term(s=0),.....]. Remove one term.
index=size(Tr_Current.x,2);
Arc_c(2*index)=[];
Arc_c(index+1)=[];
Cur_x(2*index)=[];
Cur_x(index+1)=[];
Cur_y(2*index)=[];
Cur_y(index+1)=[];
Cur_z(2*index)=[];
Cur_z(index+1)=[];

%Array of Position vectors of curent CTCR. Appended along the sections
Prv_x=[Tr_Prev.y(31,:),Tr_Prev.y(1,:),Tr_Prev.y(17,:)];
Prv_y=[Tr_Prev.y(32,:),Tr_Prev.y(2,:),Tr_Prev.y(18,:)];
Prv_z=[Tr_Prev.y(33,:),Tr_Prev.y(3,:),Tr_Prev.y(19,:)];


%Remove the repeated term. The term at the boundary is repeated twice when
%coupled. i.e., [.......term(s=l),term(s=0),.....]. Remove one term.
index=size(Tr_Prev.x,2);
Arc_p(2*index)=[];
Arc_p(index+1)=[];

Prv_x(index+1)=[];
Prv_x(2*index)=[];
Prv_y(index+1)=[];
Prv_y(2*index)=[];
Prv_z(index+1)=[];
Prv_z(2*index)=[];

%Discretize the total length of the CTCR
Cur_Arc=linspace(0,U_Current(4)+U_Current(5)+U_Current(6),Steps);
Prv_Arc=linspace(0,U_Prev(4)+U_Prev(5)+U_Prev(6),Steps);

%% Interpolate the Array of position vectors using a spline along the length of CTCR
C_x=spline(Arc_c,Cur_x,Cur_Arc);
C_y=spline(Arc_c,Cur_y,Cur_Arc);
C_z=spline(Arc_c,Cur_z,Cur_Arc);

P_x=spline(Arc_p,Prv_x,Prv_Arc);
P_y=spline(Arc_p,Prv_y,Prv_Arc);
P_z=spline(Arc_p,Prv_z,Prv_Arc);

% Evaluate the discreet area (sum of trapezoid areas) between the both 3D CTCR shapes.
% The vertices of the trapezoid are the mesh points along which CTCR is
% discretized.
earea=enclosed_area1([C_x;C_y;C_z],[P_x;P_y;P_z]); 

%% Gradient computation
% Repeat the same procedure as above for internal derivative evaluations. But, using perturbations in control
% parameters U and initial conditions
[prv3,prv2,prv1,prv_t3,prv_t2,prv_t1]=IVP_trajectory(IC_Prev',U_Prev);
l1=U_Prev(4);
l2=U_Prev(5);
l3=U_Prev(6);
l =[l3*prv_t3',l3 + l2*prv_t2', l3 + l2 + l1*prv_t1'];
B=[prv3(:,1)',prv2(:,1)',prv1(:,1)';prv3(:,2)',prv2(:,2)',prv1(:,2)';prv3(:,3)',prv2(:,3)',prv1(:,3)'];
l(size(prv_t3,1)+ size(prv_t2,1))=[];
B(:,size(prv_t3,1)+ size(prv_t2,1))=[];
l(size(prv_t3,1))=[];
B(:,size(prv_t3,1))=[];
Prv_r=spline(l,B,Prv_Arc); %Previous configuration on equally spaced given intervals
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U_Current+diff(i,:);
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC_Current',UU);
    l1=UU(4);
    l2=UU(5);
    l3=UU(6);
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];       % Arclength Steps.. Robot'd config
    Cur_Arc=linspace(0,l1+l2+l3,Steps);
    B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
    l(size(p3,1)+ size(p2,1))=[];
    B(:,size(p3,1)+ size(p2,1))=[];
    l(size(p3,1))=[];
    B(:,size(p3,1))=[];
    Cur_r=spline(l,B,Cur_Arc); %Current configuration on equally spaced intervals
    earea_f=enclosed_area1(Cur_r,Prv_r);
    Jac(i)=(earea_f-earea)/e;  
    Dbdc(i,:)=[p2(end,16)/e,p3(end,18)/e];   
end
%figure(3)
l1=U_Current(4);
l2=U_Current(5);
l3=U_Current(6);
      % Arclength Steps.. Robot'd config
Cur_Arc=linspace(0,l1+l2+l3,Steps);
ICC_diff=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e,0,0];[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e]];
for i=[1:2]
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC_Current'+ICC_diff(i,:),U_Current);%MU1
    l1=U_Current(4);
    l2=U_Current(5);
    l3=U_Current(6);
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];       % Arclength Steps.. Robot'd config
    
    B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
    l(size(p3,1)+ size(p2,1))=[];
    B(:,size(p3,1)+ size(p2,1))=[];
    l(size(p3,1))=[];
    B(:,size(p3,1))=[];
    Cur_r=spline(l,B,Cur_Arc);
    earea_f=enclosed_area1(Cur_r,Prv_r);
    Jac_y0(i)=(earea_f-earea)/e;  
    Dbdy0(i,:)=[p2(end,16)/e,p3(end,18)/e];   
    %pause(1)
end
Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));
Jac(4:6)=0;
l1=U_Current(4);
l2=U_Current(5);
l3=U_Current(6);
[cur3,cur2,cur1,cur_t3,cur_t2,cur_t1]=IVP_trajectory(IC_Current',U_Current);
B=[cur3(:,1)',cur2(:,1)',cur1(:,1)';cur3(:,2)',cur2(:,2)',cur1(:,2)';cur3(:,3)',cur2(:,3)',cur1(:,3)'];
l =[l3*cur_t3',l3 + l2*cur_t2', l3 + l2 + l1*cur_t1'];   
l(size(cur_t3,1)+ size(cur_t2,1))=[];
B(:,size(cur_t3,1)+ size(cur_t2,1))=[];
l(size(cur_t3,1))=[];
B(:,size(cur_t3,1))=[];

Cur_r=spline(l,B,Cur_Arc); %Current configuration positions at given equally spaced intervals
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];

for i=[1:6]
    UU=U_Prev+diff(i,:); %Previous configuration is varied
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC_Prev',UU);
    l1=UU(4);
    l2=UU(5);
    l3=UU(6);
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];       % Arclength Steps.. 
    Prv_Arc=linspace(0,l1+l2+l3,Steps);
    B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
    l(size(p3,1)+ size(p2,1))=[];
    B(:,size(p3,1)+ size(p2,1))=[];
    l(size(p3,1))=[];
    B(:,size(p3,1))=[];
    Prv_r=spline(l,B,Prv_Arc);
    earea_f=enclosed_area1(Cur_r,Prv_r);
    
    Jac1(i)=(earea_f-earea)/e;  
    Dbdc(i,:)=[p2(end,16)/e,p3(end,18)/e];  

end

ICC_diff=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e,0,0];[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e]];
l1=U_Prev(4);
l2=U_Prev(5);
l3=U_Prev(6);
Prv_Arc=linspace(0,l1+l2+l3,Steps);
for i=[1:2]
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC_Prev'+ICC_diff(i,:),U_Prev);
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];       % Arclength Steps.. Robot'd config
    B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
    l(size(p3,1)+ size(p2,1))=[];
    B(:,size(p3,1)+ size(p2,1))=[];
    l(size(p3,1))=[];
    B(:,size(p3,1))=[];
    Prv_r=spline(l,B,Prv_Arc);
    earea_f=enclosed_area1(Cur_r,Prv_r);
    Jac_y0(i)=(earea_f-earea)/e;  
    Dbdy0(i,:)=[p2(end,16)/e,p3(end,18)/e];  

end
%Evaluate the gradient using finite differnces (Using IVPs section in Chapter 9 of the thesis)
Jac1=Jac1 - (Jac_y0*(Dbdy0'\Dbdc'));
% For_restricting the length controls.
Jac1(4:6)=0;

function area=enclosed_area1(R_prev,R_cur) 
area=0;
for i=[1:size(R_prev,2)-1]
    (R_cur(:,i) - R_cur(:,i+1));
    a1=cross((R_cur(:,i) - R_cur(:,i+1)),(R_prev(:,i+1) - R_cur(:,i+1)))/2;
    a2=cross((R_prev(:,i) - R_prev(:,i+1)),(R_prev(:,i) - R_cur(:,i)))/2;
    area=area+norm(a1)+norm(a2);
end


