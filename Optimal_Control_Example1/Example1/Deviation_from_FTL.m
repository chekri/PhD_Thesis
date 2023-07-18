function [earea,Jac]=Deviation_from_FTL(U,Tr,UFL,TrF)
%TrF=%FTL
%Tr=%Normal Function
% perturbation for finite differences
e=5e-04; 

%Length intervals
t1=TrF{6};
t2=TrF{5};
t3=TrF{4};
l1= U(4);
l2= U(5);
l3= U(6);
B0=[0,0];

%Length controls
L1f =UFL(4)+0.99;%0.4
L2f =UFL(5);%0.75
L3f =UFL(6);

% Arclength Steps.. FTL Robot's config
Al=[L3f*t3', L3f + L2f*t2', L3f + L2f + L1f*t1'];   

% Remove the coupled term. Because the two terms correspond to the same point
N3=size(t3,1);
N2=size(t2,1); 
del_index=N3+1;
Al(del_index)=[];

XF=[(TrF{1}(:,1))',(TrF{2}(:,1))',(TrF{3}(:,1))'];
YF=[(TrF{1}(:,2))',(TrF{2}(:,2))',(TrF{3}(:,2))'];
ZF=[(TrF{1}(:,3))',(TrF{2}(:,3))',(TrF{3}(:,3))'];

XF(del_index)=[];   
YF(del_index)=[];
ZF(del_index)=[];

% Remove the coupled term. Because the two terms correspond to the same point
del_index=N3+N2;
Al(del_index)=[];
XF(del_index)=[];  
YF(del_index)=[];
ZF(del_index)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial state (state at s=0) from the solutions of CTCR-bvp solutions.
%Useful for computing gradients
IC=Tr.y(31:48,1);
%Solutions of the IVP with the prescribed control parameters adn intial
%state
[p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',U);
l1=U(4);
l2=U(5);
l3=U(6);
l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];       % Arclength Steps alson the Robot.

%Interpolation of the robots FTL's position vector in terms of its arclength.
xs=interp1(Al,XF,l,'spline');
ys=interp1(Al,YF,l,'spline');
zs=interp1(Al,ZF,l,'spline');

%FTL's position vectors along the arclength projection of the current robot's configuration. 
AA=[xs;ys;zs];
B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
%The area of the devaition from the FTL shape 
earea=enclosed_area1(AA,B,l);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computation of gradient of the objective function
%Vector of perturbtation vectors
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];

for i=[1:6]
    UU=U+diff(i,:);
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',UU);
    l1=UU(4);
    l2=UU(5);
    l3=UU(6);
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];        % Arclength Steps.. Robot'd config
    xs=interp1(Al,XF,l,'spline');
    ys=interp1(Al,YF,l,'spline');
    zs=interp1(Al,ZF,l,'spline');

    AA=[xs;ys;zs];
    %B=[[Tr.y(31,:),Tr.y(1,:),Tr.y(17,:)];[Tr.y(32,:),Tr.y(2,:),Tr.y(18,:)];[Tr.y(33,:),Tr.y(3,:),Tr.y(19,:)]];
    B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
    earea_f=enclosed_area1(AA,B,l);
    Jac(i)=(earea_f-earea)/e;  
    Dbdc(i,:)=[p2(end,16)/e,p3(end,18)/e];   
end

ICC_diff=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e,0,0];[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e]];
for i=[1:2]
    [y3,y2,y1,t3,t2,t1]=IVP_trajectory(IC'+ICC_diff(i,:),U);%MU1
    l1=U(4);
    l2=U(5);
    l3=U(6);
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];        % Arclength Steps.. Robot'd config
    xs=interp1(Al,XF,l,'spline');
    ys=interp1(Al,YF,l,'spline');
    zs=interp1(Al,ZF,l,'spline');

    AA=[xs;ys;zs];
    %B=[[Tr.y(31,:),Tr.y(1,:),Tr.y(17,:)];[Tr.y(32,:),Tr.y(2,:),Tr.y(18,:)];[Tr.y(33,:),Tr.y(3,:),Tr.y(19,:)]];
    B=[y3(:,1)',y2(:,1)',y1(:,1)';y3(:,2)',y2(:,2)',y1(:,2)';y3(:,3)',y2(:,3)',y1(:,3)'];
    earea_f=enclosed_area1(AA,B,l);
    Jac_y0(i)=(earea_f-earea)/e;  
    Dbdy0(i,:)=[y2(end,16)/e,y3(end,18)/e];   
end
%Approximation of gradient using IVP finite diferences as described in the
%thesis.
Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));
%Jac(4:6)=0;

function area=enclosed_area(Rf,R,l) 
area=trapz(vecnorm((Rf-R)));

function area=enclosed_area1(Rf,R,l) 
%Computes the distance between the sepcified position vectors and sums over
%it.
area=trapz(l,vecnorm((Rf-R)));






