function [earea,Jac]=Deviation_from_FTL(U,Tr,UFL,TrF)
% Output:   This function evaluates the deviation from the TFL objective,
%           i.e., the sum of the distances between the CTCR configuration 
%           and the FTL  shape and its its gradient with respect to control 
%           parameters using Finite difference approach ("Using IVPs" section 
%           in chapter 9 of the thesis)
%           The gradient components corresponding to Length parameters are
%           set to zero in this particular case to restrict the function's
%           dependence on them.
%Arguments: U  - Control Paramters
%           Tr - Array of position vectors of current CTCR.
%           UFL- "Follow the Leader" Control paramters
%           TrF- Array of position vectors of CTCR with the corresponding
%                UFL paramters (FTL shape) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Arclength paramters
t1=TrF{6};
t2=TrF{5};
t3=TrF{4};

%Length control paramters
l1= U(4);
l2= U(5);
l3= U(6);

%FTL length control paramters
B0=[0,0];
L1f =UFL(4);
L2f =UFL(5);
L3f =UFL(6);


%l =[l3*Tr.x,l3 + l2*Tr.x, l3 + l2 + l1*Tr.x];         % Arclength Steps.. Robot'd config
Al=[L3f*t3', L3f + L2f*t2', L3f + L2f + L1f*t1'];           % Arclength Steps.. FTL Robot's config
N3=size(t3,1);
N2=size(t2,1);
% Remove the coupled term
del_index=N3+1;
Al(del_index)=[];              % Remove the coupled term

XF=[(TrF{1}(:,1))',(TrF{2}(:,1))',(TrF{3}(:,1))'];
YF=[(TrF{1}(:,2))',(TrF{2}(:,2))',(TrF{3}(:,2))'];
ZF=[(TrF{1}(:,3))',(TrF{2}(:,3))',(TrF{3}(:,3))'];

XF(del_index)=[];   
YF(del_index)=[];
ZF(del_index)=[]; % Remove the coupled part.. repetetion

del_index=N3+N2;
Al(del_index)=[];
XF(del_index)=[];  
YF(del_index)=[];
ZF(del_index)=[];

% Perturbation parameter for finite difference computations for
% approxiamting internal derivatives
e=1e-05;

IC=Tr.y(31:48,1);
[p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',U);
l1=U(4);
l2=U(5);
l3=U(6);
% Arclength Steps.. Robot'd config
l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];  
%Interpolate CTCR position vector along the length of the CTCR
xs=interp1(Al,XF,l,'spline');
ys=interp1(Al,YF,l,'spline');
zs=interp1(Al,ZF,l,'spline');

%Position vectors along the meshed arclength for evluating the trpaezoidal
%area between the corresponding points.
AA=[xs;ys;zs];
B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
%sum of Trapezooidal areas
earea=enclosed_area1(AA,B,l);
%% Gradient evaluations. repeat the same procedure using perturbations in control parameters and intial conditions
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
    B=[y3(:,1)',y2(:,1)',y1(:,1)';y3(:,2)',y2(:,2)',y1(:,2)';y3(:,3)',y2(:,3)',y1(:,3)'];
    earea_f=enclosed_area1(AA,B,l);
    Jac_y0(i)=(earea_f-earea)/e;  
    Dbdy0(i,:)=[y2(end,16)/e,y3(end,18)/e];   
end
%Evaluation of gradients (Using IVPs section in Chapter 9 of thesis)
Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));
Jac(4:6)=0; %For_restricitng the length controls.


function area=enclosed_area1(Rf,R,l) 
%Sum area of the trapezoids from the veritces on Rf and R
area=trapz(l,vecnorm((Rf-R)));






