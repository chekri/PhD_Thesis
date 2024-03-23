function [earea,Jac]=Deviation_from_FTL(U,Tr,UFL,TrF)
%Outputs:     the deviation from FTL objective (Chapter 9 in the thesis)and its gradient
%             with respect to the control variables U
%Inputs : U  -Control Paramters
%       : Tr -The configuration of the CTCR (Its shape (X,Y,Z) in 3D space)
%       : UFL-Follow the Leader control paramters for the given problem
%       : TrF-The "Follow the Leader" configuration of the CTCR for the given UFL paramters  (Its shape (X,Y,Z) in 3D space)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Arclength variables from one-tube section, two-tube section and three-tube
%sections
t1=TrF{6};
t2=TrF{5};
t3=TrF{4};

%Length control paramters
l1= U(4);
l2= U(5);
l3= U(6);
B0=[0,0];

%FTL Length contol paramters
L1f =UFL(4);
L2f =UFL(5);
L3f =UFL(6);

%Whole CTCR's arclength  for Follow the Leader controls     
Al=[L3f*t3', L3f + L2f*t2', L3f + L2f + L1f*t1'];           % Arclength Steps.. FTL Robot's config
N3=size(t3,1);
N2=size(t2,1);

% Remove the coupled term, it appears twice
del_index=N3+1;
Al(del_index)=[];% Remove the coupled term

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


e=1e-05; %Perturbation paramter for computing internal derivatives by finite differnces.
IC=Tr.y(31:48,1);
[p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',U);
l1=U(4);
l2=U(5);
l3=U(6);
l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];        % Arclength Steps.. Robot'd config
%Interpolate the position vectors of CTCR in FTL mode as a function of
%arclength
xs=interp1(Al,XF,l,'spline');
ys=interp1(Al,YF,l,'spline');
zs=interp1(Al,ZF,l,'spline');
%Array of position vectors of CTCR in FTL mode
AA=[xs;ys;zs];
%Array of arclengths of CTCR with the given control paramters
B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
% Trapezoidal sum between the 3D FTL_CTCR curve and Current CTCR curve. 
earea=enclosed_area1(AA,B,l);

%% Gradient Computations
% Same evaluation as above but with perturbed control paramters and initial
% conditions
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];

for i=[1:6]
    UU=U+diff(i,:); %perturb the control paramter
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',UU);
    l1=UU(4);
    l2=UU(5);
    l3=UU(6);
    %CTCR arclength
    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];  
    %Interpolate the position vectors of CTCR in FTL mode as a function of arclength
    xs=interp1(Al,XF,l,'spline');
    ys=interp1(Al,YF,l,'spline');
    zs=interp1(Al,ZF,l,'spline');

    %Array of position vectors of CTCR in FTL mode
    AA=[xs;ys;zs];
    %Array of arclengths of CTCR with the control paramters
    B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
    % Trapezoidal sum between the 3D FTL_CTCR curve and Current CTCR curve. 
    earea_f=enclosed_area1(AA,B,l);
    Jac(i)=(earea_f-earea)/e;  
    Dbdc(i,:)=[p2(end,16)/e,p3(end,18)/e];   
end


ICC_diff=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e,0,0];[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e]];
for i=[1:2]
    %Evaluate the IVP with the perturbed unknown variables at s=0 
    [y3,y2,y1,t3,t2,t1]=IVP_trajectory(IC'+ICC_diff(i,:),U);%Perturb boundary condition at the end s=0
    %Length Control Paramters
    l1=U(4);
    l2=U(5);
    l3=U(6);

    l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];     
    xs=interp1(Al,XF,l,'spline');
    ys=interp1(Al,YF,l,'spline');
    zs=interp1(Al,ZF,l,'spline');

    AA=[xs;ys;zs];
    B=[y3(:,1)',y2(:,1)',y1(:,1)';y3(:,2)',y2(:,2)',y1(:,2)';y3(:,3)',y2(:,3)',y1(:,3)'];
    earea_f=enclosed_area1(AA,B,l);
    Jac_y0(i)=(earea_f-earea)/e;  
    Dbdy0(i,:)=[y2(end,16)/e,y3(end,18)/e];   
end
%Evaluate approximate gradient (see "Using IVPs" section in Chapter 9 of the thesis)
Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));


function area=enclosed_area1(Rf,R,l) 
%Sum trapezoidas area between two equally meshed 3D curves.
area=trapz(l,vecnorm((Rf-R)));






