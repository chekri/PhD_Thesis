function [Rt,Jac]=Reach_target(Rtar,U,Tr)
% Output:    This function evalauates the distance betwwen the cuurent tip 
%            position and the Target point Rtar and its gradient with
%            respect to control paramters using Finite difference approach
%            ("Using IVPs" section in chapter 9 of the thesis)
%Arguments:  Rtar- Target point vector
%            U   - Control paramters
%            Tr  - Configuration of the CTCR corresponding to U. Only its
%            boundary value is needed for computing the distance. However,
%            additional details are required for computing gradients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=1;
B0=[0,0];
xtarget=Rtar(1);
ytarget=Rtar(2);
ztarget=Rtar(3);
%Distance between Rtar and CTCR tip
Rt=norm([Tr.y(17,end),Tr.y(18,end),Tr.y(19,end)]-[xtarget,ytarget,ztarget]);

%CTCR configuration at the end s=0, for evaluating IVP
IC=Tr.y(31:48,1);

IC=Tr.y(31:48,1);
%Perturbation parameter for computing internal gradients through finite
%differences.
e=5e-04;

diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U+diff(i,:);
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',UU);    
    tar_diff=norm([p1(end,1),p1(end,2),p1(end,3)] - [xtarget,ytarget,ztarget]);
    Jac(i)=(tar_diff-Rt)/e;
    Dbdc(i,:)=[p2(end,16)/e,p3(end,18)/e];   
end

[y13,y12,y11,t30,t20,t10]=IVP_trajectory(IC'+[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e,0,0],U);%MU1
Dbdy0(1,:)=([y12(end,16),y13(end,18)]-B0 )/e;

tar_diff=norm([y11(end,1),y11(end,2),y11(end,3)] - [xtarget,ytarget,ztarget]);
Jac_y0(1)=(tar_diff-Rt)/e;

[y23,y22,y21,t30,t20,t10]=IVP_trajectory(IC'+[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,e],U);%Mu2
Dbdy0(2,:)=([y22(end,16),y23(end,18)]- B0)/e;

tar_diff=norm([y21(end,1),y21(end,2),y21(end,3)] - [xtarget,ytarget,ztarget]) ;
Jac_y0(2)=(tar_diff-Rt)/e;

%Evaluation of gradients (Using IVPs section in Chapter 9 of thesis)
Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));

