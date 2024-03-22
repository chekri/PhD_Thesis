function [Rt,Jac]=Reach_target2(Rtar,U,Tr,Load)
% Output:    This function evalauates the distance betwwen the cuurent tip 
%            position and the Target point Rtar and its gradient with
%            respect to control paramters using Finite difference approach
%            ("Using BVPs" section in chapter 9 of the thesis)
%Arguments:  Rtar- Target point vector
%            U   - Control paramters
%            Tr  - Configuration of the CTCR corresponding to U. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=1;
B0=[0,0];
xtarget=Rtar(1);
ytarget=Rtar(2);
ztarget=Rtar(3);
%Distance between Rtar and CTCR tip
Rt=norm([Tr.y(17,end),Tr.y(18,end),Tr.y(19,end)]-[xtarget,ytarget,ztarget]);

%Perturbation parameter for computing internal gradients through finite
%differences.
e=5e-05;
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U+diff(i,:);
    Tr_UU=Trajectory(UU,Load,Tr);
    Rt2=norm([Tr_UU.y(17,end),Tr_UU.y(18,end),Tr_UU.y(19,end)]-[xtarget,ytarget,ztarget]);
    Jac(i)=(Rt2-Rt)/e;   
end



