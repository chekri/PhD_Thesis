function [Rt,Jac]=Reach_target(Rtar,Tr,J,Fc)
%This "Target-pursuit" quantifies the final the distance between the end
%point and the specified target point.
%Target point
xtarget=Rtar(1);
ytarget=Rtar(2);
ztarget=Rtar(3);

%Distance between Target point and the current end point.
Rt=norm([Tr.y(17,end),Tr.y(18,end),Tr.y(19,end)]-[xtarget,ytarget,ztarget]);

Jy=zeros(size(J,2),1);
Jy(end-31)=(Tr.y(17,end)-xtarget)/Rt;
Jy(end-30)=(Tr.y(18,end)-ytarget)/Rt;
Jy(end-29)=(Tr.y(19,end)-ztarget)/Rt;

lam = linsolve(J',Jy);
Jac=-lam'*Fc;