function [Rt,Jac]=Reach_target2(Rtar,U,Tr,Load)
mu=1;
B0=[0,0];
xtarget=Rtar(1);
ytarget=Rtar(2);
ztarget=Rtar(3);

Rt=norm([Tr.y(17,end),Tr.y(18,end),Tr.y(19,end)]-[xtarget,ytarget,ztarget]);


%fig1=figure(1);
%plot3(Tr.y(17,:),Tr.y(18,:),Tr.y(19,:),'r');
%hold on;
%plot3(Tr.y(1,:),Tr.y(2,:),Tr.y(3,:),'g');
%plot3(Tr.y(31,:),Tr.y(32,:),Tr.y(33,:),'b');
%grid on;
%axis equal;

IC=Tr.y(31:48,1);
e=5e-05;
%e=5e-04;
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U+diff(i,:);
    Tr_UU=Trajectory(UU,Load,Tr);
    Rt2=norm([Tr_UU.y(17,end),Tr_UU.y(18,end),Tr_UU.y(19,end)]-[xtarget,ytarget,ztarget]);
    Jac(i)=(Rt2-Rt)/e;
    
end



