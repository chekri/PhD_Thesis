function [Rt,Jac]=Reach_target(Rtar,U,Tr)
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
e=1e-04;
%e=5e-04;
diff=[[e,0,0,0,0,0];[0,e,0,0,0,0];[0,0,e,0,0,0];[0,0,0,e,0,0];[0,0,0,0,e,0];[0,0,0,0,0,e]];
for i=[1:6]
    UU=U+diff(i,:);
    [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',UU);
    
    
   % plot3(p1(:,1),p1(:,2),p1(:,3),'r');
   % hold on;
   % plot3(p2(:,1),p2(:,2),p2(:,3),'g');
   % plot3(p3(:,1),p3(:,2),p3(:,3),'b');
   % axis equal
    
    
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

Jac=Jac - (Jac_y0*(Dbdy0'\Dbdc'));

