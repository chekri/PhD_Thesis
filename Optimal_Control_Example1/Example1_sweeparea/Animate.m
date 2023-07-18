function Animate(y5,P)
%[y5,y10]=OptimalControlwithCollocation_Jacabian()
%Rtar=[-0.20,-0.55,0.76];
%R_init= [-0.17506,0.320572,1.44679]
R_tar=[-0.2,0.44,1.1]
R_tar=[0.16,0.167,0.767]
Rtar=[0.16,0.167,0.767]
%R_tar=[-0.2,0.44,1.1]
%R_tar=[-0.053,0.175,0.933]
%R_tar=[0.851419,-0.0220887,1.25805];
%Robs=[0.07056,-0.262,0.6162];
N=floor(size(y5,2)/12)-1
sol1=Trajectory([y5(1),y5(2),y5(3),y5(4),y5(5),y5(6)],P,[]);
fig=figure(1)
plot3(sol1.y(32,:),sol1.y(33,:),sol1.y(31,:),'k','LineWidth',1);
hold on
grid on
plot3(sol1.y(2,:),sol1.y(3,:),sol1.y(1,:),'r','LineWidth',1);
plot3(sol1.y(16+2,:),sol1.y(16+3,:),sol1.y(16+1,:),'b');
axis equal
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
scatter3(sol1.y(16+2,end),sol1.y(16+3,end),sol1.y(16+1,end),'b*');
text(sol1.y(16+2,end),sol1.y(16+3,end),sol1.y(16+1,end),'\rightarrow Initial Point','FontSize',12)

scatter3(R_tar(2),R_tar(3),R_tar(1),'k*')
text(R_tar(2),R_tar(3),R_tar(1)+0.01,'\rightarrow Target','FontSize',12)
%scatter3(R_init(2),R_init(3),R_init(1),'k*')
%text(R_init(2),R_init(3),R_init(1),'\leftarrow Initial point')
%scatter3(Robs(1),Robs(2),Robs(3),'r*')
hold on
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4,0.6,0.5];
sol1=Trajectory(UFL,P,[]);
%plot3(sol1.y(32,:),sol1.y(33,:),sol1.y(31,:),'k--','LineWidth',2);
%alpha(0.4)
hold on
%plot3(sol1.y(2,:),sol1.y(3,:),sol1.y(1,:),'k--','LineWidth',2);
%alpha(0.4)
%plot3(sol1.y(16+2,:),sol1.y(16+3,:),sol1.y(16+1,:),'k--','Linewidth',2);
%alpha(0.4)



scatter3(R_tar(2),R_tar(3),R_tar(1),'k*')
text(R_tar(2),R_tar(3),R_tar(1),'\leftarrow Target','FontSize',12)
%scatter3(R_init(2),R_init(3),R_init(1),'k*')
%text(R_init(2),R_init(3),R_init(1),'\leftarrow Initial point')
%scatter3(Robs(1),Robs(2),Robs(3),'r*')
hold on
%plot3(sol1.y(16+2,:),sol1.y(16+3,:),sol1.y(16+1,:),'b','LineWidth',2);

axis equal

tipx=[];
tipy=[];
tipz=[];
theta=[];
for k=[0:N]
    12*k+1;
    sol=Trajectory([y5(12*k+1),y5(12*k+2),y5(12*k+3),y5(12*k+4),y5(12*k+5),y5(12*k+6)],P,[]);
    qend1=sol.y(20,end);
    qend2=sol.y(21,end);
    qend3=sol.y(22,end);
    qend4=sol.y(23,end);
    d31=2*(qend1*qend3 + qend2*qend4);
    d32=2*(qend2*qend3 - qend1*qend4);
    d33=-qend1*qend1 - qend2*qend2 + qend3*qend3 + qend4*qend4;
    
    if k==0
        disp('Flah')
        di31=d31;
        di32=d32;
        di33=d33;
    end
    %acos(dot([d31,d32,d33],[di31,di32,di33])/());
    theta=[theta,acosd(dot([d31,d32,d33],[di31,di32,di33])/(norm([d31,d32,d33])*norm([di31,di32,di33])))];
    %plot3(sol1.y(1,:),sol1.y(2,:),sol1.y(3,:),'k');
    hold on
    %plot3(sol1.y(16+1,:),sol1.y(16+2,:),sol1.y(16+3,:),'k');
    
    plot3(sol.y(2,:),sol.y(3,:),sol.y(1,:),'r:');
    plot3(sol.y(16+2,:),sol.y(16+3,:),sol.y(16+1,:),'b:');
    plot3(sol.y(32,:),sol.y(33,:),sol.y(31,:),'k:');
    
     if k==0
        plot3(sol.y(2,:),sol.y(3,:),sol.y(1,:),'r-','LineWidth',3);
        plot3(sol.y(16+2,:),sol.y(16+3,:),sol.y(16+1,:),'b-','LineWidth',3);
        plot3(sol.y(32,:),sol.y(33,:),sol.y(31,:),'k-','LineWidth',3);
        %plot3([sol.y(18,end),sol.y(18,end)+0.1*d32], [sol.y(19,end),sol.y(19,end)+0.1*d33],[sol.y(17,end),sol.y(17,end)+ 0.1*d31],'c-')
     end
     
      if rem(k,5)==10
         k/5
        plot3(sol.y(2,:),sol.y(3,:),sol.y(1,:),'r--','LineWidth',1);
        plot3(sol.y(16+2,:),sol.y(16+3,:),sol.y(16+1,:),'b--','LineWidth',1);
        plot3(sol.y(32,:),sol.y(33,:),sol.y(31,:),'k--','LineWidth',1);
        %plot3([sol.y(18,end),sol.y(18,end)+0.1*d32], [sol.y(19,end),sol.y(19,end)+0.1*d33],[sol.y(17,end),sol.y(17,end)+ 0.1*d31],'c-')
     end
    %plot3([sol.y(18,end),sol.y(18,end)+0.1*d32], [sol.y(19,end),sol.y(19,end)+0.1*d33],[sol.y(17,end),sol.y(17,end)+ 0.1*d31],'m:')
    e=0.01;
    %annotation('arrow',[sol.y(18,end),sol.y(19,end),sol.y(20,end)],[sol.y(18,end),sol.y(19,end),sol.y(20,end)]+ e*[d31,d32,d33])
    pause(0.1)
    %TIP
    tipx=[tipx,sol.y(16+2,end)];
    tipy=[tipy,sol.y(16+3,end)];
    tipz=[tipz,sol.y(16+1,end)];
end
%}
[d31,d32,d33]
%xlim([0,0.5])
%ylim([0,1.5])
%zlim([0,0.5])
plot3(sol.y(2,:),sol.y(3,:),sol.y(1,:),'r-','LineWidth',3);
plot3(sol.y(16+2,:),sol.y(16+3,:),sol.y(16+1,:),'b-','LineWidth',3);
plot3(sol.y(32,:),sol.y(33,:),sol.y(31,:),'k-','LineWidth',3);
plot3([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] )
patch([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] ,'k','FaceAlpha',.3)

%plot3([sol.y(18,end),sol.y(18,end)+0.1*d32], [sol.y(19,end),sol.y(19,end)+0.1*d33],[sol.y(17,end),sol.y(17,end)+ 0.1*d31],'g-')
%plot3(tipx,tipy,tipz,'g-','Linewidth',1)
xlim([-0.2,0.4])
ylim([0,1.6])
zlim([-0.40,0.25])
plot3(tipx,tipy,tipz,'m-','Linewidth',1)


fig2=figure(2)

a1=y5(12*(0:N)+1);
a2=y5(12*(0:N)+2);
a3=y5(12*(0:N)+3);
a4=y5(12*(0:N)+4);
a5=y5(12*(0:N)+5);
a6=y5(12*(0:N)+6);
subplot(2,1,1)
plot(linspace(0,1,N+1),a1,'g--');
hold on;
grid on;
plot(linspace(0,1,N+1),a2,'b--');
plot(linspace(0,1,N+1),a3,'r--');

title('Plot of rotation of tubes as a function of time t','FontSize',15)
legend('\theta^{[1]}','\theta^{[2]} ','\theta^{[3]}','FontSize',18);
xlabel('Time t','FontSize',15)
ylabel('Controls y(t)','FontSize',15)
subplot(2,1,2)
hold on;
grid on;
plot(linspace(0,1,N+1),a4,'k--');
plot(linspace(0,1,N+1),a5,'m--');
plot(linspace(0,1,N+1),a6,'c--');
xlabel('Time t','FontSize',15)
ylabel('Controls y(t)','FontSize',15)
title('Plot of lengths of tubes as a function of time t','FontSize',15)
legend('L^{[1]}','L^{[2]}','L^{[3]}','FontSize',18);

fig3=figure(3)
a1=y5(12*(0:N)+7);
a2=y5(12*(0:N)+8);
a3=y5(12*(0:N)+9);
a4=y5(12*(0:N)+10);
a5=y5(12*(0:N)+11);
a6=y5(12*(0:N)+12);
plot(linspace(0,1,N+1),a1,'g-^');
hold on;
grid on;
plot(linspace(0,1,N+1),a2,'b-^');
plot(linspace(0,1,N+1),a3,'r-^');
plot(linspace(0,1,N+1),a4,'k-^');
plot(linspace(0,1,N+1),a5,'m-^');
plot(linspace(0,1,N+1),a6,'c-^');
xlabel('Time t')
ylabel('Controls y(t)')
title('Plot of Control rates as a function of time t')

legend('\theta_{1}','\theta_{2} ','\theta_{3}','L_{1}','L_{2}','L_{3}');

fig4=figure(4)
plot(linspace(0,1,N+1),theta);
hold on
xlabel('Time t','FontSize',15)
ylabel('Orientation between prescribed orientation and the robot"s tip \Theta (in degrees) ','FontSize',15)
grid on
legend('\lambda=0','FontSize',18)
title('Angle between target orientation and tangent of the tip','FontSize',18)


