function Animationvideo_lefthanded(y5,P)


%[y5,y10]=OptimalControlwithCollocation_Jacabian()
%Rtar=[-0.20,-0.55,0.76];
%R_init= [-0.17506,0.320572,1.44679]

R_tar=[0.16,0.167,0.767]
R_tar=[0.2,0.3,1.0]
R_tar=[0.1,-0.1,1.0]
R_tar=[0.4,-0.0,1.0]
UFL=[0.5,0.5+3.1415,0.5+3.1415,0.4+2,0.6,0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RF3,RF2,RF1,t3,t2,t1]=IVP_trajectory([0,0,0,0,0,sin(UFL(1)/2),cos(UFL(1)/2),0,0,0,0,0,0,0,UFL(2)-UFL(1),0,UFL(3)-UFL(1),0],UFL);
RF3=RF3(:,1:3);
RF2=RF2(:,1:3);
RF1=RF1(:,1:3);
TrF={RF3,RF2,RF1,t3,t2,t1}

v = VideoWriter('myAnimation.mp4', 'MPEG-4');
open(v);

N=floor(size(y5,2)/12)
sol1=Trajectory([y5(1),y5(2),y5(3),y5(4),y5(5),y5(6)],P,[]);
fig=figure(1)

%set(gca,'XTickLabel',[]);
%set(gca,'YTickLabel',[]);
%set(gca,'ZTickLabel',[]);

%plot3(sol1.y(32,:),-sol1.y(33,:),sol1.y(31,:),'k--','LineWidth',2);

hold on
%plot3(sol1.y(2,:),-sol1.y(3,:),sol1.y(1,:),'k--','LineWidth',2);

%plot3(sol1.y(16+2,:),-sol1.y(16+3,:),sol1.y(16+1,:),'k--','Linewidth',2);

tipx=[];
tipy=[];
tipz=[];
theta=[];
dev=[];
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

    theta=[theta,acosd(dot([d31,d32,d33],[di31,di32,di33])/(norm([d31,d32,d33])*norm([di31,di32,di33])))];
    dev=[dev,Deviation_from_FTL([y5(12*k+1),y5(12*k+2),y5(12*k+3),y5(12*k+4),y5(12*k+5),y5(12*k+6)],sol,UFL,TrF)];

    %hold on
    
     if k==0
        plot3(sol.y(2,:),-sol.y(3,:),sol.y(1,:),'r-','LineWidth',5);
        hold on
        plot3(sol.y(16+2,:),-sol.y(16+3,:),sol.y(16+1,:),'b-','LineWidth',5);
        plot3(sol.y(32,:),-sol.y(33,:),sol.y(31,:),'k-','LineWidth',5);
        plot3([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] )
        patch([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] ,'k','FaceAlpha',.3)
        plot3([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] )
        scatter3(sol1.y(16+2,end),-sol1.y(16+3,end),sol1.y(16+1,end),'b*');
        scatter3(R_tar(2),-R_tar(3),R_tar(1),'k*')
        xlim([-0.2,1.4])
        ylim([-1.6,0])
        zlim([-0.3,1.3])
        view(60,60)

      end
     
      if rem(k,1)==0
        plot3(sol1.y(2,:),-sol1.y(3,:),sol1.y(1,:),'r-','LineWidth',5);
        hold on
        plot3(sol1.y(16+2,:),-sol1.y(16+3,:),sol1.y(16+1,:),'b-','LineWidth',5);
        plot3(sol1.y(32,:),-sol1.y(33,:),sol1.y(31,:),'k-','LineWidth',5);

        plot3(sol.y(2,:),-sol.y(3,:),sol.y(1,:),'r-','LineWidth',5);
        plot3(sol.y(16+2,:),-sol.y(16+3,:),sol.y(16+1,:),'b-','LineWidth',5);
        plot3(sol.y(32,:),-sol.y(33,:),sol.y(31,:),'k-','LineWidth',5);

        plot3([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] )
        patch([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] ,'k','FaceAlpha',.3)
        plot3([-0.1 -0.1 0.1 0.1 -0.1],[0 0 0 0 0],[-0.1 0.1 0.1 -0.1 -0.1] )
        scatter3(sol1.y(16+2,end),-sol1.y(16+3,end),sol1.y(16+1,end),'b*');
        scatter3(R_tar(2),-R_tar(3),R_tar(1),'k*')
        xlim([-0.2,1.4])
        ylim([-1.6,0])
        zlim([-0.3,1.3])
        view(-60,15)
        %axis equal
        grid on
     end
      e=0.01;
      
    tipx=[tipx,sol.y(16+2,end)];
    tipy=[tipy,sol.y(16+3,end)];
    tipz=[tipz,sol.y(16+1,end)];
    plot3(tipx,-tipy,tipz,'c-','Linewidth',3)
    drawnow
    hold off
    pause(0.5)
    frame = getframe(fig);
    writeVideo(v, frame.cdata);
end
%}
[d31,d32,d33]



z=[0,0,-0.2]
close(v)


%plot3([sol.y(18,end),sol.y(18,end)+0.1*d32], [sol.y(19,end),sol.y(19,end)+0.1*d33],[sol.y(17,end),sol.y(17,end)+ 0.1*d31],'g-')
%plot3(tipx,tipy,tipz,'g-','Linewidth',1)
%xlim([-0.2,0.4])
%ylim([0,1.8])
%zlim([-0.50,0.15])

%view(90+180,0)
%legend



