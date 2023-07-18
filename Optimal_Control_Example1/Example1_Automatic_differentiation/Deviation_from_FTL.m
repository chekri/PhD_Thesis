function [earea,Jac_area]=Deviation_from_FTL(U,Tr,UFL,TrF,J,Fc)
%TrF=%FTL
%Tr=%Normal Function
% perturbation for finite differences
%Length intervals
%tic
t1=TrF{6};
t2=TrF{5};
t3=TrF{4};
l1= U(4);
l2= U(5);
l3= U(6);
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
%figure(1)
%plot3(XF,YF,ZF,'k-o')
%hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial state (state at s=0) from the solutions of CTCR-bvp solutions.
%Useful for computing gradients
s=Tr.x;
XX=[Tr.y(31,:),Tr.y(1,:),Tr.y(17,:)];
YY=[Tr.y(32,:),Tr.y(2,:),Tr.y(18,:)];
ZZ=[Tr.y(33,:),Tr.y(3,:),Tr.y(19,:)];
l =[l3*s,l3 + l2*s, l3 + l2 + l1*s];
%{
NN=size(s,2)
XX(NN)=[];   
YY(NN)=[];
ZZ(NN)=[];
l(NN)=[];

XX(2*NN)=[];   
YY(2*NN)=[];
ZZ(2*NN)=[];
l(2*NN)=[];

plot3(XX,YY,ZZ,'r-x')
%}
%Solutions of the IVP with the prescribed control parameters and intial
%state

%[p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC',U);
%l1=U(4);
%l2=U(5);
%l3=U(6);
%l =[l3*t3',l3 + l2*t2', l3 + l2 + l1*t1'];       % Arclength Steps alson the Robot.

%Interpolation of the robots FTL's position vector in terms of its arclength.
xf=interp1(Al,XF,'spline','pp');
yf=interp1(Al,YF,'spline','pp');
zf=interp1(Al,ZF,'spline','pp');

xs=ppval(xf,l);
ys=ppval(yf,l);
zs=ppval(zf,l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xf_der=fnder(xf,1);
Yf_der=fnder(yf,1);
Zf_der=fnder(zf,1);

sec=size(l,2)/3;
lo=ones(size(l));
lo(1:sec)=l(1:sec)/l3;
grad_x3=ppval(Xf_der,l).*lo;
grad_y3=ppval(Yf_der,l).*lo;
grad_z3=ppval(Zf_der,l).*lo;
Rfu3=[grad_x3;grad_y3;grad_z3];


lo=ones(size(l));
lo((sec+1):(2*sec))=(l((sec+1):(2*sec)) - l3)/l2;
lo(1:sec)=0*l(1:sec);
grad_x2=ppval(Xf_der,l).*lo;
grad_y2=ppval(Yf_der,l).*lo;
grad_z2=ppval(Zf_der,l).*lo;
Rfu2=[grad_x2;grad_y2;grad_z2];

lo=zeros(size(l));
lo(2*sec + 1:3*sec)=(l(2*sec + 1:3*sec)-l2 -l3)/l1;
grad_x1=ppval(Xf_der,l).*lo;
grad_y1=ppval(Yf_der,l).*lo;
grad_z1=ppval(Zf_der,l).*lo;
Rfu1=[grad_x1;grad_y1;grad_z1];

%Verification of gradient
ee=0.00001;
ll3=l3+ ee ;
ll2=l2;
ll1=l1;
ll =[ll3*s,ll3 + ll2*s, ll3 + ll2 + ll1*s];
xvs=ppval(xf,ll);
yvs=ppval(yf,ll);
zvs=ppval(zf,ll);

(xvs - xs)/ee -grad_x3;
(yvs - ys)/ee-grad_y3;
(zvs - zs)/ee-grad_z3;


%{
plot3(xs,ys,zs,'k-x')
%FTL's position vectors along the arclength projection of the current robot's configuration. 
for k=1:size(l,2)
    figure(1)
    plot3([xs(k),XX(k)],[ys(k),YY(k)],[zs(k),ZZ(k)],'g')
    figure(2)
    scatter(l(k),norm([xs(k),ys(k),zs(k)]-[XX(k),YY(k),ZZ(k)]),'r')
    hold on
    grid on
end
%}
AA=[xs;ys;zs];
B=[XX;YY;ZZ];
%B=[p3(:,1)',p2(:,1)',p1(:,1)';p3(:,2)',p2(:,2)',p1(:,2)';p3(:,3)',p2(:,3)',p1(:,3)'];
%The area of the devaition from the FTL shape 
earea=enclosed_area1(AA,B,l);
[earea,Jac_area,JC]=enclosed_area(AA,B,l,Rfu1,Rfu2,Rfu3);

lam = linsolve(J',Jac_area);
Jac_area=-lam'*Fc + JC ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function area=enclosed_area1(Rf,R,l) 
%Computes the distance between the sepcified position vectors and sums over
%it.
area=trapz(l,vecnorm((Rf-R)));

function [area,Jac,JC]=enclosed_area(Rf,R,l,Rfu1,Rfu2,Rfu3) 
area=0;
Jc1=0;
Jc2=0;
Jc3=0;
Jac=zeros(48*(size(l,2)/3),1);
l3=l(size(l,2)/3);
l2=l(2*size(l,2)/3) -l3;
l1=l(end)-l2-l3;
for k=1:size(l,2)-2
    hk1=l(k+2) - l(k+1);
    hk=l(k+1) - l(k);
    area=area + hk*( norm( Rf(:,k)- R(:,k)) + norm( Rf(:,k+1)- R(:,k+1)) )/2;
    
    if k==1
        Jc1=Jc1 + hk*(Rf(:,k+1)- R(:,k+1))'*Rfu1(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1)));
        Jc2=Jc2 + hk*(Rf(:,k+1)- R(:,k+1))'*Rfu2(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1)));
        Jc3=Jc3 + hk*(Rf(:,k+1)- R(:,k+1))'*Rfu3(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1)));
    end
    
    if k>1
        Jc1=Jc1 + hk*(Rf(:,k)- R(:,k))'*Rfu1(:,k)/(2*norm( Rf(:,k)- R(:,k)))+ hk*(Rf(:,k+1)- R(:,k+1))'*Rfu1(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1))); 
        Jc2=Jc2 + hk*(Rf(:,k)- R(:,k))'*Rfu2(:,k)/(2*norm( Rf(:,k)- R(:,k))) +hk*(Rf(:,k+1)- R(:,k+1))'*Rfu2(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1)));
        Jc3=Jc3 + hk*(Rf(:,k)- R(:,k))'*Rfu3(:,k)/(2*norm( Rf(:,k)- R(:,k))) +hk*(Rf(:,k+1)- R(:,k+1))'*Rfu3(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1)));
    end
    
    if k< size(l,2)/3
        Jac(48*k+31:48*k + 33)= -(hk+hk1)*(Rf(:,k+1)- R(:,k+1))/(2*norm( Rf(:,k+1)- R(:,k+1)));
        Jc3=Jc3 + hk*( norm( Rf(:,k)- R(:,k)) + norm( Rf(:,k+1)- R(:,k+1)) )/(2*l3);
    end 

    if k> size(l,2)/3-1 && k< 2*size(l,2)/3
        j= k-size(l,2)/3;
        Jac(48*j+1:48*j+ 3)=  -(hk+hk1)*(Rf(:,k+1)- R(:,k+1))/(2*norm( Rf(:,k+1)- R(:,k+1)));
        Jc2=Jc2 + hk*( norm( Rf(:,k)- R(:,k)) + norm( Rf(:,k+1)- R(:,k+1)) )/(2*l2);
    end

    if k> 2*size(l,2)/3-1
        j=k-2*size(l,2)/3;
        Jac(48*j+17:48*j + 19)=-(hk+hk1)*(Rf(:,k+1)- R(:,k+1))/(2*norm( Rf(:,k+1)- R(:,k+1)));
        Jc1=Jc1+ hk*( norm( Rf(:,k)- R(:,k)) + norm( Rf(:,k+1)- R(:,k+1)) )/(2*l1);
    end
 
end
area=area + hk1*( norm( Rf(:,k+1)- R(:,k+1)) + norm( Rf(:,k+2)- R(:,k+2)) )/2;

Jc1=Jc1 + hk1*(Rf(:,k+1)- R(:,k+1))'*Rfu1(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1))) +  hk1*(Rf(:,k+2)- R(:,k+2))'*Rfu1(:,k+2)/(2*norm( Rf(:,k+2)- R(:,k+2)))+ hk1*( norm( Rf(:,k+1)- R(:,k+1)) + norm( Rf(:,k+2)- R(:,k+2)) )/(2*l1);
Jc2=Jc2 + hk1*(Rf(:,k+1)- R(:,k+1))'*Rfu2(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1))) +  hk1*(Rf(:,k+2)- R(:,k+2))'*Rfu2(:,k+2)/(2*norm( Rf(:,k+2)- R(:,k+2)));
Jc3=Jc3 + hk1*(Rf(:,k+1)- R(:,k+1))'*Rfu3(:,k+1)/(2*norm( Rf(:,k+1)- R(:,k+1))) +  hk1*(Rf(:,k+2)- R(:,k+2))'*Rfu3(:,k+2)/(2*norm( Rf(:,k+2)- R(:,k+2)));

Jac(48*(j+1)+17:48*(j+1) + 19)=  -hk1*(Rf(:,k+2)- R(:,k+2))/(2*norm( Rf(:,k+2)- R(:,k+2))); 
JC=[0,0,0,Jc1,Jc2,Jc3];
