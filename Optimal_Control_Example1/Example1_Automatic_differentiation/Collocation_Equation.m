function [FU,Ft] = Collocation_Equation(sol,u,Load)
nx=Load(1);
ny=Load(2);
nz=Load(3);
ang1=u(1);
ang2=u(2);
ang3=u(3);
l1= u(4);
l2= u(5);
l3= u(6);
x=sol.x;
y=sol.y;
mesh=size(x,2);
Ndim=48;
F_U=zeros(Ndim*mesh);
Y=[];
for k=1:mesh
    Y=[Y;y(:,k)];
end
Jc=zeros(Ndim*(mesh));
F=[];
YY=[];
Fo=[];
FFo=[];
Fl1=[];
Fl2=[];
Fl3=[];
for k=1:mesh-1
    h=x(k+1)-x(k);
    Jc((k-1)*Ndim + 1 : (k-1)*Ndim+ Ndim, (k-1)*Ndim + 1 : (k-1)*Ndim+ (1+1)*Ndim)=[eye(Ndim)/h,-eye(Ndim)/h];
    YK=Y((k-1)*Ndim+1:k*Ndim);
    YK1=Y((k)*Ndim+1:(k+1)*Ndim);
    [FF,FJac]=CTCR(x,(YK+YK1)/2,ang1,ang2,ang3,l1,l2,l3,[nx,ny,nz]);
    F=[F,FF];
    F_U((k-1)*Ndim + 1 : (k-1)*Ndim+ Ndim, (k-1)*Ndim + 1 : (k-1)*Ndim+ (1+1)*Ndim)=[FJac/2,FJac/2];
    Fl1=[Fl1;CTCR(x,(YK+YK1)/2,ang1,ang2,ang3,1.0,0.0,0.0,[nx,ny,nz])'];
    Fl2=[Fl2;CTCR(x,(YK+YK1)/2,ang1,ang2,ang3,0.0,1.0,0.0,[nx,ny,nz])'];
    Fl3=[Fl3;CTCR(x,(YK+YK1)/2,ang1,ang2,ang3,0.0,0.0,1.0,[nx,ny,nz])'];
    %F=[F,(CTCR(x,YK,ang1,ang2,ang3,l1,l2,l3,[nx,ny,nz]) + CTCR(x,YK1,ang1,ang2,ang3,l1,l2,l3,[nx,ny,nz]))/2.0];
    Fo=[Fo; (YK1-YK)/h- FF'];
    FFo=[FFo;FF'];
    YY=[YY;YK];
end
YY=[YY;YK1];
%Boundary Conditions
Y0=Y(1:48);
Y1=Y(end-47:end);
BC0=zeros(48);
BC1=zeros(48);

BC0(1:14,1:14)=eye(14);
BC1(1:14,31:44)=-eye(14);

BC0(15:28,17:30)=eye(14);
BC1(15:28,1:14)= -eye(14);


BC0(29:35,31:37)=eye(7);%BC
BC1(36:38,12:14)=eye(3);

BC1(39:41,20:27)=[-Y1(23+4),-Y1(22+4), Y1(21+4), Y1(20+4), Y1(23),Y1(22), -Y1(21), -Y1(20);Y1(22+4),-Y1(23+4),-Y1(20+4),Y1(21+4),-Y1(22),Y1(23),Y1(20),-Y1(21);-Y1(21+4),Y1(20+4),-Y1(23+4),Y1(22+4),Y1(21),-Y1(20),Y1(23),-Y1(22)];
BC1(42,17:30)=[2*Y1(26),2*Y1(27),2*Y1(28),Y1(20+4),Y1(21+4),Y1(22+4),Y1(23+4),Y1(20),Y1(21),Y1(22),Y1(23),2*Y1(17),2*Y1(18),2*Y1(19)];

BC0(43,45)=1;%BC
BC0(44,47)=1;%BC

BC1(45,48)=1;
BC1(46,16)=1;

BC0(47,15)=1;
BC0(48,16)=1;

BC1(47,45)=-1;
BC1(48,46)=-1;


%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC0(1:7,1:7)=eye(7);
BC1(1:7,31:37)=-eye(7);

BC1(8:14,8:14)=eye(7);
BC0(8:14,24:30)=-eye(7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC0(15,15)= 1;
BC1(15,45)=-1;
BC1(16,16)=1;

BC0(17:23,17:23)=eye(7);
BC1(17:23,1:7)=-eye(7);

%BC1(24:26,24:27)=[Y1(23),Y1(22), -Y1(21), -Y1(20);-Y1(22),Y1(23),Y1(20),-Y1(21);Y1(21),-Y1(20),Y1(23),-Y1(22)];
%BC0(27,24:30)=[Y0(20),Y0(21),Y0(22),Y0(23),2*Y0(17),2*Y0(18),2*Y0(19)];

BC1(24:26,20:27)=[-Y1(23+4),-Y1(22+4), Y1(21+4), Y1(20+4), Y1(23),Y1(22), -Y1(21), -Y1(20);Y1(22+4),-Y1(23+4),-Y1(20+4),Y1(21+4),-Y1(22),Y1(23),Y1(20),-Y1(21);-Y1(21+4),Y1(20+4),-Y1(23+4),Y1(22+4),Y1(21),-Y1(20),Y1(23),-Y1(22)];
BC0(27,17:30)=[2*Y0(26),2*Y0(27),2*Y0(28),Y0(20+4),Y0(21+4),Y0(22+4),Y0(23+4),Y0(20),Y0(21),Y0(22),Y0(23),2*Y0(17),2*Y0(18),2*Y0(19)];

BC1(28:30,28:30)=eye(3);
BC0(31:37,31:37)=eye(7);

BC1(38:40,42:44)=eye(3);
BC0(38:40,12:14)=-eye(3);

BC1(41:44,38:41)=eye(4);
BC0(41:44,8:11)=-eye(4);

BC0(45,45)=1;
BC1(46,46)=1;
BC0(46,16)=-1;
BC0(47,47)=1;
BC1(48,48)=1;
%}

BC=zeros(48,size(Y,1));
BCF=zeros(48,1);
BC(1:48,1:48)=BC0;
BC(1:48, end-47:end )=BC1;
BCF(28)=-nx;
BCF(29)=-ny;
BCF(30)=-nz;
BCF(36)= sin(ang1/2);
BCF(37)= cos(ang1/2);
BCF(45)= ang2 - ang1;
BCF(47)= ang3 - ang1;
Jc((k)*Ndim + 1 : (k)*Ndim+ Ndim, 1:48)=BC0;
Jc((k)*Ndim + 1 : (k)*Ndim+ Ndim, end-47:end)=BC1;
F=[F';-BCF];
%%%%%%%%%%%%%%%%%%%%%
FFo=[FFo;-BCF];
Fo=Jc*YY + FFo;
FU=Jc + F_U;

Fl1=[Fl1;zeros(48,1)];
Fl2=[Fl2;zeros(48,1)];
Fl3=[Fl3;zeros(48,1)];

Ft1=zeros(size(Y,1),1);
Ft2=zeros(size(Y,1),1);
Ft3=zeros(size(Y,1),1);


Ft1((mesh-1)*Ndim+34)=-cos(ang1/2)/2;
Ft1((mesh-1)*Ndim+35)= sin(ang1/2)/2;
Ft1((mesh-1)*Ndim+43)=1;
Ft1((mesh-1)*Ndim+44)=1;
Ft2((mesh-1)*Ndim+43)= -1;
Ft3((mesh-1)*Ndim+44)= -1;

Ft=[Ft1,Ft2,Ft3,Fl1,Fl2,Fl3];

%norm(Fo)
%{
fun=@(Y)Algebraic_system(Y,x,u,Load)
norm(fun(Y))

Y_i=Y;
Y(49)=Y(49)+ -0.21;

Y(48+2)=Y(48+2)+ 1.0;
Y(48*2+3)=Y(48*3+3) - 0.5
Y(48*4+3)=Y(48*3+3) - 0.5
Y(48*5+3)=Y(48*5+3) + 0.75
Y(48*7+3)=Y(48*7+3) - 0.25
Y(end-48)=Y(end-48)+0.4;

norm(Y - Y_i)
tic
xx = fsolve(fun,Y);
norm(xx - Y_i)
toc

tic
options = optimoptions('fsolve','SpecifyObjectiveGradient',true)
xxx = fsolve(fun,Y,options);
norm(xxx - Y_i)
toc
norm(xxx - xx)

tic
options = optimoptions('fsolve','SpecifyObjectiveGradient',true)
xxxy = fsolve(fun,xxx,options);
norm(xxxy - Y_i)
toc
norm(xxxy - xxx)


%tic
%A=inv(Jc - F_U)*Fc;
%A=(Jc - F_U)\Fc;
%toc
%tic
%B=sparse(Jc - F_U);
%C=inv(B)*Fc;
%J=B\Fc;
%toc

%norm(F)

%Debugging::: Verifying the correctness of Jacobian
%Write a function
function [Fo,FU]=Algebraic_system(Y,x,u,Load)
        nx=Load(1);
        ny=Load(2);
        nz=Load(3);
        ang1=u(1);
        ang2=u(2);
        ang3=u(3);
        l1= u(4);
        l2= u(5);
        l3= u(6);
        Ndim=48;
        mesh=size(Y,1)/Ndim;
        F_U=zeros(Ndim*mesh);
        Jc=zeros(Ndim*(mesh));
        F=[];
        YY=[];
        Fo=[];
        FFo=[];
        for k=1:mesh-1
            h=x(k+1)-x(k);
            YK=Y((k-1)*Ndim+1:k*Ndim);
            YK1=Y((k)*Ndim+1:(k+1)*Ndim);
            Jc((k-1)*Ndim + 1 : (k-1)*Ndim+ Ndim, (k-1)*Ndim + 1 : (k-1)*Ndim+ (1+1)*Ndim)=[eye(Ndim)/h,-eye(Ndim)/h];
            [FF,FJac]=CTCR(x,(YK+YK1)/2,ang1,ang2,ang3,l1,l2,l3,[nx,ny,nz]);
            F_U((k-1)*Ndim + 1 : (k-1)*Ndim+ Ndim, (k-1)*Ndim + 1 : (k-1)*Ndim+ (1+1)*Ndim)=[FJac/2,FJac/2];
            Fo=[Fo; (YK1-YK)/h- FF'];
            FFo=[FFo;FF'];
            YY=[YY;YK];
        end
        YY=[YY;YK1];
        %Boundary Conditions
        Y0=Y(1:48);
        Y1=Y(end-47:end);
        BC0=zeros(48);
        BC1=zeros(48);
        BC0(1:7,1:7)=eye(7);
        BC1(1:7,31:37)=-eye(7);

        BC1(8:14,8:14)=eye(7);
        BC0(8:14,24:30)=-eye(7);

        BC0(15,15)= 1;
        BC1(15,45)=-1;

        BC1(16,16)=1;

        BC0(17:23,17:23)=eye(7);
        BC1(17:23,1:7)=-eye(7);

        BC1(24:26,24:27)=[Y1(23),Y1(22), -Y1(21), -Y1(20);-Y1(22),Y1(23),Y1(20),-Y1(21);Y1(21),-Y1(20),Y1(23),-Y1(22)];
        BC0(27,24:30)=[Y0(20),Y0(21),Y0(22),Y0(23),2*Y0(17),2*Y0(18),2*Y0(19)];

        BC1(28:30,28:30)=eye(3);
        BC0(31:37,31:37)=eye(7);

        BC1(38:40,42:44)=eye(3);
        BC0(38:40,12:14)=-eye(3);

        BC1(41:44,38:41)=eye(4);
        BC0(41:44,8:11)=-eye(4);

        BC0(45,45)=1;
        BC1(46,46)=1;
        BC0(46,16)=-1;
        BC0(47,47)=1;
        BC1(48,48)=1;

        BC=zeros(48,size(Y,1));
        BCF=zeros(48,1);
        BC(1:48,1:48)=BC0;
        BC(1:48, end-47:end )=BC1;
        BCF(28)=-nx;
        BCF(29)=-ny;
        BCF(30)=-nz;
        BCF(36)=sin(ang1/2);
        BCF(37)=cos(ang1/2);
        BCF(45)=ang2 - ang1;
        BCF(47)=ang3 - ang1;
        Jc((k)*Ndim + 1 : (k)*Ndim+ Ndim, 1:48)=BC0;
        Jc((k)*Ndim + 1 : (k)*Ndim+ Ndim, end-47:end)=BC1;
        F=[F';BCF];
%%%%%%%%%%%%%%%%%%%%%
        FFo=[FFo;-BCF];
        Fo=Jc*YY + FFo;
        FU=Jc + F_U;
    end

%}

end
