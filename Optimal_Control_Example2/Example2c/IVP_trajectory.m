function [p3,p2,p1,t3,t2,t1]=IVP_trajectory(IC,C)
Th1=C(1);
Th2=C(2);
Th3=C(3);
L1=C(4);
L2=C(5);
L3=C(6);

IC3=IC;
IC3(6)=sin(Th1/2);
IC3(7)=cos(Th1/2);
IC3(15)=Th2-Th1;
IC3(17)=Th3-Th1;
[t3,p3]=ode45(@(x,y) three_tubesection(x,y,L1,L2,L3),[0 1],IC3);
sol_0=p3(end,:);
%size(p3)
IC2=sol_0(1:16);
[t2,p2]=ode45(@(x,y) two_tubesection(x,y,L1,L2,L3),[0 1],IC2);
sol_0=p2(end,:);
%size(p2)
IC1=sol_0(1:14);
[t1,p1]=ode45(@(x,y) one_tubesection(x,y,L1,L2,L3),[0 1],IC1);
%size(p1)

function F = three_tubesection(x,u,L1,L2,L3)
      uhat11=0.5;%PAR(5)%1.0
      uhat12=0.0;
      uhat13=0.0;%PAR(3)
      
      uhat21=0.8;%PAR(6)%2.0
      uhat22=0.0;
      uhat23=0.0;%PAR(3)

      uhat31=1.0;%2.0
      uhat32=0.0;
      uhat33=0.0;%PAR(3)
  
     A1=1.0;
     A2=1.2;
     A3=1.4;
     C1=1.0/1.3;
     C2=1.2/1.3;
     C3=A3/1.3;
           	
    epsilon = 0.0;
    aa1=1000.0;
    aa2=1000.0;
    aa3=1000.0;
    
    
      m3(1) = ( u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0;
      m3(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0;
      m3(3) = ( u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0;

%%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%%         the same 4 entries, but permuted and possibly negated.
%%         Here are those entries

      entry31 =  2.0*(u(4)*u(12) + u(5)*u(13) + u(6)*u(14));
      entry32 =  2.0*(u(5)*u(12) - u(4)*u(13) + u(7)*u(14));
      entry33 =  2.0*(u(6)*u(12) - u(7)*u(13) - u(4)*u(14));
      entry34 =  2.0*(u(7)*u(12) + u(6)*u(13) - u(5)*u(14));%%Done

%%       Components of D_1[q] n

      d31qn(1) =  entry31;
      d31qn(2) = -entry32;
      d31qn(3) = -entry33;
      d31qn(4) =  entry34;

%%       Components of D_2[q] n

      d32qn(1) =  entry32;
      d32qn(2) =  entry31;
      d32qn(3) = -entry34;
      d32qn(4) = -entry33;

%%       Components of D_3[q] n

      d33qn(1) =  entry33;
      d33qn(2) =  entry34;
      d33qn(3) =  entry31;
      d33qn(4) =  entry32;

%%       Components of d1

      d31(1)= u(4)*u(4) - u(5)*u(5) - u(6)*u(6) + u(7)*u(7);
      d31(2)= 2.0*u(4)*u(5) + 2.0*u(6)*u(7);
      d31(3)= 2.0*u(4)*u(6) - 2.0*u(5)*u(7);

%%       Components of d2

      d32(1)= 2.0*u(4)*u(5) - 2.0*u(6)*u(7);
      d32(2)= -u(4)*u(4) + u(5)*u(5) - u(6)*u(6) + u(7)*u(7);
      d32(3)= 2.0*u(5)*u(6) + 2.0*u(4)*u(7);

%%       Components of d3
      d33(1)= 2.0*u(4)*u(6) + 2.0*u(5)*u(7);
      d33(2)= 2.0*u(5)*u(6) - 2.0*u(4)*u(7);
      d33(3)= -u(4)*u(4) - u(5)*u(5) + u(6)*u(6) + u(7)*u(7);
%%       Components of v; note that v(3) has an extra term, due to
%%         the fact that the unstressed state is assumed to have v=[0 0 1].
      v3(1)= epsilon/aa1*(d31(1)*u(12) + d31(2)*u(13) + d31(3)*u(14));
      v3(2)= epsilon/aa2*(d32(1)*u(12) + d32(2)*u(13) + d32(3)*u(14));
      v3(3)= epsilon/aa3*(d33(1)*u(12) + d33(2)*u(13) + d33(3)*u(14)) + 1.0;


     ab1=-C2*((u(16) + u(18))/C1 - m3(3)/C1  + u(16)/C2 + uhat23- uhat13 - uhat23); 
     ab2=-C3*((u(16) + u(18))/C1 - m3(3)/C1  + u(18)/C3 + uhat33- uhat13 - uhat33);
    


    ub31=(A1*uhat11 + A2*(uhat21*cos(u(15))- uhat22*sin(u(15)) ) + A3*(uhat31*cos(u(17))- uhat32*sin(u(17))))/(A1+A2+A3);
	ub32=(A1*uhat12 + A2*(uhat21*sin(u(15))+ uhat22*cos(u(15)) ) + A3*(uhat31*sin(u(17))+ uhat32*cos(u(17))))/(A1+A2+A3);
	ub33= (ab1 + ab2 + C1*uhat13)/(C1+C2+C3);% Corrected addterm to uhat13 ------------------------


	kb31=A1+A2+A3;
	kb32=A1+A2+A3;
	kb33=C1+C2+C3; 
     
	strain3(1) = m3(1)/(A1+A2+A3) + ub31;
	strain3(2) = m3(2)/(A1+A2+A3) + ub32;
	strain3(3) = m3(3)/(C1+C2+C3) + ub33;  

	F(1) = L3*(v3(1)*d31(1) + v3(2)*d32(1) + v3(3)*d33(1));
    F(2) = L3*(v3(1)*d31(2) + v3(2)*d32(2) + v3(3)*d33(2));
	F(3) = L3*(v3(1)*d31(3) + v3(2)*d32(3) + v3(3)*d33(3));
      
	F(4) = L3*( strain3(1)*u(7)/2.0 - strain3(2)*u(6)/2.0 + strain3(3)*u(5)/2.0);
    F(5) = L3*( strain3(1)*u(6)/2.0 + strain3(2)*u(7)/2.0 - strain3(3)*u(4)/2.0);
    F(6) = L3*(-strain3(1)*u(5)/2.0 + strain3(2)*u(4)/2.0 + strain3(3)*u(7)/2.0);
    F(7) = L3*(-strain3(1)*u(4)/2.0 - strain3(2)*u(5)/2.0 - strain3(3)*u(6)/2.0);
     
	F(8) = L3*( strain3(1)*u(11)/2.0 - strain3(2)*u(10)/2.0 + strain3(3)*u(9)/2.0 - v3(1)*d31qn(1)-v3(2)*d32qn(1)-v3(3)*d33qn(1));
	F(9) = L3*( strain3(1)*u(10)/2.0 + strain3(2)*u(11)/2.0 - strain3(3)*u(8)/2.0 - v3(1)*d31qn(2)-v3(2)*d32qn(2)-v3(3)*d33qn(2));
	F(10) = L3*(-strain3(1)*u(9)/2.0 + strain3(2)*u(8)/2.0 + strain3(3)*u(11)/2.0 - v3(1)*d31qn(3)-v3(2)*d32qn(3)-v3(3)*d33qn(3));
	F(11) = L3*(-strain3(1)*u(8)/2.0 - strain3(2)*u(9)/2.0 - strain3(3)*u(10)/2.0 - v3(1)*d31qn(4)-v3(2)*d32qn(4)-v3(3)*d33qn(4));
     
	F(12) = L3*0.0;
	F(13) = L3*0.0;
	F(14) = L3*0.0;
   	F(15)= L3*((u(16) + u(18))/C1 - m3(3)/C1  + u(16)/C2 + uhat23- uhat13);  %%%Alfa2=u(15)
	F(16)= L3*(A2*uhat21*(A3*uhat31*sin(u(15)-u(17)) + A1*uhat11*sin(u(15))))/(A1+A2+A3) + L3*A2*uhat21*(m3(1)*sin(u(15)) - m3(2)*cos(u(15)))/(A1+A2+A3);
    
    F(17)= L3*((u(16) + u(18))/C1 - m3(3)/C1  + u(18)/C3 + uhat33- uhat13);  %%%Alfa3=u(17)
	F(18)= L3*(A3*uhat31*(A2*uhat21*sin(u(17)-u(15)) + A1*uhat11*sin(u(17))))/(A1+A2+A3) + A3*L3*uhat31*(m3(1)*sin(u(17)) - m3(2)*cos(u(17)))/(A1+A2+A3);
    F=F(:);
    
    function F=two_tubesection(x,u,L1,L2,L3) 
      uhat11=0.5;%PAR(5)%1.0
      uhat12=0.0;
      uhat13=0.0;%PAR(3)
      
      uhat21=0.8;%PAR(6)%2.0
      uhat22=0.0;
      uhat23=0.0;%PAR(3)

      uhat31=1.0;%2.0
      uhat32=0.0;
      uhat33=0.00;%PAR(3)
  
     A1=1.0;
     A2=1.2;
     A3=1.4;
     C1=1.0/1.3;
     C2=1.2/1.3;
     C3=A3/1.3;
           	
    epsilon = 0.0;
    aa1=1000.0;
    aa2=1000.0;
    aa3=1000.0;
    epsilon = 0.0;
    aa1=1000.0;
    aa2=1000.0;
    aa3=1000.0;
   
      m2(1) = (u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0;
      m2(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0;
      m2(3) = (u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0;

%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%         the same 4 entries, but permuted and possibly negated.
%         Here are those entries

      entry21 =  2.0*(u(4)*u(12) + u(5)*u(13) + u(6)*u(14));
      entry22 =  2.0*(u(5)*u(12) - u(4)*u(13) + u(7)*u(14));
      entry23 =  2.0*(u(6)*u(12) - u(7)*u(13) - u(4)*u(14));
      entry24 =  2.0*(u(7)*u(12) + u(6)*u(13) - u(5)*u(14));

%       Components of D_1[q] n

      d21qn(1) =  entry21;
      d21qn(2) =  -entry22;
      d21qn(3) =  -entry23;
      d21qn(4) =  entry24;

%       Components of D_2[q] n

      d22qn(1) =  entry22;
      d22qn(2) =  entry21;
      d22qn(3) =  -entry24;
      d22qn(4) =  -entry23;

%       Components of D_3[q] n

      d23qn(1) =  entry23;
      d23qn(2) =  entry24;
      d23qn(3) =  entry21;
      d23qn(4) =  entry22;

%       Components of d1

      d21(1)= u(4)*u(4)-u(5)*u(5)-u(6)*u(6)+u(7)*u(7);
      d21(2)= 2.0*u(4)*u(5)+2.0*u(6)*u(7);
      d21(3)= 2.0*u(4)*u(6)-2.0*u(5)*u(7);

%       Components of d2

      d22(1)= 2.0*u(4)*u(5)-2.0*u(6)*u(7);
      d22(2)= -u(4)*u(4)+u(5)*u(5)-u(6)*u(6)+u(7)*u(7);
      d22(3)= 2.0*u(5)*u(6)+2.0*u(4)*u(7);

%       Components of d3

      d23(1)= 2.0*u(4)*u(6)+2.0*u(5)*u(7);
      d23(2)= 2.0*u(5)*u(6)-2.0*u(4)*u(7);
      d23(3)= -u(4)*u(4)-u(5)*u(5)+u(6)*u(6)+u(7)*u(7);

%    Components of v; note that v(3) has an extra term, due to
%    the fact that the unstressed state is assumed to have v=[0 0 1].

      v2(1)= epsilon/aa1*(d21(1)*u(12)+d21(2)*u(13)+d21(3)*u(14));
      v2(2)= epsilon/aa2*(d22(1)*u(12)+d22(2)*u(13)+d22(3)*u(14));
      v2(3)= epsilon/aa3*(d23(1)*u(12)+d23(2)*u(13)+d23(3)*u(14))+ 1.0;
      
	  k21=A1+A2;
	  k22=A1+A2;
      k23=C1+C2;  
    
      ubar21=(A1*uhat11 + A2*(uhat21*cos(u(15))- uhat22*sin(u(15)) ) )/(A1+A2);
      ubar22=(A1*uhat12 + A2*(uhat21*sin(u(15))+ uhat22*cos(u(15)) ) )/(A1+A2);
      ubar23= -u(16)/C1 + m2(3)*C2/(C1*(C1+C2)) + uhat13;   % Corrected addterm to uhat13
      
      strain2(1) = m2(1)/k21 + ubar21;
      strain2(2) = m2(2)/k22 + ubar22;
      strain2(3) = m2(3)/k23 + ubar23;  %

%    2-tube section    

      F(1) = L2*(v2(1)*d21(1)+v2(2)*d22(1)+v2(3)*d23(1));
      F(2) = L2*(v2(1)*d21(2)+v2(2)*d22(2)+v2(3)*d23(2));
      F(3) = L2*(v2(1)*d21(3)+v2(2)*d22(3)+v2(3)*d23(3));
      F(4) = L2*(strain2(1)*u(7)/2.0 - strain2(2)*u(6)/2.0 + strain2(3)*u(5)/2.0);
      F(5) = L2*(strain2(1)*u(6)/2.0 + strain2(2)*u(7)/2.0 - strain2(3)*u(4)/2.0);
      F(6) = L2*(-strain2(1)*u(5)/2.0 + strain2(2)*u(4)/2.0 + strain2(3)*u(7)/2.0);
      F(7) = L2*(-strain2(1)*u(4)/2.0 - strain2(2)*u(5)/2.0 - strain2(3)*u(6)/2.0);
	  F(8) =L2*(strain2(1)*u(11)/2.0 - strain2(2)*u(10)/2.0 + strain2(3)*u(9)/2.0 -v2(1)*d21qn(1)-v2(2)*d22qn(1)-v2(3)*d23qn(1));
	  F(9) =L2*(strain2(1)*u(10)/2.0 + strain2(2)*u(11)/2.0 - strain2(3)*u(8)/2.0 -v2(1)*d21qn(2)-v2(2)*d22qn(2)-v2(3)*d23qn(2));
	  F(10) =L2*(-strain2(1)*u(9)/2.0 + strain2(2)*u(8)/2.0 + strain2(3)*u(11)/2.0-v2(1)*d21qn(3)-v2(2)*d22qn(3)-v2(3)*d23qn(3));
	  F(11) =L2*(-strain2(1)*u(8)/2.0 - strain2(2)*u(9)/2.0 - strain2(3)*u(10)/2.0-v2(1)*d21qn(4)-v2(2)*d22qn(4)-v2(3)*d23qn(4));
      F(12) = L2*0.0;
      F(13) = L2*0.0;
      F(14) = L2*0.0;
      F(15)=  L2*( ( (C1+C2)*u(16)/(C1*C2) ) - (m2(3)/C1) + (uhat23-uhat13) );   % Alpha
      ater1=( m2(1)*sin(u(15)) - m2(2)*cos(u(15)) )*A2*uhat21/(A1+A2); %Corrected -m2(2) to +ve
      ater2=((uhat21*uhat11 + uhat12*uhat22)*sin(u(15)) + (uhat11*uhat22 - uhat12*uhat21)*cos(u(15)))*A1*A2/(A1+A2);
      F(16) = L2*(ater1+ater2);  
      F=F(:);

      
function F=one_tubesection(x,u,L1,L2,L3)
      uhat11=0.5;%PAR(5)%1.0
      uhat12=0.0;
      uhat13=0.0;%PAR(3)
      
  
     A1=1.0;
     A2=1.2;
     A3=1.4;
     C1=1.0/1.3;
     C2=1.2/1.3;
     C3=A3/1.3;
      
      ubar11=uhat11;
      ubar12=uhat12;
      ubar13=uhat13;
      
      k11=A1;
      k12=A1;
      k13=C1;
    
      epsilon = 0.0;
    aa1=1000.0;
    aa2=1000.0;
    aa3=1000.0;
      m1(1) = (u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0;
      m1(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0;
      m1(3) = (u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0;

%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%         the same 4 entries, but permuted and possibly negated.
%         Here are those entries

      entry11 =  2.0*(u(4)*u(12) + u(5)*u(13) + u(6)*u(14));
      entry12 =  2.0*(u(5)*u(12) - u(4)*u(13) + u(7)*u(14));
      entry13 =  2.0*(u(6)*u(12) - u(7)*u(13) - u(4)*u(14));
      entry14 =  2.0*(u(7)*u(12) + u(6)*u(13) - u(5)*u(14));

%       Components of D_1[q] n

      d11qn(1) =  entry11;
      d11qn(2) =  -entry12;
      d11qn(3) =  -entry13;
      d11qn(4) =  entry14;

%       Components of D_2[q] n

      d12qn(1) =  entry12;
      d12qn(2) =  entry11;
      d12qn(3) =  -entry14;
      d12qn(4) =  -entry13;

%       Components of D_3[q] n

      d13qn(1) =  entry13;
      d13qn(2) =  entry14;
      d13qn(3) =  entry11;
      d13qn(4) =  entry12;

%       Components of d1

      d11(1)= u(4)*u(4)-u(5)*u(5)-u(6)*u(6)+u(7)*u(7);
      d11(2)= 2.0*u(4)*u(5)+2.0*u(6)*u(7);
      d11(3)= 2.0*u(4)*u(6)-2.0*u(5)*u(7);

%       Components of d2

      d12(1)= 2.0*u(4)*u(5)-2.0*u(6)*u(7);
      d12(2)= -u(4)*u(4)+u(5)*u(5)-u(6)*u(6)+u(7)*u(7);
      d12(3)= 2.0*u(5)*u(6)+2.0*u(4)*u(7);

%       Components of d3
      d13(1)= 2.0*u(4)*u(6)+2.0*u(5)*u(7);
      d13(2)= 2.0*u(5)*u(6)-2.0*u(4)*u(7);
      d13(3)= -u(4)*u(4)-u(5)*u(5)+u(6)*u(6)+u(7)*u(7);
%       Components of v; note that v(3) has an extra term, due to
%         the fact that the unstressed state is assumed to have v=[0 0 1].
      v1(1)= epsilon/aa1*(d11(1)*u(12)+d11(2)*u(13)+d11(3)*u(14));
      v1(2)= epsilon/aa2*(d12(1)*u(12)+d12(2)*u(13)+d12(3)*u(14));
      v1(3)= epsilon/aa3*(d13(1)*u(12)+d13(2)*u(13)+d13(3)*u(14))+ 1.0;
%       Calculate components of the strain
      strain1(1) = m1(1)/k11 + uhat11;
      strain1(2) = m1(2)/k12 + uhat12;
      strain1(3) = m1(3)/k13 + uhat13;
      
     % strain1(1) = strain2(1)
     % strain1(2) = strain2(2)
     % strain1(3) = strain2(3)
%       Now calculate the right hand sides
%       Additional terms      
     
      F(1) = L1*(v1(1)*d11(1)+v1(2)*d12(1)+v1(3)*d13(1));
      F(2) = L1*(v1(1)*d11(2)+v1(2)*d12(2)+v1(3)*d13(2));
      F(3) = L1*(v1(1)*d11(3)+v1(2)*d12(3)+v1(3)*d13(3));
      F(4) = L1*(strain1(1)*u(7)/2.0 - strain1(2)*u(6)/2.0 + strain1(3)*u(5)/2.0);
      F(5) = L1*(strain1(1)*u(6)/2.0 + strain1(2)*u(7)/2.0 - strain1(3)*u(4)/2.0);
      F(6) = L1*(-strain1(1)*u(5)/2.0 + strain1(2)*u(4)/2.0 + strain1(3)*u(7)/2.0);
      F(7) = L1*(-strain1(1)*u(4)/2.0 - strain1(2)*u(5)/2.0 - strain1(3)*u(6)/2.0);
      F(8) = L1*(strain1(1)*u(11)/2.0- strain1(2)*u(10)/2.0 + strain1(3)*u(9)/2.0 -v1(1)*d11qn(1)-v1(2)*d12qn(1)-v1(3)*d13qn(1));
      F(9) = L1*(strain1(1)*u(10)/2.0+ strain1(2)*u(11)/2.0-strain1(3)*u(8)/2.0 -v1(1)*d11qn(2)-v1(2)*d12qn(2)-v1(3)*d13qn(2));
      F(10) = L1*(-strain1(1)*u(9)/2.0+ strain1(2)*u(8)/2.0+strain1(3)*u(11)/2.0 -v1(1)*d11qn(3)-v1(2)*d12qn(3)-v1(3)*d13qn(3));
      F(11) = L1*(-strain1(1)*u(8)/2.0- strain1(2)*u(9)/2.0-strain1(3)*u(10)/2.0 -v1(1)*d11qn(4)-v1(2)*d12qn(4)-v1(3)*d13qn(4));
      F(12) = L1*0.0;
      F(13) = L1*0.0;
      F(14) = L1*0.0;
       F=F(:);
        
    
  