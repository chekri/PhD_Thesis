function XYZ=Trajectory(u,Load,prevsol)
Th1=u(1);
Th2=u(2);
Th3=u(3);
L1= u(4);
L2= u(5);
L3= u(6);
%If inital guess is not given, start from a stright solution
if size(prevsol,1)==0
    prevsol = bvpinit(linspace(0,1,10),@mat4init);
    PrvSol = bvp4c(@(x,y)robot2rods(x,y,0,0,0,0.4,0.5,0.6,[0,0,0]), @(x,y)bcs(x,y,0,3.1415,0,0.4,0.5,0.6,Load), prevsol);
    prevsol=PrvSol;
else 
    prevsol=bvpinit(prevsol,[0 1]);
end
opts = bvpset('RelTol',1e-05,'AbsTol',1e-05);

XYZ=bvp4c(@(x,y)robot2rods(x,y,Th1,Th2,Th3,L1,L2,L3,Load), @(x,y)bcs(x,y,Th1,Th2,Th3,L1,L2,L3,Load), prevsol, opts);
      
 function yinit = mat4init(s) % initial guess function
    L1=0.0001;
    L2=0.0001;
    L3=0.0001;
    yinit=zeros(1,48);
    yinit(3) = L3 + L2*s;
    yinit(7) = 1.0;
    yinit(19) = L3 + L2 + L1*s;
    yinit(23) = 1.0;
    
    yinit(33) = L3*s;
    yinit(37) = 1.0;
        
function F = robot2rods(s,u,th1,th2,th3,ul1,ul2,ul3,Load)
 
      uhat11=0.5;%PAR(5)%1.0
      uhat12=0.0;
      uhat13=0.0;%PAR(3)
      
      uhat21=0.8;%PAR(6)%2.0
      uhat22=0.0;
      uhat23=0.0;%PAR(3)

      uhat31=1.0;%2.0
      uhat32=0.0;
      uhat33=0.0;%PAR(3)
  
 %     uhat11=0.5;%PAR(5)%1.0
      uhat12=0.0;
      uhat13=0.0;%PAR(3)
      
 %     uhat21=1.5;%PAR(6)%2.0
      uhat22=0.0;
      uhat23=0.0;%PAR(3)

 %     uhat31=2.0;%2.0
      uhat32=0.0;
      uhat33=0.0;%PAR(3)
  
      
      
      
     A1=1.0;
     A2=1.2;
     A3=1.4;
     C1=A1/1.3;
     C2=A2/1.3;
     C3=A3/1.3;
     
     L1=ul1;%PAR(1)
     L2=ul2;%PAR(4)
     L3=ul3;%PAR(10)
      	
   
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

%       Components of D_3[q] n

      d23qn(1) =  entry23;
      d23qn(2) =  entry24;
      d23qn(3) =  entry21;
      d23qn(4) =  entry22;


%       Components of d3

      d23(1)= 2.0*u(4)*u(6)+2.0*u(5)*u(7);
      d23(2)= 2.0*u(5)*u(6)-2.0*u(4)*u(7);
      d23(3)= -u(4)*u(4)-u(5)*u(5)+u(6)*u(6)+u(7)*u(7);

      
	  k21=A1+A2;
	  k22=A1+A2;
      k23=C1+C2;  
    
      ubar21=(A1*uhat11 + A2*(uhat21*cos(u(15)) - uhat22*sin(u(15)) ) )/(A1+A2);
      ubar22=(A1*uhat12 + A2*(uhat21*sin(u(15)) + uhat22*cos(u(15)) ) )/(A1+A2);
      ubar23= -u(16)/C1 + m2(3)*C2/(C1*(C1+C2)) + uhat13;   % Corrected addterm to uhat13
      
      strain2(1) = m2(1)/k21 + ubar21;
      strain2(2) = m2(2)/k22 + ubar22;
      strain2(3) = m2(3)/k23 + ubar23;  %

%    2-tube section    

      F(1) = L2*(d23(1));
      F(2) = L2*(d23(2));
      F(3) = L2*(d23(3));
      F(4) = L2*(strain2(1)*u(7)/2.0 - strain2(2)*u(6)/2.0 + strain2(3)*u(5)/2.0);
      F(5) = L2*(strain2(1)*u(6)/2.0 + strain2(2)*u(7)/2.0 - strain2(3)*u(4)/2.0);
      F(6) = L2*(-strain2(1)*u(5)/2.0 + strain2(2)*u(4)/2.0 + strain2(3)*u(7)/2.0);
      F(7) = L2*(-strain2(1)*u(4)/2.0 - strain2(2)*u(5)/2.0 - strain2(3)*u(6)/2.0);
	  F(8) =L2*(strain2(1)*u(11)/2.0 - strain2(2)*u(10)/2.0 + strain2(3)*u(9)/2.0 -d23qn(1));
	  F(9) =L2*(strain2(1)*u(10)/2.0 + strain2(2)*u(11)/2.0 - strain2(3)*u(8)/2.0 -d23qn(2));
	  F(10) =L2*(-strain2(1)*u(9)/2.0 + strain2(2)*u(8)/2.0 + strain2(3)*u(11)/2.0-d23qn(3));
	  F(11) =L2*(-strain2(1)*u(8)/2.0 - strain2(2)*u(9)/2.0 - strain2(3)*u(10)/2.0-d23qn(4));
      F(12) = L2*0.0;
      F(13) = L2*0.0;
      F(14) = L2*0.0;
      F(15)=  L2*( ( (C1+C2)*u(16)/(C1*C2) ) - (m2(3)/C1) + (uhat23-uhat13) );   % Alpha
      ater1=( m2(1)*sin(u(15)) - m2(2)*cos(u(15)) )*A2*uhat21/(A1+A2); %Corrected -m2(2) to +ve
      ater2=((uhat21*uhat11 + uhat12*uhat22)*sin(u(15)) + (uhat11*uhat22 - uhat12*uhat21)*cos(u(15)))*A1*A2/(A1+A2);
      F(16) = L2*(ater1+ater2);  
   
      
      k11=A1;
      k12=A1;
      k13=C1;
     
      m1(1) = (u(23)*u(24) + u(22)*u(25) - u(21)*u(26) - u(20)*u(27))/2.0;
      m1(2) = (-u(22)*u(24) + u(23)*u(25) + u(20)*u(26) - u(21)*u(27))/2.0;
      m1(3) = (u(21)*u(24) - u(20)*u(25) + u(23)*u(26) - u(22)*u(27))/2.0;

%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%         the same 4 entries, but permuted and possibly negated.
%         Here are those entries

      entry11 =  2.0*(u(20)*u(28) + u(21)*u(29) + u(22)*u(30));
      entry12 =  2.0*(u(21)*u(28) - u(20)*u(29) + u(23)*u(30));
      entry13 =  2.0*(u(22)*u(28) - u(23)*u(29) - u(20)*u(30));
      entry14 =  2.0*(u(23)*u(28) + u(22)*u(29) - u(21)*u(30));


%       Components of D_3[q] n

      d13qn(1) =  entry13;
      d13qn(2) =  entry14;
      d13qn(3) =  entry11;
      d13qn(4) =  entry12;


%       Components of d3
      d13(1)= 2.0*u(20)*u(22)+2.0*u(21)*u(23);
      d13(2)= 2.0*u(21)*u(22)-2.0*u(20)*u(23);
      d13(3)= -u(20)*u(20)-u(21)*u(21)+u(22)*u(22)+u(23)*u(23);

%       Calculate components of the strain
      strain1(1) = m1(1)/k11 + uhat11;
      strain1(2) = m1(2)/k12 + uhat12;
      strain1(3) = m1(3)/k13 + uhat13;
      
%       Now calculate the right hand sides
%       Additional terms      
     
      F(17) = L1*(d13(1));
      F(18) = L1*(d13(2));
      F(19) = L1*(d13(3));
      F(20) = L1*(strain1(1)*u(23)/2.0 - strain1(2)*u(22)/2.0 + strain1(3)*u(21)/2.0);
      F(21) = L1*(strain1(1)*u(22)/2.0 + strain1(2)*u(23)/2.0 - strain1(3)*u(20)/2.0);
      F(22) = L1*(-strain1(1)*u(21)/2.0 + strain1(2)*u(20)/2.0 + strain1(3)*u(23)/2.0);
      F(23) = L1*(-strain1(1)*u(20)/2.0 - strain1(2)*u(21)/2.0 - strain1(3)*u(22)/2.0);
      F(24) = L1*(strain1(1)*u(27)/2.0- strain1(2)*u(26)/2.0 + strain1(3)*u(25)/2.0 -d13qn(1));
      F(25) = L1*(strain1(1)*u(26)/2.0+ strain1(2)*u(27)/2.0-strain1(3)*u(24)/2.0 -d13qn(2));
      F(26) = L1*(-strain1(1)*u(25)/2.0+ strain1(2)*u(24)/2.0+strain1(3)*u(27)/2.0 -d13qn(3));
      F(27) = L1*(-strain1(1)*u(24)/2.0- strain1(2)*u(25)/2.0-strain1(3)*u(26)/2.0 -d13qn(4));
      F(28) = L1*0.0;
      F(29) = L1*0.0;
      F(30) = L1*0.0;
    
      m3(1) = ( u(23+14)*u(24+14) + u(22+14)*u(25+14) - u(21+14)*u(26+14) - u(20+14)*u(27+14))/2.0;
      m3(2) = (-u(22+14)*u(24+14) + u(23+14)*u(25+14) + u(20+14)*u(26+14) - u(21+14)*u(27+14))/2.0;
      m3(3) = ( u(21+14)*u(24+14) - u(20+14)*u(25+14) + u(23+14)*u(26+14) - u(22+14)*u(27+14))/2.0;

%%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%%         the same 4 entries, but permuted and possibly negated.
%%         Here are those entries

      entry31 =  2.0*(u(20+14)*u(28+14) + u(21+14)*u(29+14) + u(22+14)*u(30+14));
      entry32 =  2.0*(u(21+14)*u(28+14) - u(20+14)*u(29+14) + u(23+14)*u(30+14));
      entry33 =  2.0*(u(22+14)*u(28+14) - u(23+14)*u(29+14) - u(20+14)*u(30+14));
      entry34 =  2.0*(u(23+14)*u(28+14) + u(22+14)*u(29+14) - u(21+14)*u(30+14));%%Done

%%       Components of D_3[q] n

      d33qn(1) =  entry33;
      d33qn(2) =  entry34;
      d33qn(3) =  entry31;
      d33qn(4) =  entry32;


%%       Components of d3
      d33(1)= 2.0*u(20+14)*u(22+14) + 2.0*u(21+14)*u(23+14);
      d33(2)= 2.0*u(21+14)*u(22+14) - 2.0*u(20+14)*u(23+14);
      d33(3)= -u(20+14)*u(20+14) - u(21+14)*u(21+14) + u(22+14)*u(22+14) + u(23+14)*u(23+14);
%%       Components of v; note that v(3) has an extra term, due to
%%         the fact that the unstressed state is assumed to have v=[0 0 1].

     ab1=-C2*((u(46) + u(48))/C1 - m3(3)/C1  + u(46)/C2 + uhat23- uhat13 - uhat23); 
     ab2=-C3*((u(46) + u(48))/C1 - m3(3)/C1  + u(48)/C3 + uhat33- uhat13 - uhat33);
    

    ub31=(A1*uhat11 + A2*(uhat21*cos(u(45))- uhat22*sin(u(45)) ) + A3*(uhat31*cos(u(47))- uhat32*sin(u(47))))/(A1+A2+A3);
	ub32=(A1*uhat12 + A2*(uhat21*sin(u(45))+ uhat22*cos(u(45)) ) + A3*(uhat31*sin(u(47))+ uhat32*cos(u(47))))/(A1+A2+A3);
	ub33= (ab1 + ab2 + C1*uhat13)/(C1+C2+C3);

     
	strain3(1) = m3(1)/(A1+A2+A3) + ub31;
	strain3(2) = m3(2)/(A1+A2+A3) + ub32;
	%strain3(3) = m3(3)/(C1+C2+C3) + ub33;  
	strain3(3) = (m3(3)- u(46) -u(48)) /C1 + uhat13;  

	F(31) = L3*(d33(1));
    F(32) = L3*(d33(2));
	F(33) = L3*(d33(3));
      
	F(34) = L3*( strain3(1)*u(23+14)/2.0 - strain3(2)*u(22+14)/2.0 + strain3(3)*u(21+14)/2.0);
    F(35) = L3*( strain3(1)*u(22+14)/2.0 + strain3(2)*u(23+14)/2.0 - strain3(3)*u(20+14)/2.0);
    F(36) = L3*(-strain3(1)*u(21+14)/2.0 + strain3(2)*u(20+14)/2.0 + strain3(3)*u(23+14)/2.0);
    F(37) = L3*(-strain3(1)*u(20+14)/2.0 - strain3(2)*u(21+14)/2.0 - strain3(3)*u(22+14)/2.0);
     
	F(38) = L3*( strain3(1)*u(27+14)/2.0 - strain3(2)*u(26+14)/2.0 + strain3(3)*u(25+14)/2.0 -d33qn(1));
	F(39) = L3*( strain3(1)*u(26+14)/2.0 + strain3(2)*u(27+14)/2.0 - strain3(3)*u(24+14)/2.0 -d33qn(2));
	F(40) = L3*(-strain3(1)*u(25+14)/2.0 + strain3(2)*u(24+14)/2.0 + strain3(3)*u(27+14)/2.0 -d33qn(3));
	F(41) = L3*(-strain3(1)*u(24+14)/2.0 - strain3(2)*u(25+14)/2.0 - strain3(3)*u(26+14)/2.0 -d33qn(4));
     
	F(42) = L3*0.0;
	F(43) = L3*0.0;
	F(44) = L3*0.0;
   	F(45)= L3*((u(46) + u(48))/C1 - m3(3)/C1  + u(46)/C2 + uhat23- uhat13);  %%%Alfa2=u(45)
	F(46)= L3*(A2*uhat21*(A3*uhat31*sin(u(45)-u(47)) + A1*uhat11*sin(u(45))))/(A1+A2+A3) + L3*A2*uhat21*(m3(1)*sin(u(45)) - m3(2)*cos(u(45)))/(A1+A2+A3);
    F(47)= L3*((u(46) + u(48))/C1 - m3(3)/C1  + u(48)/C3 + uhat33- uhat13);  %%%Alfa3=u(47)
	F(48)= L3*(A3*uhat31*(A2*uhat21*sin(u(47)-u(45)) + A1*uhat11*sin(u(47))))/(A1+A2+A3) + A3*L3*uhat31*(m3(1)*sin(u(47)) - m3(2)*cos(u(47)))/(A1+A2+A3);
      

function FB = bcs(u0,u1,Th1,Th2,Th3,L1,L2,L3,Load) 
      nz=Load(1);%PAR(7)*0.0
      nx=Load(2);%PAR(7)*sin(3.1415/4)
      ny=Load(3);%PAR(7)*cos(3.1415/4)
      ang1=Th1;%PAR()
      ang2=Th2;%PAR(2)
      ang3=Th3;%PAR(12)
      
      e1=0.0; 	
   	  e2=0.0;
   	  e3=0.0;
   
     % q1,q2,q3,q4, u20,u21,u22,u23   
   	  d111= u1(20)*u1(20) - u1(21)*u1(21) - u1(22)*u1(22) + u1(23)*u1(23);
      d112= 2.0*u1(20)*u1(21) + 2.0*u1(22)*u1(23);
      d113= 2.0*u1(20)*u1(22) - 2.0*u1(21)*u1(23);
   
   	  d121= 2.0*u1(20)*u1(21) - 2.0*u1(22)*u1(23);
      d122= -u1(20)*u1(20) + u1(21)*u1(21)- u1(22)*u1(22) + u1(23)*u1(23);
      d123= 2.0*u1(21)*u1(22) + 2.0*u1(20)*u1(23);

      d131= 2.0*u1(20)*u1(22) + 2.0*u1(21)*u1(23);
      d132= 2.0*u1(21)*u1(22) - 2.0*u1(20)*u1(23);
      d133= -u1(20)*u1(20) - u1(21)*u1(21) + u1(22)*u1(22) + u1(23)*u1(23);
   
     % r X N
      Ndir1=nx*d111 + ny*d112 + nz*d113;  
      Ndir2=nx*d121 + ny*d122 + nz*d123; 
      Ndir3=nx*d131 + ny*d132 + nz*d133;  
      %Matching Conditions
           
      FB(1)=u0(1) - u1(31);
      FB(2)=u0(2) - u1(32);
      FB(3)=u0(3) - u1(33);
	  	  
      FB(4)=u0(4) - u1(34);
      FB(5)=u0(5) - u1(35); 
      FB(6)=u0(6) - u1(36); 
      FB(7)=u0(7) - u1(37); % Dirichlet BCs coupled
      
      FB(8) = u1(8) -u0(24);
      FB(9) = u1(9) -u0(25); 
      FB(10)=u1(10) -u0(26); 
      FB(11)=u1(11) -u0(27);
	       
     
      FB(12)=u1(12) - u0(28); 
      FB(13)=u1(13) - u0(29);
      FB(14)=u1(14) - u0(30); %Neumann Coupled between the 2-tube and 1-tube
      
        
      FB(15)=u0(15) - u1(45); %Dirichlet Condition on the end of tube twist
      FB(16)=u1(16);         %Neumann Condition on the end of tube twist
      
      
      FB(17)=u0(17) - u1(1);
      FB(18)=u0(18) - u1(2);
      FB(19)=u0(19) - u1(3);
      
      
      FB(20)=u0(20) - u1(4);
      FB(21)=u0(21) - u1(5);
      FB(22)=u0(22) - u1(6);
      FB(23)=u0(23) - u1(7); %Dirichlet BCs coupled at the boundary between 2 tube section and 1 tube

     
      FB(24)= ( u1(16+7)*u1(16+8) + u1(16+6)*u1(16+9) - u1(16+5)*u1(16+10) - u1(16+4)*u1(16+11))/2.0 +(e2*Ndir3 - e3*Ndir2);
      FB(25)= (-u1(16+6)*u1(16+8) + u1(16+7)*u1(16+9) + u1(16+4)*u1(16+10) - u1(16+5)*u1(16+11))/2.0 +(e3*Ndir1 - e1*Ndir3);%
      FB(26)= ( u1(16+5)*u1(16+8) - u1(16+4)*u1(16+9) + u1(16+7)*u1(16+10) - u1(16+6)*u1(16+11))/2.0 +(e1*Ndir2 - e2*Ndir1);
	     
	  FB(27)=  u0(24)*u0(20) + u0(25)*u0(21)+ u0(26)*u0(22)+ u0(27)*u0(23) + 2.0*(u0(17)*u0(28) + u0(18)*u0(29) + u0(19)*u0(30)); %Neumann BCs on the tip
      
      FB(28)=u1(28) + nx;
      FB(29)=u1(29) + ny;
      FB(30)=u1(30) + nz;%Neumann BCs on the tip
                  
      FB(31)=u0(31); 
      FB(32)=u0(32);
      FB(33)=u0(33);
      
      FB(34)=u0(34);
      FB(35)=u0(35);
      FB(36)=u0(36) - sin(ang1/2);
      FB(37)=u0(37) - cos(ang1/2);  % Dirichlet BCs coupled 
   
       
      FB(38)=u1(42) - u0(12);
      FB(39)=u1(43) - u0(13);
      FB(40)=u1(44) - u0(14);
      
      FB(41)=u1(38) - u0(8); %mb31 - mb21
      FB(42)=u1(39) - u0(9); %mb32 - mb22
      FB(43)=u1(40) - u0(10);%mb33 - mb23
      FB(44)=u1(41) - u0(11);%Neumann Conditions Coupled
	
      FB(45)=u0(45) - (ang2 - ang1); 
      FB(46)=u1(46) - u0(16);
      FB(47)=u0(47) - (ang3 - ang1);% Dirichlet BCs coupled 
  
      FB(48)=u1(48);