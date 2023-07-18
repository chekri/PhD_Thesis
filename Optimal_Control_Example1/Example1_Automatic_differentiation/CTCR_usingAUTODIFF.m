function [JC,JC1] = CTCR_usingAUTODIFF(y,u)
      syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 u31 u32 u33 u34 u35 u36 u37 u38 u39 u40 u41 u42 u43 u44 u45 u46 u47 u48
      syms A1 A2 A3 C1 C2 C3 uhat11 uhat12 uhat13 uhat21 uhat22 uhat23 uhat31 uhat32 uhat33 L1 L2 L3
   	
      m2(1) = (u7*u8 + u6*u9 - u5*u10 - u4*u11)/2.0;
      m2(2) = (-u6*u8 + u7*u9 + u4*u10 - u5*u11)/2.0;
      m2(3) = (u5*u8 - u4*u9 + u7*u10 - u6*u11)/2.0;

%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%         the same 4 entries, but permuted and possibly negated.
%         Here are those entries

      entry21 =  2.0*(u4*u12 + u5*u13 + u6*u14);
      entry22 =  2.0*(u5*u12 - u4*u13 + u7*u14);
      entry23 =  2.0*(u6*u12 - u7*u13 - u4*u14);
      entry24 =  2.0*(u7*u12 + u6*u13 - u5*u14);

%       Components of D_3[q] n

      d23qn(1) =  entry23;
      d23qn(2) =  entry24;
      d23qn(3) =  entry21;
      d23qn(4) =  entry22;

%       Components of d1

      d21(1)= u4*u4-u5*u5-u6*u6+u7*u7;
      d21(2)= 2.0*u4*u5+2.0*u6*u7;
      d21(3)= 2.0*u4*u6-2.0*u5*u7;

%       Components of d2

      d22(1)= 2.0*u4*u5-2.0*u6*u7;
      d22(2)= -u4*u4+u5*u5-u6*u6+u7*u7;
      d22(3)= 2.0*u5*u6+2.0*u4*u7;

%       Components of d3

      d23(1)= 2.0*u4*u6+2.0*u5*u7;
      d23(2)= 2.0*u5*u6-2.0*u4*u7;
      d23(3)= -u4*u4-u5*u5+u6*u6+u7*u7;

      
	  k21=A1+A2;
	  k22=A1+A2;
      k23=C1+C2;  
    
      ubar21=(A1*uhat11 + A2*(uhat21*cos(u15) - uhat22*sin(u15) ) )/(A1+A2);
      ubar22=(A1*uhat12 + A2*(uhat21*sin(u15) + uhat22*cos(u15) ) )/(A1+A2);
      ubar23= -u16/C1 + m2(3)*C2/(C1*(C1+C2)) + uhat13;   % Corrected addterm to uhat13
      
      strain2(1) = m2(1)/k21 + ubar21;
      strain2(2) = m2(2)/k22 + ubar22;
      strain2(3) = m2(3)/k23 + ubar23;  %

%    2-tube section    

      F(1) = L2*(d23(1));
      F(2) = L2*(d23(2));
      F(3) = L2*(d23(3));
      F(4) = L2*(strain2(1)*u7/2.0 - strain2(2)*u6/2.0 + strain2(3)*u5/2.0);
      F(5) = L2*(strain2(1)*u6/2.0 + strain2(2)*u7/2.0 - strain2(3)*u4/2.0);
      F(6) = L2*(-strain2(1)*u5/2.0 + strain2(2)*u4/2.0 + strain2(3)*u7/2.0);
      F(7) = L2*(-strain2(1)*u4/2.0 - strain2(2)*u5/2.0 - strain2(3)*u6/2.0);
	  F(8) =L2*(strain2(1)*u11/2.0 - strain2(2)*u10/2.0 + strain2(3)*u9/2.0 -d23qn(1));
	  F(9) =L2*(strain2(1)*u10/2.0 + strain2(2)*u11/2.0 - strain2(3)*u8/2.0 -d23qn(2));
	  F(10) =L2*(-strain2(1)*u9/2.0 + strain2(2)*u8/2.0 + strain2(3)*u11/2.0-d23qn(3));
	  F(11) =L2*(-strain2(1)*u8/2.0 - strain2(2)*u9/2.0 - strain2(3)*u10/2.0-d23qn(4));
      F(12) = L2*0.0;
      F(13) = L2*0.0;
      F(14) = L2*0.0;
      F(15)=  L2*( ( (C1+C2)*u16/(C1*C2) ) - (m2(3)/C1) + (uhat23-uhat13) );   % Alpha
      ater1=( m2(1)*sin(u15) - m2(2)*cos(u15) )*A2*uhat21/(A1+A2); %Corrected -m2(2) to +ve
      ater2=((uhat21*uhat11 + uhat12*uhat22)*sin(u15) + (uhat11*uhat22 - uhat12*uhat21)*cos(u15))*A1*A2/(A1+A2);
      F(16) = L2*(ater1+ater2);  
   
      
      k11=A1;
      k12=A1;
      k13=C1;
     
      m1(1) = (u23*u24 + u22*u25 - u21*u26 - u20*u27)/2.0;
      m1(2) = (-u22*u24 + u23*u25 + u20*u26 - u21*u27)/2.0;
      m1(3) = (u21*u24 - u20*u25 + u23*u26 - u22*u27)/2.0;

%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%         the same 4 entries, but permuted and possibly negated.
%         Here are those entries

      entry11 =  2.0*(u20*u28 + u21*u29 + u22*u30);
      entry12 =  2.0*(u21*u28 - u20*u29 + u23*u30);
      entry13 =  2.0*(u22*u28 - u23*u29 - u20*u30);
      entry14 =  2.0*(u23*u28 + u22*u29 - u21*u30);

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

      d11(1)= u20*u20-u21*u21-u22*u22+u23*u23;
      d11(2)= 2.0*u20*u21+2.0*u22*u23;
      d11(3)= 2.0*u20*u22-2.0*u21*u23;

%       Components of d2

      d12(1)= 2.0*u20*u21-2.0*u22*u23;
      d12(2)= -u20*u20+u21*u21-u22*u22+u23*u23;
      d12(3)= 2.0*u21*u22+2.0*u20*u23;

%       Components of d3
      d13(1)= 2.0*u20*u22+2.0*u21*u23;
      d13(2)= 2.0*u21*u22-2.0*u20*u23;
      d13(3)= -u20*u20-u21*u21+u22*u22+u23*u23;

%       Calculate components of the strain
      strain1(1) = m1(1)/k11 + uhat11;
      strain1(2) = m1(2)/k12 + uhat12;
      strain1(3) = m1(3)/k13 + uhat13;
      
%       Now calculate the right hand sides
%       Additional terms      
     
      F(17) = L1*(d13(1));
      F(18) = L1*(d13(2));
      F(19) = L1*(d13(3));
      F(20) = L1*(strain1(1)*u23/2.0 - strain1(2)*u22/2.0 + strain1(3)*u21/2.0);
      F(21) = L1*(strain1(1)*u22/2.0 + strain1(2)*u23/2.0 - strain1(3)*u20/2.0);
      F(22) = L1*(-strain1(1)*u21/2.0 + strain1(2)*u20/2.0 + strain1(3)*u23/2.0);
      F(23) = L1*(-strain1(1)*u20/2.0 - strain1(2)*u21/2.0 - strain1(3)*u22/2.0);
      F(24) = L1*(strain1(1)*u27/2.0 - strain1(2)*u26/2.0 + strain1(3)*u25/2.0 -d13qn(1));
      F(25) = L1*(strain1(1)*u26/2.0 + strain1(2)*u27/2.0 - strain1(3)*u24/2.0 -d13qn(2));
      F(26) = L1*(-strain1(1)*u25/2.0+ strain1(2)*u24/2.0 + strain1(3)*u27/2.0 -d13qn(3));
      F(27) = L1*(-strain1(1)*u24/2.0- strain1(2)*u25/2.0 -strain1(3)*u26/2.0 -d13qn(4));
      F(28) = L1*0.0;
      F(29) = L1*0.0;
      F(30) = L1*0.0;
    
      
      
      
      m3(1) = ( u37*u38 + u36*u39 - u35*u40 - u34*u41)/2.0;
      m3(2) = (-u36*u38 + u37*u39 + u34*u40 - u35*u41)/2.0;
      m3(3) = ( u35*u38 - u34*u39 + u37*u40 - u36*u41)/2.0;

%%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
%%         the same 4 entries, but permuted and possibly negated.
%%         Here are those entries

      entry31 =  2.0*(u34*u42 + u35*u43 + u36*u44);
      entry32 =  2.0*(u35*u42 - u34*u43 + u37*u44);
      entry33 =  2.0*(u36*u42 - u37*u43 - u34*u44);
      entry34 =  2.0*(u37*u42 + u36*u43 - u35*u44);

%%       Components of D_3[q] n

      d33qn(1) =  entry33;
      d33qn(2) =  entry34;
      d33qn(3) =  entry31;
      d33qn(4) =  entry32;


%%       Components of d3
      d33(1)= 2.0*u34*u36 + 2.0*u35*u37;
      d33(2)= 2.0*u35*u36 - 2.0*u34*u37;
      d33(3)= -u34*u34 - u35*u35 + u36*u36 + u37*u37;
%%       Components of v; note that v(3) has an extra term, due to
%%         the fact that the unstressed state is assumed to have v=[0 0 1].

     ab1=-C2*((u46 + u48)/C1 - m3(3)/C1  + u46/C2 + uhat23- uhat13 - uhat23); 
     ab2=-C3*((u46 + u48)/C1 - m3(3)/C1  + u48/C3 + uhat33- uhat13 - uhat33);
    

    ub31=(A1*uhat11 + A2*(uhat21*cos(u45)- uhat22*sin(u45) ) + A3*(uhat31*cos(u47)- uhat32*sin(u47)))/(A1+A2+A3);
	ub32=(A1*uhat12 + A2*(uhat21*sin(u45)+ uhat22*cos(u45) ) + A3*(uhat31*sin(u47)+ uhat32*cos(u47)))/(A1+A2+A3);
	ub33= (ab1 + ab2 + C1*uhat13)/(C1+C2+C3);

     
	strain3(1) = m3(1)/(A1+A2+A3) + ub31;
	strain3(2) = m3(2)/(A1+A2+A3) + ub32;
	strain3(3) = (m3(3)- u46 -u48) /C1 + uhat13;  

	F(31) = L3*(d33(1));
    F(32) = L3*(d33(2));
	F(33) = L3*(d33(3));
      
	F(34) = L3*( strain3(1)*u37/2.0 - strain3(2)*u36/2.0 + strain3(3)*u35/2.0);
    F(35) = L3*( strain3(1)*u36/2.0 + strain3(2)*u37/2.0 - strain3(3)*u34/2.0);
    F(36) = L3*(-strain3(1)*u35/2.0 + strain3(2)*u34/2.0 + strain3(3)*u37/2.0);
    F(37) = L3*(-strain3(1)*u34/2.0 - strain3(2)*u35/2.0 - strain3(3)*u36/2.0);
     
	F(38) = L3*( strain3(1)*u41/2.0 - strain3(2)*u40/2.0 + strain3(3)*u39/2.0 -d33qn(1));
	F(39) = L3*( strain3(1)*u40/2.0 + strain3(2)*u41/2.0 - strain3(3)*u38/2.0 -d33qn(2));
	F(40) = L3*(-strain3(1)*u39/2.0 + strain3(2)*u38/2.0 + strain3(3)*u41/2.0 -d33qn(3));
	F(41) = L3*(-strain3(1)*u38/2.0 - strain3(2)*u39/2.0 - strain3(3)*u40/2.0 -d33qn(4));
     
	F(42) = L3*0.0;
	F(43) = L3*0.0;
	F(44) = L3*0.0;
   	F(45)= L3*((u46 + u48)/C1 - m3(3)/C1  + u46/C2 + uhat23- uhat13);  %%%Alfa2=u(45)
	F(46)= L3*(A2*uhat21*(A3*uhat31*sin(u45-u47) + A1*uhat11*sin(u45)))/(A1+A2+A3) + L3*A2*uhat21*(m3(1)*sin(u45) - m3(2)*cos(u45))/(A1+A2+A3);
    F(47)= L3*((u46 + u48)/C1 - m3(3)/C1  + u48/C3 + uhat33- uhat13);  %%%Alfa3=u(47)
	F(48)= L3*(A3*uhat31*(A2*uhat21*sin(u47-u45) + A1*uhat11*sin(u47)))/(A1+A2+A3) + A3*L3*uhat31*(m3(1)*sin(u47) - m3(2)*cos(u47))/(A1+A2+A3);
      

  %  [dydu1,dydu2,dydu3,dydu4,dydu5,dydu6,dydu7] = dlgradient(F(1),u1,u2,u3,u4,u5,u6,u7);
  %  DFDU1=[dydu1,dydu2,dydu3,dydu4,dydu5,dydu6,dydu7];
  

  JC=    jacobian(F,[u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44,u45,u46,u47,u48]);
  %JC1=subs(JC,u1,y(1),u2,y(2),u3,y(3),u4,y(4),u5,y(5),u6,y(6),u7,y(7),u8,y(8),u9,y(9),u10,y(10),u11,y(11),u12,y(12),u13,y(13),u14,y(14),u15,y(15),u16,y(16),u17,y(17),u18,y(18),u19,y(19),u20,y(20),u21,y(21),u22,y(22),u23,y(23),u24,y(24),u25,y(25),u26,y(26),u27,y(27),u28,y(28),u29,y(29),u30,y(30),u31,y(31),u32,y(32),u33,y(33),u34,y(34),u35,y(35),u36,y(36),u37,y(37),u38,y(38),u39,y(39),u40,y(40),u41,y(41),u42,y(42),u43,y(43),u44,y(44),u45,y(45),u46,y(46),u47,y(47),u48,y(48))
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
     C1=A1/1.3;
     C2=A2/1.3;
     C3=A3/1.3;
  
     L1=u(4);%PAR(1)
     L2=u(5);%PAR(4)
     L3=u(6);%PAR(10)
 
   JC1=subs(JC);
%  [u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44,u45,u46,u47,u48]=y';
  JC1=subs(JC1,[u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44,u45,u46,u47,u48],vpa(y'));
end

