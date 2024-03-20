
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION entry11,entry12,entry13,entry14,k21,k22,k23,k11,k12,k13
      DOUBLE PRECISION ubar11,ubar12,ubar13,A1,A2,C1,C2,addterm21
      DOUBLE PRECISION pi,L1,L2,L3,ab1,ab2,ab3
      DOUBLE PRECISION kb31,kb32,kb33,ub31,ub32,ub33,entry31,entry32,entry33,entry34,A3,C3
      DOUBLE PRECISION m1(3),strain1(3),d11(3),d12(3),d13(3),v1(3),d11qn(4),d12qn(4),d13qn(4)
      DOUBLE PRECISION m2(3),strain2(3),d21(3),d22(3),d23(3),v2(3),d21qn(4),d22qn(4),d23qn(4)
      DOUBLE PRECISION m3(3),strain3(3),d31(3),d32(3),d33(3),v3(3),d31qn(4),d32qn(4),d33qn(4)
      DOUBLE PRECISION uhat11,uhat12,uhat13,uhat21,uhat22,uhat23,uhat31,uhat32,uhat33,ubar21,ubar22,ubar23
      DOUBLE PRECISION ater1,ater2,entry21,entry22,entry23,entry24


! Intrinsic curvatures of the constitute tubes 
      uhat11=PAR(5)
      uhat12=0.0d0
      uhat13=PAR(3)
      
      uhat21=PAR(6)
      uhat22=0.0d0
      uhat23=PAR(3)

      uhat31=PAR(11)
      uhat32=0.0d0
      uhat33=PAR(3)
      
    !Bending and Torsion stiffnesses of the constitute tubes from ()
     A1=0.03487
     A2=0.09422
     A3=3.4
     C1=A1/1.3
     C2=A2/1.3
     C3=A3/1.3
     
   !Lengths of the CTCR sections  
     L2=PAR(1)
     L1=PAR(4)
     L3=PAR(10)
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Governing equations corresponding to 3-tube section   

      k21=A1+A2
      k22=A1+A2
      k23=C1+C2     


      m2(1) = (u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0d0
      m2(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0d0
      m2(3) = (u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0d0

!       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
!         the same 4 entries, but permuted and possibly negated.
!         Here are those entries

      entry21 =  2.0D0*(u(4)*u(12) + u(5)*u(13) + U(6)*u(14))
      entry22 =  2.0D0*(u(5)*u(12) - u(4)*u(13) + U(7)*u(14))
      entry23 =  2.0D0*(u(6)*u(12) - u(7)*u(13) - U(4)*u(14))
      entry24 =  2.0D0*(u(7)*u(12) + u(6)*u(13) - U(5)*u(14))


!       Components of D_3[q] n

      d23qn(1) =  entry23
      d23qn(2) =  entry24
      d23qn(3) =  entry21
      d23qn(4) =  entry22


!       Components of d3

      d23(1)= 2.0d0*u(4)*u(6)+2.0d0*u(5)*u(7)
      d23(2)= 2.0d0*u(5)*u(6)-2.0d0*u(4)*u(7)
      d23(3)= -u(4)*u(4)-u(5)*u(5)+u(6)*u(6)+u(7)*u(7)

      ubar21=(A1*uhat11 + A2*(uhat21*COS(u(15))- uhat22*SIN(u(15)) ) )/(A1+A2)
      ubar22=(A1*uhat12 + A2*(uhat21*SIN(u(15))+ uhat22*COS(u(15)) ) )/(A1+A2)
      ubar23= -U(16)/C1 + m2(3)*C2/(C1*(C1+C2)) + uhat13   
      
      strain2(1) = m2(1)/k21 + ubar21
      strain2(2) = m2(2)/k22 + ubar22
      strain2(3) = m2(3)/k23 + ubar23  


      F(1) = L2*(d23(1))
      F(2) = L2*(d23(2))
      F(3) = L2*(d23(3))
      F(4) = L2*(strain2(1)*u(7)/2.0d0 - strain2(2)*u(6)/2.0d0 + strain2(3)*u(5)/2.0d0)
      F(5) = L2*(strain2(1)*u(6)/2.0d0 + strain2(2)*u(7)/2.0d0 - strain2(3)*u(4)/2.0d0)
      F(6) = L2*(-strain2(1)*u(5)/2.0d0 + strain2(2)*u(4)/2.0d0 + strain2(3)*u(7)/2.0d0)
      F(7) = L2*(-strain2(1)*u(4)/2.0d0 - strain2(2)*u(5)/2.0d0 - strain2(3)*u(6)/2.0d0)
	F(8) =L2*(strain2(1)*u(11)/2.0d0 - strain2(2)*u(10)/2.0d0 + strain2(3)*u(9)/2.0d0 -d23qn(1))
	F(9) =L2*(strain2(1)*u(10)/2.0d0 + strain2(2)*u(11)/2.0d0 - strain2(3)*u(8)/2.0d0 -d23qn(2))
	F(10) =L2*(-strain2(1)*u(9)/2.0d0 + strain2(2)*u(8)/2.0d0 + strain2(3)*u(11)/2.0d0-d23qn(3))
	F(11) =L2*(-strain2(1)*u(8)/2.0d0 - strain2(2)*u(9)/2.0d0 - strain2(3)*u(10)/2.0d0-d23qn(4))
      F(12) = L2*0.0D0
      F(13) = L2*0.0D0
      F(14) = L2*0.0D0
      F(15)=  L2*( ( (C1+C2)*u(16)/(C1*C2) ) - (m2(3)/C1) + (uhat23-uhat13) )   ! Alpha
      ater1=( m2(1)*sin(u(15)) - m2(2)*cos(u(15)) )*A2*uhat21/(A1+A2) !Corrected -m2(2) to +ve
      ater2=((uhat21*uhat11 + uhat12*uhat22)*sin(u(15)) + (uhat11*uhat22 - uhat12*uhat21)*cos(u(15)))*A1*A2/(A1+A2)
      F(16) = L2*(ater1+ater2)  
   
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Governing equations corresponding to 1-tube section         
      ubar11=uhat11
      ubar12=uhat12
      ubar13=uhat13
      
      k11=A1
      k12=A1
      k13=C1
     
      m1(1) = (u(23)*u(24) + u(22)*u(25) - u(21)*u(26) - u(20)*u(27))/2.0d0
      m1(2) = (-u(22)*u(24) + u(23)*u(25) + u(20)*u(26) - u(21)*u(27))/2.0d0
      m1(3) = (u(21)*u(24) - u(20)*u(25) + u(23)*u(26) - u(22)*u(27))/2.0d0

!       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
!         the same 4 entries, but permuted and possibly negated.
!         Here are those entries

      entry11 =  2.0D0*(u(20)*u(28) + u(21)*u(29) + u(22)*u(30))
      entry12 =  2.0D0*(u(21)*u(28) - u(20)*u(29) + u(23)*u(30))
      entry13 =  2.0D0*(u(22)*u(28) - u(23)*u(29) - u(20)*u(30))
      entry14 =  2.0D0*(u(23)*u(28) + u(22)*u(29) - u(21)*u(30))


      d13qn(1) =  entry13
      d13qn(2) =  entry14
      d13qn(3) =  entry11
      d13qn(4) =  entry12


!       Components of d3
      d13(1)= 2.0d0*u(20)*u(22)+2.0d0*u(21)*u(23)
      d13(2)= 2.0d0*u(21)*u(22)-2.0d0*u(20)*u(23)
      d13(3)= -u(20)*u(20)-u(21)*u(21)+u(22)*u(22)+u(23)*u(23)

!       Calculate components of the strain
      strain1(1) = m1(1)/k11 + uhat11
      strain1(2) = m1(2)/k12 + uhat12
      strain1(3) = m1(3)/k13 + uhat13
      
 
!       Now calculate the right hand sides
!       Additional terms      
     
      F(17) = L1*(d13(1))
      F(18) = L1*(d13(2))
      F(19) = L1*(d13(3))
      F(20) = L1*(strain1(1)*u(23)/2.0d0 - strain1(2)*u(22)/2.0d0 + strain1(3)*u(21)/2.0d0)
      F(21) = L1*(strain1(1)*u(22)/2.0d0 + strain1(2)*u(23)/2.0d0 - strain1(3)*u(20)/2.0d0)
      F(22) = L1*(-strain1(1)*u(21)/2.0d0 + strain1(2)*u(20)/2.0d0 + strain1(3)*u(23)/2.0d0)
      F(23) = L1*(-strain1(1)*u(20)/2.0d0 - strain1(2)*u(21)/2.0d0 - strain1(3)*u(22)/2.0d0)
F(24)=L1*(strain1(1)*u(27)/2.0d0- strain1(2)*u(26)/2.0d0 + strain1(3)*u(25)/2.0d0 -d13qn(1))
F(25)=L1*(strain1(1)*u(26)/2.0d0+ strain1(2)*u(27)/2.0d0-strain1(3)*u(24)/2.0d0 -d13qn(2))
F(26)=L1*(-strain1(1)*u(25)/2.0d0+ strain1(2)*u(24)/2.0d0+strain1(3)*u(27)/2.0d0 -d13qn(3))
F(27)=L1*(-strain1(1)*u(24)/2.0d0- strain1(2)*u(25)/2.0d0-strain1(3)*u(26)/2.0d0 -d13qn(4) )
      F(28) = L1*0.0D0
      F(29) = L1*0.0D0
      F(30) = L1*0.0D0
  
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Governing equations corresponding to 3-tube section     
      
      m3(1) = ( u(23+14)*u(24+14) + u(22+14)*u(25+14) - u(21+14)*u(26+14) - u(20+14)*u(27+14))/2.0;
      m3(2) = (-u(22+14)*u(24+14) + u(23+14)*u(25+14) + u(20+14)*u(26+14) - u(21+14)*u(27+14))/2.0;
      m3(3) = ( u(21+14)*u(24+14) - u(20+14)*u(25+14) + u(23+14)*u(26+14) - u(22+14)*u(27+14))/2.0;

!%       The 4-vectors D_1[q] n, D_2[q] n, and D_3[q] n contain
!%         the same 4 entries, but permuted and possibly negated.
!%         Here are those entries

      entry31 =  2.0*(u(20+14)*u(28+14) + u(21+14)*u(29+14) + u(22+14)*u(30+14));
      entry32 =  2.0*(u(21+14)*u(28+14) - u(20+14)*u(29+14) + u(23+14)*u(30+14));
      entry33 =  2.0*(u(22+14)*u(28+14) - u(23+14)*u(29+14) - u(20+14)*u(30+14));
      entry34 =  2.0*(u(23+14)*u(28+14) + u(22+14)*u(29+14) - u(21+14)*u(30+14));


!%       Components of D_3[q] n

      d33qn(1) =  entry33;
      d33qn(2) =  entry34;
      d33qn(3) =  entry31;
      d33qn(4) =  entry32;



!%       Components of d3
      d33(1)= 2.0*u(20+14)*u(22+14) + 2.0*u(21+14)*u(23+14);
      d33(2)= 2.0*u(21+14)*u(22+14) - 2.0*u(20+14)*u(23+14);
      d33(3)= -u(20+14)*u(20+14) - u(21+14)*u(21+14) + u(22+14)*u(22+14) + u(23+14)*u(23+14);
      


     ab1=-C2*((U(47) + U(48))/C1 - m3(3)/C1  + U(47)/C2 + uhat23- uhat13 - uhat23) 
     ab2=-C3*((U(47) + U(48))/C1 - m3(3)/C1  + U(48)/C3 + uhat33- uhat13 - uhat33)
    


        ub31=(A1*uhat11 + A2*(uhat21*COS(u(45))- uhat22*SIN(u(45)) ) + A3*(uhat31*COS(u(46))- uhat32*SIN(u(46))))/(A1+A2+A3)
	ub32=(A1*uhat12 + A2*(uhat21*SIN(u(45))+ uhat22*COS(u(45)) ) + A3*(uhat31*SIN(u(46))+ uhat32*COS(u(46))))/(A1+A2+A3)
	ub33= (ab1 + ab2 + C1*uhat13)/(C1+C2+C3)! Corrected addterm to uhat13 ------------------------


	kb31=A1+A2+A3;
	kb32=A1+A2+A3;
	kb33=C1+C2+C3; 
     
	strain3(1) = m3(1)/(A1+A2+A3) + ub31
	strain3(2) = m3(2)/(A1+A2+A3) + ub32
	strain3(3) = m3(3)/(C1+C2+C3) + ub33  

	F(31) = L3*(d33(1));
        F(32) = L3*(d33(2));
	F(33) = L3*(d33(3));
      
	F(34) = L3*( strain3(1)*u(23+14)/2.0 - strain3(2)*u(22+14)/2.0 + strain3(3)*u(21+14)/2.0);
    F(35) = L3*( strain3(1)*u(22+14)/2.0 + strain3(2)*u(23+14)/2.0 - strain3(3)*u(20+14)/2.0);
    F(36) = L3*(-strain3(1)*u(21+14)/2.0 + strain3(2)*u(20+14)/2.0 + strain3(3)*u(23+14)/2.0);
    F(37) = L3*(-strain3(1)*u(20+14)/2.0 - strain3(2)*u(21+14)/2.0 - strain3(3)*u(22+14)/2.0);
     
	F(38) = L3*( strain3(1)*u(27+14)/2.0 - strain3(2)*u(26+14)/2.0 + strain3(3)*u(25+14)/2.0 &
        - d33qn(1));
	F(39) = L3*( strain3(1)*u(26+14)/2.0 + strain3(2)*u(27+14)/2.0 - strain3(3)*u(24+14)/2.0 &
        - d33qn(2));
	F(40) = L3*(-strain3(1)*u(25+14)/2.0 + strain3(2)*u(24+14)/2.0 + strain3(3)*u(27+14)/2.0 &
        - d33qn(3));
	F(41) = L3*(-strain3(1)*u(24+14)/2.0 - strain3(2)*u(25+14)/2.0 - strain3(3)*u(26+14)/2.0 &
        - d33qn(4));
     
	F(42) = L3*0.0;
	F(43) = L3*0.0;
	F(44) = L3*0.0;
      

	F(45)= L3*((U(47) + U(48))/C1 - m3(3)/C1  + U(47)/C2 + uhat23- uhat13)  !!%Alfa2=U(45)
	F(46)= L3*((U(47) + U(48))/C1 - m3(3)/C1  + U(48)/C3 + uhat33- uhat13)  !!%Alfa3=U(46)


	F(47)=L3*(A2*uhat21*(A3*uhat31*SIN(U(45)-U(46)) + A1*uhat11*SIN(U(45))))/(A1+A2+A3) &
       + L3*A2*uhat21*(m3(1)*SIN(U(45)) - m3(2)*COS(U(45)))/(A1+A2+A3)
	F(48)=L3*(A3*uhat31*(A2*uhat21*SIN(U(46)-U(45)) + A1*uhat11*SIN(U(46))))/(A1+A2+A3) & 
       + A3*L3*uhat31*(m3(1)*SIN(U(46)) - m3(2)*COS(U(46)))/(A1+A2+A3)


 	
      RETURN
      END

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION k1,k2,k3,ubar1,ubar2,pi
      pi = 4.d0*ATAN(1.d0)
      PAR(1)=0.0001d0 		
      PAR(2)=0.0d0 		
      PAR(3)=0.0d0
      PAR(4)= 0.0001d0 
      PAR(5)= 0.0d0
	  PAR(6)= 0.0d0
	  PAR(7)= 0.0d0
	  PAR(8)= 0.0d0
	  PAR(9)= 0.0d0
	  PAR(10)=0.0001d0
	  PAR(11)=0.0d0 
      PAR(12)=0.0d0  
      
      U= 0.0d0
      U(3)=PAR(10) + PAR(1)*T   
      U(7)=1.0d0

U(19)=PAR(1) + PAR(10) + PAR(4)*T   
      U(23)=1.0d0
      
      U(33)=PAR(10)*T
      U(37)=1.0d0
!
      RETURN
      END
!                                                                      
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
      DOUBLE PRECISION pi,nx,ny,nz,e1,d111,d112,d113,d121,d122,d123,d131,d132,d133,e2,e3,Ndir1,Ndir2,Ndir3,Del1,Del2,Del3,bas
      DOUBLE PRECISION mb21,mb31,mb22,mb32,mb23,mb33,intb2,intb3,ang1,ang2,ang3
      INTEGER I1,I

      pi = 4.d0*ATAN(1.d0)
      
      nz=PAR(7)*0.0d0
      ny=PAR(7)*SIN(PAR(8))
      nx=PAR(7)*COS(PAR(8))
      ang1=PAR(9)
      ang2=PAR(2)
      ang3=PAR(12)
      
      
   	  e1=0.0d0 	
   	  e2=0.0d0
   	  e3=0.0d0
   
     ! q1,q2,q3,q4, u20,u21,u22,u23   
      d111= u1(20)*u1(20) - u1(21)*u1(21) - u1(22)*u1(22) + u1(23)*u1(23)
      d112= 2.0d0*u1(20)*u1(21) + 2.0d0*u1(22)*u1(23)
      d113= 2.0d0*u1(20)*u1(22) - 2.0d0*u1(21)*u1(23)
   
   	  d121= 2.0d0*u1(20)*u1(21) - 2.0d0*u1(22)*u1(23)
      d122= -u1(20)*u1(20) + u1(21)*u1(21)- u1(22)*u1(22) + u1(23)*u1(23)
      d123= 2.0d0*u1(21)*u1(22) + 2.0d0*u1(20)*u1(23)

      d131= 2.0d0*u1(20)*u1(22) + 2.0d0*u1(21)*u1(23)
      d132= 2.0d0*u1(21)*u1(22) - 2.0d0*u1(20)*u1(23)
      d133= -u1(20)*u1(20) - u1(21)*u1(21) + u1(22)*u1(22) + u1(23)*u1(23)
   
     ! r X N
      Ndir1=nx*d111 + ny*d112 + nz*d113  
      Ndir2=nx*d121 + ny*d122 + nz*d123  
      Ndir3=nx*d131 + ny*d132 + nz*d133  
      
      Del1=e1*d111 + e2*d121 + e3*d131
      Del2=e1*d112 + e2*d122 + e3*d132
      Del3=e1*d113 + e2*d123 + e3*d133 
      
     !Boundary Conditions
           
      FB(1)=U0(1) - U1(31);
      FB(2)=U0(2) - U1(32);
      FB(3)=U0(3) - U1(33);
	  	  
      FB(4)=U0(4) - U1(34);
      FB(5)=U0(5) - U1(35); 
      FB(6)=U0(6) - U1(36); 
      FB(7)=U0(7) - U1(37); ! Dirichlet BCs coupled
      
      FB(8)= U1(8) -U0(24);
      FB(9)= U1(9) -U0(25); 
      FB(10)=U1(10)-U0(26); 
      FB(11)=U1(11)-U0(27);
	       
     
      FB(12)=U1(12) - U0(28); 
      FB(13)=U1(13) - U0(29);
      FB(14)=U1(14) - U0(30); !Neumann Coupled between the 2-tube and 1-tube
      
        
      FB(15)=U0(15) - U1(45) !Dirichlet Condition on the end of tube twist
      FB(16)=U1(16);         !Neumann Condition on the end of tube twist
      
      
      FB(17)=U0(17) - U1(1);
      FB(18)=U0(18) - U1(2);
      FB(19)=U0(19) - U1(3);
      
      FB(20)=U0(20) - U1(4);
      FB(21)=U0(21) - U1(5);
      FB(22)=U0(22) - U1(6);
      FB(23)=U0(23) - U1(7); !Dirichlet BCs coupled at the boundary between 2 tube section and 1 tube
      
      FB(24)= ( U1(16+7)*U1(16+8) + U1(16+6)*U1(16+9) - U1(16+5)*U1(16+10) - U1(16+4)*U1(16+11))/2.0d0 +(e2*Ndir3 - e3*Ndir2)
      FB(25)= (-U1(16+6)*U1(16+8) + U1(16+7)*U1(16+9) + U1(16+4)*U1(16+10) - U1(16+5)*U1(16+11))/2.0d0 +(e3*Ndir1 - e1*Ndir3)!
      FB(26)= ( U1(16+5)*U1(16+8) - U1(16+4)*U1(16+9) + U1(16+7)*U1(16+10) - U1(16+6)*U1(16+11))/2.0d0 +(e1*Ndir2 - e2*Ndir1)
	    
      FB(27)=  U0(24)*U0(20) + U0(25)*U0(21)+ U0(26)*U0(22)+ U0(27)*U0(23) + 2.0d0*(U0(17)*U0(28) + U0(18)*U0(29) + U0(19)*U0(30)) !Neumann BCs on the tip
      
      FB(28)=U1(28) + nx;
      FB(29)=U1(29) + ny;
      FB(30)=U1(30) + nz;!Neumann BCs on the tip
                       
      FB(31)=U0(31); 
      FB(32)=U0(32);
      FB(33)=U0(33);
      
      FB(34)=U0(34);
      FB(35)=U0(35);
      FB(36)=U0(36) - SIN(ang1/2);
      FB(37)=U0(37) - COS(ang1/2);  ! Dirichlet BCs coupled 
   
       
      FB(38)=U1(42) - U0(12);
      FB(39)=U1(43) - U0(13);
      FB(40)=U1(44) - U0(14);
      
      FB(41)=U1(38) - U0(8);
      FB(42)=U1(39) - U0(9); 
      FB(43)=U1(40) - U0(10);
      FB(44)=U1(41) - U0(11);!Neumann Conditions Coupled
	
      FB(45)=U0(45) - ang2; 
      FB(46)=U0(46) - ang3;! Dirichlet BCs coupled 
      FB(47)=U1(47) - U0(16);
      FB(48)=U1(48);
      
      RETURN
      END
!                                                                      

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
