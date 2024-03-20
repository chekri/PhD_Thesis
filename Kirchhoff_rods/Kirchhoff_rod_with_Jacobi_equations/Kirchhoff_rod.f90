
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION epsilon,A1,C1,uhat1,uhat2,uhat3,pi,L1
      DOUBLE PRECISION entry1,entry2,entry3,entry4
      DOUBLE PRECISION m(3),d3qn(4),strain(3),d3(3)
      DOUBLE PRECISION w201,w202,w203,w204,w205,w206,w207,w208,w209,w210,w211,w212,w213,w214,w215,w216,w217,w218
      DOUBLE PRECISION w221,w222,w223,w224,w225,w226,w227,w228,w231,w232,w233,w234,w235,w236,w237,w238
      DOUBLE PRECISION w241,w242,w243,w244,w245,w246,w247,w248,w251,w252,w253,w254,w255,w256,w257,w258
      DOUBLE PRECISION w261,w262,w263,w264,w265,w266,w267,w268,w271,w272,w273,w274,w275,w276,w277,w278
      INTEGER I,I1
!       Right-hand side from Eqs. (4.2.9-12) in Dichmann, Li, and Maddocks,
!       "Hamiltonian Formulations and Symmetries in Rod Mechanics", in the form to be fed to AUTO.

!       PARAMETERS: 
!    	1: Intrinsic curvature of the rod 
!	2: Torsion component of the rod  
!	3: Length of the rod. 
!	4: Tip Load
!	5: Angle of the clamped end 
!	6: Load orientation
!	7: Basis (To provide initial guess for solutions, we perform continuation from 0 vector to Identity (e) vector at s=l for Jacobi IVPs)

      !Daroboux vector of the elastic rod
      uhat1=PAR(1)
      uhat2=0.0d0
      uhat3=PAR(2)
      
      ! Specify the bending stiffness A1 and Torsion stiffness C1
      A1=1.0
      C1=1.0/1.3
     
      ! Length of the elastic rod
      L1=PAR(3)
      
     ! Calculate three components of moment in local coordinates   
      m(1) = ( u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0d0
      m(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0d0
      m(3) = ( u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0d0

!       The 4 entries in D_3[q] n 
      entry1 =  2.0D0*(u(4)*u(12) + u(5)*u(13) + U(6)*u(14))
      entry2 =  2.0D0*(u(5)*u(12) - u(4)*u(13) + U(7)*u(14))
      entry3 =  2.0D0*(u(6)*u(12) - u(7)*u(13) - U(4)*u(14))
      entry4 =  2.0D0*(u(7)*u(12) + u(6)*u(13) - U(5)*u(14))


!       Components of D_3[q] n

      d3qn(1) =  entry3
      d3qn(2) =  entry4
      d3qn(3) =  entry1
      d3qn(4) =  entry2

!       Components of d3

      d3(1)= 2.0d0*u(4)*u(6)+2.0d0*u(5)*u(7)
      d3(2)= 2.0d0*u(5)*u(6)-2.0d0*u(4)*u(7)
      d3(3)= -u(4)*u(4)-u(5)*u(5)+u(6)*u(6)+u(7)*u(7)

!     Strain components
      strain(1) = m(1)/A1 + uhat1
      strain(2) = m(2)/A1 + uhat2
      strain(3) = m(3)/C1 + uhat3  

!     Hamiltonian from of the rod equilibria      

      F(1) = L1*(d3(1))
      F(2) = L1*(d3(2))
      F(3) = L1*(d3(3))
      F(4) = L1*(strain(1)*u(7)/2.0d0 - strain(2)*u(6)/2.0d0 + strain(3)*u(5)/2.0d0)
      F(5) = L1*(strain(1)*u(6)/2.0d0 + strain(2)*u(7)/2.0d0 - strain(3)*u(4)/2.0d0)
      F(6) = L1*(-strain(1)*u(5)/2.0d0 + strain(2)*u(4)/2.0d0 + strain(3)*u(7)/2.0d0)
      F(7) = L1*(-strain(1)*u(4)/2.0d0 - strain(2)*u(5)/2.0d0 - strain(3)*u(6)/2.0d0)
      F(8) = L1*(strain(1)*u(11)/2.0d0 - strain(2)*u(10)/2.0d0 + strain(3)*u(9)/2.0d0 -d3qn(1))
      F(9) = L1*(strain(1)*u(10)/2.0d0 + strain(2)*u(11)/2.0d0 - strain(3)*u(8)/2.0d0 -d3qn(2)) 
      F(10) =L1*(-strain(1)*u(9)/2.0d0 + strain(2)*u(8)/2.0d0 + strain(3)*u(11)/2.0d0-d3qn(3))
      F(11) =L1*(-strain(1)*u(8)/2.0d0 - strain(2)*u(9)/2.0d0 - strain(3)*u(10)/2.0d0-d3qn(4))
      F(12) =L1*0.0D0
      F(13) =L1*0.0D0
      F(14) =L1*0.0D0

   ! Jacobi equations


 	w201= -U(7)*U(11)/(4.0d0*A1) - U(6)*U(10)/(4.0d0*A1)- U(5)*U(9)/(4.0d0*C1)  + 0.0d0                  
	w202= -U(7)*U(10)/(4.0d0*A1) + U(6)*U(11)/(4.0d0*A1)+ U(5)*U(8)/(4.0d0*C1)  + (strain(3)/2.0d0)
	w203=  U(7)*U(9)/(4.0d0*A1)  + U(6)*U(8)/(4.0d0*A1) - U(5)*U(11)/(4.0d0*C1) - (strain(2)/2.0d0)
	w204=  U(7)*U(8)/(4.0d0*A1)  - U(6)*U(9)/(4.0d0*A1) + U(5)*U(10)/(4.0d0*C1) + (strain(1)/2.0d0)

	w205=  U(7)*U(7)/(4.0d0*A1) + U(6)*U(6)/(4.0d0*A1) + U(5)*U(5)/(4.0d0*C1)
	w206=  U(6)*U(7)/(4.0d0*A1) - U(7)*U(6)/(4.0d0*A1) - U(4)*U(5)/(4.0d0*C1)
	w207= -U(5)*U(7)/(4.0d0*A1) - U(4)*U(6)/(4.0d0*A1) + U(7)*U(5)/(4.0d0*C1)
	w208= -U(4)*U(7)/(4.0d0*A1) + U(5)*U(6)/(4.0d0*A1) - U(6)*U(5)/(4.0d0*C1)

	w211= -U(6)*U(11)/(4.0d0*A1)+ U(7)*U(10)/(4.0d0*A1)+ U(4)*U(9)/(4.0d0*C1) - strain(3)/2.0d0
	w212= -U(6)*U(10)/(4.0d0*A1)- U(7)*U(11)/(4.0d0*A1)- U(4)*U(8)/(4.0d0*C1) + 0.0d0
	w213=  U(6)*U(9)/(4.0d0*A1) - U(7)*U(8)/(4.0d0*A1) + U(4)*U(11)/(4.0d0*C1) + strain(1)/2.0d0
	w214=  U(6)*U(8)/(4.0d0*A1) + U(7)*U(9)/(4.0d0*A1) - U(4)*U(10)/(4.0d0*C1) + strain(2)/2.0d0
	
	w215=  U(7)*U(6)/(4.0d0*A1) - U(6)*U(7)/(4.0d0*A1) - U(5)*U(4)/(4.0d0*C1)
	w216=  U(6)*U(6)/(4.0d0*A1) + U(7)*U(7)/(4.0d0*A1) + U(4)*U(4)/(4.0d0*C1)
	w217= -U(5)*U(6)/(4.0d0*A1) + U(4)*U(7)/(4.0d0*A1) - U(7)*U(4)/(4.0d0*C1)
	w218= -U(4)*U(6)/(4.0d0*A1) - U(5)*U(7)/(4.0d0*A1) + U(6)*U(4)/(4.0d0*C1)

	w221= +U(5)*U(11)/(4.0d0*A1)+ U(4)*U(10)/(4.0d0*A1) - U(7)*U(9)/(4.0d0*C1)  + strain(2)/2.0d0
	w222= +U(5)*U(10)/(4.0d0*A1)- U(4)*U(11)/(4.0d0*A1) + U(7)*U(8)/(4.0d0*C1)  - strain(1)/2.0d0
	w223= -U(5)*U(9)/(4.0d0*A1) - U(4)*U(8)/(4.0d0*A1)  - U(7)*U(11)/(4.0d0*C1) + 0.0d0
	w224= -U(5)*U(8)/(4.0d0*A1) + U(4)*U(9)/(4.0d0*A1)  + U(7)*U(10)/(4.0d0*C1) + strain(3)/2.0d0
	
	w225= -U(7)*U(5)/(4.0d0*A1) - U(6)*U(4)/(4.0d0*A1) + U(5)*U(7)/(4.0d0*C1)
	w226= -U(6)*U(5)/(4.0d0*A1) + U(7)*U(4)/(4.0d0*A1) - U(4)*U(7)/(4.0d0*C1)
	w227= U(5)*U(5)/(4.0d0*A1) + U(4)*U(4)/(4.0d0*A1)  + U(7)*U(7)/(4.0d0*C1)
	w228= U(4)*U(5)/(4.0d0*A1) - U(5)*U(4)/(4.0d0*A1)  - U(6)*U(7)/(4.0d0*C1)
		
	w231= +U(4)*U(11)/(4.0d0*A1)- U(5)*U(10)/(4.0d0*A1)+ U(6)*U(9)/(4.0d0*C1)  - strain(1)/2.0d0
	w232= +U(4)*U(10)/(4.0d0*A1)+ U(5)*U(11)/(4.0d0*A1)- U(6)*U(8)/(4.0d0*C1)  - strain(2)/2.0d0	
	w233= -U(4)*U(9)/(4.0d0*A1) + U(5)*U(8)/(4.0d0*A1) + U(6)*U(11)/(4.0d0*C1) - strain(3)/2.0d0
	w234= -U(4)*U(8)/(4.0d0*A1) - U(5)*U(9)/(4.0d0*A1) - U(6)*U(10)/(4.0d0*C1) + 0.0d0
	
	w235= -U(7)*U(4)/(4.0d0*A1) + U(6)*U(5)/(4.0d0*A1) - U(5)*U(6)/(4.0d0*C1)
	w236= -U(6)*U(4)/(4.0d0*A1) - U(7)*U(5)/(4.0d0*A1) + U(4)*U(6)/(4.0d0*C1)
	w237=  U(5)*U(4)/(4.0d0*A1) - U(4)*U(5)/(4.0d0*A1) - U(7)*U(6)/(4.0d0*C1)
	w238=  U(4)*U(4)/(4.0d0*A1) + U(5)*U(5)/(4.0d0*A1) + U(6)*U(6)/(4.0d0*C1)
     

	w241= -U(11)*U(11)/(4.0d0*A1) - U(10)*U(10)/(4.0d0*A1) - U(9)*U(9)/(4.0d0*C1) + 2.0d0*U(14)
	w242= -U(11)*U(10)/(4.0d0*A1) + U(10)*U(11)/(4.0d0*A1) + U(9)*U(8)/(4.0d0*C1)
	w243=  U(11)*U(9)/(4.0d0*A1)  + U(10)*U(8)/(4.0d0*A1)  - U(9)*U(11)/(4.0d0*C1) - 2.0d0*U(12)
	w244=  U(11)*U(8)/(4.0d0*A1)  - U(10)*U(9)/(4.0d0*A1)  + U(9)*U(10)/(4.0d0*C1) + 2.0d0*U(13)
	
	w245=  U(7)*U(11)/(4.0d0*A1) + U(6)*U(10)/(4.0d0*A1) + U(5)*U(9)/(4.0d0*C1) + 0.0d0
	w246=  U(6)*U(11)/(4.0d0*A1) - U(7)*U(10)/(4.0d0*A1) - U(4)*U(9)/(4.0d0*C1) + strain(3)/2.0d0
	w247= -U(5)*U(11)/(4.0d0*A1) - U(4)*U(10)/(4.0d0*A1) + U(7)*U(9)/(4.0d0*C1) - strain(2)/2.0d0
	w248= -U(4)*U(11)/(4.0d0*A1) + U(5)*U(10)/(4.0d0*A1) - U(6)*U(9)/(4.0d0*C1) + strain(1)/2.0d0
	

	
	w251= -U(10)*U(11)/(4.0d0*A1) + U(11)*U(10)/(4.0d0*A1) + U(8)*U(9)/(4.0d0*C1)
	w252= -U(10)*U(10)/(4.0d0*A1) - U(11)*U(11)/(4.0d0*A1) - U(8)*U(8)/(4.0d0*C1)  + 2.0D0*U(14) 
	w253=  U(10)*U(9)/(4.0d0*A1)  - U(11)*U(8)/(4.0d0*A1)  + U(8)*U(11)/(4.0d0*C1) - 2.0D0*U(13)  	
	w254=  U(10)*U(8)/(4.0d0*A1)  + U(11)*U(9)/(4.0d0*A1)  - U(8)*U(10)/(4.0d0*C1) - 2.0D0*U(12) 

	w255=  U(7)*U(10)/(4.0d0*A1) - U(6)*U(11)/(4.0d0*A1) - U(5)*U(8)/(4.0d0*C1) - strain(3)/2.0d0
	w256=  U(6)*U(10)/(4.0d0*A1) + U(7)*U(11)/(4.0d0*A1) + U(4)*U(8)/(4.0d0*C1) + 0.0d0
	w257= -U(5)*U(10)/(4.0d0*A1) + U(4)*U(11)/(4.0d0*A1) - U(7)*U(8)/(4.0d0*C1) + strain(1)/2.0d0
	w258= -U(4)*U(10)/(4.0d0*A1) - U(5)*U(11)/(4.0d0*A1) + U(6)*U(8)/(4.0d0*C1) + strain(2)/2.0d0
	

	w261= +U(9)*U(11)/(4.0d0*A1) + U(8)*U(10)/(4.0d0*A1) - U(11)*U(9)/(4.0d0*C1) - 2.0D0*U(12) 
	w262= +U(9)*U(10)/(4.0d0*A1) - U(8)*U(11)/(4.0d0*A1) + U(11)*U(8)/(4.0d0*C1) - 2.0D0*U(13) 
	w263= -U(9)*U(9)/(4.0d0*A1)  - U(8)*U(8)/(4.0d0*A1)  - U(11)*U(11)/(4.0d0*C1)- 2.0D0*U(14)  
	w264= -U(9)*U(8)/(4.0d0*A1)  + U(8)*U(9)/(4.0d0*A1)  + U(11)*U(10)/(4.0d0*C1) 
	
	w265= -U(7)*U(9)/(4.0d0*A1) - U(6)*U(8)/(4.0d0*A1) + U(5)*U(11)/(4.0d0*C1) + strain(2)/2.0d0
	w266= -U(6)*U(9)/(4.0d0*A1) + U(7)*U(8)/(4.0d0*A1) - U(4)*U(11)/(4.0d0*C1) - strain(1)/2.0d0
	w267=  U(5)*U(9)/(4.0d0*A1) + U(4)*U(8)/(4.0d0*A1) + U(7)*U(11)/(4.0d0*C1) + 0.0d0
	w268=  U(4)*U(9)/(4.0d0*A1) - U(5)*U(8)/(4.0d0*A1) - U(6)*U(11)/(4.0d0*C1) + strain(3)/2.0d0
	
	w271=  U(8)*U(11)/(4.0d0*A1) - U(9)*U(10)/(4.0d0*A1) + U(10)*U(9)/(4.0d0*C1) + 2.0D0*U(13) 
	w272=  U(8)*U(10)/(4.0d0*A1) + U(9)*U(11)/(4.0d0*A1) - U(10)*U(8)/(4.0d0*C1) - 2.0D0*U(12) 
	w273= -U(8)*U(9)/(4.0d0*A1)  + U(9)*U(8)/(4.0d0*A1)  + U(10)*U(11)/(4.0d0*C1)
	w274= -U(8)*U(8)/(4.0d0*A1)  - U(9)*U(9)/(4.0d0*A1)  - U(10)*U(10)/(4.0d0*C1)- 2.0D0*U(14)  	
	
	w275= -U(7)*U(8)/(4.0d0*A1) + U(6)*U(9)/(4.0d0*A1) - U(5)*U(10)/(4.0d0*C1) - strain(1)/2.0d0
	w276= -U(6)*U(8)/(4.0d0*A1) - U(7)*U(9)/(4.0d0*A1) + U(4)*U(10)/(4.0d0*C1) - strain(2)/2.0d0
	w277=  U(5)*U(8)/(4.0d0*A1) - U(4)*U(9)/(4.0d0*A1) - U(7)*U(10)/(4.0d0*C1) - strain(3)/2.0d0
	w278=  U(4)*U(8)/(4.0d0*A1) + U(5)*U(9)/(4.0d0*A1) + U(6)*U(10)/(4.0d0*C1) - 0.0d0


	DO I=1,6
  	I1=14 + 14*(I-1)
	F(I1+1)=L1*(  2.0d0*U(6)*U(I1+4) + 2.0d0*u(7)*U(I1+5)  + 2.0d0*U(4)*U(I1+6) + 2.0d0*U(5)*U(I1+7))
 	F(I1+2)=L1*( -2.0d0*u(7)*U(I1+4) + 2.0d0*U(6)*U(I1+5)  + 2.0d0*U(5)*U(I1+6) - 2.0d0*U(4)*U(I1+7))	
        F(I1+3)=L1*( -2.0d0*U(4)*U(I1+4) - 2.0d0*U(5)*U(I1+5)  + 2.0d0*U(6)*U(I1+6) + 2.0d0*u(7)*U(I1+7))
 

        F(I1+4)=L1*(w201*U(I1+4) + w202*U(I1+5) + w203*U(I1+6) + w204*U(I1+7) + w205*U(I1+8) + w206*U(I1+9) + w207*U(I1+10) &
            + w208*U(I1+11))
        F(I1+5)=L1*(w211*U(I1+4) + w212*U(I1+5) + w213*U(I1+6) + w214*U(I1+7) + w215*U(I1+8) + w216*U(I1+9) + w217*U(I1+10) &
            + w218*U(I1+11))
        F(I1+6)=L1*(w221*U(I1+4) + w222*U(I1+5) + w223*U(I1+6) + w224*U(I1+7) + w225*U(I1+8) + w226*U(I1+9) + w227*U(I1+10) &
            + w228*U(I1+11)) 
        F(I1+7)=L1*(w231*U(I1+4) + w232*U(I1+5) + w233*U(I1+6) + w234*U(I1+7) + w235*U(I1+8) + w236*U(I1+9) + w237*U(I1+10) &
            + w238*U(I1+11)) 

	    F(I1+8)= L1*(w241*U(I1+4) + w242*U(I1+5) + w243*U(I1+6) + w244*U(I1+7) + w245*U(I1+8) + w246*U(I1+9) + w247*U(I1+10)  &
	        + w248*U(I1+11) + 2.0d0*( -U(6)*U(I1+12) + U(7)*U(I1+13) + U(4)*U(I1+14)))		     
        F(I1+9)= L1*(w251*U(I1+4) + w252*U(I1+5) + w253*U(I1+6) + w254*U(I1+7) + w255*U(I1+8) + w256*U(I1+9) + w257*U(I1+10)  &
            + w258*U(I1+11) + 2.0d0*( -U(7)*U(I1+12) - U(6)*U(I1+13) + U(5)*U(I1+14))) 
        F(I1+10)= L1*(w261*U(I1+4) + w262*U(I1+5) + w263*U(I1+6) + w264*U(I1+7) + w265*U(I1+8) + w266*U(I1+9) + w267*U(I1+10) &
            + w268*U(I1+11) + 2.0d0*( -U(4)*U(I1+12) - U(5)*U(I1+13) - U(6)*U(I1+14)))  
        F(I1+11)= L1*(w271*U(I1+4) + w272*U(I1+5) + w273*U(I1+6) + w274*U(I1+7) + w275*U(I1+8) + w276*U(I1+9) + w277*U(I1+10) &
            + w278*U(I1+11) + 2.0d0*( -U(5)*U(I1+12) + U(4)*U(I1+13) - U(7)*U(I1+14))) 
  
        F(I1+12)=L1*0.0d0
        F(I1+13)=L1*0.0d0
        F(I1+14)=L1*0.0d0	 
    
    ENDDO
 	
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
      PAR(1)=0.0d0 		!uhat11
      PAR(2)=0.0d0 		!uhat13
      PAR(3)=0.00001d0 
      PAR(4)=0.0d0 
      PAR(5)=0.0d0 
      PAR(6)=0.0d0 
      PAR(7)=0.0d0  
      U= 0.0d0
      
	 
      U(3)=PAR(1)*T! SIN(ubar1*T)/ubar1
      U(7)=1.0d0 !COS(ubar1*T/2.0d0)

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
      DOUBLE PRECISION pi,nx,ny,nz,e1,d111,d112,d113,d121,d122,d123,d131,d132,d133,e2,e3,Ndir1,Ndir2,Ndir3,bas
      DOUBLE PRECISION DNdir1,DNdir2,DNdir3,Del1,Del2,Del3
      INTEGER I1,I

      pi = 4.d0*ATAN(1.d0)
      bas=PAR(7)
      
      
      nz=PAR(4)*0.0d0
      nx=PAR(4)*SIN(PAR(6))
      ny=PAR(4)*COS(PAR(6))
     
      ! Arm of the load [e1,e2,e3]
       e1=0.0d0
       e2=PAR(7)
       e3=0.00d0
   	  
   
    ! Director frame at the tip 
     ! Components of d1
      d111= u1(4)*u1(4) - u1(5)*u1(5) - u1(6)*u1(6) + u1(7)*u1(7)
      d112= 2.0d0*u1(4)*u1(5) + 2.0d0*u1(6)*u1(7)
      d113= 2.0d0*u1(4)*u1(6) - 2.0d0*u1(5)*u1(7)
   
      ! Components of d2
      d121= 2.0d0*u1(4)*u1(5) - 2.0d0*u1(6)*u1(7)
      d122= -u1(4)*u1(4) + u1(5)*u1(5)- u1(6)*u1(6) + u1(7)*u1(7)
      d123= 2.0d0*u1(5)*u1(6) + 2.0d0*u1(4)*u1(7)

      ! Components of d3
      d131= 2.0d0*u1(4)*u1(6) + 2.0d0*u1(5)*u1(7)
      d132= 2.0d0*u1(5)*u1(6) - 2.0d0*u1(4)*u1(7)
      d133= -u1(4)*u1(4) - u1(5)*u1(5) + u1(6)*u1(6) + u1(7)*u1(7)
   
     ! r X N
      Ndir1=nx*d111 + ny*d112 + nz*d113  
      Ndir2=nx*d121 + ny*d122 + nz*d123  
      Ndir3=nx*d131 + ny*d132 + nz*d133  
      
      Del1=e1*d111 + e2*d121 + e3*d131
      Del2=e1*d112 + e2*d122 + e3*d132
      Del3=e1*d113 + e2*d123 + e3*d133 
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Dirichlet Boundary conditions for position vector r and quaternion q at clamped end s=0
      FB(1)=U0(1)
      FB(2)=U0(2)
      FB(3)=U0(3)
	  	  
      FB(4)=U0(4) 
      FB(5)=U0(5) 
      FB(6)=U0(6) - SIN(-PAR(5)/2)
      FB(7)=U0(7) - COS(-PAR(5)/2)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      ! Natural boundary conditions at the free end s=L
      ! Moment at the tip of the rod   


      FB(8)= ( U1(7)*U1(8) + U1(6)*U1(9) - U1(5)*U1(10) - U1(4)*U1(11))/2.0d0 +(e2*Ndir3 - e3*Ndir2)
      FB(9)= (-U1(6)*U1(8) + U1(7)*U1(9) + U1(4)*U1(10) - U1(5)*U1(11))/2.0d0 +(e3*Ndir1 - e1*Ndir3)!
      FB(10)=( U1(5)*U1(8) - U1(4)*U1(9) + U1(7)*U1(10) -U1(6)*U1(11))/2.0d0  +(e1*Ndir2 - e2*Ndir1)
            
      !Integral mu.q + 2 r.n =0  (Dichmann, Li, and Maddocks)
      FB(11)=  U1(8)*U1(4) + U1(9)*U1(5)+ U1(10)*U1(6)+ U1(11)*U1(7) + 2.0d0*(U1(1)*U1(12) + U1(2)*U1(13) + U1(3)*U1(14)) 
            
     ! Force at the tip of the rod
      FB(12)=U1(12) + nx
      FB(13)=U1(13) + ny
      FB(14)=U1(14) + nz
        
        DO I=1,6
	    I1=14 + 14*(I-1)
      	
   
        ! -------------------------------->
	FB(I1+1)= U1(I1+1) !r ! dx=0
        FB(I1+2)= U1(I1+2)    ! dy=0
        FB(I1+3)= U1(I1+3)    ! dz=0
        
        FB(I1+4)= U1(I1+4) !q ! dq1=0
        FB(I1+5)= U1(I1+5)    ! dq2=0
        FB(I1+6)= U1(I1+6)    ! dq3=0
        FB(I1+7)= U1(I1+7)    ! dq4=0
        
        DNdir1 = 2.0d0*nx*(U1(4)*U1(I1+4) - U1(5)*U1(I1+5) - U1(6)*U1(I1+6) + U1(7)*U1(I1+7)) &
                +2.0d0*ny*(U1(4)*U1(I1+5) + U1(I1+4)*U1(5) + U1(6)*U1(I1+7) + U1(I1+6)*U1(7)) &
                +2.0d0*nz*(U1(4)*U1(I1+6) + U1(I1+4)*U1(6) - U1(5)*U1(I1+7) - U1(I1+5)*U1(7))
        
        DNdir2 = 2.0d0*nx*(U1(4)*U1(I1+5) + U1(I1+4)*U1(5) - U1(6)*U1(I1+7) - U1(I1+6)*U1(7)) &
                +2.0d0*ny*(-U1(4)*U1(I1+4) + U1(5)*U1(I1+5) - U1(6)*U1(I1+6) + U1(7)*U1(I1+7)) &
                +2.0d0*nz*(U1(5)*U1(I1+6) + U1(I1+5)*U1(6) + U1(4)*U1(I1+7) + U1(I1+4)*U1(7))
          
               
        DNdir3 = 2.0d0*nx*(U1(4)*U1(I1+6) + U1(I1+4)*U1(6) + U1(5)*U1(I1+7) + U1(I1+5)*U1(7)) &
                +2.0d0*ny*(U1(5)*U1(I1+6) + U1(I1+5)*U1(6) - U1(4)*U1(I1+7) - U1(I1+4)*U1(7)) &
                +2.0d0*nz*(-U1(4)*U1(I1+4) - U1(5)*U1(I1+5) + U1(6)*U1(I1+6) + U1(7)*U1(I1+7)) 
        

	FB(I1+8)= (U1(7)*U1(I1+8) + U1(6)*U1(I1+9) - U1(5)*U1(I1+10)  - U1(4)*U1(I1+11)&
	         + U1(I1+7)*U1(8) + U1(I1+6)*U1(9) - U1(I1+5)*U1(10) - U1(I1+4)*U1(11))/2.0d0 + (e2*DNdir3 - e3*DNdir2)
	    
	FB(I1+9)=(-U1(6)*U1(I1+8) + U1(7)*U1(I1+9) + U1(4)*U1(I1+10) - U1(5)*U1(I1+11) &
	         - U1(I1+6)*U1(8) + U1(I1+7)*U1(9) + U1(I1+4)*U1(10) - U1(I1+5)*U1(11) )/2.0d0 + (e3*DNdir1 - e1*DNdir3)!
      	
	FB(I1+10)=(U1(5)*U1(I1+8)- U1(4)*U1(I1+9) + U1(7)*U1(I1+10) - U1(6)*U1(I1+11)  &
             + U1(I1+5)*U1(8) - U1(I1+4)*U1(9) + U1(I1+7)*U1(10) - U1(I1+6)*U1(11) )/2.0d0 + (e1*DNdir2 - e2*DNdir1)
	    


        
      	FB(I1+11)=U1(4)*U1(I1+8) + U1(5)*U1(I1+9) + U1(6)*U1(I1+10)+ U1(7)*U1(I1+11)                   &
      	          +U1(I1+4)*U1(8) + U1(I1+5)*U1(9) + U1(I1+6)*U1(10) + U1(I1+7)*U1(11)                 &
      	+ 2.0d0*(U1(1)*U1(I1+12) + U1(2)*U1(I1+13) + U1(3)*U1(I1+14) + U1(I1+1)*U1(12) + U1(I1+2)*U1(13) + U1(I1+3)*U1(14))  ! d(mu.q + 2r.n)=0
      
        
        FB(I1+12)= U1(I1+12) !x dn1=0
        FB(I1+13)= U1(I1+13) !x dn2=0
        FB(I1+14)= U1(I1+14) !x dn3=0
        ENDDO
      
      I1=14
      FB(I1+1)=U1(I1+1) - bas
      I1=28
      FB(I1+2)=U1(I1+2) - bas
      I1=42
      FB(I1+3)=U1(I1+3) - bas
      
      I1=56
      FB(I1+4)=U1(I1+4) - bas*U1(7)
      FB(I1+5)=U1(I1+5) - bas*U1(6)
      FB(I1+6)=U1(I1+6) + bas*U1(5)
      FB(I1+7)=U1(I1+7) + bas*U1(4)
      
      I1=70
      FB(I1+4)=U1(I1+4) + bas*U1(6)
      FB(I1+5)=U1(I1+5) - bas*U1(7)
      FB(I1+6)=U1(I1+6) - bas*U1(4)
      FB(I1+7)=U1(I1+7) + bas*U1(5)
      
      I1=84
      FB(I1+4)=U1(I1+4) - bas*U1(5)
      FB(I1+5)=U1(I1+5) + bas*U1(4)
      FB(I1+6)=U1(I1+6) - bas*U1(7)
      FB(I1+7)=U1(I1+7) + bas*U1(6)
      
      
      RETURN
      END
!                                                                      

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
