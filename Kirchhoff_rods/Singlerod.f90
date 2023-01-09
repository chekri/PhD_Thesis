
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

!       Right-hand side from Eqs. (4.2.9-12) in Dichmann, Li, and Maddocks,
!       "Hamiltonian Formulations and Symmetries in Rod Mechanics", in the form to be fed to AUTO.

!       PARAMETERS: 
!    	1: Intrinsic curvature of the rod 
!	2: Torsion component of the rod  
!	3: Length of the rod. 
!	4: Tip Load
!	5: Angle of the clamped end 
!	6: Load orientation
!	7: Length of arm    

      !Daroboux vector of the elastic rod
      uhat1=PAR(1)!1.0d0
      uhat2=0.0d0
      uhat3=PAR(2)
      
      ! Specify the bending stiffness A1 and Torsion stiffness C1
      A1=1.0
      C1=1.0/1.3

     
      ! Length of the elastic rod
      L1=PAR(3)
      
     ! Calculate three components of moment in local coordnates   
      m(1) = (u(7)*u(8) + u(6)*u(9) - u(5)*u(10) - u(4)*u(11))/2.0d0
      m(2) = (-u(6)*u(8) + u(7)*u(9) + u(4)*u(10) - u(5)*u(11))/2.0d0
      m(3) = (u(5)*u(8) - u(4)*u(9) + u(7)*u(10) - u(6)*u(11))/2.0d0

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
      strain(3) = m(3)/C1 + uhat3  !

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

 
 	
      RETURN
      END

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION pi
      pi = 4.d0*ATAN(1.d0)
      PAR(1)=0.0d0 		!uhat11
      PAR(2)=0.0d0 		!uhat13
      PAR(3)=0.00001d0          ! Length of the tube should be non-zero 
      PAR(4)=0.0d0 
      PAR(5)=0.0d0 
      PAR(6)=0.0d0 
      PAR(7)=0.00d0
	  
      U= 0.0d0
      
      ! Straight solution as an intial guess
      U(3)=PAR(1)*T
 	
      ! Fixed orientation as an intial guess
      U(7)=1.0d0 

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
      DOUBLE PRECISION pi,nx,ny,nz,e1,e2,e3,d111,d112,d113,d121,d122,d123,d131,d132,d133,Ndir1,Ndir2,Ndir3,Del1,Del2,Del3


      pi = 4.d0*ATAN(1.d0)
      
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
          
      RETURN
      END
!                                                                      

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
