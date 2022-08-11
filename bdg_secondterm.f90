!take the secondary  interaction  paring term into cosideration.
!so the delta terms  should change a little.
!and the explanation of the consistent equation should also change in order to obtain  the  superconductor parameter.
Program Bdg_second_square
      Implicit None
      Integer :: ngrid =5
      Integer :: N
      Double Precision :: miu = 2.0d0, v = 4.0d0, eps = 10.0e-7, T = 10.0e-22
	
      N = 2 * ngrid * ngrid
      Call calculate_pairterm (ngrid, N, miu, v, eps, T)
	
End Program
!
!
Subroutine calculate_pairterm (ng, N, miu, v, eps, T)
      Implicit None
!
      Integer :: N, ng, ndot
      Integer, Parameter :: LWORK = 1000 !LWORK>=MAX(2*N-1)
      Double Precision :: miu, v, eps, T, e
      Integer :: INFO, i, j, k
	  integer::u,d,r,l
	  double precision,parameter::pi=3.1415926
      Complex * 16 :: c
	
      Complex * 16, Dimension (N, N) :: H
	
      Double Precision, Dimension (N) :: W
	
      
	
      Complex * 16, Dimension (LWORK) :: WORK
!
      Double Precision, Dimension (3*N-2) :: RWORK
	
      Complex * 16, Dimension (ng*ng, ng*ng) :: delta0, delta1, delta
!
      ndot = ng * ng
	
	
	
	  !use  the  randomfunction to the pairing  term.
	  !SO we intially give  the  paring  term the  same complex  number.
      delta0 = 0.0d0
	  !x direction
      Do j = 0, ng - 1
         Do i = 1, ng - 1
            Call cplex_creation (c)
            delta0 (i+j*ng, j*ng+i+1)= c
			delta0 (j*ng+i+1,i+j*ng)=c
         End Do
      End Do
	 
     !y direction
      Do i = 1, ng
         Do j = 0, ng - 2
            Call cplex_creation (c)
            delta0 (i+ng*j, i+ng+ng*j) = c
			delta0 (i+ng+ng*j, i+ng*j) = c
         End Do
      End Do
	 !x,y boundery
      Do j = 0, ng - 1
         Call cplex_creation (c)
         delta0 (1+j*ng, ng+j*ng) = c
		 delta0(ng+j*ng,1+j*ng)=c
      End Do
      Do j = 1, ng
         Call cplex_creation (c)
         delta0 (j, ng*(ng-1)+j) = c
		 delta0 (ng*(ng-1)+j,j) = c
      End Do
	
	
	
20    H = 0.0d0
!
	  !reconfined the pairing term(secondary)
	  !x  direction
      Do i = 1, ndot
         Do j = 1, ndot
            H (i, ndot+j) = delta0 (i, j)
         End Do
      End Do
	
!
     !~~~~~~~~~~~~~~~~~~~~~~~~!
     !the none paring_term	
      Do i = 1, ndot
         H (i, i) = - miu
      End Do
	
	 !x  direction
      Do j = 0, ng - 1
         Do i = 1, ng - 1
            H (i+j*ng, j*ng+i+1) = - 1.0d0
         End Do
      End Do
	 !ydirection
      Do i = 1, ng
         Do j = 0, ng - 2
            H (i+ng*j, i+ng+ng*j) = - 1.0d0
         End Do
      End Do
	 !x,y boundery
      Do j = 0, ng - 1
         H (1+j*ng, ng+j*ng) = - 1.0d0
      End Do
      Do j = 1, ng
         H (j, ng*(ng-1)+j) = - 1.0d0
      End Do
	
	!construct the -h  matricx
      Do i = ndot + 1, N
         H (i, i) = miu
      End Do
	 !X-DIRECTION
      Do j = 0, ng - 1
         Do i = 1, ng - 1
            H (ndot+i+j*ng, ndot+j*ng+i+1) = 1.0d0
         End Do
      End Do
	 !ydirection
      Do i = 1, ng
         Do j = 0, ng - 2
            H (i+ng*j+ndot, i+ng+ng*j+ndot) = 1.0d0
         End Do
      End Do
     !x,y boundery
      Do j = 0, ng - 1
         H (ndot+1+j*ng, ndot+ng+j*ng) = 1.0d0
      End Do
      Do j = 1, ng
         H (j+ndot, ndot+ng*(ng-1)+j) = 1.0d0
      End Do
!
	
      Call zheev ('V', 'U', N, H, N, W, WORK, LWORK, RWORK, INFO)
	
      delta1 = 0.0d0
	  !x direction
      Do j = 0, ng - 1
         Do i = 1, ng - 1
            Do k = 1, N
               delta1 (i+j*ng, j*ng+i+1) = delta1 (i+j*ng, j*ng+i+1) + &
              & (H(i+j*ng, k)*conjg(H(j*ng+i+1+ndot, k))+H(j*ng+i+1,k)*conjg(H(i+j*ng+ndot, k))) * Tanh (W(k)/(2.0d0*T))
			  
			  delta1 (j*ng+i+1, i+j*ng) = delta1 (j*ng+i+1, i+j*ng) + &
              & (H(j*ng+i+1, k)*conjg(H(i+j*ng+ndot, k))+H(i+j*ng,k)*conjg(H(j*ng+i+1+ndot, k))) * Tanh (W(k)/(2.0d0*T))
            End Do
         End Do
      End Do
	
     !y direction
      Do i = 1, ng
         Do j = 0, ng - 2
            Do k = 1, N
               delta1 (i+ng*j, i+ng+ng*j) = delta1 (i+ng*j, i+ng+ng*j) &
              & + (H(i+ng*j, k)*conjg(H(i+ng+ng*j+ndot, &
              & k))+H(i+ng+ng*j, k)*conjg(H(i+ng*j+ndot, k))) * Tanh &
              & (W(k)/(2.0d0*T))
			  
			   delta1 (i+ng+ng*j, i+ng*j) = delta1 (i+ng+ng*j, i+ng*j) &
              & + (H(i+ng+ng*j, k)*conjg(H(i+ng*j+ndot, &
              & k))+H(i+ng*j, k)*conjg(H(i+ng+ng*j+ndot, k))) * Tanh &
              & (W(k)/(2.0d0*T))
            End Do
         End Do
      End Do
	
	
	 !x-boundery
      Do j = 0, ng - 1
         Do k = 1, N
            delta1 (1+j*ng, ng+j*ng) = delta1 (1+j*ng, ng+j*ng) + &
           & (H(1+j*ng, k)*conjg(H(ng+j*ng+ndot, k))+H(ng+j*ng, &
           & k)*conjg(H(1+j*ng+ndot, k))) * Tanh (W(k)/(2.0d0*T))
		   
		    delta1 (ng+j*ng,1+j*ng) = delta1 ( ng+j*ng,1+j*ng) + &
           & (H(ng+j*ng, k)*conjg(H(1+j*ng+ndot, k))+H(1+j*ng, &
           & k)*conjg(H(ng+j*ng+ndot, k))) * Tanh (W(k)/(2.0d0*T))
		   
         End Do
      End Do
	  !y-boundery
      Do j = 1, ng
         Do k = 1, N
            delta1 (j, ng*(ng-1)+j) = delta1 (j, ng*(ng-1)+j) + (H(j, &
           & k)*conjg(H(ng*(ng-1)+j+ndot, k))+H(ng*(ng-1)+j, &
           & k)*conjg(H(j+ndot, k))) * Tanh (W(k)/(2.0d0*T))
		   
		   delta1 (ng*(ng-1)+j, j) = delta1 (ng*(ng-1)+j, j) + (H(ng*(ng-1)+j, &
           & k)*conjg(H(j+ndot, k))+H(j, &
           & k)*conjg(H(ng*(ng-1)+j+ndot, k))) * Tanh (W(k)/(2.0d0*T))
         End Do
      End Do
	   delta1=delta1*v/4.0d0
       delta = delta1 - delta0
       delta0= delta1
	
      e = 0.0d0
	
      Do i = 1, ndot
         Do j = 1, ndot
            e = e + Abs (delta(i, j))
         End Do
      End Do
      e = e / (2*ndot)
      Print *, "the error:"
      Print *, e
      If (e > eps) Then
         Go To 20
      Else
         
     Open (10, File='mo_u=2.0_v=4.0.txt')
	 Open (30, File='fujiao_u=2.0_v=4.0.txt')
	 open (40, File='fushu_u=2.0_v=4.0.txt')
	
	
	 !~~~~~~output the value of the paring_term named delta0~~~~~~!
   do i=1,ndot
   call  position_find(i,ng,u,d,r,l)
   write(40,*),delta0(i,u),delta0(i,d),delta0(i,l),delta0(i,r)
   write(10,*),abs(delta0(i,u)),abs(delta0(i,d)),abs(delta0(i,l)),abs(delta0(i,r))
   write(30,*),acos(real(delta0(i,u))/abs(delta0(i,u)))/pi,acos(real(delta0(i,d))/abs(delta0(i,d)))/pi,&
   &acos(real(delta0(i,l))/abs(delta0(i,l)))/pi,acos(real(delta0(i,r))/abs(delta0(i,r)))/pi
   print *,'change another form of output:'
   print *,u,d,l,r
   print *,delta0(i,u),delta0(i,d),delta0(i,l),delta0(i,r)
   end do
	
    End If
    Close (10)
	close (30)
	close (40)
  End Subroutine
!
!
Subroutine cplex_creation (c)
      Implicit None
      Complex * 16 :: c
      Double Precision :: r1, r2
      Call random_seed ()
      Call random_number (r1)
      Call random_number (r2)
	  c = cmplx (r1, r2)
End Subroutine


subroutine position_find(num,ng,u,d,r,l)
implicit none
integer::num,ng,u,d,r,l
   if(num==1) then
      u=num+ng*(ng-1)
	  d=num+ng
      r=num+1
	  l=ng*1
	  else if(num==ng) then
	   u=num+ng*(ng-1)
	   d=num+ng
	   r=1
	   l=num-1
	   
	   else if(num==1+ng*(ng-1)) then
	    u=num-ng
		d=1
		r=num+1
		l=ng*ng
		
		else if(num==ng*ng)then
		u=num-ng
		d=ng
		r=num-ng+1
		l=num-1
		
		else if(num>1.and.num<ng) then
		u=num+ng*(ng-1)
		d=num+ng
		r=num+1
		l=num-1
		
		else if(num>1+ng*(ng-1).and.num<ng*ng)  then
		u=num-ng
		d=num-ng*(ng-1)
		r=num+1
		l=num-1
		
		else if(mod(num,ng)==1) then
		u=num-ng
		d=num+ng
		r=num+1
		l=num+ng-1
		
		else if(mod(num,ng)==0)then
		u=num-ng
		d=num+ng
		r=num-ng+1
		l=num-1
		
		else  
		u=num-ng
		d=num+ng
		r=num+1 
		l=num-1
		end if
		
end subroutine 

