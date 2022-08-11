! to solve the BDG eqaution
! in a square lattice,N X N grids
! the first nearest neighbors
! the on site singltet pairing term
! parameter: t=1.0,u=2.0,v=1.0


Program Bdg_square
      Implicit None
      Integer  :: ngrid = 10
	  integer::N
      double precision :: u = 2.0d0, v =1.0d0,eps=10.0e-15,T=10.0e-22
	  
	  N=2*ngrid*ngrid
      call calculate_pairterm (ngrid,N, u, v,eps,T)
	  
End Program
!
!

Subroutine calculate_pairterm (ng,N, u, v,eps,T)
      Implicit None
   
      Integer :: N,ng,ndot
      Integer,parameter :: LWORK = 1000   !LWORK>=MAX(2*N-1)
      double precision :: u, v,r1,r2,eps,T,e,z0
      Integer :: INFO, i, j
	
      Complex * 16, Dimension (N, N) :: H
	
      Double Precision, Dimension (N) :: W
	  
	  Double Precision, Dimension (ng) :: x,y
	
      Complex * 16, Dimension (LWORK) :: WORK

      Double Precision, Dimension (3*N-2) :: RWORK
	  
	  complex  *16, Dimension (ng*ng) :: delta0,delta1,delta
       
	  ndot=ng*ng
	  
	  !in  order to write in the  form of  a  vector,so  i  add the  value of  the z0  axis.
	  z0=0.0d0
	  
	  !use  the  randomfunction to the pairing  term.
	  !SO we intially give  the  paring  term the  same complex  number.
	  call  random_seed()
      call  random_number(r1)
	  call  random_number(r2)
      delta0=cmplx(r1,r2)
	 
	 
 20   H=0.0d0
 
	  do i=1,ndot
	  H(i,ndot+i)=delta0(i)
	  end do
	  
      Do i = 1, ndot
         H (i, i) = -u
      End Do
	
	 !x  direction
	  do j=0,ng-1
      do i=1,ng-1
         H (i+j*ng, j*ng+i+1) = -1.0d0
      end do
	  end do
	  
     !ydirection 
	  do i=1,ng
	  do j=0,ng-2
	   H(i+ng*j,i+ng+ng*j)=-1.0d0
	  end do
   	  end do
	 
	 
	 !x,y boundery
	  do j=0,ng-1
	    H(1+j*ng,ng+j*ng)=-1.0d0
	  end do
	  do j=1,ng
	    H(j,ng*(ng-1)+j)=-1.0d0
	  end do
	
	!construct the -h  matricx
      Do i = ndot + 1, N
         H (i, i)= u
      End Do
	  
	!X-DIRECTION
      do j=0,ng-1
      do i=1,ng-1
         H (ndot+i+j*ng, ndot+j*ng+i+1) = 1.0d0
      end Do
	  end do
	  
    !ydirection 
	  do i=1,ng
	  do j=0,ng-2
	    H(i+ng*j+ndot,i+ng+ng*j+ndot)=1.0d0
	  end do
   	  end do
	 
	 
	!x,y boundery
	  do j=0,ng-1
	   H(ndot+1+j*ng,ndot+ng+j*ng)=1.0d0
	  end do
	  do j=1,ng
	   H(j+ndot,ndot+ng*(ng-1)+j)=1.0d0
	  end do 
    
	 
	 Call zheev ('V', 'U', N, H, N, W, WORK, LWORK, RWORK, INFO)
	
	 delta1=0.0d0
	 do  j=1,ndot
	 do  i=1,N
	 delta1(j)=delta1(j)+H(j,i)*conjg(H(j+ndot,i))*tanh(W(i)/(2.0d0*T))
	 end do
	 delta1(j)=delta1(j)*v/2.0d0
	 end do
	 
	 delta=delta1-delta0
	 delta0=delta1
      
	 
	 e=0.0d0
	
	 do i=1,ndot
	   e=e+abs(delta(i))
	 end do
	 print *,"the error:"
	 print *,e/ndot
	 if (e/ndot>eps)  then
	 goto 20 
     else 
     open(10,file='pairing_u=2.0_v=1.0.txt') 
	!open(40,file='delta_xyz.txt')
	 do i=1,ndot
	! write(40,200),real(delta0(i)),aimag(delta0(i)),z0
     write(10,*),abs(delta0(i)),acos(real(delta0(i))/abs(delta0(i)))
     end do	 

     end if
   
	 
	! open(30,file='xyz0_position.txt')
	! do i=1,ng
	! x(i)=0.0d0+(i-1)*2.0d0/(ng-1)
	! do j=1,ng
	!y(j)=0.0d0+(j-1)*2.0d0/(ng-1)
	!write(30,200),x(i),y(j),z0
	!end do
	!end do
!100  format(f15.10,f15.10) 
!200  format(3f12.6)
     
	 close(10) 
	 close(30)
	 close(40)
End Subroutine


