!***************************************************************
!OH-[H2O] Potential Energy Surface, version BHB-1, 2003-10-15
!   Made by Bastiaan J. Braams, Xinchuan Huang, Joel M. Bowman
!   (Intrinsic version: PES-2)
!
!reference:
!   Huang, X.; Braams, B.J.; Carter, S.; Bowman, J.M.  
!     J. Am. Chem. Soc. (2004),  126(16),  5042-5043
!
!   CCSD(T)/aug-cc-pVTZ based, without any further adjustments
!
!notes:
!   Requiring h3o2.pes2.ifchcm.dat present in same directory
!
!   Fortran 90 compiler required 
!   (Testified with pgf90 v.4.x - v.5.x)   
!   
!   The prepot() should be called before first time calling 
!   getpot(V,cart_in) 
!
!   Cart_in(1:5,3) is the cartesian x,y,z- coordinates for five
!   atoms in order of H O H O H. (in bohr)
!
!   Returned V is in a.u., related to H3O2- minimum.
!
!All rights reserved. Contact bowman@euch4e.chem.emory.edu for
!details or any latest upgrades.
!***************************************************************

        subroutine prepot() 

!        implicit none
        implicit double precision (a-h,o-z)
        implicit integer (i-n) 

        double precision V, cart_in(5,3), v0, coef(0:2020)

        integer i,j,k,i1,j1,k1

        common/NCOE/m,mr
        common/h3o2coef/dc0,dc1,dc2,dw0,dw1,dw2,coef

        m=1976 ; mr=15

        open(20,file='h3o2.pes2.ifchcm.dat',status='old')
        read(20,*)
        read(20,*)
        read(20,*)
        read(20,*)dc0
        read(20,*)dc1
        read(20,*)dc2
        read(20,*)dw0
        read(20,*)dw1
        read(20,*)dw2
        read(20,*)
        read(20,*)(coef(i1),i1=0,m+3*mr-1)
!        write(*,*)(coef(i1),i1=0,m+3*mr-1)
        close(20)

        return
        end 

!*******************************************************************


        subroutine getpot(V,cart_in)

        implicit none

        integer m,mr
        double precision dc0,dc1,dc2,dw0,dw1,dw2,coef(0:2020)

        common/NCOE/m,mr
        common/h3o2coef/dc0,dc1,dc2,dw0,dw1,dw2,coef

        double precision, intent(out):: V
        double precision, intent(in) :: cart_in(5,3)

        double precision cart0(5,3),cart1(5,3)
        integer i,j,k,l,i1,j1,k1,l1,i2,j2,k2,l2
        

        double precision rvec(0:3),d0(0:4,0:4),r0(0:4,0:4),
     $  vec(0:m+3*mr-1)
        double precision xnuc(0:2,0:4)

         do j=1,3
          xnuc(j-1,0)=cart_in(2,j)
          xnuc(j-1,1)=cart_in(4,j)
          xnuc(j-1,2)=cart_in(1,j)
          xnuc(j-1,3)=cart_in(3,j)
          xnuc(j-1,4)=cart_in(5,j)
        end do

        call getd0 (xnuc, dc0, dc1, dc2, dw0, dw1, dw2, d0, r0)
        call getvec (m, d0, vec(0:m-1))
        call getrvec (4, r0, rvec)
        do l = 0, mr-1
          do k = 0, 2
            vec(m+3*l+k) = rvec(k+1)*vec(l)
          enddo
        enddo
        V = dot_product(coef,vec)

! set C1 min as potential zero-point
        V=(V+70.5692461013296d0)/219474.63067d0
 
! set C2 linear SP as potential zero-point
!        V=(V+2.031866418016736d0)/219474.63067d0

        return
        end subroutine getpot 
!********************************************************

      subroutine getd0 (xn, dc0, dc1, dc2, dw0, dw1, dw2, d0, r0)
      implicit none
      double precision xn(0:2,0:4), dc0, dc1, dc2, dw0, dw1, dw2,
     $  d0(0:4,0:4),  r0(0:4,0:4)
      integer i, j
      double precision t0
! Note: H3O2-.  Nuclei 0..4; i in 0..1 for the O, i in 2..4 for the H.
      do i = 0, 4
       d0(i,i) = 0.d0
       r0(i,i) = 0.d0
      enddo
      t0 = dsqrt(dot_product(xn(0:2,1)-xn(0:2,0),xn(0:2,1)-xn(0:2,0)))
      d0(0,1) = (dlog(t0)-dc0)/dw0
      d0(1,0) = d0(0,1)
      r0(0,1) = t0
      r0(1,0) = t0
      do j = 2, 4
       do i = 0, 1
        t0 = dsqrt(dot_product(xn(0:2,j)-xn(0:2,i),xn(0:2,j)-xn(0:2,i)))
        d0(i,j) = (dlog(t0)-dc1)/dw1
        d0(j,i) = d0(i,j)
        r0(i,j) = t0
        r0(j,i) = t0
       enddo
      enddo
      do i = 2, 4
       do j = i+1, 4
        t0 = dsqrt(dot_product(xn(0:2,j)-xn(0:2,i),xn(0:2,j)-xn(0:2,i)))
        d0(i,j) = (dlog(t0)-dc2)/dw2
        d0(j,i) = d0(i,j)
        r0(i,j) = t0
        r0(j,i) = t0
       enddo
      enddo
 
!      write(2001,*)t0,d0(1,0),r0(3,4),d0(1,2)

      return
      end
!*****************************************************************

      subroutine getvec (m, d, vec)
      implicit none
      integer m
      double precision d(0:4,0:4), vec(0:m-1)
      integer j0, j1, j2, j3, j4, i, j, k, k3, k4, k5, k6, k7
      double precision x(0:2), y(0:4), z(0:7), u(0:6), v(0:1), t0
      double precision d2(0:4,0:4), d3(0:4,0:4),d4(0:4,0:4),d5(0:4,0:4)
      double precision her2, her3, her4, her5 , her6, her7
      her2(t0) = (4*t0**2-2)/dsqrt(dble(8*2))
      her3(t0) = (8*t0**2-12)*t0/dsqrt(dble(16*6))
      her4(t0) = ((16*t0**2-48)*t0**2+12)/dsqrt(dble(32*24))
      her5(t0) = ((32*t0**2-160)*t0**2+120)*t0/dsqrt(dble(64*120))
      her6(t0) = (((64*t0**2-480)*t0**2+720)*t0**2-120)/dsqrt(
     $dble(128*720))
      her7(t0) = (((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/dsqrt(
     $dble(256*5040))
!-----------------------------------------------------------------------
! version for H3O2-
!-----------------------------------------------------------------------
! Test for compatibility
      if (.not.(m.eq.1.or.m.eq.4.or.m.eq.15.or.m.eq.48.or. 
     $   m.eq.139.or.m.eq.364.or.m.eq.869.or.m.eq.1976)) then
       stop 'getvec - wrong dimension'
      endif
! Preparation
      do j0 = 0, 4
       do j1 = j0+1, 4
        d2(j0,j1) = her2(d(j0,j1))
        d2(j1,j0) = d2(j0,j1)
        d3(j0,j1) = her3(d(j0,j1))
        d3(j1,j0) = d3(j0,j1)
        d4(j0,j1) = her4(d(j0,j1))
        d4(j1,j0) = d4(j0,j1)
        d5(j0,j1) = her5(d(j0,j1))
        d5(j1,j0) = d5(j0,j1)
       enddo
      enddo
      x = 0.d0 ; y = 0.d0 ; z = 0.d0 ; u = 0.d0 ; v = 0.d0
! Generators using two Oxygens
      do j0 = 0, 1
       do j1 = 0, 1
       if (j1.ne.j0) then
        x(0) = x(0)+d(j0,j1)/2
        do j2 = 2, 4
         y(0) = y(0)+d(j0,j2)*d(j1,j2)/6
         z(0) = z(0)+d2(j0,j2)*d(j1,j2)/6
         u(0) = u(0)+d3(j0,j2)*d(j1,j2)/6
         do j3 = 2, 4
         if (j3.ne.j2) then
          z(1) = z(1)+d(j0,j2)*d(j1,j2)*d(j2,j3)/12
          z(2) = z(2)+d(j0,j2)*d(j1,j3)*d(j2,j3)/12
          do j4 = 2, 4
          if (j4.ne.j3.and.j4.ne.j2) then
           u(6) = u(6)+d(j0,j2)*d(j0,j3)*d(j1,j4)*d(j2,j4)/12
          endif
          enddo
         endif
         enddo
        enddo
       endif
       enddo
      enddo
! Generators using one Oxygen
      do j0 = 0, 1
       do j1 = 2, 4
        x(1) = x(1)+d(j0,j1)/6
        y(1) = y(1)+d2(j0,j1)/6
        z(3) = z(3)+d3(j0,j1)/6
        u(1) = u(1)+d4(j0,j1)/6
        v(0) = v(0)+d5(j0,j1)/6
        do j2 = 2, 4
        if (j2.ne.j1) then
         y(2) = y(2)+d(j0,j1)*d(j1,j2)/12
         y(3) = y(3)+d(j0,j1)*d(j0,j2)/12
         z(4) = z(4)+d2(j0,j1)*d(j0,j2)/12
         z(5) = z(5)+d2(j0,j1)*d(j1,j2)/12
         u(2) = u(2)+d3(j0,j1)*d(j0,j2)/12
         u(3) = u(3)+d3(j0,j1)*d(j1,j2)/12
         v(1) = v(1)+d4(j0,j1)*d(j1,j2)/12
         do j3 = 2, 4
         if (j3.ne.j2.and.j3.ne.j1) then
          z(6) = z(6)+d(j0,j1)*d(j1,j2)*d(j1,j3)/12
          u(4) = u(4)+d2(j0,j1)*d(j1,j2)*d(j1,j3)/12
          u(5) = u(5)+d(j0,j1)*d(j1,j2)*d(j1,j3)*d(j0,j2)/12
         endif
         enddo
        endif
        enddo
       enddo
      enddo
! Generators using no Oxygen
      do j0 = 2, 4
       do j1 = 2, 4
       if (j1.ne.j0) then
        x(2) = x(2)+d(j0,j1)/6
        y(4) = y(4)+d2(j0,j1)/6
        z(7) = z(7)+d3(j0,j1)/6
       endif
       enddo
      enddo
! Compute vec(0:m-1)
! constant term
      vec(0) = 1
! first degree terms
      if (4.le.m) then
       vec(1:3) = x(0:2)
      endif
! second degree terms
      if (15.le.m) then
       vec(4) = her2(x(0))
       vec(5) = x(0)*x(1)
       vec(6) = x(0)*x(2)
       vec(7) = her2(x(1))
       vec(8) = x(1)*x(2)
       vec(9) = her2(x(2))
       vec(10:14) = y(0:4)
      endif
! third degree terms
      if (48.le.m) then
       k3 = 15
       vec(k3) = her3(x(0))
       vec(k3+1) = her2(x(0))*x(1)
       vec(k3+2) = her2(x(0))*x(2)
       vec(k3+3) = x(0)*her2(x(1))
       vec(k3+4) = x(0)*x(1)*x(2)
       vec(k3+5) = x(0)*her2(x(2))
       vec(k3+6:k3+10) = x(0)*y(0:4)
       vec(k3+11) = her3(x(1))
       vec(k3+12) = her2(x(1))*x(2)
       vec(k3+13) = x(1)*her2(x(2))
       vec(k3+14:k3+18) = x(1)*y(0:4)
       vec(k3+19) = her3(x(2))
       vec(k3+20:k3+24) = x(2)*y(0:4)
       vec(k3+25:k3+32) = z(0:7)
      endif
! fourth degree terms
      if (139.le.m) then
       k4 = 48
       vec(k4) = her4(x(0))
       vec(k4+1) = her3(x(0))*x(1)
       vec(k4+2) = her3(x(0))*x(2)
       vec(k4+3) = her2(x(0))*her2(x(1))
       vec(k4+4) = her2(x(0))*x(1)*x(2)
       vec(k4+5) = her2(x(0))*her2(x(2))
       vec(k4+6:k4+10) = her2(x(0))*y(0:4)
       vec(k4+11) = x(0)*her3(x(1))
       vec(k4+12) = x(0)*her2(x(1))*x(2)
       vec(k4+13) = x(0)*x(1)*her2(x(2))
       vec(k4+14:k4+18) = x(0)*x(1)*y(0:4)
       vec(k4+19) = x(0)*her3(x(2))
       vec(k4+20:k4+24) = x(0)*x(2)*y(0:4)
       vec(k4+25:k4+32) = x(0)*z(0:7)
       vec(k4+33) = her4(x(1))
       vec(k4+34) = her3(x(1))*x(2)
       vec(k4+35) = her2(x(1))*her2(x(2))
       vec(k4+36:k4+40) = her2(x(1))*y(0:4)
       vec(k4+41) = x(1)*her3(x(2))
       vec(k4+42:k4+46) = x(1)*x(2)*y(0:4)
       vec(k4+47:k4+54) = x(1)*z(0:7)
       vec(k4+55) = her4(x(2))
       vec(k4+56:k4+60) = her2(x(2))*y(0:4)
       vec(k4+61:k4+68) = x(2)*z(0:7)
       vec(k4+69) = her2(y(0))
       vec(k4+70) = y(0)*y(1)
       vec(k4+71) = y(0)*y(2)
       vec(k4+72) = y(0)*y(3)
       vec(k4+73) = y(0)*y(4)
       vec(k4+74) = her2(y(1))
       vec(k4+75) = y(1)*y(2)
       vec(k4+76) = y(1)*y(3)
       vec(k4+77) = y(1)*y(4)
       vec(k4+78) = her2(y(2))
       vec(k4+79) = y(2)*y(3)
       vec(k4+80) = y(2)*y(4)
       vec(k4+81) = her2(y(3))
       vec(k4+82) = y(3)*y(4)
       vec(k4+83) = her2(y(4))
       vec(k4+84:k4+90) = u(0:6)
      endif
! fifth degree terms
      if (364.le.m) then
       k5=139
       vec(k5:k5+90) = x(0)*vec(k4:k4+90)
       vec(k5+91:k5+148) = x(1)*vec(k4+33:k4+90)
       vec(k5+149:k5+184) = x(2)*vec(k4+55:k4+90)
       vec(k5+185:k5+192) = y(0)*vec(k3+25:k3+32)
       vec(k5+193:k5+200) = y(1)*vec(k3+25:k3+32)
       vec(k5+201:k5+208) = y(2)*vec(k3+25:k3+32)
       vec(k5+209:k5+216) = y(3)*vec(k3+25:k3+32)
       vec(k5+217) = y(4)*vec(k3+25)
       vec(k5+218) = y(4)*vec(k3+26)
       vec(k5+219) = y(4)*vec(k3+27)
       vec(k5+220) = y(4)*vec(k3+29)
       vec(k5+221) = y(4)*vec(k3+31)
       vec(k5+222) = y(4)*vec(k3+32)
       vec(k5+223) = v(0)
       vec(k5+224) = v(1)
      endif
! sixth degree terms
! (according to Kemper there should be 525 of them; here I'm computing
! 505 terms from the lower-degree ones, but only 523 of them are
! independent)
      if (869.le.m) then
       k6 = 364
       vec(k6:k6+224) = x(0)*vec(k5:k5+224)
       vec(k6+225:k6+358) = x(1)*vec(k5+91:k5+224)
       vec(k6+359:k6+434) = x(2)*vec(k5+149:k5+224)
       vec(k6+435:k6+454) = y(0)*vec(k4+69:k4+88)  !!!
       vec(k6+455:k6+473) = y(1)*vec(k4+74:k4+92)  !!!
       vec(k6+474:k6+486) = y(2)*vec(k4+78:k4+90)
       vec(k6+487:k6+496) = y(3)*vec(k4+81:k4+90)
       vec(k6+497:k6+504) = y(4)*vec(k4+83:k4+90)
      endif
! seventh degree terms
! (neither complete nor independent)
      if (1976.le.m) then
       k7 = 869
       vec(k7:k7+504) = x(0)*vec(k6:k6+504)
       vec(k7+505:k7+784) = x(1)*vec(k6+225:k6+504)
       vec(k7+785:k7+930) = x(2)*vec(k6+359:k6+504)
       vec(k7+931:k7+970) = y(0)*vec(k5+185:k5+224)
       vec(k7+971:k7+1002) = y(1)*vec(k5+193:k5+224)
       vec(k7+1003:k7+1026) = y(2)*vec(k5+201:k5+224)
       vec(k7+1027:k7+1042) = y(3)*vec(k5+209:k5+224)
       vec(k7+1043:k7+1050) = y(4)*vec(k5+217:k5+224)
       vec(k7+1051:k7+1057) = z(0)*vec(k4+84:k4+90)
       vec(k7+1058:k7+1064) = z(1)*vec(k4+84:k4+90)
       vec(k7+1065:k7+1071) = z(2)*vec(k4+84:k4+90)
       vec(k7+1072:k7+1078) = z(3)*vec(k4+84:k4+90)
       vec(k7+1079:k7+1085) = z(4)*vec(k4+84:k4+90)
       vec(k7+1086:k7+1092) = z(5)*vec(k4+84:k4+90)
       vec(k7+1093:k7+1099) = z(6)*vec(k4+84:k4+90)
       vec(k7+1100:k7+1106) = z(7)*vec(k4+84:k4+90)
      endif
!      do i = 1,m
!         write(50,*)i,vec(i)
!      enddo
      return
      end
!*****************************************************

      subroutine getrvec (m, r, vec)
      implicit none
      integer m
      double precision r(0:4,0:4), vec(0:m-1)
      integer j0, j1, j2, i, j, k
      double precision x(0:2), y(0:4), r1(0:4,0:4), r2(0:4,0:4)
!-----------------------------------------------------------------------
! Test for compatibility
      if (.not.(m.eq.1.or.m.eq.4.or.m.eq.15)) then
       stop 'getrvec - wrong dimension'
      endif
! Computation
      x = 0.d0 ; y = 0.d0
      do i = 0, 4
       do j = 0, 4
        if (i.eq.j) then
         r1(i,j) = 0.d0
        else
         r1(i,j) = dexp(-r(i,j))/r(i,j)
        endif
       enddo
      enddo
      do j0 = 0, 1
       do j1 = 0, 1
       if (j1.ne.j0) then
        x(0) = x(0)+r1(j0,j1)/2
        do j2 = 2, 4
         y(0) = y(0)+r1(j0,j2)*r1(j1,j2)/6
        enddo
       endif
       enddo
      enddo
      do j0 = 0, 1
       do j1 = 2, 4
        x(1) = x(1)+r1(j0,j1)/6
        y(1) = y(1)+r2(j0,j1)/6
        do j2 = 2, 4
        if (j2.ne.j1) then
         y(2) = y(2)+r1(j0,j1)*r1(j1,j2)/12
         y(3) = y(3)+r1(j0,j1)*r1(j0,j2)/12
        endif
        enddo
       enddo
      enddo
! Generators using no Oxygen
      do j0 = 2, 4
       do j1 = 2, 4
       if (j1.ne.j0) then
        x(2) = x(2)+r1(j0,j1)/6
        y(4) = y(4)+r2(j0,j1)/6
       endif
       enddo
      enddo
      vec(0) = 1
      if (4.le.m) then
       vec(1:3) = x(0:2)
      endif
      if (15.le.m) then
       vec(4) = x(0)**2
       vec(5) = x(0)*x(1)
       vec(6) = x(0)*x(2)
       vec(7) = x(1)**2
       vec(8) = x(1)*x(2)
       vec(9) = x(2)**2
       vec(10:14) = y(0:4)
      endif
      return
      end
!***************************************************

