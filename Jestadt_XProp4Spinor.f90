!...Jestadt Six Spinor QLA
!...1D X-prop algorithm, y-polarization only

    Program Jestadt 
    implicit none
    double precision, parameter  ::  pi = acos(-1.0)
    double precision, parameter  ::  epslon = 0.3  
    double precision, parameter :: n1 = 1.0000d0
    double precision, parameter :: n2 = 2.5000d0
    integer, parameter :: L =6000
    integer, parameter :: L1= 3000, L2 = 4000
    
    double complex,  parameter  ::  icc = dcmplx(0., 1.0d0)
    double precision, dimension(0:L) ::   c11, s11

    integer  ::  time, tstep,ntime, output, iy , yinit, yfinal, ymid, ict,i,k,kk,ii,n,m, outputc
        
    double precision ::  waveno, ampl, theta, E1,E2,B1,B2, vwinh, shrp, slpe
    double precision ::  Poynt1, Debh1, Poynt2, Debh2 ,Poynt3, Debh3
    double precision, dimension(0:L)  ::  Exf, Eyf, Ezf, Bxf, Byf , Bzf, nindex
    !...double precision, dimension(0:L1)  ::  Ey1, Bz1, Ey1a, Bz1a, Ey2a, Bz2a    !...double precision, dimension(L1:L2)  ::  Ey2, Bz2
    !...double precision, dimension(L2:L)  ::  Ey3, Bz3
    double complex,dimension(0:L)  ::  inhomog,cinhomo, sinhomo
    double complex,dimension(0:L) :: Fxrsp, Fyrsp, Fzrsp, q0, q1, q2, q3, q4, q5 
    double complex,dimension(0:L) :: Fxrsm, Fyrsm, Fzrsm
    double complex,dimension(0:L) :: q1t,q2t,q4t,q5t
    ntime = 16000 ; output = 0  ;  outputc = 0
    tstep = 200


 !   do i = 0, L
 1   nindex(i) = 1.00d0
 !   enddo
  
  !  ...  refractive index of medium :  sharper refractrive index change                                                  
  shrp = 0.22 !  the smaller this is, the wider the refractive index part                                               
    shrp = 1.00
    do i = 0, L2
 !    nindex(i) = 0.5*(n1+n2) - 0.5*(n1-n2)*Tanh( 0.5*shrp*(i-0.5*L2))                                                  
    nindex(i) = 0.5*(n1+n2) - 0.5*(n1-n2)*Tanh( 0.5*shrp*(i-L1))
       enddo
    do i = L2 - 100, L
      nindex(i) = 0.5*(n1+n2) + 0.5*(n1-n2)*Tanh( 0.5*shrp*(i-L2))
        enddo
  
    do i = 1, L  
    write(9,*) i,nindex(i)
    enddo
    
    
        Exf = 0. ;  Eyf = 0. ; Ezf = 0. ; Bxf = 0. ; Byf = 0. ; Bzf = 0.
        !...Ey1 = 0.  ;  Ey2 = 0.  ;  Bz1 = 0.  ; Bz2 = 0.
        q1 = 0.; q2 =0.; q4 = 0.; q5 = 0.
        
        ampl = 0.01d0
        yinit = 2300 ; yfinal = 2300;  ymid = (yfinal+ yinit)/2
  !  ... different polarization  (-Ez, By)      
        do iy = 0, L1
        Ezf(iy) =- ampl*exp(-epslon*epslon*(iy-ymid)*(iy-ymid)/1500.) 
        enddo
        Byf = - nindex*Ezf 
        
        
        q0 = (nindex*Exf + icc*Bxf)   !  John uses different normalization 
        q1 = (nindex*Eyf + icc*Byf)    ! ...  factor  sqrt(0.5)
        q2 =(nindex*Ezf + icc*Bzf)
           q3 = (nindex*Exf - icc*Bxf)
           q4 = (nindex*Eyf - icc*Byf)
           q5 = (nindex*Ezf - icc*Bzf)
           
        do m = 0,L
          theta = - epslon*0.25/nindex(m)    
          c11(m) = cos(theta)   ;   s11(m) = sin(theta)
        enddo  
           
        inhomog =  + icc *epslon *0.25*( Cshift(nindex, +1) - Cshift(nindex, -1))/(nindex*nindex)
        cinhomo = cos(inhomog)  ;  sinhomo = sin(inhomog)
    
    do time = 0, ntime
        
        q1t = c11*q1 - icc*s11*q2
        q2t = c11*q2 - icc*s11*q1 
        q4t = c11*q4 + icc*s11*q5
        q5t = icc*s11*q4 + c11*q5
        q1 = Cshift(q1t, +1) ; q4 = Cshift(q4t, +1)
        q2 = q2t ; q5 = q5t
        
        
        q1t = c11*q1 + icc*s11*q2
        q2t = icc*s11*q1 + c11*q2 
        q4t = c11*q4 - icc*s11*q5
        q5t = c11*q5 - icc*s11*q4 
        q1 = Cshift(q1t, -1) ; q4 = Cshift(q4t, -1)
        q2 = q2t; q5 = q5t 
        
        
        q1t = c11*q1 - icc*s11*q2
        q2t = c11*q2 - icc*s11*q1
        q4t = c11*q4 + icc*s11*q5
        q5t = icc*s11*q4 + c11*q5
        q2 = Cshift(q2t, -1) ; q5 = Cshift(q5t, -1)
        q1 = q1t; q4 = q4t 
        
        
        q1t = c11*q1 + icc*s11*q2
        q2t = icc*s11*q1 + c11*q2
        q4t = c11*q4 - icc*s11*q5
        q5t = c11*q5 - icc*s11*q4 
        q2 = Cshift(q2t, +1) ; q5 = Cshift(q5t, +1)
        q1 = q1t; q4 = q4t 
        
        
        q1t = c11*q1 + icc*s11*q2
        q2t = icc*s11*q1 + c11*q2
        q4t = c11*q4 - icc*s11*q5
        q5t = c11*q5 - icc*s11*q4 
        q1 = Cshift(q1t, -1) ; q4 = Cshift(q4t, -1)
        q2 = q2t ; q5 = q5t
        
        
        q1t = c11*q1 - icc*s11*q2
        q2t = c11*q2 - icc*s11*q1 
        q4t = c11*q4 + icc*s11*q5
        q5t = icc*s11*q4 + c11*q5 
        q1 = Cshift(q1t,+1) ; q4 = Cshift(q4t, +1)
        q2 = q2t ; q5 = q5t
        
        
        q1t = c11*q1 + icc*s11*q2 
        q2t = icc*s11*q1 + c11*q2 
        q4t = c11*q4 - icc*s11*q5
        q5t = c11*q5 - icc*s11*q4 
        q2 = Cshift(q2t, +1) ; q5 = Cshift(q5t, +1)
        q1 = q1t ; q4 = q4t;
        
        
        q1t = c11*q1 - icc*s11*q2
        q2t = c11*q2 - icc*s11*q1  
        q4t = c11*q4 + icc*s11*q5 
        q5t = icc*s11*q4 + c11*q5
        q2 = Cshift(q2t, -1) ; q5 = Cshift(q5t, -1)
        q1 = q1t ; q4 = q4t;
        
        
        q1t = cinhomo*q1 - sinhomo*q2
        q2t = sinhomo*q1 + cinhomo*q2
        q4t = cinhomo*q4 + sinhomo*q5
        q5t = -sinhomo*q4 + cinhomo*q5
        
        q1 = cinhomo*q1t - sinhomo*q5t
        q2 = cinhomo*q2t + sinhomo*q4t
        q4 = sinhomo*q2t + cinhomo*q4t
        q5 = -sinhomo*q1t + cinhomo*q5t
        
        do iy = 0, L
!        Eyf(iy) = 0.5*Real(q1(iy) + q4(iy))/nindex(iy)   !  John has extra factor  sqrt(2)
!        Bzf(iy) = + 0.5*AIMAG(q2(iy) - q5(iy))   !   -ve  ---- changed it to plus
       enddo
        
   Ezf = 0.5* (q2 + q5)/nindex
   Byf = 0.5 * AIMAG(q1 - q4)

        if(time == output)then
            do iy = 0, L
            write(100000+time,*)iy,Byf(iy)
            write(600000+time,*)iy,Ezf(iy)
            enddo
            output = output + tstep
        endif
        
 !       if(time == outputc) then
    !...Ex1a = Abs(Ex1);  By1a = Abs(By1)
    !...Ex2a = Abs(Ex2);  By2a = Abs(By2)
          
    !...write(30,*) time,  maxval(Ex1a), maxval(By1a), maxval(Ex2a), maxval(By2a)
    !...write(31,*) time,  maxloc(Ex1a), maxloc(By1a), maxloc(Ex2a), maxloc(By2a)
 !       write(32,*) time, Poynt1, Poynt2, Abs(Poynt1) + Poynt2
 !       write(33,*)  time, Debh1, Debh2, Debh3, Debh1 + Debh2 + Debh3
 !       write(34,*) time, Abs(Poynt1) + Abs(Poynt2) + Abs(Poynt3)
    !...write(34,*) time, Poynt1 + Poynt2 + Poynt3
 !       outputc = outputc + tstep/10
 !       endif
        
    enddo
    write(6,*)  "finished run", n1, n2, epslon   
      end  program Jestadt
      
           
           
           
