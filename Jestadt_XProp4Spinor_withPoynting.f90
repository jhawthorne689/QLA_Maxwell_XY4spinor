!...1D Y-Prop Jestadt Algorithm (4 spinor for 1D)
!...Z-polarization

    Program Jestadt_YProp
    implicit none
    double precision, parameter  ::  pi = acos(-1.0)
    double precision, parameter  ::  epslon = 0.3  
    double precision, parameter :: n1 = 1.0000d0
    double precision, parameter :: n2 = 2.0000d0
    integer, parameter :: L =6000
    integer, parameter :: L1= 3000, L2 = 3500
    
    double complex,  parameter  ::  icc = dcmplx(0., 1.0d0)
    double precision, dimension(0:L) ::   c11, s11, c11v, s11v
    
    integer  ::  time, tstep,ntime, output, iy , yinit, yfinal, ymid, ict,i,k,kk,ii,n,m, outputc
    
    double precision ::  ampl, theta, shrp, slpe !,waveno
    double precision ::  theta_v
    
    double precision :: Poynt, Poynt_st, Poyntr
    double precision, dimension(0:L)  ::  Ex, Ey, Ez, Bx, By, Bz, nindex, vac
    double precision, dimension(0:L)  ::  Exv, Eyv, Ezv, Bxv, Byv, Bzv
    double complex,dimension(0:L)  ::  inhomog,cinhomo, sinhomo
    double complex,dimension(0:L)  ::  inhomog_v,cinhomo_v, sinhomo_v
    double complex,dimension(0:L) :: q0, q1, q2, q3, q4, q5
    double complex,dimension(0:L) :: q0v, q1v, q2v, q3v, q4v, q5v
    double complex,dimension(0:L) :: q0t,q2t,q3t,q5t
    double complex,dimension(0:L) :: q0tv,q2tv,q3tv,q5tv
    
    ntime = 10000 ; output = 0  ;  outputc = 0
    tstep = 200
    
  
     shrp = 1.00

     do i = 0, L
        vac(i) = 1.00
     end do

     do i = 0, L2
        nindex(i) = 0.5*(n1+n2) - 0.5*(n1-n2)*Tanh( 0.5*shrp*(i-L1))
     enddo
    
     do i = L2 - 100, L
        nindex(i) = 0.5*(n1+n2) + 0.5*(n1-n2)*Tanh( 0.5*shrp*(i-L2))
     enddo
    
     do i = 0, L  
     write(9,*) i,nindex(i)
     enddo
     
     Ex = 0. ;  Ey = 0. ; Ez = 0. ; Bx = 0. ; By = 0. ; Bz = 0.
     Exv = 0. ;  Eyv = 0. ; Ezv = 0. ; Bxv = 0. ; Byv = 0. ; Bzv = 0.
     q0 = 0.; q2 =0.; q3 = 0.; q5 = 0.
     q0v = 0.; q2v = 0.; q3v = 0.; q5v = 0.
     
     ampl = 0.01d0
     yinit = 2300 ; yfinal = 2300;  ymid = (yfinal+ yinit)/2
     
     do iy = 0, L1
        Ez(iy) = ampl*exp(-epslon*epslon*(iy-ymid)*(iy-ymid)/1500.) 
     enddo
     Bx = nindex*Ez !...Z-polarization
     Ezv = Ez
     Bxv = Ezv

      q0 = (nindex*Ex + icc*Bx)   
      q1 = (nindex*Ey + icc*By)    
      q2 = (nindex*Ez + icc*Bz)
      q3 = (nindex*Ex - icc*Bx)
      q4 = (nindex*Ey - icc*By)
      q5 = (nindex*Ez - icc*Bz)
      
     !vaccuum field qubit-field components
      q0v = (Exv + icc*Bxv)   
      q1v = (Eyv + icc*Byv)    
      q2v = (Ezv + icc*Bzv)
      q3v = (Exv - icc*Bxv)
      q4v = (Eyv - icc*Byv)
      q5v = (Ezv - icc*Bzv)
      
      do m = 0,L
          theta =  epslon*0.25/nindex(m) 
          theta_v = epslon*0.25
          
          c11(m) = cos(theta)   ;   s11(m) = sin(theta)
          c11v(m) = cos(theta_v) ; s11v(m) = sin(theta_v)
          
      enddo
      
      inhomog =  + icc *epslon *0.25*( Cshift(nindex, +1) - Cshift(nindex, -1))/(nindex*nindex)
      cinhomo = cos(inhomog)  ;  sinhomo = sin(inhomog)
      
      inhomog_v = + icc *epslon *0.25*(Cshift(vac, +1) - Cshift(vac, -1)) !vac exists because of this line
      cinhomo_v = cos(inhomog_v) ; sinhomo_v = sin(inhomog_v)
      
      
      !start of main do loop
      do time = 0, ntime
     
      q0t = c11*q0 - icc*s11*q2
      q2t = c11*q2 - icc*s11*q0
      q3t = c11*q3 + icc*s11*q5             !...s03py.Cy
      q5t = c11*q5 + icc*s11*q3
      q0 = Cshift(q0t, +1); q3 = Cshift(q3t, +1)
      q2 = q2t; q5 = q5t
      
      
      q0t = c11*q0 + icc*s11*q2
      q2t = c11*q2 + icc*s11*q0
      q3t = c11*q3 - icc*s11*q5             !...s03my.Cya
      q5t = c11*q5 - icc*s11*q3
      q0 = Cshift(q0t, -1); q3 = Cshift(q3t, -1)
      q2 = q2t; q5 = q5t
      
      q0t = c11*q0 - icc*s11*q2
      q2t = c11*q2 - icc*s11*q0
      q3t = c11*q3 + icc*s11*q5             !...s25my.Cy
      q5t = c11*q5 + icc*s11*q3
      q2 = Cshift(q2t, -1); q5 = Cshift(q5t, -1)
      q0 = q0t; q3 = q3t
      
      
      q0t = c11*q0 + icc*s11*q2
      q2t = c11*q2 + icc*s11*q0
      q3t = c11*q3 - icc*s11*q5             !...s25py.Cya
      q5t = c11*q5 - icc*s11*q3
      q2 = Cshift(q2t, +1); q5 = Cshift(q5t, +1)
      q0 = q0t; q3 = q3t;
      
      
      q0t = c11*q0 + icc*s11*q2
      q2t = c11*q2 + icc*s11*q0
      q3t = c11*q3 - icc*s11*q5             !...s03my.Cya
      q5t = c11*q5 - icc*s11*q3
      q0 = Cshift(q0t, -1); q3 = Cshift(q3t, -1)
      q2 = q2t; q5 = q5t
      
      
      q0t = c11*q0 - icc*s11*q2
      q2t = c11*q2 - icc*s11*q0
      q3t = c11*q3 + icc*s11*q5             !...s03py.Cy
      q5t = c11*q5 + icc*s11*q3
      q0 = Cshift(q0t, +1); q3 = Cshift(q3t, +1)
      q2 = q2t; q5 = q5t
      
      
      q0t = c11*q0 + icc*s11*q2
      q2t = c11*q2 + icc*s11*q0
      q3t = c11*q3 - icc*s11*q5             !...s25py.Cya
      q5t = c11*q5 - icc*s11*q3
      q2 = Cshift(q2t, +1); q5 = Cshift(q5t, +1)
      q0 = q0t; q3 = q3t
      
      
      q0t = c11*q0 - icc*s11*q2
      q2t = c11*q2 - icc*s11*q0
      q3t = c11*q3 + icc*s11*q5             !...s25my.Cy
      q5t = c11*q5 + icc*s11*q3
      q2 = Cshift(q2t, -1); q5 = Cshift(q5t, -1)
      q0 = q0t; q3 = q3t
      
      
      q0t = cinhomo*q0 + sinhomo*q2
      q2t = -sinhomo*q0 + cinhomo*q2
      q3t = cinhomo*q3 - sinhomo*q5
      q5t = sinhomo*q3 + cinhomo*q5
      
      q0 = cinhomo*q0t + sinhomo*q5t
      q2 = cinhomo*q2t - sinhomo*q3t
      q3 = -sinhomo*q2t + cinhomo*q3t
      q5 = sinhomo*q0t + cinhomo*q5t
      
!      ...end of QLA portion of do-loop
!      
!      
!      ...start of vaccuum portion calculation
      
      if (time == 1400) then
            Ezv = Ez
            Bxv = Bx
      end if
      
      q0tv = c11v*q0v - icc*s11v*q2v
      q2tv = c11v*q2v - icc*s11v*q0v
      q3tv = c11v*q3v + icc*s11v*q5v             !...s03py.Cy
      q5tv = c11v*q5v + icc*s11v*q3v
      q0v = Cshift(q0tv, +1); q3v = Cshift(q3tv, +1)
      q2v = q2tv; q5v = q5tv
      
      
      q0tv = c11v*q0v + icc*s11v*q2v
      q2tv = c11v*q2v + icc*s11v*q0v
      q3tv = c11v*q3v - icc*s11v*q5v             !...s03my.Cya
      q5tv = c11v*q5v - icc*s11v*q3v
      q0v = Cshift(q0tv, -1); q3v = Cshift(q3tv, -1)
      q2v = q2tv; q5v = q5tv
      
      q0tv = c11v*q0v - icc*s11v*q2v
      q2tv = c11v*q2v - icc*s11v*q0v
      q3tv = c11v*q3v + icc*s11v*q5v             !...s25my.Cy
      q5tv = c11v*q5v + icc*s11v*q3v
      q2v = Cshift(q2tv, -1); q5v = Cshift(q5tv, -1)
      q0v = q0tv; q3v = q3tv
      
      
      q0tv = c11v*q0v + icc*s11v*q2v
      q2tv = c11v*q2v + icc*s11v*q0v
      q3tv = c11v*q3v - icc*s11v*q5v             !...s25py.Cya
      q5tv = c11v*q5v - icc*s11v*q3v
      q2v = Cshift(q2tv, +1); q5v = Cshift(q5tv, +1)
      q0v = q0tv; q3v = q3tv;
      
      
      q0tv = c11v*q0v + icc*s11v*q2v
      q2tv = c11v*q2v + icc*s11v*q0v
      q3tv = c11v*q3v - icc*s11v*q5v             !...s03my.Cya
      q5tv = c11v*q5v - icc*s11v*q3v
      q0v = Cshift(q0tv, -1); q3v = Cshift(q3tv, -1)
      q2v = q2tv; q5v = q5tv
      
      
      q0tv = c11v*q0v - icc*s11v*q2v
      q2tv = c11v*q2v - icc*s11v*q0v
      q3tv = c11v*q3v + icc*s11v*q5v             !...s03py.Cy
      q5tv = c11v*q5v + icc*s11v*q3v
      q0v = Cshift(q0tv, +1); q3v = Cshift(q3tv, +1)
      q2v = q2tv; q5v = q5tv
      
      
      q0tv = c11v*q0v + icc*s11v*q2v
      q2tv = c11v*q2v + icc*s11v*q0v
      q3tv = c11v*q3v - icc*s11v*q5v             !...s25py.Cya
      q5tv = c11v*q5v - icc*s11v*q3v
      q2v = Cshift(q2tv, +1); q5v = Cshift(q5tv, +1)
      q0v = q0tv; q3v = q3tv
      
      
      q0tv = c11v*q0v - icc*s11v*q2v
      q2tv = c11v*q2v - icc*s11v*q0v
      q3tv = c11v*q3v + icc*s11v*q5v             !...s25my.Cy
      q5tv = c11v*q5v + icc*s11v*q3v
      q2v = Cshift(q2tv, -1); q5v = Cshift(q5tv, -1)
      q0v = q0tv; q3v = q3tv
      
      
      q0tv = cinhomo_v*q0v + sinhomo_v*q2v
      q2tv = -sinhomo_v*q0v + cinhomo_v*q2v
      q3tv = cinhomo_v*q3v - sinhomo_v*q5v
      q5tv = sinhomo_v*q3v + cinhomo_v*q5v
      
      q0v = cinhomo_v*q0tv + sinhomo_v*q5tv
      q2v = cinhomo_v*q2tv - sinhomo_v*q3tv
      q3v = -sinhomo_v*q2tv + cinhomo_v*q3tv
      q5v = sinhomo_v*q0tv + cinhomo_v*q5tv
      
      
      
      Ez = 0.5*(q2 + q5)/nindex
      Bx = 0.5* AIMAG(q0 - q3) !this calculates the fields components from qubits for QLA
      
      Ezv = 0.5*(q2v + q5v)  !we don't have to calc this, no output necessary but for checking correct
      Bxv = 0.5* AIMAG(q0v - q3v)
        
    
      
      if (time == outputc) then
         Poynt = 0.0
         Poynt_st = 0.0
         Poyntr = 0.0
      
        if (time >= 1400 .and. time <= 3400) then

            do iy = 0, L1
                Poyntr = Poyntr + ABS(Ezv(iy)- Ez(iy))*ABS(Bxv(iy)- Bx(iy))
                ! Ez field component switches sign when pulse hits front of dielectric
                ! but keeping the abs on both will make it easier to extend later
            end do
            
            do iy = 0, L1
                Poyntr = Poyntr + ABS(Ezv(iy)*Bxv(iy))
            end do
            
            !print *, time, Poyntr
            Poynt = Poyntr
            
            do iy = L1, L
                Poynt = Poynt + ABS(Ez(iy)*Bx(iy))
            end do

        else
            do iy = 0,L
                Poynt = Poynt + ABS( Ez(iy)*Bx(iy) )
            end do
        end if
      
!      do iy = 0, L
!            Poynt_st = Poynt_st + Ez(iy)*Bx(iy)
!      end do
      
        !print *, time, Poynt
        write(34,*) time, Poynt 
        !write(35,*) time, Poynt_st
        outputc = outputc + tstep/10
      
      endif
      
          
      if(time == output)then
            do iy = 0, L
            write(100000+time,*)iy,Bx(iy)
            write(200000+time,*)iy,Ez(iy)
!            write(800000+time,*)iy,Bxv(iy)
!            write(900000+time,*)iy,Ezv(iy)
            enddo
            output = output + tstep
      endif
      
      enddo
           
    write(6,*)  "finished run", n1, n2, epslon 
    end Program Jestadt_YProp
