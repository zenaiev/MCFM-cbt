      subroutine dgen2m(r,dp,wt2,wt3,*)
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'mxdim.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'breit.f'
      integer j,nu

      double precision r(mxdim),p(mxpart,4),dp(mxpart,4),xx(2),mass
      double precision sqrts,ymax,yave,ydif,xjac,y3,y4,phi,wt0,wt2,w3,vs
      double precision vsqmin,vsqmax,pt,s34,xmin,rtshat,udif,trmass,beta
      double precision dydif,dyave,dy3,dy4,dxx(2),dtrmass,dbeta,drtshat
      double precision dpt,wt3
      character*6 mscase
      common/mscase/mscase
      common/energy/sqrts
      parameter(wt0=1d0/16d0/pi)
      common/x1x2/xx

      if(mscase .eq. 'numeric') then
        write(6,*) 'Error: cannot calculate analytic derivatives with '
        write(6,*) '   mscase = ',mscase
        stop         
      endif

      do j=1,mxpart     
      do nu=1,4     
      dp(j,nu)=0d0
      enddo     
      enddo   

      wt2=0d0
      
      if (n3 .eq. 1) then
c--- generate s34 according to a Breit-Wigner, for gg->H 
        call breitw(r(3),wsqmin,wsqmax,mass3,width3,s34,w3)
      if     (hdecaymode .eq. 'bqba') then
        mass=mb
      elseif (hdecaymode .eq. 'tlta') then
        mass=mtau
      else
        write(6,*) 'Unanticipated hdecaymode in gen2m: ',hdecaymode
        stop
      endif
      else
c--- no resonance, for tt~,bb~,cc~
        mass=mass2
        vsqmin=4d0*mass**2
        vsqmax=sqrts**2
        xmin=vsqmin/vsqmax
        s34=(vsqmax-vsqmin)*r(3)+vsqmin
        w3=(vsqmax-vsqmin)
      endif

      rtshat=dsqrt(s34)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(1)-1d0)
c----udif=tanh(ydif)
      beta=dsqrt(1d0-4d0*mass**2/s34)

      if((mscase .eq. 'pt4') .or. (mscase .eq. 'm34')) then
         pt=dsqrt(s34/4d0 - mass**2)*r(2)
         dpt=-mass*r(2)**2*r(3)/pt
         udif=dsqrt((sqrts**2-4d0*mass**2)*(1d0-r(2)**2)*r(3)
     .        /s34)
         ydif=half*dlog((1d0+udif)/(1d0-udif))
         trmass=rtshat/(2d0*dcosh(ydif))
         dtrmass=trmass*(4d0*mass*(1-r(2)**2*r(3))
     .        /((sqrts**2-4d0*mass**2)*r(2)**2*r(3)+4d0*mass**2))
         w3=w3*r(2)/dsqrt(1d0-r(2)**2)
      elseif(mscase .eq. 'y4') then
         udif=beta*(two*r(2)-1d0)
         ydif=half*dlog((1d0+udif)/(1d0-udif))
         trmass=rtshat/(2d0*dcosh(ydif))
         dtrmass=trmass*(mass-4d0*mass*(1d0-r(2))*r(2)*r(3))/
     .        (mass**2+(sqrts**2-4d0*mass**2)*(1d0-r(2))*r(2)*r(3))
         pt=dsqrt(trmass**2-mass**2)
         dpt=-4d0*mass*(1d0-r(2))*r(2)*r(3)/pt
      endif

      xjac=four*ymax*beta
          
      y3=yave+ydif
      y4=yave-ydif

      dydif=-udif/(1d0-udif**2)*4d0*mass*sqrts**2
     .     /((sqrts**2-4d0*mass**2)*s34)

      dyave=4d0*mass*(1-2d0*r(1))*(1d0-r(3))/s34
      dy3=dyave+dydif
      dy4=dyave-dydif

      xjac=xjac*w3
      phi=2d0*pi*r(4)

      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)
      dxx(1) = xx(1)*8d0*mass*(1-r(1))*(1-r(3))/s34
      dxx(2) = xx(2)*8d0*mass*(1-r(3))*r(1)/s34

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
        write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1),xx(2)  
        return 1 
      endif

      dp(1,4)=-0.5*dxx(1)*sqrts
      dp(1,1)=0d0
      dp(1,2)=0d0
      dp(1,3)=-0.5d0*dxx(1)*sqrts

      dp(2,4)=-0.5d0*dxx(2)*sqrts
      dp(2,1)=0d0
      dp(2,2)=0d0
      dp(2,3)=+0.5d0*dxx(2)*sqrts
      
      dp(3,4)=+dtrmass*dcosh(y3)+trmass*dy3*dsinh(y3)
      dp(3,1)=+dpt*dsin(phi)
      dp(3,2)=+dpt*dcos(phi)
      dp(3,3)=+dtrmass*dsinh(y3)+trmass*dy3*dcosh(y3)

      dp(4,4)=+dtrmass*dcosh(y4)+trmass*dy4*dsinh(y4)
      dp(4,1)=-dpt*dsin(phi)
      dp(4,2)=-dpt*dcos(phi)
      dp(4,3)=+dtrmass*dsinh(y4)+trmass*dy4*dcosh(y4)

      if((mscase .eq. 'pt4') .or. (mscase .eq. 'm34')) then
         wt2=-4d0*wt0*(4d0*mass*(1d0-r(3))/s34*beta*w3
     &        +4d0*mass*sqrts**2*r(2)*beta/s34*ymax/dsqrt(1d0-r(2)**2)
     &        +8d0*mass*ymax*beta*r(2)/dsqrt(1d0-r(2)**2))/sqrts**2
      else if(mscase .eq. 'y4') then
         wt2=-16d0*wt0*beta*mass*((1d0-r(3))*w3
     &        +(sqrts**2+2d0*s34)*ymax)/s34/sqrts**2
      endif
         wt3=-4d0*wt0*(beta*w3/rtshat -
     &        4d0*mass**2/dsqrt((sqrts**2-4d0*mass**2)*r(3))/
     &        s34*ymax*w3)/sqrts**2

      return

      end
