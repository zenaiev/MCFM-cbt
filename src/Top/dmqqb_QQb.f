      subroutine dmqqb_QQb(r,p,dmp,dmsq,dmsqmtt,dmsqp4) 
      implicit none

************************************************************************
*     Author: Matthew Dowling                                          *
*     September 2013.                                                  *
*     calculate the derivative required                                *
*     for ttbar crsec in the msbar scheme                              *
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'msq_cs.f'
      include 'mxdim.f'
      include 'breit.f'
      include 'msbar.f'
      include 'process.f'

      integer j,k,cs
      logical first
      double precision dmp(mxpart,4),p(mxpart,4)
      double precision wtqqb,wtgg,t1,t2,ro
      double precision wtqqbmtt,wtggmtt,wtqqbp4,wtggp4
      double precision dt1,dt2,dro,ds(mxpart,mxpart)
      double precision dromtt,dt2mtt,dt1mtt,mtt
      double precision dt2pt,dt1pt,pt_ms
      double precision dt2p4,dt1p4
      double precision r(mxdim)
      double precision sqrts
      common/energy/sqrts
      data first/.true./
      save first

      if (case .ne. 'tt_tot' .and. case .ne. 'bb_tot' .and. 
     $    case .ne. 'cc_tot') then
        write(6,*) 'Only implemented for tt(bb,cc)bar cross section'
        call exit(-1)
      endif

      if (first) then
      first=.false.
      write(6,*) 'Heavy Quark mass:',mass2
      write(6,*) 'msbar:',msbar
      endif 

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      dmsq(j,k)=0d0
      do cs=0,2
      msq_cs(cs,j,k)=0d0
      enddo
      enddo
      enddo
      call dotem(4,p,s)

      ds(1,3)=2d0*(dmp(1,4)*p(3,4)-dmp(1,3)*p(3,3)-dmp(1,2)*p(3,2)
     . -dmp(1,1)*p(3,1) + p(1,4)*dmp(3,4)-p(1,3)*dmp(3,3)
     . -p(1,2)*dmp(3,2)-p(1,1)*dmp(3,1))

      ds(2,3)=2d0*(dmp(2,4)*p(3,4)-dmp(2,3)*p(3,3)-dmp(2,2)*p(3,2)
     . -dmp(2,1)*p(3,1) + p(2,4)*dmp(3,4)-p(2,3)*dmp(3,3)
     . -p(2,2)*dmp(3,2)-p(2,1)*dmp(3,1))

      ds(1,2)=2d0*(dmp(1,4)*p(2,4)-dmp(1,3)*p(2,3)-dmp(1,2)*p(2,2)
     . -dmp(1,1)*p(2,1) + p(1,4)*dmp(2,4)-p(1,3)*dmp(2,3)
     . -p(1,2)*dmp(2,2)-p(1,1)*dmp(2,1))

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)

      dt1=-ds(1,3)/s(1,2)-t1/s(1,2)*ds(1,2)
      dt2=-ds(2,3)/s(1,2)-t2/s(1,2)*ds(1,2)

      mtt=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
      pt_ms=(mtt**2-4d0*mass2**2)*r(2)*(1d0-r(2))
      
      dt1mtt=-4d0*t1/mtt/(4d0-mtt/mass2**2*(mtt-
     .     dsqrt((mtt**2-4d0*mass2**2)/(1d0-r(2)**2))))

      dt2mtt=-4d0*t2/mtt/(4d0-mtt/mass2**2*(mtt+
     .     dsqrt((mtt**2-4d0*mass2**2)/(1d0-r(2)**2))))

      dt1p4=r(2)/dsqrt((1-r(2)**2)*((sqrts**2-4d0*mass2**2)
     .        *r(3)+4d0*mass2**2))

      dt2p4=-dt1p4
      
c      dt1udif = -t1/(1d0+dsqrt(1d0-4d0*mass2**2/mtt**2)*(1d0-2d0*r(2)))
c      dt2udif = t2/(1d0-dsqrt(1d0-4d0*mass2**2/mtt**2)*(1d0-2d0*r(2)))

      ro=4d0*mass2**2/s(1,2)
      dro=ro*(2d0/mass2 - ds(1,2)/s(1,2))
      dromtt=-ro*2d0/mtt
      
      wtqqb=gsq**2*4d0/9d0*(2d0*t1*dt1 + 2d0*t2*dt2 + dro/2d0)
      wtgg=gsq**2*(-1d0/6d0/(t1*t2)*(dt1/t1+dt2/t2)
     .   *(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))
     .   +(1d0/6d0/t1/t2-3d0/8d0)
     .   *((2d0*t1+ro**2/(t1**2*t2*4d0))*dt1
     .   +(2d0*t2+ro**2/(t1*t2**2*4d0))*dt2
     .   +(dro-ro/2d0/(t1*t2)*dro)))

      wtqqbmtt=gsq**2*4d0/9d0*(2d0*t1*dt1mtt+2d0*t2*dt2mtt+dromtt/2d0)
      wtggmtt=gsq**2*(-1d0/6d0/(t1*t2)*(dt1mtt/t1+dt2mtt/t2)
     .   *(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))
     .   +(1d0/6d0/t1/t2-3d0/8d0)
     .   *((2d0*t1+ro**2/(t1**2*t2*4d0))*dt1mtt
     .   +(2d0*t2+ro**2/(t1*t2**2*4d0))*dt2mtt
     .   +(dromtt-ro/2d0/(t1*t2)*dromtt)))

      wtqqbp4=gsq**2*4d0/9d0*(2d0*t1*dt1p4+2d0*t2*dt2p4)
      wtggp4=gsq**2*(-1d0/6d0/(t1*t2)*(dt1p4/t1+dt2p4/t2)
     .   *(t1**2+t2**2+ro-0.25d0*ro**2/(t1*t2))
     .   +(1d0/6d0/t1/t2-3d0/8d0)
     .   *((2d0*t1+ro**2/(t1**2*t2*4d0))*dt1p4
     .   +(2d0*t2+ro**2/(t1*t2**2*4d0))*dt2p4))

      do j=-nf,nf
      k=-j
      if ((j .eq. 0) .and. (k .eq. 0)) then
         dmsq(j,k)=wtgg
         dmsqmtt(j,k)=wtggmtt
         dmsqp4(j,k)=wtggp4
      elseif ((j .gt. 0) .and. (k.lt.0)) then
         dmsq(j,k)=wtqqb
         dmsqmtt(j,k)=wtqqbmtt
         dmsqp4(j,k)=wtqqbp4
      elseif ((j .lt. 0) .and. (k.gt.0)) then
         dmsq(j,k)=wtqqb
         dmsqmtt(j,k)=wtqqbmtt
         dmsqp4(j,k)=wtqqbp4
      endif
      enddo
      return
      end
