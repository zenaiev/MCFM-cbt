      double precision function msbint(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
      include 'stopscales.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'dm_params.f'
      include 'qcdcouple.f'
      include 'msbar.f'
c --- DSW. To store flavour information :
      include 'nflav.f'
c      include 'b0.f'
c --- DSW.
c --- to access HQ masses for numerical derivative calculation
      include 'breit.f'
      integer pflav,pbarflav
c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk,ii
      double precision dfx1(-nf:nf),dfx2(-nf:nf),dmp(mxpart,4),
     . dmpswt,msb_surf_max,msb_surf_min
      double precision r(mxdim),W,sqrts,xmsq,val,val2,ptmp,
     . fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart,
     . fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     . fxb1(-nf:nf),fxb2(-nf:nf)
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      double precision xx(2),flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat
      double precision xmsq_bypart(-1:1,-1:1),lord_bypart(-1:1,-1:1)
      logical bin,first,includedipole
c      double precision gx1(-nf:nf),gx2(-nf:nf)
c      integer idum
c      COMMON/ranno/idum
      character*30 runstring
      double precision b1scale,q2scale,q1scale,b2scale
      external qg_tbq,BSYqqb_QQbdk_gvec,qqb_QQbdk,qg_tbqdk,qg_tbqdk_gvec
      common/runstring/runstring
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/bypart/lord_bypart
      common/bqscale/b1scale,q2scale,q1scale,b2scale
      data p/56*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      external qq_tchan_ztq,qq_tchan_ztq_mad
      external qq_tchan_htq,qq_tchan_htq_mad,qq_tchan_htq_amp
      double precision masssave,dm,lord,dlord,deriv,l
      double precision lowint
      double precision coef
      logical dynamicscalesave

      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      msbint=0d0
      
      if (mscase .eq. 'numeric') then
        W=sqrts**2
        p(:,:)=0d0
        nvec=npart+2
        call gen2m(r,p,pswt,*999)
        ! need scale to calculate alpha_s
        ! also need to prevent varying mass in dynamic scale
        if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
        dynamicscalesave=dynamicscale
        dynamicscale=.false.
        !scale=rscalestart
        !facscale=fscalestart
        masssave=mass2
        dm=mass2*1.0d-4
        l=4d0/3d0
        !l=4d0/3d0+dlog(scale**2/mass2**2)
        coef=-1.0d0*as/pi*l*masssave/dm
        lord=lowint(r,coef*wgt)
        mass2=mass2+dm
        coef=as/pi*l*masssave/dm
        dlord=lowint(r,coef*wgt)
        deriv=(dlord-lord)/dm
        mass2=masssave
        msbint=as/pi*l*mass2*deriv
        dynamicscale=dynamicscalesave
        return
      endif

c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0d0
      W=sqrts**2
      p(:,:)=0d0

      if ( ((case .eq. 'tt_tot') .or. (case .eq. 'bb_tot') .or. 
     $ (case .eq. 'cc_tot')) .and. (msbar)) then
         npart=2
         call gen2m_ms(r,p,pswt,*999)
         call dgen2m(r,dmp,dmpswt,dmpswtmtt,*999)
      endif

      nvec=npart+2
      call dotem(nvec,p,s)

c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
         
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)


      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts
      
      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      
     
      call qqb_QQb(p,msq)
      call dmqqb_QQb(r,p,dmp,dmsq,dmsqmtt,dmsqp4)

      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0d0
      enddo
      enddo 

      xmsq=0d0
      dxmsq=0d0
      dxmsqmtt=0d0
      dxmsqp4=0d0
      dxmsqy4=0d0

      currentPDF=0

      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
cmd- for tt_tot in msbar scheme
      mt=mass2
      dflux = -flux*8d0*mt*(1d0-r(3))/((W-4d0*mt**2)*r(3)+4d0*mt**2)
      dfluxmtt = -2d0*flux/dsqrt(4d0*mt**2+(W-4d0*mt**2)*r(3))

 777  continue
      xmsq=0d0
      if (PDFerrors) then
         call InitPDF(currentPDF)
      endif

      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)
      call dfdist(ih1,xx(1),facscale,dfx1)
      call dfdist(ih2,xx(2),facscale,dfx2)

      
      do j=-nflav,nflav
      do k=-nflav,nflav    

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (omitgg) then 
      if ((j.eq.0) .and. (k.eq.0)) goto 20
      endif

cmd- calculate deriavtive of hadronic cross-section
c         xmsq = xmsq+msq(j,k)
         xmsqjk = fx1(j)*fx2(k)*msq(j,k)
c         dxmsq = dxmsq+dmsq(j,k)
         dxmsqjk = dfx1(j)*fx2(k)*msq(j,k)
     .   *(xx(1)*8d0*mt*(1d0-r(1))*(1d0-r(3))
     .   /((W-4d0*mt**2)*r(3)+4d0*mt**2))
     .   + fx1(j)*dfx2(k)*msq(j,k)
     .   *(xx(2)*8d0*mt*r(1)*(1d0-r(3))
     .   /((W-4d0*mt**2)*r(3)+4d0*mt**2))
     .   + fx1(j)*fx2(k)*dmsq(j,k)

c            dxmsqmtt=dxmsqmtt+dmsqmtt(j,k)
            dxmsqmttjk=
     .   dfx1(j)*fx2(k)*msq(j,k)*(2d0*(1d0-r(1))*xx(1)/
     .     dsqrt(4d0*mt**2+(W-4d0*mt**2)*r(3)))
     .+  fx1(j)*dfx2(k)*msq(j,k)*2d0*r(1)*xx(2)/
     .     dsqrt(4d0*mt**2+(W-4d0*mt**2)*r(3))
     .+  fx1(j)*fx2(k)*dmsqmtt(j,k)

c            dxmsqp4=dxmsqp4+dmsqp4(j,k)
         dxmsqp4jk=
     .   fx1(j)*fx2(k)*dmsqp4(j,k)

c         dxmsqy4=0d0
         dxmsqy4jk=
     .        dfx1(j)*fx2(k)*xx(1)*msq(j,k)
     .       -fx1(j)*dfx2(k)*xx(2)*msq(j,k)

         xmsq=xmsq+xmsqjk
         dxmsq=dxmsq+dxmsqjk
         dxmsqmtt=dxmsqmtt+dxmsqmttjk
         dxmsqp4=dxmsqp4+dxmsqp4jk
         dxmsqy4=dxmsqy4+dxmsqy4jk

      if     (j .gt. 0) then
        sgnj=+1
      elseif (j .lt. 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k .gt. 0) then
        sgnk=+1
      elseif (k .lt. 0) then
        sgnk=-1
      else
        sgnk=0
      endif

cmd   NOTE: This is not properly implemented yet.
cmd   It does not include the derivative term so I multiply by 0d0
      if (currentPDF .eq. 0) then
        xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+xmsqjk*0d0
      endif
      
 20   continue
      enddo
      enddo

      if (mscase .eq. 'm34') then
cmd-derivative of jacobian w.r.t. mtt
         msjacob = flux*pswt*xmsq*(8d0*mt/(W-4d0*mt**2))
cmd-partial derivative of lomsq w.r.t. mtt
     .    - (dfluxmtt*pswt*xmsq
     .    +  flux*dmpswtmtt*xmsq
     .    +  flux*pswt*dxmsqmtt)*
     .       4d0*mt*(1d0-r(3))/dsqrt(4d0*mt**2+(W-4d0*mt**2)*r(3))
      elseif (mscase .eq. 'pt4') then
cmd-derivative of jacobian w.r.t. pt4
         msjacob = flux*pswt*xmsq*(4d0*mt/(W-4d0*mt**2))
cmd-partial derivative of xmsq w.r.t. pt4
     .    - (flux*pswt*dxmsqp4)*
     .       (2d0*mt*r(2)*r(3)/dsqrt((W-4d0*mt**2)*r(3)))
      elseif (mscase .eq. 'y4') then
         msjacob =
     .        -(flux*pswt*dxmsqy4)*(
     .        (4d0*mt*(1d0-2d0*r(1))*(1d0-r(3))/((W-4d0*mt**2)*r(3)
     .        +4d0*mt**2))-mt*W*(1d0-2d0*r(2))*dsqrt((W-4d0*mt**2)*r(3)/
     .        ((W-4d0*mt**2)*r(3)+4d0*mt**2))/((W-4d0*mt**2)*(
     .        (W-4d0*mt**2)*(1d0-r(2))*r(2)*r(3)+mt**2)))
cmd-derivative of jacobian w.r.t. y4
     .       + flux*pswt*xmsq*(
     .        8d0*mt*(1d0-r(3))/(((W-4d0*mt**2)*r(3)+4d0*mt**2)*
     .        dlog(W/((W-4d0*mt**2)*r(3)+4d0*mt**2))))
      endif

      if (currentPDF .eq. 0) then
        msbint=
!     .      as/pi*(4d0/3d0+dlog(scale**2/mt**2))*mt*
     .      as/pi*(4d0/3d0)*mt*
     .      ((dflux*pswt*xmsq
     .    + flux*dmpswt*xmsq
     .    + flux*pswt*dxmsq)
     .     + msjacob
     .    )/BrnRat
      endif
c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        PDFwgt(currentPDF)=
!     .      as/pi*(4d0/3d0+dlog(scale**2/mt**2))*mass2*
     .      as/pi*(4d0/3d0)*mass2*
     .      ((dflux*pswt*xmsq
     .    + flux*dmpswt*xmsq
     .    + flux*pswt*dxmsq)
     .     + msjacob
     .    )/BrnRat*wgt/itmx

        PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .     +PDFwgt(currentPDF)
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
      endif    

cmd   This is also not implemented properly. 
      if (creatent) then
         
        wt_gg=xmsq_bypart(0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     .        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
      endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=lord_bypart(j,k)+
     .       wgt*flux*pswt*xmsq_bypart(j,k)/BrnRat
      enddo
      enddo

      val=msbint*wgt
      val2=val**2

c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
        wtmax=dabs(val)
      endif

c      if(rescale) then 
c         call rescale_pjet(pjet)
c      endif

      if (bin) then
c ---   DSW. If the user has not selected to generate
c ---   events, still call nplotter here in order to
c ---   fill histograms/ntuples with weighted events :
        if (.not.evtgen) then
          call nplotter(pjet,val,val2,0)
c--- POWHEG-style output if requested
        if (writepwg) then
            call pwhgplotter(p,pjet,val,0)
        endif
        endif
      endif
c --- Check weights :
      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
        wtabs = dabs(val)
        if (ran2() .lt. (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
          if (wtabs.lt.wtmax) then
            newwt = 1d0
          else
            newwt = wtabs/wtmax
          endif
          if (newwt .gt. 1.0d0) then
            write(6,*) 'WARNING : lowint : event with |weight| > 1.',
     +            ' |weight| = ',newwt
          endif
c ---     just in case the weight was negative :
          newwt = newwt*dsign(1d0,val)
c          call nplotter(pjet,newwt,newwt,0)
c ---     DSW. If I'm storing the event, I need to make a decision
c ---     about the flavours :
          call decide_flavour(pflav,pbarflav)
          call storeevent(pjet,newwt,pflav,pbarflav)
        endif
      endif

      return
      
 999  continue
      ntotzero=ntotzero+1
      
      return
      end
