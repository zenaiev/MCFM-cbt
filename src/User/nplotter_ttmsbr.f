      subroutine nplotter_ttmsbr(p,wt,wt2,switch,hq)
      implicit none
      include 'vegas_common.f'
      include 'bbproc.f'
      include 'clustering.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'process.f'
      include 'removebr.f'
      include 'nodecay.f'
      include 'useet.f'
      include 'plabel.f'
      include 'outputflags.f'
      include 'msbar.f'
cz Add single-top b fraction, Z. Sullivan 1/25/05
      double precision bwgt
      common/btagging/ bwgt

cz //

      integer n,switch,i5,i6,i7,nu,nplotmax,hq
      character tag*4
      double precision
     & ETARAP
     & ,PT
     & ,WT
     & ,WT2
     & ,YRAP
     & ,DPHI_LL
     & ,M_LL
     & ,MTRANS
     & ,SCUT1
     & ,SCUT2
     & ,yraptwo
      integer jet(mxpart),jetstart,ibbar,iz,izj,iztmp,
     . inotb,ilight1,ilight2
      integer iy
cz //

      double precision wtbbar,wtnotb,wtlight1,wtlight2
      double precision m34,etmiss,misset,m345678,
     . p(mxpart,4),fphi,HT,etcharm,deltaeta,cosdeltaphi
      double precision eta3,eta4,eta5,eta6,eta7,eta8,eta9,eta10,
     . eta34,eta56,y34
      double precision r45,r56,m345,m348,m567,m678
      double precision pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,
     . pt34,pt56,oldpt(5:7)
      double precision pt345,eta345,y345,y3,y4
      double precision bclustmass,etvec(4),tmp5(4),tmp6(4),tmp7(4)
      integer ssi(4),tmpi,j,k
      double precision ssd(4),tmpd
      integer nproc,eventpart,ib1,ib2,nqcdjets,nqcdstart
      logical first,jetmerge
      logical jetevent
      character*2 ptet
      character(len=1024) :: histname
      common/nplotmax/nplotmax
      common/nproc/nproc
      common/nqcdjets/nqcdjets,nqcdstart
      common/jetmerge/jetmerge
      common/hwwvars/dphi_ll,m_ll,mtrans,scut1,scut2
      data first/.true./
      save first
      

c--- Set up string for pt or Et
      if (useEt) then
        ptet='Et'
      else
        ptet='pt'
      endif
      
      if (first) then
        tag='book'
c--- ensure we initialize all possible histograms
        eventpart=npart+3
        eta4=1d3
        pt4=0d0
        goto 99
      else
        tag='plot'
      endif

      eventpart=npart-switch+2

      m34=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)

      y3=yrap(3,p)
      pt3=pt(3,p)        
      y4=yrap(4,p)
      pt4=pt(4,p)
      !do iy=1,4
      !  p34(1,iy)=p(3,iy)+p(4,iy)
      !enddo
      !y34=yrap(1,p34)
      y34=yraptwo(3,4,p)
      !print*,y34,m34
      !print*,p34
      !print*,p(3,1),p(3,2),p(3,3),p(3,4)
      !print*,p(4,1),p(4,2),p(4,3),p(4,4)
        
   99 continue

      n=nextnplot 
c hq=1 ttbar, hq=2 bbbar, hq=3 ccbar      
      if ( mscase .eq. 'y4' ) then
         if(hq.eq.1) then
          call bookplot(n,tag,'y4',y4,wt,wt2,-4d0,4d0,1.0d0,'lin')
         elseif(hq.eq.2) then
          call bookplot(n,tag,'y4',y4,wt,wt2,-8d0,8d0,1.0d0,'lin')
         elseif(hq.eq.3) then
          call bookplot(n,tag,'y4',y4,wt,wt2,-10d0,10d0,1.0d0,'lin')
         endif
         n=n+1
      elseif ( mscase .eq. 'pt4' ) then
         if(hq.eq.1) then
          call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,600d0,50d0,'lin')
         elseif(hq.eq.2) then
          call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,30d0,2.0d0,'lin')
         elseif(hq.eq.3) then
          call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,12d0,0.5d0,'lin')
         endif
         n=n+1
      elseif ( mscase .eq. 'm34' ) then
         if(hq.eq.1) then
          call bookplot(n,tag,'m34',m34,wt,wt2,0d0,1200d0,100d0,'lin')
         elseif(hq.eq.2) then
          call bookplot(n,tag,'m34',m34,wt,wt2,0d0,80d0,4.0d0,'lin')
         elseif(hq.eq.3) then
          !call bookplot(n,tag,'m34',m34,wt,wt2,0d0,15d0,0.5d0,'lin')
          call bookplot(n,tag,'m34',m34,wt,wt2,0d0,40d0,1.0d0,'lin')
         endif
         n=n+1 
      elseif ( mscase .eq. 'numeric' ) then
         if(hq.eq.1) then
          call bookplot(n,tag,'y4',y4,wt,wt2,-4d0,4d0,1.0d0,'lin')
          n=n+1 
          call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,600d0,50d0,'lin')
          n=n+1 
          call bookplot(n,tag,'m34',m34,wt,wt2,0d0,1200d0,100d0,'lin')
          n=n+1 
          ! double-differential pT-y
          do iy=1,6
            if((tag.eq.'book') .or. 
     . (y3.ge.((iy-1)*0.5d0).and.y3.le.(iy*0.5d0)))then
              write (histname, "(A7,I1)") "pTybin", iy
              call bookplot(n,tag,trim(histname),pt3,wt,wt2,
     .                      0.0d0,600.0d0,100.0d0,'lin')
            endif
            n=n+1
          enddo
         elseif(hq.eq.2) then
          call bookplot(n,tag,'y4',y4,wt,wt2,-8d0,8d0,1.0d0,'lin')
          n=n+1 
          call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,30d0,2.0d0,'lin')
          n=n+1 
          call bookplot(n,tag,'m34',m34,wt,wt2,0d0,80d0,4.0d0,'lin')
          n=n+1 
          ! double-differential pT-y
          do iy=1,5
            if((tag.eq.'book') .or. 
     . ((y3-2.0d0).ge.((iy-1)*0.5d0).and.(y3-2.0d0).le.(iy*0.5d0)))then
              write (histname, "(A7,I1)") "pTybin", iy
              call bookplot(n,tag,trim(histname),pt3,wt,wt2,
     .                      0.0d0,20.0d0,1.0d0,'lin')
            endif
            n=n+1
          enddo
         elseif(hq.eq.3) then
          call bookplot(n,tag,'y4',y4,wt,wt2,-10d0,10d0,1.0d0,'lin')
          n=n+1 
          call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,12d0,0.5d0,'lin')
          n=n+1 
          call bookplot(n,tag,'m34',m34,wt,wt2,0d0,15d0,0.5d0,'lin')
          n=n+1 
          ! double-differential pT-y
          do iy=1,10
            if((tag.eq.'book') .or. 
     . (y3.ge.((iy-1)*0.5d0).and.y3.le.(iy*0.5d0)))then
              write (histname, "(A7,I0.2)") "pTybin", iy
              call bookplot(n,tag,trim(histname),pt3,wt,wt2,
     .                      0.0d0,15.0d0,1.0d0,'lin')
            endif
            n=n+1
          enddo
          ! double-differential M-y
          do iy=1,10
            if((tag.eq.'book') .or. 
     . (y34.ge.((iy-1)*0.5d0).and.y34.le.(iy*0.5d0)))then
              write (histname, "(A6,I0.2)") "MYbin", iy
              call bookplot(n,tag,trim(histname),m34,wt,wt2,
     .                      0.0d0,30.0d0,2.0d0,'lin')
            endif
            n=n+1
          enddo
         endif
      else
        write(6,*) 'mscase not recognized'
        write(6,*) '   mscase = ',mscase
        stop         
      endif
            
      n=n-1
      
c--- ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
cz      
      bwgt=0d0  ! for safety
cz //

      return 
      end 
