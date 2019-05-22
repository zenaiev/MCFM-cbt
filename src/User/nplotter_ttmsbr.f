      subroutine nplotter_ttmsbr(p,wt,wt2,switch)
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

      integer n,switch,i5,i6,i7,nu,nplotmax
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
      integer jet(mxpart),jetstart,ibbar,iz,izj,iztmp,
     . inotb,ilight1,ilight2
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

      y4=yrap(4,p)
      pt4=pt(4,p)        
        
   99 continue

      n=nextnplot 
      if ( mscase .eq. 'y4' ) then
         call bookplot(n,tag,'y4',y4,wt,wt2,-4d0,4d0,0.2d0,'lin')
         n=n+1
      elseif ( mscase .eq. 'pt4' ) then
         call bookplot(n,tag,ptet//'4',pt4,wt,wt2,25d0,50d0,0.5d0,'lin')
         n=n+1
         call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,150d0,5d0,'log')
         n=n+1
         call bookplot(n,tag,ptet//'4',pt4,wt,wt2,0d0,600d0,20d0,'log')
         n=n+1
      elseif ( mscase .eq. 'm34' ) then
         call bookplot(n,tag,'m34',m34,wt,wt2,0d0,200d0,5d0,'lin')
         n=n+1 
         call bookplot(n,tag,'m34',m34,wt,wt2,0d0,500d0,10d0,'log')
         n=n+1
      elseif ( mscase .eq. 'numeric' ) then
         stop
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
