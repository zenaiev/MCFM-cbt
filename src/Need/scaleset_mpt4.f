      subroutine scaleset_mpt4(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3 and 4
      implicit none
      include 'constants.f'
      include 'process.f'
      double precision p(mxpart,4),mu0
      integer i

      if((case .eq. 'W_only') .or.
     &   (case .eq. 'Z_only') .or.
     &   (case .eq. 'W_1jet') .or.
     &   (case .eq. 'W_2jet') .or.
     &   (case .eq. 'W_3jet') .or.
     &   (case .eq. 'Z_1jet') .or.
     &   (case .eq. 'Z_2jet') .or.
     &   (case .eq. 'Z_3jet') .or.
     &   (case .eq. 'Wbbbar') .or.
     &   (case .eq. 'Wbbmas') .or.
     &   (case .eq. 'Zbbbar') .or.
     &   (case .eq. 'gamgam') .or.
     &   (case .eq. 'ggfus0') .or.
     &   (case .eq. 'ggfus1') .or.
     &   (case .eq. 'ggfus2') .or.
     &   (case .eq. 'ggfus3') .or.
     &   (case .eq. 'gagajj') .or.
     &   (case .eq. 'qg_tbq') .or.
     &   (case .eq. 'tt_tot') .or.
     &   (case .eq. 'bb_tot') .or.
     &   (case .eq. 'cc_tot') .or.
     &   (case .eq. 'httjet') .or.
     &   (case .eq. 'Higaga') .or.
     &   (case .eq. 'Hgagaj') .or.
     &   (case .eq. 'qq_Hqq') .or.
     &   (case .eq. 'dm_jet') .or.
     &   (case .eq. 'dm_gam') .or.
     &   (case .eq. 'dm2jet') .or.
     &   (case .eq. 'dm_gaj') .or.
     &   (case .eq. 'dmgamj') .or.
     &   (case .eq. 'qq_Hgg')) then
!        mu0=(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
!     &     -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2       
        mu0=(4.0*(p(4,4)**2-p(4,3)**2-p(4,2)**2-p(4,1)**2)
     &      +1.0*(p(4,2)**2+p(4,1)**2))
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*) 'dynamicscale m(34) not supported for this process.'
        stop
      endif
!      print*,'p3 = ',p(3,1),p(3,2),p(3,3),p(3,4),' mu = ',mu0
!      print*,'p4 = ',p(4,1),p(4,2),p(4,3),p(4,4),' mu = ',mu0
!      print*,'m34 = ',dsqrt((p(3,4)+p(4,4))**2-(p(3,3)+p(4,3))**2-
!     & (p(3,2)+p(4,2))**2-(p(3,1)+p(4,1))**2)
      
      return
      end
      
