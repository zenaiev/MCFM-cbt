*****************
* LHAPDF version*
*****************
      subroutine dfdist(ih,x,xmu,fx)
      implicit none
      double precision fx(-5:5),x,xmu,fPDFm(-6:6),fPDFp(-6:6)
      double precision xp,xm,del
      integer Iprtn,ih,Irt
c---  ih1=+1 proton 
c---  ih1=-1 pbar 
      
      del = 0.001d0*x
      xp = x+del
      xm = x-del
C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
 
      call evolvePDF(xp,xmu,fPDFp)
      call evolvePDF(xm,xmu,fPDFm)

      if (ih.eq.1) then
        do Iprtn=-5,5
          fx(+Iprtn)=(fPDFp(+Iprtn)/xp - fPDFm(+Iprtn)/xm)/(2d0*del)
        enddo
      elseif(ih.eq.-1) then
        do Iprtn=-5,5
          fx(+Iprtn)=(fPDFp(-Iprtn)/xp - fPDFm(-Iprtn)/xm)/(2d0*del)
        enddo
      endif
                     
      return
      end

  

