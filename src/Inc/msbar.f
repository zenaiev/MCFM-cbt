      logical msbar
      double precision dpt_ms,dflux,dpswt,dxmsq,dxx(2)      
      double precision dshat,dmsq(-nf:nf,-nf:nf),dfluxmtt
      double precision wt3,dmpswtmtt,dmsqmtt(-nf:nf,-nf:nf)
      double precision dxmsqmtt,dxmsqpt,msjacob
      double precision dxmsqjk,dxmsqmttjk,dxmsqp4jk,dxmsqy4jk
      double precision dmsqpt(-nf:nf,-nf:nf)
      double precision dmsqp4(-nf:nf,-nf:nf),dmpswtp4
      double precision dxmsqp4,dxmsqy4
      character*16 mscase
      common/mscase/mscase
      common/msbar/msbar
      common/dmsbar/dpt_ms,dflux,dpswt,dxmsq,dxx,
     . dshat,dfluxmtt,wt3,dxmsqmtt,dxmsqp4,msjacob,dxmsqy4
c      common/dmsq/dmsq

