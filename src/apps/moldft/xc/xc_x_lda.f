c:X_LDAsubrstart

c    Generated: Wed Jan 29 09:08:46 GMT 2003

      subroutine uks_x_lda
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     P.A.M. Dirac
c     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      vrhoa(i) = 0.D0
      vrhob(i) = -0.12407009817988D1*t1
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = -0.12407009817988D1*t4
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      vrhoa(i) = 0.D0
      vrhob(i) = -0.12407009817988D1*t1
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t5 = t1**2
      v2rhob2(i) = -0.4135669939329333D0/t5
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = 0.D0
      t5 = t1**2
      v2rhoa2(i) = -0.4135669939329333D0/t5
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = -0.12407009817988D1*t4
      t9 = t1**2
      v2rhoa2(i) = -0.4135669939329333D0/t9
      t12 = t4**2
      v2rhob2(i) = -0.4135669939329333D0/t12
      v2rhoab(i) = 0.D0
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_lda
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     P.A.M. Dirac
c     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      t5 = t1**2
      v2rhoa2(i) = -0.6564966812284644D0/t5
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:X_LDAsubrend
