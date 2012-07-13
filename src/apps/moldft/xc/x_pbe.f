c:X_PBEsubrstart
c    Generated: Wed Sep  3 12:48:32 GMT 2003
      subroutine uks_x_pbe
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
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
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t4 = rhob**2
      t5 = t2**2
      zk(i) = -0.9305257363491D0*t2*rhob*(0.1804D1-0.804D0/(1.D0
     #+0.449276922095889D-2*sigmabb/t5/t4))
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t4 = rhoa**2
      t5 = t2**2
      zk(i) = -0.9305257363491D0*t2*rhoa*(0.1804D1-0.804D0/(1.D0
     #+0.449276922095889D-2*sigmaaa/t5/t4))
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t6 = rhoa**2
      t7 = t4**2
      t18 = rhob**(1.D0/3.D0)
      t20 = rhob**2
      t21 = t18**2
      zk(i) = -0.9305257363491D0*t4*rhoa*(0.1804D1-0.804D0/(1.D0
     #+0.449276922095889D-2*sigmaaa/t7/t6))-0.9305257363491D0*t18*rhob
     #*(0.1804D1-0.804D0/(1.D0+0.449276922095889D-2*sigmabb/t21/t20))
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
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmabb/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      vrhoa(i) = 0.D0
      t20 = t10**2
      t21 = 1/t20
      vrhob(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmabb
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t3*t21
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmaaa/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmaaa
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.3361232459613792D-2/t3*t21
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t12 = 1.D0+0.449276922095889D-2*sigmaaa/t7/t6
      t15 = 0.1804D1-0.804D0/t12
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t20 = rhob**2
      t21 = t18**2
      t26 = 1.D0+0.449276922095889D-2*sigmabb/t21/t20
      t29 = 0.1804D1-0.804D0/t26
      zk(i) = -0.9305257363491D0*t5*t15-0.9305257363491D0*t19*t29
      t36 = t12**2
      t37 = 1/t36
      vrhoa(i) = -0.12407009817988D1*t4*t15+0.8963286558970112D-2/t4
     #/t6*t37*sigmaaa
      t45 = t26**2
      t46 = 1/t45
      vrhob(i) = -0.12407009817988D1*t18*t29+0.8963286558970112D-2/t18
     #/t20*t46*sigmabb
      vsigmaaa(i) = -0.3361232459613792D-2/t5*t37
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t19*t46
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
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
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmabb/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      vrhoa(i) = 0.D0
      t20 = t10**2
      t21 = 1/t20
      vrhob(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmabb
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t3*t21
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigmabb**2
      v2rhob2(i) = -0.4135669939329333D0/t5*t13-0.8963286558970112D-2
     #/t2/t4/rhob*t21*sigmabb+0.2147732158441357D-3/t37/t4*t41*t43
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.3020248347808158D-4/t37*t41
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmaaa/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmaaa
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.3361232459613792D-2/t3*t21
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigmaaa**2
      v2rhoa2(i) = -0.4135669939329333D0/t5*t13-0.8963286558970112D-2
     #/t2/t4/rhoa*t21*sigmaaa+0.2147732158441357D-3/t37/t4*t41*t43
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2sigmaaa2(i) = 0.3020248347808158D-4/t37*t41
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t12 = 1.D0+0.449276922095889D-2*sigmaaa/t7/t6
      t15 = 0.1804D1-0.804D0/t12
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t20 = rhob**2
      t21 = t18**2
      t26 = 1.D0+0.449276922095889D-2*sigmabb/t21/t20
      t29 = 0.1804D1-0.804D0/t26
      zk(i) = -0.9305257363491D0*t5*t15-0.9305257363491D0*t19*t29
      t36 = t12**2
      t37 = 1/t36
      t38 = 1/t4/t6*t37
      vrhoa(i) = -0.12407009817988D1*t4*t15+0.8963286558970112D-2*t38
     #*sigmaaa
      t45 = t26**2
      t46 = 1/t45
      t47 = 1/t18/t20*t46
      vrhob(i) = -0.12407009817988D1*t18*t29+0.8963286558970112D-2*t47
     #*sigmabb
      vsigmaaa(i) = -0.3361232459613792D-2/t5*t37
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t19*t46
      t65 = t6**2
      t69 = 1/t36/t12
      t71 = sigmaaa**2
      v2rhoa2(i) = -0.4135669939329333D0/t7*t15-0.8963286558970112D-2
     #/t4/t6/rhoa*t37*sigmaaa+0.2147732158441357D-3/t65/t6*t69*t71
      t83 = t20**2
      t87 = 1/t45/t26
      t89 = sigmabb**2
      v2rhob2(i) = -0.4135669939329333D0/t21*t29-0.8963286558970112D-2
     #/t18/t20/rhob*t46*sigmabb+0.2147732158441357D-3/t83/t20*t87*t89
      v2rhoab(i) = 0.D0
      v2rhoasigmaaa(i) = 0.4481643279485056D-2*t38
     #-0.8053995594155087D-4/t65/rhoa*t69*sigmaaa
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      v2rhobsigmabb(i) = 0.4481643279485056D-2*t47
     #-0.8053995594155087D-4/t83/rhob*t87*sigmabb
      v2sigmaaa2(i) = 0.3020248347808158D-4/t65*t69
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.3020248347808158D-4/t83*t87
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_pbe
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
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
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t4 = rho**2
      t5 = t2**2
      zk(i) = -0.7385587663820224D0*t2*rho*(0.1804D1-0.804D0/(1.D0
     #+0.7131826587600489D-2*sigma/t5/t4))
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t10 = 1.D0+0.7131826587600489D-2*sigma/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.7385587663820224D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.9847450218426965D0*t2*t13+0.1129303341188623D-1/t2
     #/t4*t21*sigma
      vsigmaaa(i) = -0.1693955011782934D-1/t3*t21
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t10 = 1.D0+0.7131826587600489D-2*sigma/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.7385587663820224D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      t22 = 1/t2/t4*t21
      vrhoa(i) = -0.9847450218426965D0*t2*t13+0.1129303341188623D-1
     #*t22*sigma
      vsigmaaa(i) = -0.1693955011782934D-1/t3*t21
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigma**2
      v2rhoa2(i) = -0.6564966812284644D0/t5*t13-0.2258606682377246D-1
     #/t2/t4/rho*t21*sigma+0.8590928633765426D-3/t37/t4*t41*t43
      v2rhoasigmaaa(i) = 0.2258606682377246D-1*t22
     #-0.644319647532407D-3/t37/rho*t41*sigma
      v2sigmaaa2(i) = 0.9664794712986105D-3/t37*t41
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
c:X_PBEsubrend
