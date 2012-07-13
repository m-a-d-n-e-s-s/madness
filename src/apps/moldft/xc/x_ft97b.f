c:X_FT97Bsubrstart
c    Generated: Thu Jan 30 10:36:52 GMT 2003
      subroutine uks_x_ft97b
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
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
      t8 = 0.2913644D-2+0.9474169D-3*sigmabb/(0.6255746320201D7+sigmabb)
      t10 = rhob**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t17 = dlog(sigmabb*t13+dsqrt(1+sigmabb**2*t13**2))
      t18 = t17**2
      t23 = dsqrt(1.D0+9.D0*t14*sigmabb*t13*t18)
      zk(i) = -0.9305257363491D0*t2*rhob*(1.D0+0.1074661302677646D1*t8
     #*sigmabb*t13/t23)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t8 = 0.2913644D-2+0.9474169D-3*sigmaaa/(0.6255746320201D7+sigmaaa)
      t10 = rhoa**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t17 = dlog(sigmaaa*t13+dsqrt(1+sigmaaa**2*t13**2))
      t18 = t17**2
      t23 = dsqrt(1.D0+9.D0*t14*sigmaaa*t13*t18)
      zk(i) = -0.9305257363491D0*t2*rhoa*(1.D0+0.1074661302677646D1*t8
     #*sigmaaa*t13/t23)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t10 = 0.2913644D-2+0.9474169D-3*sigmaaa/(0.6255746320201D7
     #+sigmaaa)
      t12 = rhoa**2
      t13 = t4**2
      t15 = 1/t13/t12
      t16 = t10**2
      t19 = dlog(sigmaaa*t15+dsqrt(1+sigmaaa**2*t15**2))
      t20 = t19**2
      t25 = dsqrt(1.D0+9.D0*t16*sigmaaa*t15*t20)
      t33 = rhob**(1.D0/3.D0)
      t39 = 0.2913644D-2+0.9474169D-3*sigmabb/(0.6255746320201D7
     #+sigmabb)
      t41 = rhob**2
      t42 = t33**2
      t44 = 1/t42/t41
      t45 = t39**2
      t48 = dlog(sigmabb*t44+dsqrt(1+sigmabb**2*t44**2))
      t49 = t48**2
      t54 = dsqrt(1.D0+9.D0*t45*sigmabb*t44*t49)
      zk(i) = -0.9305257363491D0*t4*rhoa*(1.D0+0.1074661302677646D1
     #*t10*sigmaaa*t15/t25)-0.9305257363491D0*t33*rhob*(1.D0
     #+0.1074661302677646D1*t39*sigmabb*t44/t54)
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
      t4 = 0.6255746320201D7+sigmabb
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmabb*t5
      t9 = t8*sigmabb
      t10 = rhob**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmabb
      t17 = dlog(sigmabb*t13+dsqrt(1+sigmabb**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      vrhoa(i) = 0.D0
      t35 = 1/t11/t10/rhob
      t41 = t13/t23/t22
      t45 = sigmabb**2
      t47 = t10**2
      t54 = 1/t2/t47/rhob
      t57 = dsqrt(1.D0+t45*t54)
      t58 = 1/t57
      vrhob(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*(
     #-0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9*t41*(
     #-24.D0*t15*t35*t18-48.D0*t14*t45/t2/t47/t10*t17*t58))
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t70 = t4**2
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmabb/t70
      vsigmabb(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t74
     #*sigmabb*t25+0.1074661302677646D1*t8*t13*t24
     #-0.5373306513388232D0*t9*t41*(18.D0*t9*t19*t74+9.D0*t14*t13*t18
     #+18.D0*t15*t54*t17*t58))
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = 0.6255746320201D7+sigmaaa
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmaaa*t5
      t9 = t8*sigmaaa
      t10 = rhoa**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmaaa
      t17 = dlog(sigmaaa*t13+dsqrt(1+sigmaaa**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      t35 = 1/t11/t10/rhoa
      t41 = t13/t23/t22
      t45 = sigmaaa**2
      t47 = t10**2
      t54 = 1/t2/t47/rhoa
      t57 = dsqrt(1.D0+t45*t54)
      t58 = 1/t57
      vrhoa(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*(
     #-0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9*t41*(
     #-24.D0*t15*t35*t18-48.D0*t14*t45/t2/t47/t10*t17*t58))
      vrhob(i) = 0.D0
      t70 = t4**2
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmaaa/t70
      vsigmaaa(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t74
     #*sigmaaa*t25+0.1074661302677646D1*t8*t13*t24
     #-0.5373306513388232D0*t9*t41*(18.D0*t9*t19*t74+9.D0*t14*t13*t18
     #+18.D0*t15*t54*t17*t58))
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = 0.6255746320201D7+sigmaaa
      t7 = 1/t6
      t10 = 0.2913644D-2+0.9474169D-3*sigmaaa*t7
      t11 = t10*sigmaaa
      t12 = rhoa**2
      t13 = t4**2
      t15 = 1/t13/t12
      t16 = t10**2
      t17 = t16*sigmaaa
      t19 = dlog(sigmaaa*t15+dsqrt(1+sigmaaa**2*t15**2))
      t20 = t19**2
      t21 = t15*t20
      t24 = 1.D0+9.D0*t17*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t15*t26
      t30 = 1.D0+0.1074661302677646D1*t11*t27
      t33 = rhob**(1.D0/3.D0)
      t34 = t33*rhob
      t35 = 0.6255746320201D7+sigmabb
      t36 = 1/t35
      t39 = 0.2913644D-2+0.9474169D-3*sigmabb*t36
      t40 = t39*sigmabb
      t41 = rhob**2
      t42 = t33**2
      t44 = 1/t42/t41
      t45 = t39**2
      t46 = t45*sigmabb
      t48 = dlog(sigmabb*t44+dsqrt(1+sigmabb**2*t44**2))
      t49 = t48**2
      t50 = t44*t49
      t53 = 1.D0+9.D0*t46*t50
      t54 = dsqrt(t53)
      t55 = 1/t54
      t56 = t44*t55
      t59 = 1.D0+0.1074661302677646D1*t40*t56
      zk(i) = -0.9305257363491D0*t5*t30-0.9305257363491D0*t34*t59
      t66 = 1/t13/t12/rhoa
      t72 = t15/t25/t24
      t76 = sigmaaa**2
      t78 = t12**2
      t85 = 1/t4/t78/rhoa
      t88 = dsqrt(1.D0+t76*t85)
      t89 = 1/t88
      vrhoa(i) = -0.12407009817988D1*t4*t30-0.9305257363491D0*t5*(
     #-0.2865763473807057D1*t11*t66*t26-0.5373306513388232D0*t11*t72*(
     #-24.D0*t17*t66*t20-48.D0*t16*t76/t4/t78/t12*t19*t89))
      t104 = 1/t42/t41/rhob
      t110 = t44/t54/t53
      t114 = sigmabb**2
      t116 = t41**2
      t123 = 1/t33/t116/rhob
      t126 = dsqrt(1.D0+t114*t123)
      t127 = 1/t126
      vrhob(i) = -0.12407009817988D1*t33*t59-0.9305257363491D0*t34*(
     #-0.2865763473807057D1*t40*t104*t55-0.5373306513388232D0*t40*t110
     #*(-24.D0*t46*t104*t49-48.D0*t45*t114/t33/t116/t41*t48*t127))
      t139 = t6**2
      t143 = 0.9474169D-3*t7-0.9474169D-3*sigmaaa/t139
      vsigmaaa(i) = -0.9305257363491D0*t5*(0.1074661302677646D1*t143
     #*sigmaaa*t27+0.1074661302677646D1*t10*t15*t26
     #-0.5373306513388232D0*t11*t72*(18.D0*t11*t21*t143+9.D0*t16*t15
     #*t20+18.D0*t17*t85*t19*t89))
      vsigmaab(i) = 0.D0
      t168 = t35**2
      t172 = 0.9474169D-3*t36-0.9474169D-3*sigmabb/t168
      vsigmabb(i) = -0.9305257363491D0*t34*(0.1074661302677646D1*t172
     #*sigmabb*t56+0.1074661302677646D1*t39*t44*t55
     #-0.5373306513388232D0*t40*t110*(18.D0*t40*t50*t172+9.D0*t45*t44
     #*t49+18.D0*t46*t123*t48*t127))
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
      t4 = 0.6255746320201D7+sigmabb
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmabb*t5
      t9 = t8*sigmabb
      t10 = rhob**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmabb
      t17 = dlog(sigmabb*t13+dsqrt(1+sigmabb**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      vrhoa(i) = 0.D0
      t33 = t10*rhob
      t35 = 1/t11/t33
      t40 = 1/t23/t22
      t41 = t13*t40
      t45 = sigmabb**2
      t46 = t14*t45
      t47 = t10**2
      t54 = 1/t2/t47/rhob
      t56 = 1.D0+t45*t54
      t57 = dsqrt(t56)
      t58 = 1/t57
      t62 = -24.D0*t15*t35*t18-48.D0*t46/t2/t47/t10*t17*t58
      t66 = -0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9
     #*t41*t62
      vrhob(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*t66
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t70 = t4**2
      t71 = 1/t70
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmabb*t71
      t75 = t74*sigmabb
      t78 = t8*t13
      t91 = 18.D0*t9*t19*t74+9.D0*t14*t13*t18+18.D0*t15*t54*t17*t58
      t92 = t41*t91
      vsigmabb(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t75
     #*t25+0.1074661302677646D1*t78*t24-0.5373306513388232D0*t9*t92)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t104 = 1/t11/t47
      t112 = t22**2
      t115 = t13/t23/t112
      t116 = t62**2
      t132 = t47**2
      t133 = t132*t10
      t135 = 1/t56
      t139 = t45**2
      t146 = 1/t57/t56
      v2rhob2(i) = -0.4135669939329333D0/t11*t28-0.24814019635976D1*t2
     #*t66-0.9305257363491D0*t3*(0.1050779940395921D2*t9*t104*t24
     #+0.2865763473807057D1*t9*t35*t40*t62+0.8059959770082348D0*t9
     #*t115*t116-0.5373306513388232D0*t9*t41*(88.D0*t15*t104*t18
     #+432.D0*t46/t2/t47/t33*t17*t58+128.D0*t14*t45*sigmabb/t133*t135
     #-128.D0*t14*t139/t11/t132/t47*t17*t146))
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t162 = -0.18948338D-2*t71+0.18948338D-2*sigmabb/t70/t4
      t174 = t91**2
      t178 = t74**2
      v2sigmabb2(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t162
     #*sigmabb*t25+0.2149322605355293D1*t74*t13*t24
     #-0.1074661302677646D1*t75*t92-0.1074661302677646D1*t78*t40*t91
     #+0.8059959770082348D0*t9*t115*t174-0.5373306513388232D0*t9*t41*
     #(18.D0*t178*sigmabb*t19+36.D0*t78*t18*t74+72.D0*t9*t54*t17*t74
     #*t58+18.D0*t9*t19*t162+36.D0*t14*t54*t17*t58+18.D0*t15/t132*t135
     #-18.D0*t46/t11/t133*t17*t146))
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = 0.6255746320201D7+sigmaaa
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmaaa*t5
      t9 = t8*sigmaaa
      t10 = rhoa**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmaaa
      t17 = dlog(sigmaaa*t13+dsqrt(1+sigmaaa**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      t33 = t10*rhoa
      t35 = 1/t11/t33
      t40 = 1/t23/t22
      t41 = t13*t40
      t45 = sigmaaa**2
      t46 = t14*t45
      t47 = t10**2
      t54 = 1/t2/t47/rhoa
      t56 = 1.D0+t45*t54
      t57 = dsqrt(t56)
      t58 = 1/t57
      t62 = -24.D0*t15*t35*t18-48.D0*t46/t2/t47/t10*t17*t58
      t66 = -0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9
     #*t41*t62
      vrhoa(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*t66
      vrhob(i) = 0.D0
      t70 = t4**2
      t71 = 1/t70
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmaaa*t71
      t75 = t74*sigmaaa
      t78 = t8*t13
      t91 = 18.D0*t9*t19*t74+9.D0*t14*t13*t18+18.D0*t15*t54*t17*t58
      t92 = t41*t91
      vsigmaaa(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t75
     #*t25+0.1074661302677646D1*t78*t24-0.5373306513388232D0*t9*t92)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t104 = 1/t11/t47
      t112 = t22**2
      t115 = t13/t23/t112
      t116 = t62**2
      t132 = t47**2
      t133 = t132*t10
      t135 = 1/t56
      t139 = t45**2
      t146 = 1/t57/t56
      v2rhoa2(i) = -0.4135669939329333D0/t11*t28-0.24814019635976D1*t2
     #*t66-0.9305257363491D0*t3*(0.1050779940395921D2*t9*t104*t24
     #+0.2865763473807057D1*t9*t35*t40*t62+0.8059959770082348D0*t9
     #*t115*t116-0.5373306513388232D0*t9*t41*(88.D0*t15*t104*t18
     #+432.D0*t46/t2/t47/t33*t17*t58+128.D0*t14*t45*sigmaaa/t133*t135
     #-128.D0*t14*t139/t11/t132/t47*t17*t146))
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t162 = -0.18948338D-2*t71+0.18948338D-2*sigmaaa/t70/t4
      t174 = t91**2
      t178 = t74**2
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t162
     #*sigmaaa*t25+0.2149322605355293D1*t74*t13*t24
     #-0.1074661302677646D1*t75*t92-0.1074661302677646D1*t78*t40*t91
     #+0.8059959770082348D0*t9*t115*t174-0.5373306513388232D0*t9*t41*
     #(18.D0*t178*sigmaaa*t19+36.D0*t78*t18*t74+72.D0*t9*t54*t17*t74
     #*t58+18.D0*t9*t19*t162+36.D0*t14*t54*t17*t58+18.D0*t15/t132*t135
     #-18.D0*t46/t11/t133*t17*t146))
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
      t6 = 0.6255746320201D7+sigmaaa
      t7 = 1/t6
      t10 = 0.2913644D-2+0.9474169D-3*sigmaaa*t7
      t11 = t10*sigmaaa
      t12 = rhoa**2
      t13 = t4**2
      t15 = 1/t13/t12
      t16 = t10**2
      t17 = t16*sigmaaa
      t19 = dlog(sigmaaa*t15+dsqrt(1+sigmaaa**2*t15**2))
      t20 = t19**2
      t21 = t15*t20
      t24 = 1.D0+9.D0*t17*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t15*t26
      t30 = 1.D0+0.1074661302677646D1*t11*t27
      t33 = rhob**(1.D0/3.D0)
      t34 = t33*rhob
      t35 = 0.6255746320201D7+sigmabb
      t36 = 1/t35
      t39 = 0.2913644D-2+0.9474169D-3*sigmabb*t36
      t40 = t39*sigmabb
      t41 = rhob**2
      t42 = t33**2
      t44 = 1/t42/t41
      t45 = t39**2
      t46 = t45*sigmabb
      t48 = dlog(sigmabb*t44+dsqrt(1+sigmabb**2*t44**2))
      t49 = t48**2
      t50 = t44*t49
      t53 = 1.D0+9.D0*t46*t50
      t54 = dsqrt(t53)
      t55 = 1/t54
      t56 = t44*t55
      t59 = 1.D0+0.1074661302677646D1*t40*t56
      zk(i) = -0.9305257363491D0*t5*t30-0.9305257363491D0*t34*t59
      t64 = t12*rhoa
      t66 = 1/t13/t64
      t67 = t66*t26
      t71 = 1/t25/t24
      t72 = t15*t71
      t73 = t66*t20
      t76 = sigmaaa**2
      t77 = t16*t76
      t78 = t12**2
      t81 = 1/t4/t78/t12
      t85 = 1/t4/t78/rhoa
      t87 = 1.D0+t76*t85
      t88 = dsqrt(t87)
      t89 = 1/t88
      t90 = t81*t19*t89
      t93 = -24.D0*t17*t73-48.D0*t77*t90
      t94 = t72*t93
      t97 = -0.2865763473807057D1*t11*t67-0.5373306513388232D0*t11*t94
      vrhoa(i) = -0.12407009817988D1*t4*t30-0.9305257363491D0*t5*t97
      t102 = t41*rhob
      t104 = 1/t42/t102
      t105 = t104*t55
      t109 = 1/t54/t53
      t110 = t44*t109
      t111 = t104*t49
      t114 = sigmabb**2
      t115 = t45*t114
      t116 = t41**2
      t119 = 1/t33/t116/t41
      t123 = 1/t33/t116/rhob
      t125 = 1.D0+t114*t123
      t126 = dsqrt(t125)
      t127 = 1/t126
      t128 = t119*t48*t127
      t131 = -24.D0*t46*t111-48.D0*t115*t128
      t132 = t110*t131
      t135 = -0.2865763473807057D1*t40*t105-0.5373306513388232D0*t40
     #*t132
      vrhob(i) = -0.12407009817988D1*t33*t59-0.9305257363491D0*t34*t135
      t139 = t6**2
      t140 = 1/t139
      t143 = 0.9474169D-3*t7-0.9474169D-3*sigmaaa*t140
      t144 = t143*sigmaaa
      t147 = t10*t15
      t160 = 18.D0*t11*t21*t143+9.D0*t16*t15*t20+18.D0*t17*t85*t19*t89
      t161 = t72*t160
      t164 = 0.1074661302677646D1*t144*t27+0.1074661302677646D1*t147
     #*t26-0.5373306513388232D0*t11*t161
      vsigmaaa(i) = -0.9305257363491D0*t5*t164
      vsigmaab(i) = 0.D0
      t168 = t35**2
      t169 = 1/t168
      t172 = 0.9474169D-3*t36-0.9474169D-3*sigmabb*t169
      t173 = t172*sigmabb
      t176 = t39*t44
      t189 = 18.D0*t40*t50*t172+9.D0*t45*t44*t49+18.D0*t46*t123*t48*t127
      t190 = t110*t189
      t193 = 0.1074661302677646D1*t173*t56+0.1074661302677646D1*t176
     #*t55-0.5373306513388232D0*t40*t190
      vsigmabb(i) = -0.9305257363491D0*t34*t193
      t202 = 1/t13/t78
      t206 = t66*t71
      t210 = t24**2
      t212 = 1/t25/t210
      t213 = t15*t212
      t214 = t93**2
      t229 = t16*t76*sigmaaa
      t230 = t78**2
      t231 = t230*t12
      t233 = 1/t87
      t237 = t76**2
      t244 = 1/t88/t87
      v2rhoa2(i) = -0.4135669939329333D0/t13*t30-0.24814019635976D1*t4
     #*t97-0.9305257363491D0*t5*(0.1050779940395921D2*t11*t202*t26
     #+0.2865763473807057D1*t11*t206*t93+0.8059959770082348D0*t11*t213
     #*t214-0.5373306513388232D0*t11*t72*(88.D0*t17*t202*t20+432.D0
     #*t77/t4/t78/t64*t19*t89+128.D0*t229/t231*t233-128.D0*t16*t237
     #/t13/t230/t78*t19*t244))
      t261 = 1/t42/t116
      t265 = t104*t109
      t269 = t53**2
      t271 = 1/t54/t269
      t272 = t44*t271
      t273 = t131**2
      t288 = t45*t114*sigmabb
      t289 = t116**2
      t290 = t289*t41
      t292 = 1/t125
      t296 = t114**2
      t303 = 1/t126/t125
      v2rhob2(i) = -0.4135669939329333D0/t42*t59-0.24814019635976D1
     #*t33*t135-0.9305257363491D0*t34*(0.1050779940395921D2*t40*t261
     #*t55+0.2865763473807057D1*t40*t265*t131+0.8059959770082348D0*t40
     #*t272*t273-0.5373306513388232D0*t40*t110*(88.D0*t46*t261*t49
     #+432.D0*t115/t33/t116/t102*t48*t127+128.D0*t288/t290*t292-128.D0
     #*t45*t296/t42/t289/t116*t48*t303))
      v2rhoab(i) = 0.D0
      t344 = t19*t89
      t345 = t344*t143
      v2rhoasigmaaa(i) = -0.12407009817988D1*t4*t164-0.9305257363491D0
     #*t5*(-0.2865763473807057D1*t144*t67-0.2865763473807057D1*t10*t66
     #*t26+0.1432881736903529D1*t11*t206*t160-0.5373306513388232D0
     #*t144*t94-0.5373306513388232D0*t147*t71*t93+0.8059959770082348D0
     #*t11*t15*t212*t93*t160-0.5373306513388232D0*t11*t72*(-48.D0*t11
     #*t73*t143-24.D0*t16*t66*t20-144.D0*t17*t90-96.D0*t10*t76*t81
     #*t345-48.D0*t77/t230/rhoa*t233+48.D0*t229/t13/t230/t64*t19*t244))
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      t397 = t48*t127
      t398 = t397*t172
      v2rhobsigmabb(i) = -0.12407009817988D1*t33*t193
     #-0.9305257363491D0*t34*(-0.2865763473807057D1*t173*t105
     #-0.2865763473807057D1*t39*t104*t55+0.1432881736903529D1*t40*t265
     #*t189-0.5373306513388232D0*t173*t132-0.5373306513388232D0*t176
     #*t109*t131+0.8059959770082348D0*t40*t44*t271*t131*t189
     #-0.5373306513388232D0*t40*t110*(-48.D0*t40*t111*t172-24.D0*t45
     #*t104*t49-144.D0*t46*t128-96.D0*t39*t114*t119*t398-48.D0*t115
     #/t289/rhob*t292+48.D0*t288/t42/t289/t102*t48*t303))
      t425 = -0.18948338D-2*t140+0.18948338D-2*sigmaaa/t139/t6
      t437 = t160**2
      t441 = t143**2
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(0.1074661302677646D1*t425
     #*sigmaaa*t27+0.2149322605355293D1*t143*t15*t26
     #-0.1074661302677646D1*t144*t161-0.1074661302677646D1*t147*t71
     #*t160+0.8059959770082348D0*t11*t213*t437-0.5373306513388232D0
     #*t11*t72*(18.D0*t441*sigmaaa*t21+36.D0*t147*t20*t143+72.D0*t11
     #*t85*t345+18.D0*t11*t21*t425+36.D0*t16*t85*t344+18.D0*t17/t230
     #*t233-18.D0*t77/t13/t231*t19*t244))
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t479 = -0.18948338D-2*t169+0.18948338D-2*sigmabb/t168/t35
      t491 = t189**2
      t495 = t172**2
      v2sigmabb2(i) = -0.9305257363491D0*t34*(0.1074661302677646D1
     #*t479*sigmabb*t56+0.2149322605355293D1*t172*t44*t55
     #-0.1074661302677646D1*t173*t190-0.1074661302677646D1*t176*t109
     #*t189+0.8059959770082348D0*t40*t272*t491-0.5373306513388232D0
     #*t40*t110*(18.D0*t495*sigmabb*t50+36.D0*t176*t49*t172+72.D0*t40
     #*t123*t398+18.D0*t40*t50*t479+36.D0*t45*t123*t397+18.D0*t46/t289
     #*t292-18.D0*t115/t42/t290*t48*t303))
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
      
      
      subroutine rks_x_ft97b
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
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
      t9 = 0.2913644D-2+0.236854225D-3*sigma/(0.6255746320201D7+0.25D0
     #*sigma)
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t9**2
      t19 = dlog(0.1587401051968199D1*sigma*t14+dsqrt(1
     #+0.2519842099789745D1*sigma**2*t14**2))
      t20 = t19**2
      t25 = dsqrt(1.D0+0.142866094677138D2*t15*sigma*t14*t20)
      zk(i) = -0.7385587663820224D0*t2*rho*(1.D0+0.1705918482380012D1
     #*t9*sigma*t14/t25)
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
      t5 = 0.6255746320201D7+0.25D0*sigma
      t6 = 1/t5
      t9 = 0.2913644D-2+0.236854225D-3*sigma*t6
      t10 = t9*sigma
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t9**2
      t16 = t15*sigma
      t19 = dlog(0.1587401051968199D1*sigma*t14+dsqrt(1
     #+0.2519842099789745D1*sigma**2*t14**2))
      t20 = t19**2
      t21 = t14*t20
      t24 = 1.D0+0.142866094677138D2*t16*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t14*t26
      t30 = 1.D0+0.1705918482380012D1*t10*t27
      zk(i) = -0.7385587663820224D0*t3*t30
      t37 = 1/t12/t11/rho
      t43 = t14/t25/t24
      t47 = sigma**2
      t49 = t11**2
      t56 = 1/t2/t49/rho
      t60 = dsqrt(1.D0+0.2519842099789746D1*t47*t56)
      t61 = 1/t60
      vrhoa(i) = -0.9847450218426965D0*t2*t30-0.3692793831910112D0*t3*
     #(-0.9098231906026728D1*t10*t37*t26-0.8529592411900058D0*t10*t43*
     #(-0.7619525049447357D2*t16*t37*t20-0.2419048415798156D3*t15*t47
     #/t2/t49/t11*t19*t61))
      t73 = t5**2
      t77 = 0.9474169D-3*t6-0.236854225D-3*sigma/t73
      vsigmaaa(i) = -0.7385587663820224D0*t3*(0.1705918482380012D1*t77
     #*sigma*t27+0.6823673929520046D1*t9*t14*t26-0.8529592411900058D0
     #*t10*t43*(0.2857321893542759D2*t10*t21*t77+0.5714643787085518D2
     #*t15*t14*t20+0.1814286311848617D3*t16*t56*t19*t61))
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
      t5 = 0.6255746320201D7+0.25D0*sigma
      t6 = 1/t5
      t9 = 0.2913644D-2+0.236854225D-3*sigma*t6
      t10 = t9*sigma
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t9**2
      t16 = t15*sigma
      t19 = dlog(0.1587401051968199D1*sigma*t14+dsqrt(1
     #+0.2519842099789745D1*sigma**2*t14**2))
      t20 = t19**2
      t21 = t14*t20
      t24 = 1.D0+0.142866094677138D2*t16*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t14*t26
      t30 = 1.D0+0.1705918482380012D1*t10*t27
      zk(i) = -0.7385587663820224D0*t3*t30
      t35 = t11*rho
      t37 = 1/t12/t35
      t38 = t37*t26
      t42 = 1/t25/t24
      t43 = t14*t42
      t44 = t37*t20
      t47 = sigma**2
      t48 = t15*t47
      t49 = t11**2
      t52 = 1/t2/t49/t11
      t56 = 1/t2/t49/rho
      t59 = 1.D0+0.2519842099789746D1*t47*t56
      t60 = dsqrt(t59)
      t61 = 1/t60
      t62 = t52*t19*t61
      t65 = -0.7619525049447357D2*t16*t44-0.2419048415798156D3*t48*t62
      t66 = t43*t65
      t69 = -0.9098231906026728D1*t10*t38-0.8529592411900058D0*t10*t66
      vrhoa(i) = -0.9847450218426965D0*t2*t30-0.3692793831910112D0*t3
     #*t69
      t73 = t5**2
      t74 = 1/t73
      t77 = 0.9474169D-3*t6-0.236854225D-3*sigma*t74
      t78 = t77*sigma
      t81 = t9*t14
      t94 = 0.2857321893542759D2*t10*t21*t77+0.5714643787085518D2*t15
     #*t14*t20+0.1814286311848617D3*t16*t56*t19*t61
      t95 = t43*t94
      t98 = 0.1705918482380012D1*t78*t27+0.6823673929520046D1*t81*t26
     #-0.8529592411900058D0*t10*t95
      vsigmaaa(i) = -0.7385587663820224D0*t3*t98
      t107 = 1/t12/t49
      t111 = t37*t42
      t115 = t24**2
      t117 = 1/t25/t115
      t118 = t14*t117
      t119 = t65**2
      t134 = t15*t47*sigma
      t135 = t49**2
      t136 = t135*t11
      t138 = 1/t59
      t142 = t47**2
      t149 = 1/t60/t59
      v2rhoa2(i) = -0.6564966812284644D0/t12*t30-0.1969490043685393D1
     #*t2*t69-0.3692793831910112D0*t3*(0.6672036731086267D2*t10*t107
     #*t26+0.9098231906026728D1*t10*t111*t65+0.1279438861785009D1*t10
     #*t118*t119-0.8529592411900058D0*t10*t43*(0.5587651702928062D3
     #*t16*t107*t20+0.4354287148436682D4*t48/t2/t49/t35*t19*t61
     #+2048.D0*t134/t136*t138-0.3250997354430873D4*t15*t142/t12/t135
     #/t49*t19*t149))
      t190 = t19*t61
      t191 = t190*t77
      v2rhoasigmaaa(i) = -0.9847450218426965D0*t2*t98
     #-0.3692793831910112D0*t3*(-0.9098231906026728D1*t78*t38
     #-0.3639292762410691D2*t9*t37*t26+0.4549115953013364D1*t10*t111
     #*t94-0.8529592411900058D0*t78*t66-0.3411836964760023D1*t81*t42
     #*t65+0.1279438861785009D1*t10*t14*t117*t65*t94
     #-0.8529592411900058D0*t10*t43*(-0.1523905009889471D3*t10*t44*t77
     #-0.3047810019778943D3*t15*t37*t20-0.2902858098957788D4*t16*t62
     #-0.4838096831596313D3*t9*t47*t52*t191-1536.D0*t48/t135/rho*t138
     #+0.2438248015823154D4*t134/t12/t135/t35*t19*t149))
      t218 = -0.18948338D-2*t74+0.47370845D-3*sigma/t73/t5
      t230 = t94**2
      t234 = t77**2
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(0.1705918482380012D1
     #*t218*sigma*t27+0.1364734785904009D2*t77*t14*t26
     #-0.1705918482380012D1*t78*t95-0.6823673929520046D1*t81*t42*t94
     #+0.1279438861785009D1*t10*t118*t230-0.8529592411900058D0*t10*t43
     #*(0.2857321893542759D2*t234*sigma*t21+0.2285857514834207D3*t81
     #*t20*t77+0.7257145247394469D3*t10*t56*t191+0.2857321893542759D2
     #*t10*t21*t218+0.1451429049478894D4*t15*t56*t190+1152.D0*t16/t135
     #*t138-0.1828686011867366D4*t48/t12/t136*t19*t149))
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
c:X_FT97Bsubrend
