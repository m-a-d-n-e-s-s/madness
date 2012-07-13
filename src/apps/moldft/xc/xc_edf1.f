c:XC_EDF1subrstart
c    Generated: Wed Jul 23 17:14:37 GMT 2003
      subroutine uks_xc_edf1
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     R.D. Adamson, P.M.W. Gill, and J.A. Pople
c     Empirical density functionals
c     Chem. Phys. Lett. 284 (1998) 6-11
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
      t2 = rhob**2
      t3 = rhob**(1.D0/3.D0)
      t4 = t3**2
      t7 = sigmabb/t4/t2
      t8 = dsqrt(sigmabb)
      t9 = t3*rhob
      t11 = t8/t9
      t12 = dlog(t11+dsqrt(1+t11**2))
      t13 = t11*t12
      zk(i) = (-0.9593273689405774D0-0.3640595D-1*t7/(1.D0+0.21D-1*t13
     #)+0.35481306D-1*t7/(1.D0+0.252D-1*t13))*t9
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**2
      t3 = rhoa**(1.D0/3.D0)
      t4 = t3**2
      t7 = sigmaaa/t4/t2
      t8 = dsqrt(sigmaaa)
      t9 = t3*rhoa
      t11 = t8/t9
      t12 = dlog(t11+dsqrt(1+t11**2))
      t13 = t11*t12
      zk(i) = (-0.9593273689405774D0-0.3640595D-1*t7/(1.D0+0.21D-1*t13
     #)+0.35481306D-1*t7/(1.D0+0.252D-1*t13))*t9
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**2
      t5 = rhoa**(1.D0/3.D0)
      t6 = t5**2
      t7 = t6*t4
      t9 = sigmaaa/t7
      t10 = dsqrt(sigmaaa)
      t11 = t5*rhoa
      t13 = t10/t11
      t14 = dlog(t13+dsqrt(1+t13**2))
      t15 = t13*t14
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t31 = t30*t28
      t33 = sigmabb/t31
      t34 = dsqrt(sigmabb)
      t35 = t29*rhob
      t37 = t34/t35
      t38 = dlog(t37+dsqrt(1+t37**2))
      t39 = t37*t38
      t52 = rho**(1.D0/3.D0)
      t53 = 1/t52
      t56 = 1/(1.D0+0.3505D0*t53)
      t58 = 1/rho
      t59 = rhob*t58
      t62 = 0.25D0*t53
      t63 = dexp(-t62)
      t65 = rho**2
      t67 = t52**2
      t74 = t53*t56
      t96 = 0.6666666666666667D0*t65
      zk(i) = (-0.9593273689405774D0-0.3640595D-1*t9/(1.D0+0.21D-1*t15
     #)+0.35481306D-1*t9/(1.D0+0.252D-1*t15))*t11+(
     #-0.9593273689405774D0-0.3640595D-1*t33/(1.D0+0.21D-1*t39)
     #+0.35481306D-1*t33/(1.D0+0.252D-1*t39))*t35-0.22D0*t56*rhoa*t59
     #-0.869D-2*t63*t56/t67/t65/rho*(rhoa*rhob*(0.3646239897876478D2
     #*t7+0.3646239897876478D2*t31+(0.2611111111111111D1
     #-0.9722222222222222D-1*t53-0.1363055555555556D0*t74)*sigma-1.D0*
     #(0.25D1-0.1388888888888889D-1*t53-0.1947222222222222D-1*t74)*
     #(sigmaaa+sigmabb)-0.1111111111111111D0*(t62+0.3505D0*t74-11.D0)*
     #(rhoa*t58*sigmaaa+t59*sigmabb))-0.6666666666666667D0*t65*sigma+
     #(t96-1.D0*t4)*sigmabb+(t96-1.D0*t28)*sigmaaa)
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
      t2 = rhob**2
      t3 = rhob**(1.D0/3.D0)
      t4 = t3**2
      t6 = 1/t4/t2
      t7 = sigmabb*t6
      t8 = dsqrt(sigmabb)
      t9 = t3*rhob
      t10 = 1/t9
      t11 = t8*t10
      t12 = dlog(t11+dsqrt(1+t11**2))
      t13 = t11*t12
      t15 = 1.D0+0.21D-1*t13
      t16 = 1/t15
      t20 = 1.D0+0.252D-1*t13
      t21 = 1/t20
      t24 = -0.9593273689405774D0-0.3640595D-1*t7*t16+0.35481306D-1*t7
     #*t21
      zk(i) = t24*t9
      vrhoa(i) = 0.D0
      t28 = sigmabb/t4/t2/rhob
      t31 = t15**2
      t32 = 1/t31
      t36 = t8/t3/t2*t12
      t39 = dsqrt(1.D0+t7)
      t40 = 1/t39
      t41 = t28*t40
      t49 = t20**2
      t50 = 1/t49
      vrhob(i) = (0.9708253333333333D-1*t28*t16+0.3640595D-1*t7*t32*(
     #-0.28D-1*t36-0.28D-1*t41)-0.94616816D-1*t28*t21-0.35481306D-1*t7
     #*t50*(-0.336D-1*t36-0.336D-1*t41))*t9+0.1333333333333333D1*t24*t3
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t65 = 1/t8*t10*t12
      t67 = t6*t40
      vsigmabb(i) = (-0.3640595D-1*t6*t16+0.3640595D-1*t7*t32*
     #(0.105D-1*t65+0.105D-1*t67)+0.35481306D-1*t6*t21-0.35481306D-1
     #*t7*t50*(0.126D-1*t65+0.126D-1*t67))*t9
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**2
      t3 = rhoa**(1.D0/3.D0)
      t4 = t3**2
      t6 = 1/t4/t2
      t7 = sigmaaa*t6
      t8 = dsqrt(sigmaaa)
      t9 = t3*rhoa
      t10 = 1/t9
      t11 = t8*t10
      t12 = dlog(t11+dsqrt(1+t11**2))
      t13 = t11*t12
      t15 = 1.D0+0.21D-1*t13
      t16 = 1/t15
      t20 = 1.D0+0.252D-1*t13
      t21 = 1/t20
      t24 = -0.9593273689405774D0-0.3640595D-1*t7*t16+0.35481306D-1*t7
     #*t21
      zk(i) = t24*t9
      t28 = sigmaaa/t4/t2/rhoa
      t31 = t15**2
      t32 = 1/t31
      t36 = t8/t3/t2*t12
      t39 = dsqrt(1.D0+t7)
      t40 = 1/t39
      t41 = t28*t40
      t49 = t20**2
      t50 = 1/t49
      vrhoa(i) = (0.9708253333333333D-1*t28*t16+0.3640595D-1*t7*t32*(
     #-0.28D-1*t36-0.28D-1*t41)-0.94616816D-1*t28*t21-0.35481306D-1*t7
     #*t50*(-0.336D-1*t36-0.336D-1*t41))*t9+0.1333333333333333D1*t24*t3
      vrhob(i) = 0.D0
      t65 = 1/t8*t10*t12
      t67 = t6*t40
      vsigmaaa(i) = (-0.3640595D-1*t6*t16+0.3640595D-1*t7*t32*
     #(0.105D-1*t65+0.105D-1*t67)+0.35481306D-1*t6*t21-0.35481306D-1
     #*t7*t50*(0.126D-1*t65+0.126D-1*t67))*t9
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**2
      t5 = rhoa**(1.D0/3.D0)
      t6 = t5**2
      t7 = t6*t4
      t8 = 1/t7
      t9 = sigmaaa*t8
      t10 = dsqrt(sigmaaa)
      t11 = t5*rhoa
      t12 = 1/t11
      t13 = t10*t12
      t14 = dlog(t13+dsqrt(1+t13**2))
      t15 = t13*t14
      t17 = 1.D0+0.21D-1*t15
      t18 = 1/t17
      t22 = 1.D0+0.252D-1*t15
      t23 = 1/t22
      t26 = -0.9593273689405774D0-0.3640595D-1*t9*t18+0.35481306D-1*t9
     #*t23
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t31 = t30*t28
      t32 = 1/t31
      t33 = sigmabb*t32
      t34 = dsqrt(sigmabb)
      t35 = t29*rhob
      t36 = 1/t35
      t37 = t34*t36
      t38 = dlog(t37+dsqrt(1+t37**2))
      t39 = t37*t38
      t41 = 1.D0+0.21D-1*t39
      t42 = 1/t41
      t46 = 1.D0+0.252D-1*t39
      t47 = 1/t46
      t50 = -0.9593273689405774D0-0.3640595D-1*t33*t42+0.35481306D-1
     #*t33*t47
      t52 = rho**(1.D0/3.D0)
      t53 = 1/t52
      t55 = 1.D0+0.3505D0*t53
      t56 = 1/t55
      t57 = t56*rhoa
      t58 = 1/rho
      t59 = rhob*t58
      t62 = 0.25D0*t53
      t63 = dexp(-t62)
      t64 = t63*t56
      t65 = rho**2
      t67 = t52**2
      t69 = 1/t67/t65/rho
      t70 = rhoa*rhob
      t74 = t53*t56
      t76 = 0.2611111111111111D1-0.9722222222222222D-1*t53
     #-0.1363055555555556D0*t74
      t81 = sigmaaa+sigmabb
      t85 = t62+0.3505D0*t74-11.D0
      t89 = rhoa*t58*sigmaaa+t59*sigmabb
      t92 = 0.3646239897876478D2*t7+0.3646239897876478D2*t31+t76*sigma
     #-1.D0*(0.25D1-0.1388888888888889D-1*t53-0.1947222222222222D-1
     #*t74)*t81-0.1111111111111111D0*t85*t89
      t96 = 0.6666666666666667D0*t65
      t97 = 1.D0*t4
      t100 = 1.D0*t28
      t103 = t70*t92-0.6666666666666667D0*t65*sigma+(t96-t97)*sigmabb+
     #(t96-t100)*sigmaaa
      zk(i) = t26*t11+t50*t35-0.22D0*t57*t59-0.869D-2*t64*t69*t103
      t110 = sigmaaa/t6/t4/rhoa
      t113 = t17**2
      t114 = 1/t113
      t118 = t10/t5/t4*t14
      t121 = dsqrt(1.D0+t9)
      t122 = 1/t121
      t123 = t110*t122
      t131 = t22**2
      t132 = 1/t131
      t149 = t85*t58
      t160 = t55**2
      t161 = 1/t160
      t167 = 0.2570333333333333D-1*t161*rhoa*rhob/t52/t65
      t168 = 1/t65
      t169 = rhob*t168
      t171 = 0.22D0*t57*t169
      t172 = t65**2
      t174 = 1/t172/rho
      t178 = 0.7241666666666667D-3*t174*t63*t56*t103
      t182 = 0.1015281666666667D-2*t63*t161*t174*t103
      t187 = 0.3186333333333333D-1*t64/t67/t172*t103
      t189 = 1/t52/rho
      t191 = t189*t56
      t195 = 1/t67/rho*t161
      t230 = 0.869D-2*t64*t69*(t70*((0.3240740740740741D-1*t189
     #+0.4543518518518519D-1*t191-0.1592503240740741D-1*t195)*sigma
     #-1.D0*(0.462962962962963D-2*t189+0.6490740740740741D-2*t191
     #-0.227500462962963D-2*t195)*t81-0.1111111111111111D0*(
     #-0.8333333333333333D-1*t189-0.1168333333333333D0*t191
     #+0.4095008333333333D-1*t195)*t89-0.1111111111111111D0*t85*(-1.D0
     #*rhoa*t168*sigmaaa-1.D0*t169*sigmabb))-0.1333333333333333D1*rho
     #*sigma+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho
     #*sigmaaa)
      vrhoa(i) = (0.9708253333333333D-1*t110*t18+0.3640595D-1*t9*t114*
     #(-0.28D-1*t118-0.28D-1*t123)-0.94616816D-1*t110*t23
     #-0.35481306D-1*t9*t132*(-0.336D-1*t118-0.336D-1*t123))*t11
     #+0.1333333333333333D1*t26*t5-0.22D0*t56*rhob*t58-0.869D-2*t64
     #*t69*(rhob*t92+t70*(0.9723306394337274D2*t6*rhoa
     #-0.1111111111111111D0*t149*sigmaaa)-2.D0*rhoa*sigmabb)-t167+t171
     #-t178-t182+t187-t230
      t234 = sigmabb/t30/t28/rhob
      t237 = t41**2
      t238 = 1/t237
      t242 = t34/t29/t28*t38
      t245 = dsqrt(1.D0+t33)
      t246 = 1/t245
      t247 = t234*t246
      t255 = t46**2
      t256 = 1/t255
      vrhob(i) = (0.9708253333333333D-1*t234*t42+0.3640595D-1*t33*t238
     #*(-0.28D-1*t242-0.28D-1*t247)-0.94616816D-1*t234*t47
     #-0.35481306D-1*t33*t256*(-0.336D-1*t242-0.336D-1*t247))*t35
     #+0.1333333333333333D1*t50*t29-0.22D0*t57*t58-0.869D-2*t64*t69*
     #(rhoa*t92+t70*(0.9723306394337274D2*t30*rhob
     #-0.1111111111111111D0*t149*sigmabb)-2.D0*rhob*sigmaaa)-t167+t171
     #-t178-t182+t187-t230
      t286 = 1/t10*t12*t14
      t288 = t8*t122
      t304 = 0.1388888888888889D-1*t53
      t305 = 0.1947222222222222D-1*t74
      t319 = t64*t69*(t70*t76-0.6666666666666667D0*t65)
      t320 = 0.869D-2*t319
      vsigmaaa(i) = (-0.3640595D-1*t8*t18+0.3640595D-1*t9*t114*
     #(0.105D-1*t286+0.105D-1*t288)+0.35481306D-1*t8*t23-0.35481306D-1
     #*t9*t132*(0.126D-1*t286+0.126D-1*t288))*t11-0.869D-2*t64*t69*
     #(t70*(-0.25D1+t304+t305-0.1111111111111111D0*t85*rhoa*t58)+t96
     #-t100)-t320
      vsigmaab(i) = -0.1738D-1*t319
      t326 = 1/t34*t36*t38
      t328 = t32*t246
      vsigmabb(i) = (-0.3640595D-1*t32*t42+0.3640595D-1*t33*t238*
     #(0.105D-1*t326+0.105D-1*t328)+0.35481306D-1*t32*t47
     #-0.35481306D-1*t33*t256*(0.126D-1*t326+0.126D-1*t328))*t35
     #-0.869D-2*t64*t69*(t70*(-0.25D1+t304+t305-0.1111111111111111D0
     #*t85*rhob*t58)+t96-t97)-t320
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
      t2 = rhob**2
      t3 = rhob**(1.D0/3.D0)
      t4 = t3**2
      t6 = 1/t4/t2
      t7 = sigmabb*t6
      t8 = dsqrt(sigmabb)
      t9 = t3*rhob
      t10 = 1/t9
      t11 = t8*t10
      t12 = dlog(t11+dsqrt(1+t11**2))
      t13 = t11*t12
      t15 = 1.D0+0.21D-1*t13
      t16 = 1/t15
      t20 = 1.D0+0.252D-1*t13
      t21 = 1/t20
      t24 = -0.9593273689405774D0-0.3640595D-1*t7*t16+0.35481306D-1*t7
     #*t21
      zk(i) = t24*t9
      vrhoa(i) = 0.D0
      t25 = t2*rhob
      t28 = sigmabb/t4/t25
      t31 = t15**2
      t32 = 1/t31
      t36 = t8/t3/t2*t12
      t38 = 1.D0+t7
      t39 = dsqrt(t38)
      t40 = 1/t39
      t41 = t28*t40
      t43 = -0.28D-1*t36-0.28D-1*t41
      t44 = t32*t43
      t49 = t20**2
      t50 = 1/t49
      t53 = -0.336D-1*t36-0.336D-1*t41
      t54 = t50*t53
      t57 = 0.9708253333333333D-1*t28*t16+0.3640595D-1*t7*t44
     #-0.94616816D-1*t28*t21-0.35481306D-1*t7*t54
      vrhob(i) = t57*t9+0.1333333333333333D1*t24*t3
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t65 = 1/t8*t10*t12
      t67 = t6*t40
      t69 = 0.105D-1*t65+0.105D-1*t67
      t77 = 0.126D-1*t65+0.126D-1*t67
      vsigmabb(i) = (-0.3640595D-1*t6*t16+0.3640595D-1*t7*t32*t69
     #+0.35481306D-1*t6*t21-0.35481306D-1*t7*t50*t77)*t9
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t82 = t2**2
      t85 = sigmabb/t4/t82
      t91 = 1/t31/t15
      t92 = t43**2
      t99 = t8/t3/t25*t12
      t101 = t85*t40
      t103 = sigmabb**2
      t109 = 1/t39/t38
      t110 = t103/t3/t82/t25*t109
      t121 = 1/t49/t20
      t122 = t53**2
      v2rhob2(i) = (-0.3559692888888889D0*t85*t16-0.1941650666666667D0
     #*t28*t44-0.728119D-1*t7*t91*t92+0.3640595D-1*t7*t32*
     #(0.6533333333333333D-1*t99+0.14D0*t101-0.3733333333333333D-1
     #*t110)+0.3469283253333333D0*t85*t21+0.189233632D0*t28*t54
     #+0.70962612D-1*t7*t121*t122-0.35481306D-1*t7*t50*(0.784D-1*t99
     #+0.168D0*t101-0.448D-1*t110))*t9+0.2666666666666667D1*t57*t3
     #+0.4444444444444444D0*t24/t4
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t143 = t69**2
      t150 = 1/t8/sigmabb*t10*t12
      t154 = 1/sigmabb*t6*t40
      t159 = 1/t3/t82/rhob*t109
      t168 = t77**2
      v2sigmabb2(i) = (0.728119D-1*t6*t32*t69-0.728119D-1*t7*t91*t143
     #+0.3640595D-1*t7*t32*(-0.525D-2*t150+0.525D-2*t154-0.525D-2*t159
     #)-0.70962612D-1*t6*t50*t77+0.70962612D-1*t7*t121*t168
     #-0.35481306D-1*t7*t50*(-0.63D-2*t150+0.63D-2*t154-0.63D-2*t159))
     #*t9
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**2
      t3 = rhoa**(1.D0/3.D0)
      t4 = t3**2
      t6 = 1/t4/t2
      t7 = sigmaaa*t6
      t8 = dsqrt(sigmaaa)
      t9 = t3*rhoa
      t10 = 1/t9
      t11 = t8*t10
      t12 = dlog(t11+dsqrt(1+t11**2))
      t13 = t11*t12
      t15 = 1.D0+0.21D-1*t13
      t16 = 1/t15
      t20 = 1.D0+0.252D-1*t13
      t21 = 1/t20
      t24 = -0.9593273689405774D0-0.3640595D-1*t7*t16+0.35481306D-1*t7
     #*t21
      zk(i) = t24*t9
      t25 = t2*rhoa
      t28 = sigmaaa/t4/t25
      t31 = t15**2
      t32 = 1/t31
      t36 = t8/t3/t2*t12
      t38 = 1.D0+t7
      t39 = dsqrt(t38)
      t40 = 1/t39
      t41 = t28*t40
      t43 = -0.28D-1*t36-0.28D-1*t41
      t44 = t32*t43
      t49 = t20**2
      t50 = 1/t49
      t53 = -0.336D-1*t36-0.336D-1*t41
      t54 = t50*t53
      t57 = 0.9708253333333333D-1*t28*t16+0.3640595D-1*t7*t44
     #-0.94616816D-1*t28*t21-0.35481306D-1*t7*t54
      vrhoa(i) = t57*t9+0.1333333333333333D1*t24*t3
      vrhob(i) = 0.D0
      t65 = 1/t8*t10*t12
      t67 = t6*t40
      t69 = 0.105D-1*t65+0.105D-1*t67
      t77 = 0.126D-1*t65+0.126D-1*t67
      vsigmaaa(i) = (-0.3640595D-1*t6*t16+0.3640595D-1*t7*t32*t69
     #+0.35481306D-1*t6*t21-0.35481306D-1*t7*t50*t77)*t9
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t82 = t2**2
      t85 = sigmaaa/t4/t82
      t91 = 1/t31/t15
      t92 = t43**2
      t99 = t8/t3/t25*t12
      t101 = t85*t40
      t103 = sigmaaa**2
      t109 = 1/t39/t38
      t110 = t103/t3/t82/t25*t109
      t121 = 1/t49/t20
      t122 = t53**2
      v2rhoa2(i) = (-0.3559692888888889D0*t85*t16-0.1941650666666667D0
     #*t28*t44-0.728119D-1*t7*t91*t92+0.3640595D-1*t7*t32*
     #(0.6533333333333333D-1*t99+0.14D0*t101-0.3733333333333333D-1
     #*t110)+0.3469283253333333D0*t85*t21+0.189233632D0*t28*t54
     #+0.70962612D-1*t7*t121*t122-0.35481306D-1*t7*t50*(0.784D-1*t99
     #+0.168D0*t101-0.448D-1*t110))*t9+0.2666666666666667D1*t57*t3
     #+0.4444444444444444D0*t24/t4
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t143 = t69**2
      t150 = 1/t8/sigmaaa*t10*t12
      t154 = 1/sigmaaa*t6*t40
      t159 = 1/t3/t82/rhoa*t109
      t168 = t77**2
      v2sigmaaa2(i) = (0.728119D-1*t6*t32*t69-0.728119D-1*t7*t91*t143
     #+0.3640595D-1*t7*t32*(-0.525D-2*t150+0.525D-2*t154-0.525D-2*t159
     #)-0.70962612D-1*t6*t50*t77+0.70962612D-1*t7*t121*t168
     #-0.35481306D-1*t7*t50*(-0.63D-2*t150+0.63D-2*t154-0.63D-2*t159))
     #*t9
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
      t4 = rhoa**2
      t5 = rhoa**(1.D0/3.D0)
      t6 = t5**2
      t7 = t6*t4
      t8 = 1/t7
      t9 = sigmaaa*t8
      t10 = dsqrt(sigmaaa)
      t11 = t5*rhoa
      t12 = 1/t11
      t13 = t10*t12
      t14 = dlog(t13+dsqrt(1+t13**2))
      t15 = t13*t14
      t17 = 1.D0+0.21D-1*t15
      t18 = 1/t17
      t22 = 1.D0+0.252D-1*t15
      t23 = 1/t22
      t26 = -0.9593273689405774D0-0.3640595D-1*t9*t18+0.35481306D-1*t9
     #*t23
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t31 = t30*t28
      t32 = 1/t31
      t33 = sigmabb*t32
      t34 = dsqrt(sigmabb)
      t35 = t29*rhob
      t36 = 1/t35
      t37 = t34*t36
      t38 = dlog(t37+dsqrt(1+t37**2))
      t39 = t37*t38
      t41 = 1.D0+0.21D-1*t39
      t42 = 1/t41
      t46 = 1.D0+0.252D-1*t39
      t47 = 1/t46
      t50 = -0.9593273689405774D0-0.3640595D-1*t33*t42+0.35481306D-1
     #*t33*t47
      t52 = rho**(1.D0/3.D0)
      t53 = 1/t52
      t55 = 1.D0+0.3505D0*t53
      t56 = 1/t55
      t57 = t56*rhoa
      t58 = 1/rho
      t59 = rhob*t58
      t62 = 0.25D0*t53
      t63 = dexp(-t62)
      t64 = t63*t56
      t65 = rho**2
      t66 = t65*rho
      t67 = t52**2
      t69 = 1/t67/t66
      t70 = rhoa*rhob
      t71 = 0.3646239897876478D2*t7
      t72 = 0.3646239897876478D2*t31
      t74 = t53*t56
      t76 = 0.2611111111111111D1-0.9722222222222222D-1*t53
     #-0.1363055555555556D0*t74
      t77 = t76*sigma
      t81 = sigmaaa+sigmabb
      t83 = 1.D0*(0.25D1-0.1388888888888889D-1*t53
     #-0.1947222222222222D-1*t74)*t81
      t85 = t62+0.3505D0*t74-11.D0
      t89 = rhoa*t58*sigmaaa+t59*sigmabb
      t91 = 0.1111111111111111D0*t85*t89
      t92 = t71+t72+t77-t83-t91
      t96 = 0.6666666666666667D0*t65
      t97 = 1.D0*t4
      t100 = 1.D0*t28
      t103 = t70*t92-0.6666666666666667D0*t65*sigma+(t96-t97)*sigmabb+
     #(t96-t100)*sigmaaa
      zk(i) = t26*t11+t50*t35-0.22D0*t57*t59-0.869D-2*t64*t69*t103
      t107 = t4*rhoa
      t109 = 1/t6/t107
      t110 = sigmaaa*t109
      t113 = t17**2
      t114 = 1/t113
      t116 = 1/t5/t4
      t118 = t10*t116*t14
      t120 = 1.D0+t9
      t121 = dsqrt(t120)
      t122 = 1/t121
      t123 = t110*t122
      t125 = -0.28D-1*t118-0.28D-1*t123
      t126 = t114*t125
      t131 = t22**2
      t132 = 1/t131
      t135 = -0.336D-1*t118-0.336D-1*t123
      t136 = t132*t135
      t139 = 0.9708253333333333D-1*t110*t18+0.3640595D-1*t9*t126
     #-0.94616816D-1*t110*t23-0.35481306D-1*t9*t136
      t143 = t56*rhob
      t147 = t6*rhoa
      t149 = t85*t58
      t152 = 0.9723306394337274D2*t147-0.1111111111111111D0*t149*sigmaaa
      t156 = rhob*t92+t70*t152-2.D0*rhoa*sigmabb
      t160 = t55**2
      t161 = 1/t160
      t162 = t161*rhoa
      t164 = 1/t52/t65
      t167 = 0.2570333333333333D-1*t162*rhob*t164
      t168 = 1/t65
      t169 = rhob*t168
      t171 = 0.22D0*t57*t169
      t172 = t65**2
      t173 = t172*rho
      t174 = 1/t173
      t175 = t174*t63
      t176 = t56*t103
      t178 = 0.7241666666666667D-3*t175*t176
      t179 = t63*t161
      t182 = 0.1015281666666667D-2*t179*t174*t103
      t184 = 1/t67/t172
      t187 = 0.3186333333333333D-1*t64*t184*t103
      t189 = 1/t52/rho
      t191 = t189*t56
      t195 = 1/t67/rho*t161
      t197 = 0.3240740740740741D-1*t189+0.4543518518518519D-1*t191
     #-0.1592503240740741D-1*t195
      t208 = -0.8333333333333333D-1*t189-0.1168333333333333D0*t191
     #+0.4095008333333333D-1*t195
      t216 = -1.D0*rhoa*t168*sigmaaa-1.D0*t169*sigmabb
      t219 = t197*sigma-1.D0*(0.462962962962963D-2*t189
     #+0.6490740740740741D-2*t191-0.227500462962963D-2*t195)*t81
     #-0.1111111111111111D0*t208*t89-0.1111111111111111D0*t85*t216
      t227 = t70*t219-0.1333333333333333D1*rho*sigma
     #+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho*sigmaaa
      t230 = 0.869D-2*t64*t69*t227
      vrhoa(i) = t139*t11+0.1333333333333333D1*t26*t5-0.22D0*t143*t58
     #-0.869D-2*t64*t69*t156-t167+t171-t178-t182+t187-t230
      t231 = t28*rhob
      t233 = 1/t30/t231
      t234 = sigmabb*t233
      t237 = t41**2
      t238 = 1/t237
      t240 = 1/t29/t28
      t242 = t34*t240*t38
      t244 = 1.D0+t33
      t245 = dsqrt(t244)
      t246 = 1/t245
      t247 = t234*t246
      t249 = -0.28D-1*t242-0.28D-1*t247
      t250 = t238*t249
      t255 = t46**2
      t256 = 1/t255
      t259 = -0.336D-1*t242-0.336D-1*t247
      t260 = t256*t259
      t263 = 0.9708253333333333D-1*t234*t42+0.3640595D-1*t33*t250
     #-0.94616816D-1*t234*t47-0.35481306D-1*t33*t260
      t270 = t30*rhob
      t274 = 0.9723306394337274D2*t270-0.1111111111111111D0*t149*sigmabb
      t278 = rhoa*t92+t70*t274-2.D0*rhob*sigmaaa
      vrhob(i) = t263*t35+0.1333333333333333D1*t50*t29-0.22D0*t57*t58
     #-0.869D-2*t64*t69*t278-t167+t171-t178-t182+t187-t230
      t284 = 1/t10
      t286 = t284*t12*t14
      t288 = t8*t122
      t290 = 0.105D-1*t286+0.105D-1*t288
      t291 = t114*t290
      t298 = 0.126D-1*t286+0.126D-1*t288
      t299 = t132*t298
      t302 = -0.3640595D-1*t8*t18+0.3640595D-1*t9*t291+0.35481306D-1
     #*t8*t23-0.35481306D-1*t9*t299
      t304 = 0.1388888888888889D-1*t53
      t305 = 0.1947222222222222D-1*t74
      t306 = t85*rhoa
      t309 = -0.25D1+t304+t305-0.1111111111111111D0*t306*t58
      t311 = t70*t309+t96-t100
      t317 = t70*t76-0.6666666666666667D0*t65
      t319 = t64*t69*t317
      t320 = 0.869D-2*t319
      vsigmaaa(i) = t302*t11-0.869D-2*t64*t69*t311-t320
      vsigmaab(i) = -0.1738D-1*t319
      t324 = 1/t34
      t326 = t324*t36*t38
      t328 = t32*t246
      t330 = 0.105D-1*t326+0.105D-1*t328
      t331 = t238*t330
      t338 = 0.126D-1*t326+0.126D-1*t328
      t339 = t256*t338
      t342 = -0.3640595D-1*t32*t42+0.3640595D-1*t33*t331+0.35481306D-1
     #*t32*t47-0.35481306D-1*t33*t339
      t344 = t85*rhob
      t347 = -0.25D1+t304+t305-0.1111111111111111D0*t344*t58
      t349 = t70*t347+t96-t97
      vsigmabb(i) = t342*t35-0.869D-2*t64*t69*t349-t320
      t362 = t4**2
      t365 = sigmaaa/t6/t362
      t371 = 1/t113/t17
      t372 = t125**2
      t379 = t10/t5/t107*t14
      t381 = t365*t122
      t383 = sigmaaa**2
      t389 = 1/t121/t120
      t390 = t383/t5/t362/t107*t389
      t401 = 1/t131/t22
      t402 = t135**2
      t416 = 1/t160/t55
      t418 = rhob*t69
      t420 = 0.6006012222222222D-2*t416*rhoa*t418
      t423 = 0.1448333333333333D-2*t175*t56*t227
      t424 = t172*t65
      t426 = 1/t52/t424
      t427 = t426*t63
      t429 = 0.6034722222222222D-4*t427*t176
      t432 = 0.1692136111111111D-3*t427*t161*t103
      t435 = 0.2030563333333333D-2*t179*t174*t227
      t439 = 0.2372374827777778D-3*t63*t416*t426*t103
      t442 = 0.6372666666666667D-1*t64*t184*t227
      t444 = t164*t56
      t448 = 1/t67/t65*t161
      t450 = 1/t66
      t451 = t450*t416
      t474 = rhob*t450
      t488 = 0.869D-2*t64*t69*(t70*((-0.4320987654320988D-1*t164
     #-0.6058024691358025D-1*t444+0.3185006481481481D-1*t448
     #-0.3721149239197531D-2*t451)*sigma-1.D0*(-0.617283950617284D-2
     #*t164-0.8654320987654321D-2*t444+0.4550009259259259D-2*t448
     #-0.5315927484567901D-3*t451)*t81-0.1111111111111111D0*
     #(0.1111111111111111D0*t164+0.1557777777777778D0*t444
     #-0.8190016666666667D-1*t448+0.9568669472222222D-2*t451)*t89
     #-0.2222222222222222D0*t208*t216-0.1111111111111111D0*t85*(2.D0
     #*rhoa*t450*sigmaaa+2.D0*t474*sigmabb))-0.1333333333333333D1
     #*sigma+0.1333333333333333D1*sigmabb+0.1333333333333333D1*sigmaaa)
      t490 = t208*t58
      t493 = t85*t168
      t500 = t64*t69*(rhob*t219+t70*(-0.1111111111111111D0*t490
     #*sigmaaa+0.1111111111111111D0*t493*sigmaaa))
      t502 = -0.869D-2*t64*t69*(2.D0*rhob*t152+0.1620551065722879D3
     #*t147*rhob-2.D0*sigmabb)+(-0.3559692888888889D0*t365*t18
     #-0.1941650666666667D0*t110*t126-0.728119D-1*t9*t371*t372
     #+0.3640595D-1*t9*t114*(0.6533333333333333D-1*t379+0.14D0*t381
     #-0.3733333333333333D-1*t390)+0.3469283253333333D0*t365*t23
     #+0.189233632D0*t110*t136+0.70962612D-1*t9*t401*t402
     #-0.35481306D-1*t9*t132*(0.784D-1*t379+0.168D0*t381-0.448D-1*t390
     #))*t11-t420-t423-t429-t432-t435-t439+t442-t488-0.1738D-1*t500
      t504 = t179*t174*t156
      t507 = t64*t184*t156
      t510 = t161*rhob*t164
      t512 = t143*t168
      t515 = t175*t56*t156
      t522 = 1/t424
      t525 = 0.8799107777777778D-2*t179*t522*t103
      t527 = 0.44D0*t57*t474
      t530 = 0.6276111111111111D-2*t522*t63*t176
      t535 = 0.8567777777777778D-1*t162*rhob/t52/t66
      t540 = 0.1486955555555556D0*t64/t67/t173*t103
      t541 = -0.2030563333333333D-2*t504+0.6372666666666667D-1*t507
     #-0.5140666666666667D-1*t510+0.44D0*t512-0.1448333333333333D-2
     #*t515+0.2666666666666667D1*t139*t5+0.4444444444444444D0*t26/t6
     #+t525-t527+t530+t535-t540
      v2rhoa2(i) = t502+t541
      t542 = -t420-t423-t429-t432-t435-t439+t442-t488+t525-t527+t530
      t557 = t162*t164
      t559 = t57*t168
      t562 = t175*t56*t278
      t565 = t179*t174*t278
      t567 = t28**2
      t570 = sigmabb/t30/t567
      t576 = 1/t237/t41
      t577 = t249**2
      t584 = t34/t29/t231*t38
      t586 = t570*t246
      t588 = sigmabb**2
      t594 = 1/t245/t244
      t595 = t588/t29/t567/t231*t594
      t606 = 1/t255/t46
      t607 = t259**2
      t621 = t64*t184*t278
      t632 = t64*t69*(rhoa*t219+t70*(-0.1111111111111111D0*t490
     #*sigmabb+0.1111111111111111D0*t493*sigmabb))
      t634 = t535+0.2666666666666667D1*t263*t29+0.4444444444444444D0
     #*t50/t30-0.869D-2*t64*t69*(2.D0*rhoa*t274+0.1620551065722879D3
     #*rhoa*t270-2.D0*sigmaaa)-0.5140666666666667D-1*t557+0.44D0*t559
     #-0.1448333333333333D-2*t562-0.2030563333333333D-2*t565+(
     #-0.3559692888888889D0*t570*t42-0.1941650666666667D0*t234*t250
     #-0.728119D-1*t33*t576*t577+0.3640595D-1*t33*t238*
     #(0.6533333333333333D-1*t584+0.14D0*t586-0.3733333333333333D-1
     #*t595)+0.3469283253333333D0*t570*t47+0.189233632D0*t234*t260
     #+0.70962612D-1*t33*t606*t607-0.35481306D-1*t33*t256*(0.784D-1
     #*t584+0.168D0*t586-0.448D-1*t595))*t35-t540
     #+0.6372666666666667D-1*t621-0.1738D-1*t632
      v2rhob2(i) = t542+t634
      t640 = -t420-t423-t435-t439+t442-t429-t432-t488-0.869D-2*t500
     #-0.1015281666666667D-2*t504+0.3186333333333333D-1*t507
     #-0.7241666666666667D-3*t515-0.2570333333333333D-1*t510
      t656 = 0.22D0*t512+t525-t527+t530+t535-t540
     #-0.2570333333333333D-1*t557+0.22D0*t559-0.7241666666666667D-3
     #*t562-0.1015281666666667D-2*t565+0.3186333333333333D-1*t621
     #-0.869D-2*t632-0.869D-2*t64*t69*(t71+t72+t77-t83-t91+rhob*t274
     #+rhoa*t152)-0.22D0*t56*t58
      v2rhoab(i) = t640+t656
      t661 = t8*t114
      t669 = t284*t116*t14
      t671 = t109*t122
      t677 = sigmaaa/t5/t362/t4*t389
      t687 = t8*t132
      t707 = 0.1111111111111111D0*t70*t149
      t714 = 0.7241666666666667D-3*t175*t56*t311
      t717 = 0.1015281666666667D-2*t179*t174*t311
      t720 = 0.3186333333333333D-1*t64*t184*t311
      t721 = 0.462962962962963D-2*t189
      t722 = 0.6490740740740741D-2*t191
      t723 = 0.227500462962963D-2*t195
      t731 = 0.1333333333333333D1*rho
      t735 = 0.869D-2*t64*t69*(t70*(-t721-t722+t723
     #-0.1111111111111111D0*t208*rhoa*t58+0.1111111111111111D0*t306
     #*t168)+t731)
      t737 = t64*t418*t76
      t738 = 0.869D-2*t737
      t740 = t175*t56*t317
      t741 = 0.7241666666666667D-3*t740
      t743 = t179*t174*t317
      t744 = 0.1015281666666667D-2*t743
      t746 = t64*t184*t317
      t747 = 0.3186333333333333D-1*t746
      t752 = t64*t69*(t70*t197-0.1333333333333333D1*rho)
      t753 = 0.869D-2*t752
      v2rhoasigmaaa(i) = (0.9708253333333333D-1*t109*t18
     #-0.9708253333333333D-1*t110*t291+0.3640595D-1*t661*t125
     #-0.728119D-1*t9*t371*t125*t290+0.3640595D-1*t9*t114*(-0.14D-1
     #*t669-0.42D-1*t671+0.14D-1*t677)-0.94616816D-1*t109*t23
     #+0.94616816D-1*t110*t299-0.35481306D-1*t687*t135+0.70962612D-1
     #*t9*t401*t135*t298-0.35481306D-1*t9*t132*(-0.168D-1*t669
     #-0.504D-1*t671+0.168D-1*t677))*t11+0.1333333333333333D1*t302*t5
     #-0.869D-2*t64*t69*(rhob*t309-t707)-t714-t717+t720-t735-t738-t741
     #-t744+t747-t753
      t755 = 0.1448333333333333D-2*t740
      t756 = 0.2030563333333333D-2*t743
      t757 = 0.6372666666666667D-1*t746
      t758 = 0.1738D-1*t752
      v2rhoasigmaab(i) = -0.1738D-1*t737-t755-t756+t757-t758
      t767 = 0.7241666666666667D-3*t175*t56*t349
      t770 = 0.1015281666666667D-2*t179*t174*t349
      t773 = 0.3186333333333333D-1*t64*t184*t349
      t784 = 0.869D-2*t64*t69*(t70*(-t721-t722+t723
     #-0.1111111111111111D0*t208*rhob*t58+0.1111111111111111D0*t344
     #*t168)+t731)
      v2rhoasigmabb(i) = -0.869D-2*t64*t69*(rhob*t347-2.D0*rhoa)-t767
     #-t770+t773-t784-t738-t741-t744+t747-t753
      t793 = t64*t69*rhoa*t76
      t794 = 0.869D-2*t793
      v2rhobsigmaaa(i) = -0.869D-2*t64*t69*(rhoa*t309-2.D0*rhob)-t714
     #-t717+t720-t735-t794-t741-t744+t747-t753
      v2rhobsigmaab(i) = -0.1738D-1*t793-t755-t756+t757-t758
      t800 = t32*t238
      t808 = t324*t240*t38
      t810 = t233*t246
      t816 = sigmabb/t29/t567/t28*t594
      t826 = t32*t256
      v2rhobsigmabb(i) = (0.9708253333333333D-1*t233*t42
     #-0.9708253333333333D-1*t234*t331+0.3640595D-1*t800*t249
     #-0.728119D-1*t33*t576*t249*t330+0.3640595D-1*t33*t238*(-0.14D-1
     #*t808-0.42D-1*t810+0.14D-1*t816)-0.94616816D-1*t233*t47
     #+0.94616816D-1*t234*t339-0.35481306D-1*t826*t259+0.70962612D-1
     #*t33*t606*t259*t338-0.35481306D-1*t33*t256*(-0.168D-1*t808
     #-0.504D-1*t810+0.168D-1*t816))*t35+0.1333333333333333D1*t342*t29
     #-0.869D-2*t64*t69*(rhoa*t347-t707)-t767-t770+t773-t784-t794-t741
     #-t744+t747-t753
      t851 = t290**2
      t858 = 1/t10/sigmaaa*t12*t14
      t862 = 1/sigmaaa*t8*t122
      t867 = 1/t5/t362/rhoa*t389
      t875 = t298**2
      v2sigmaaa2(i) = (0.728119D-1*t661*t290-0.728119D-1*t9*t371*t851
     #+0.3640595D-1*t9*t114*(-0.525D-2*t858+0.525D-2*t862-0.525D-2
     #*t867)-0.70962612D-1*t687*t298+0.70962612D-1*t9*t401*t875
     #-0.35481306D-1*t9*t132*(-0.63D-2*t858+0.63D-2*t862-0.63D-2*t867)
     #)*t11
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t889 = t330**2
      t896 = 1/t34/sigmabb*t36*t38
      t900 = 1/sigmabb*t32*t246
      t905 = 1/t29/t567/rhob*t594
      t913 = t338**2
      v2sigmabb2(i) = (0.728119D-1*t800*t330-0.728119D-1*t33*t576*t889
     #+0.3640595D-1*t33*t238*(-0.525D-2*t896+0.525D-2*t900-0.525D-2
     #*t905)-0.70962612D-1*t826*t338+0.70962612D-1*t33*t606*t913
     #-0.35481306D-1*t33*t256*(-0.63D-2*t896+0.63D-2*t900-0.63D-2*t905
     #))*t35
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
      
      
      subroutine rks_xc_edf1
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     R.D. Adamson, P.M.W. Gill, and J.A. Pople
c     Empirical density functionals
c     Chem. Phys. Lett. 284 (1998) 6-11
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
      t2 = rho**2
      t3 = rho**(1.D0/3.D0)
      t4 = t3**2
      t5 = t4*t2
      t7 = sigma/t5
      t8 = dsqrt(sigma)
      t9 = t3*rho
      t11 = t8/t9
      t13 = dlog(0.1259921049894873D1*t11+dsqrt(1+0.1587401051968199D1
     #*t11**2))
      t14 = t11*t13
      t28 = 1/t3
      t31 = 1/(1.D0+0.3505D0*t28)
      t34 = 0.25D0*t28
      t35 = dexp(-t34)
      t42 = t28*t31
      zk(i) = 0.7937005259840997D0*(-0.9593273689405774D0
     #-0.5779084332790167D-1*t7/(1.D0+0.2645834204779234D-1*t14)
     #+0.5632306246960559D-1*t7/(1.D0+0.317500104573508D-1*t14))*t9
     #-0.55D-1*t31*rho-0.869D-2*t35*t31/t4/t2/rho*(0.25D0*t2*
     #(0.1148493600075277D2*t5+(0.2611111111111111D1
     #-0.9722222222222222D-1*t28-0.1363055555555556D0*t42)*sigma-0.5D0
     #*(0.25D1-0.1388888888888889D-1*t28-0.1947222222222222D-1*t42)
     #*sigma-0.2777777777777778D-1*(t34+0.3505D0*t42-11.D0)*sigma)
     #-0.4583333333333333D0*t2*sigma)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**2
      t3 = rho**(1.D0/3.D0)
      t4 = t3**2
      t5 = t4*t2
      t6 = 1/t5
      t7 = sigma*t6
      t8 = dsqrt(sigma)
      t9 = t3*rho
      t10 = 1/t9
      t11 = t8*t10
      t13 = dlog(0.1259921049894873D1*t11+dsqrt(1+0.1587401051968199D1
     #*t11**2))
      t14 = t11*t13
      t16 = 1.D0+0.2645834204779234D-1*t14
      t17 = 1/t16
      t21 = 1.D0+0.317500104573508D-1*t14
      t22 = 1/t21
      t25 = -0.9593273689405774D0-0.5779084332790167D-1*t7*t17
     #+0.5632306246960559D-1*t7*t22
      t28 = 1/t3
      t30 = 1.D0+0.3505D0*t28
      t31 = 1/t30
      t34 = 0.25D0*t28
      t35 = dexp(-t34)
      t36 = t35*t31
      t39 = 1/t4/t2/rho
      t42 = t28*t31
      t44 = 0.2611111111111111D1-0.9722222222222222D-1*t28
     #-0.1363055555555556D0*t42
      t52 = t34+0.3505D0*t42-11.D0
      t55 = 0.1148493600075277D2*t5+t44*sigma-0.5D0*(0.25D1
     #-0.1388888888888889D-1*t28-0.1947222222222222D-1*t42)*sigma
     #-0.2777777777777778D-1*t52*sigma
      t60 = 0.25D0*t2*t55-0.4583333333333333D0*t2*sigma
      zk(i) = 0.7937005259840997D0*t25*t9-0.55D-1*t31*rho-0.869D-2*t36
     #*t39*t60
      t64 = sigma*t39
      t67 = t16**2
      t68 = 1/t67
      t72 = t8/t3/t2*t13
      t76 = dsqrt(1.D0+0.1587401051968199D1*t7)
      t77 = 1/t76
      t78 = t64*t77
      t86 = t21**2
      t87 = 1/t86
      t102 = t4*rho
      t106 = t52/rho*sigma
      t111 = rho*sigma
      t117 = t30**2
      t118 = 1/t117
      t121 = t2**2
      t123 = 1/t121/rho
      t138 = t10*t31
      t140 = 1/t102
      t141 = t140*t118
      s1 = 0.3968502629920499D0*(0.3082178310821422D0*t64*t17
     #+0.5779084332790167D-1*t7*t68*(-0.705555787941129D-1*t72
     #-0.8889445891021917D-1*t78)-0.3003896665045631D0*t64*t22
     #-0.5632306246960559D-1*t7*t87*(-0.8466669455293548D-1*t72
     #-0.106673350692263D0*t78))*t9+0.10582673679788D1*t25*t3-0.55D-1
     #*t31-0.869D-2*t36*t39*(0.5D0*rho*t55+0.25D0*t2*
     #(0.3062649600200738D2*t102-0.2777777777777778D-1*t106)-0.25D0
     #*t111)
      vrhoa(i) = s1-0.6425833333333333D-2*t118*t28
     #-0.7241666666666667D-3*t123*t35*t31*t60-0.1015281666666667D-2
     #*t35*t118*t123*t60+0.3186333333333333D-1*t36/t4/t121*t60
     #-0.869D-2*t36*t39*(0.25D0*t2*((0.3240740740740741D-1*t10
     #+0.4543518518518519D-1*t138-0.1592503240740741D-1*t141)*sigma
     #-0.5D0*(0.462962962962963D-2*t10+0.6490740740740741D-2*t138
     #-0.227500462962963D-2*t141)*sigma-0.2777777777777778D-1*(
     #-0.8333333333333333D-1*t10-0.1168333333333333D0*t138
     #+0.4095008333333333D-1*t141)*sigma+0.2777777777777778D-1*t106)
     #-0.6666666666666667D0*t111)
      t170 = 1/t8*t10*t13
      t172 = t6*t77
      vsigmaaa(i) = 0.7937005259840997D0*(-0.2311633733116067D0*t6*t17
     #+0.5779084332790167D-1*t7*t68*(0.5291668409558467D-1*t170
     #+0.6667084418266438D-1*t172)+0.2252922498784224D0*t6*t22
     #-0.5632306246960559D-1*t7*t87*(0.6350002091470161D-1*t170
     #+0.8000501301919725D-1*t172))*t9+0.9655555555555556D-3*t36*t140
     #-0.3476D-1*t36*t39*(0.25D0*t2*t44-0.6666666666666667D0*t2)
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
      t2 = rho**2
      t3 = rho**(1.D0/3.D0)
      t4 = t3**2
      t5 = t4*t2
      t6 = 1/t5
      t7 = sigma*t6
      t8 = dsqrt(sigma)
      t9 = t3*rho
      t10 = 1/t9
      t11 = t8*t10
      t13 = dlog(0.1259921049894873D1*t11+dsqrt(1+0.1587401051968199D1
     #*t11**2))
      t14 = t11*t13
      t16 = 1.D0+0.2645834204779234D-1*t14
      t17 = 1/t16
      t21 = 1.D0+0.317500104573508D-1*t14
      t22 = 1/t21
      t25 = -0.9593273689405774D0-0.5779084332790167D-1*t7*t17
     #+0.5632306246960559D-1*t7*t22
      t28 = 1/t3
      t30 = 1.D0+0.3505D0*t28
      t31 = 1/t30
      t34 = 0.25D0*t28
      t35 = dexp(-t34)
      t36 = t35*t31
      t37 = t2*rho
      t39 = 1/t4/t37
      t40 = 0.1148493600075277D2*t5
      t42 = t28*t31
      t44 = 0.2611111111111111D1-0.9722222222222222D-1*t28
     #-0.1363055555555556D0*t42
      t45 = t44*sigma
      t50 = 0.5D0*(0.25D1-0.1388888888888889D-1*t28
     #-0.1947222222222222D-1*t42)*sigma
      t52 = t34+0.3505D0*t42-11.D0
      t54 = 0.2777777777777778D-1*t52*sigma
      t55 = t40+t45-t50-t54
      t60 = 0.25D0*t2*t55-0.4583333333333333D0*t2*sigma
      zk(i) = 0.7937005259840997D0*t25*t9-0.55D-1*t31*rho-0.869D-2*t36
     #*t39*t60
      t64 = sigma*t39
      t67 = t16**2
      t68 = 1/t67
      t70 = 1/t3/t2
      t72 = t8*t70*t13
      t75 = 1.D0+0.1587401051968199D1*t7
      t76 = dsqrt(t75)
      t77 = 1/t76
      t78 = t64*t77
      t80 = -0.705555787941129D-1*t72-0.8889445891021917D-1*t78
      t81 = t68*t80
      t86 = t21**2
      t87 = 1/t86
      t90 = -0.8466669455293548D-1*t72-0.106673350692263D0*t78
      t91 = t87*t90
      t94 = 0.3082178310821422D0*t64*t17+0.5779084332790167D-1*t7*t81
     #-0.3003896665045631D0*t64*t22-0.5632306246960559D-1*t7*t91
      t102 = t4*rho
      t104 = 1/rho
      t105 = t52*t104
      t106 = t105*sigma
      t108 = 0.3062649600200738D2*t102-0.2777777777777778D-1*t106
      t111 = rho*sigma
      t113 = 0.5D0*rho*t55+0.25D0*t2*t108-0.25D0*t111
      t117 = t30**2
      t118 = 1/t117
      t121 = t2**2
      t122 = t121*rho
      t123 = 1/t122
      t124 = t123*t35
      t125 = t31*t60
      t128 = t35*t118
      t133 = 1/t4/t121
      t138 = t10*t31
      t140 = 1/t102
      t141 = t140*t118
      t143 = 0.3240740740740741D-1*t10+0.4543518518518519D-1*t138
     #-0.1592503240740741D-1*t141
      t154 = -0.8333333333333333D-1*t10-0.1168333333333333D0*t138
     #+0.4095008333333333D-1*t141
      t158 = t143*sigma-0.5D0*(0.462962962962963D-2*t10
     #+0.6490740740740741D-2*t138-0.227500462962963D-2*t141)*sigma
     #-0.2777777777777778D-1*t154*sigma+0.2777777777777778D-1*t106
      t162 = 0.25D0*t2*t158-0.6666666666666667D0*t111
      vrhoa(i) = 0.3968502629920499D0*t94*t9+0.10582673679788D1*t25*t3
     #-0.55D-1*t31-0.869D-2*t36*t39*t113-0.6425833333333333D-2*t118
     #*t28-0.7241666666666667D-3*t124*t125-0.1015281666666667D-2*t128
     #*t123*t60+0.3186333333333333D-1*t36*t133*t60-0.869D-2*t36*t39*t162
      t168 = 1/t8
      t170 = t168*t10*t13
      t172 = t6*t77
      t174 = 0.5291668409558467D-1*t170+0.6667084418266438D-1*t172
      t175 = t68*t174
      t182 = 0.6350002091470161D-1*t170+0.8000501301919725D-1*t172
      t183 = t87*t182
      t186 = -0.2311633733116067D0*t6*t17+0.5779084332790167D-1*t7
     #*t175+0.2252922498784224D0*t6*t22-0.5632306246960559D-1*t7*t183
      t194 = 0.25D0*t2*t44-0.6666666666666667D0*t2
      vsigmaaa(i) = 0.7937005259840997D0*t186*t9+0.9655555555555556D-3
     #*t36*t140-0.3476D-1*t36*t39*t194
      t205 = t70*t31
      t207 = t6*t118
      t209 = 1/t37
      t211 = 1/t117/t30
      t212 = t209*t211
      t231 = t154*t104*sigma
      t235 = t52/t2*sigma
      t261 = sigma*t133
      t267 = 1/t67/t16
      t268 = t80**2
      t275 = t8/t3/t37*t13
      t277 = t261*t77
      t279 = sigma**2
      t285 = 1/t76/t75
      t286 = t279/t3/t121/t37*t285
      t297 = 1/t86/t21
      t298 = t90**2
      t321 = t121*t2
      t323 = 1/t3/t321
      t333 = 1/t321
      t345 = rho*t108
      t359 = t323*t35
      t365 = -0.4744749655555556D-3*t35*t211*t323*t60
     #+0.1274533333333333D0*t36*t133*t113-0.2896666666666667D-2*t124
     #*t31*t113+0.1759821555555556D-1*t128*t333*t60
     #+0.1255222222222222D-1*t333*t35*t125-0.2973911111111111D0*t36/t4
     #/t122*t60-0.869D-2*t36*t39*(t40+t45-t50-t54+t345)-0.869D-2*t36
     #*t39*(t345+0.2552208000167282D2*t5-0.5D0*sigma)
     #-0.2896666666666667D-2*t124*t31*t162-0.1206944444444444D-3*t359
     #*t125-0.3384272222222222D-3*t359*t118*t60
      s1 = -0.4061126666666667D-2*t128*t123*t162+0.1274533333333333D0
     #*t36*t133*t162-0.1738D-1*t36*t39*(0.25D0*t2*((
     #-0.4320987654320988D-1*t70-0.6058024691358025D-1*t205
     #+0.3185006481481481D-1*t207-0.3721149239197531D-2*t212)*sigma
     #-0.5D0*(-0.617283950617284D-2*t70-0.8654320987654321D-2*t205
     #+0.4550009259259259D-2*t207-0.5315927484567901D-3*t212)*sigma
     #-0.2777777777777778D-1*(0.1111111111111111D0*t70
     #+0.1557777777777778D0*t205-0.8190016666666667D-1*t207
     #+0.9568669472222222D-2*t212)*sigma+0.5555555555555556D-1*t231
     #-0.5555555555555556D-1*t235)-0.6666666666666667D0*sigma)
     #-0.3476D-1*t36*t39*(0.5D0*rho*t158+0.25D0*t2*(
     #-0.2777777777777778D-1*t231+0.2777777777777778D-1*t235))
     #-0.4061126666666667D-2*t128*t123*t113
      v2rhoa2(i) = s1-0.3003006111111111D-2*t211*t140
     #+0.3968502629920499D0*(-0.2260264094602376D1*t261*t17
     #-0.6164356621642845D0*t64*t81-0.1155816866558033D0*t7*t267*t268
     #+0.5779084332790167D-1*t7*t68*(0.3292593677058602D0*t275
     #+0.8889445891021917D0*t277-0.3762964202352688D0*t286)
     #+0.2202857554366796D1*t261*t22+0.6007793330091263D0*t64*t91
     #+0.1126461249392112D0*t7*t297*t298-0.5632306246960559D-1*t7*t87*
     #(0.3951112412470322D0*t275+0.106673350692263D1*t277
     #-0.4515557042823225D0*t286))*t9+0.2116534735957599D1*t94*t3
     #+0.7055115786525331D0*t25/t4-0.8567777777777778D-2*t118*t10+t365
      t370 = t6*t68
      t378 = t168*t70*t13
      t380 = t39*t77
      t383 = sigma*t323*t285
      t393 = t6*t87
      s1 = 0.3968502629920499D0*(0.1232871324328569D1*t39*t17
     #-0.3082178310821422D0*t64*t175+0.2311633733116067D0*t370*t80
     #-0.1155816866558033D0*t7*t267*t80*t174+0.5779084332790167D-1*t7
     #*t68*(-0.1411111575882258D0*t378-0.533366753461315D0*t380
     #+0.2822223151764516D0*t383)-0.1201558666018253D1*t39*t22
     #+0.3003896665045631D0*t64*t183-0.2252922498784224D0*t393*t90
     #+0.1126461249392112D0*t7*t297*t90*t182-0.5632306246960559D-1*t7
     #*t87*(-0.169333389105871D0*t378-0.640040104153578D0*t380
     #+0.3386667782117419D0*t383))*t9+0.10582673679788D1*t186*t3
     #-0.869D-2*t36*t39*(-0.9444444444444444D0*rho
     #-0.2777777777777778D-1*rho*t52)+0.8046296296296296D-4*t209*t35
     #*t31+0.1128090740740741D-3*t128*t209+0.1335685185185185D-1*t36*t6
      v2rhoasigmaaa(i) = s1-0.1738D-1*t36*t39*(0.25D0*t2*(-0.1D-22*t10
     #-0.1D-22*t138+0.5555555555555556D-1*t105)+0.1333333333333333D1
     #*rho)-0.1738D-1*t36*t6*t44-0.2896666666666667D-2*t124*t31*t194
     #-0.4061126666666667D-2*t128*t123*t194+0.1274533333333333D0*t36
     #*t133*t194-0.3476D-1*t36*t39*(0.25D0*t2*t143
     #-0.1333333333333333D1*rho)
      t458 = t174**2
      t465 = 1/t8/sigma*t10*t13
      t469 = 1/sigma*t6*t77
      t473 = 1/t3/t122*t285
      t481 = t182**2
      v2sigmaaa2(i) = 0.7937005259840997D0*(0.4623267466232134D0*t370
     #*t174-0.1155816866558033D0*t7*t267*t458+0.5779084332790167D-1*t7
     #*t68*(-0.1058333681911693D0*t465+0.1333416883653288D0*t469
     #-0.2116667363823387D0*t473)-0.4505844997568447D0*t393*t182
     #+0.1126461249392112D0*t7*t297*t481-0.5632306246960559D-1*t7*t87*
     #(-0.1270000418294032D0*t465+0.1600100260383945D0*t469
     #-0.2540000836588064D0*t473))*t9
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
c:XC_EDF1subrend
