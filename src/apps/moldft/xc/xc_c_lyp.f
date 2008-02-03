c:C_LYPsubrstart

c    Generated: Wed Mar 10 10:16:01 GMT 2004

      subroutine uks_c_lyp
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     C. Lee, W. Yang, and R.G. Parr
c     Development of the Colle-Salvetti correlation-energy formula into
c     a functional of the electron density
c     Phys. Rev. B37 (1988) 785-789
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
      zk(i) = 0.D0
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      zk(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rho**(1.D0/3.D0)
      t5 = 1/t4
      t8 = 1/(1.D0+0.349D0*t5)
      t10 = 1/rho
      t11 = rhob*t10
      t14 = 0.2533D0*t5
      t15 = dexp(-t14)
      t17 = rho**2
      t19 = t4**2
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t34 = t5*t8
      t56 = 0.6666666666666667D0*t17
      zk(i) = -0.19672D0*t8*rhoa*t11-0.649176D-2*t15*t8/t19/t17/rho*
     #(rhoa*rhob*(0.3646239897876478D2*t25*t23+0.3646239897876478D2
     #*t30*t28+(0.2611111111111111D1-0.9850555555555556D-1*t5
     #-0.1357222222222222D0*t34)*sigma-1.D0*(0.25D1
     #-0.1407222222222222D-1*t5-0.1938888888888889D-1*t34)*(sigmaaa
     #+sigmabb)-0.1111111111111111D0*(t14+0.349D0*t34-11.D0)*(rhoa*t10
     #*sigmaaa+t11*sigmabb))-0.6666666666666667D0*t17*sigma+(t56-1.D0
     #*t23)*sigmabb+(t56-1.D0*t28)*sigmaaa)
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
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rho**(1.D0/3.D0)
      t5 = 1/t4
      t7 = 1.D0+0.349D0*t5
      t8 = 1/t7
      t9 = t8*rhoa
      t10 = 1/rho
      t11 = rhob*t10
      t14 = 0.2533D0*t5
      t15 = dexp(-t14)
      t16 = t15*t8
      t17 = rho**2
      t19 = t4**2
      t21 = 1/t19/t17/rho
      t22 = rhoa*rhob
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t34 = t5*t8
      t36 = 0.2611111111111111D1-0.9850555555555556D-1*t5
     #-0.1357222222222222D0*t34
      t41 = sigmaaa+sigmabb
      t45 = t14+0.349D0*t34-11.D0
      t49 = rhoa*t10*sigmaaa+t11*sigmabb
      t52 = 0.3646239897876478D2*t25*t23+0.3646239897876478D2*t30*t28
     #+t36*sigma-1.D0*(0.25D1-0.1407222222222222D-1*t5
     #-0.1938888888888889D-1*t34)*t41-0.1111111111111111D0*t45*t49
      t56 = 0.6666666666666667D0*t17
      t57 = 1.D0*t23
      t60 = 1.D0*t28
      t63 = t22*t52-0.6666666666666667D0*t17*sigma+(t56-t57)*sigmabb+
     #(t56-t60)*sigmaaa
      zk(i) = -0.19672D0*t9*t11-0.649176D-2*t16*t21*t63
      t73 = t45*t10
      t84 = t7**2
      t85 = 1/t84
      t91 = 0.2288509333333333D-1*t85*rhoa*rhob/t4/t17
      t92 = 1/t17
      t93 = rhob*t92
      t95 = 0.19672D0*t9*t93
      t96 = t17**2
      t98 = 1/t96/rho
      t102 = 0.548120936D-3*t98*t15*t8*t63
      t106 = 0.75520808D-3*t15*t85*t98*t63
      t111 = 0.2380312D-1*t16/t19/t96*t63
      t113 = 1/t4/rho
      t115 = t113*t8
      t119 = 1/t19/rho*t85
      t154 = 0.649176D-2*t16*t21*(t22*((0.3283518518518519D-1*t113
     #+0.4524074074074074D-1*t115-0.1578901851851852D-1*t119)*sigma
     #-1.D0*(0.4690740740740741D-2*t113+0.6462962962962963D-2*t115
     #-0.2255574074074074D-2*t119)*t41-0.1111111111111111D0*(
     #-0.8443333333333333D-1*t113-0.1163333333333333D0*t115
     #+0.4060033333333333D-1*t119)*t49-0.1111111111111111D0*t45*(-1.D0
     #*rhoa*t92*sigmaaa-1.D0*t93*sigmabb))-0.1333333333333333D1*rho
     #*sigma+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho
     #*sigmaaa)
      vrhoa(i) = -0.19672D0*t8*rhob*t10-0.649176D-2*t16*t21*(rhob*t52
     #+t22*(0.9723306394337274D2*t25*rhoa-0.1111111111111111D0*t73
     #*sigmaaa)-2.D0*rhoa*sigmabb)-t91+t95-t102-t106+t111-t154
      vrhob(i) = -0.19672D0*t9*t10-0.649176D-2*t16*t21*(rhoa*t52+t22*
     #(0.9723306394337274D2*t30*rhob-0.1111111111111111D0*t73*sigmabb)
     #-2.D0*rhob*sigmaaa)-t91+t95-t102-t106+t111-t154
      t170 = 0.1407222222222222D-1*t5
      t171 = 0.1938888888888889D-1*t34
      t185 = t16*t21*(t22*t36-0.6666666666666667D0*t17)
      t186 = 0.649176D-2*t185
      vsigmaaa(i) = -0.649176D-2*t16*t21*(t22*(-0.25D1+t170+t171
     #-0.1111111111111111D0*t45*rhoa*t10)+t56-t60)-t186
      vsigmaab(i) = -0.1298352D-1*t185
      vsigmabb(i) = -0.649176D-2*t16*t21*(t22*(-0.25D1+t170+t171
     #-0.1111111111111111D0*t45*rhob*t10)+t56-t57)-t186
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
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2rhob2(i) = 0.D0
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      v2rhoa2(i) = 0.D0
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2sigmaaa2(i) = 0.D0
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
      t4 = rho**(1.D0/3.D0)
      t5 = 1/t4
      t7 = 1.D0+0.349D0*t5
      t8 = 1/t7
      t9 = t8*rhoa
      t10 = 1/rho
      t11 = rhob*t10
      t14 = 0.2533D0*t5
      t15 = dexp(-t14)
      t16 = t15*t8
      t17 = rho**2
      t18 = t17*rho
      t19 = t4**2
      t21 = 1/t19/t18
      t22 = rhoa*rhob
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t27 = 0.3646239897876478D2*t25*t23
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t32 = 0.3646239897876478D2*t30*t28
      t34 = t5*t8
      t36 = 0.2611111111111111D1-0.9850555555555556D-1*t5
     #-0.1357222222222222D0*t34
      t37 = t36*sigma
      t41 = sigmaaa+sigmabb
      t43 = 1.D0*(0.25D1-0.1407222222222222D-1*t5
     #-0.1938888888888889D-1*t34)*t41
      t45 = t14+0.349D0*t34-11.D0
      t49 = rhoa*t10*sigmaaa+t11*sigmabb
      t51 = 0.1111111111111111D0*t45*t49
      t52 = t27+t32+t37-t43-t51
      t56 = 0.6666666666666667D0*t17
      t57 = 1.D0*t23
      t60 = 1.D0*t28
      t63 = t22*t52-0.6666666666666667D0*t17*sigma+(t56-t57)*sigmabb+
     #(t56-t60)*sigmaaa
      zk(i) = -0.19672D0*t9*t11-0.649176D-2*t16*t21*t63
      t67 = t8*rhob
      t71 = t25*rhoa
      t73 = t45*t10
      t76 = 0.9723306394337274D2*t71-0.1111111111111111D0*t73*sigmaaa
      t80 = rhob*t52+t22*t76-2.D0*rhoa*sigmabb
      t84 = t7**2
      t85 = 1/t84
      t86 = t85*rhoa
      t88 = 1/t4/t17
      t91 = 0.2288509333333333D-1*t86*rhob*t88
      t92 = 1/t17
      t93 = rhob*t92
      t95 = 0.19672D0*t9*t93
      t96 = t17**2
      t97 = t96*rho
      t98 = 1/t97
      t99 = t98*t15
      t100 = t8*t63
      t102 = 0.548120936D-3*t99*t100
      t103 = t15*t85
      t106 = 0.75520808D-3*t103*t98*t63
      t108 = 1/t19/t96
      t111 = 0.2380312D-1*t16*t108*t63
      t113 = 1/t4/rho
      t115 = t113*t8
      t119 = 1/t19/rho*t85
      t121 = 0.3283518518518519D-1*t113+0.4524074074074074D-1*t115
     #-0.1578901851851852D-1*t119
      t132 = -0.8443333333333333D-1*t113-0.1163333333333333D0*t115
     #+0.4060033333333333D-1*t119
      t140 = -1.D0*rhoa*t92*sigmaaa-1.D0*t93*sigmabb
      t143 = t121*sigma-1.D0*(0.4690740740740741D-2*t113
     #+0.6462962962962963D-2*t115-0.2255574074074074D-2*t119)*t41
     #-0.1111111111111111D0*t132*t49-0.1111111111111111D0*t45*t140
      t151 = t22*t143-0.1333333333333333D1*rho*sigma
     #+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho*sigmaaa
      t154 = 0.649176D-2*t16*t21*t151
      vrhoa(i) = -0.19672D0*t67*t10-0.649176D-2*t16*t21*t80-t91+t95
     #-t102-t106+t111-t154
      t158 = t30*rhob
      t162 = 0.9723306394337274D2*t158-0.1111111111111111D0*t73*sigmabb
      t166 = rhoa*t52+t22*t162-2.D0*rhob*sigmaaa
      vrhob(i) = -0.19672D0*t9*t10-0.649176D-2*t16*t21*t166-t91+t95
     #-t102-t106+t111-t154
      t170 = 0.1407222222222222D-1*t5
      t171 = 0.1938888888888889D-1*t34
      t172 = t45*rhoa
      t175 = -0.25D1+t170+t171-0.1111111111111111D0*t172*t10
      t177 = t22*t175+t56-t60
      t183 = t22*t36-0.6666666666666667D0*t17
      t185 = t16*t21*t183
      t186 = 0.649176D-2*t185
      vsigmaaa(i) = -0.649176D-2*t16*t21*t177-t186
      vsigmaab(i) = -0.1298352D-1*t185
      t188 = t45*rhob
      t191 = -0.25D1+t170+t171-0.1111111111111111D0*t188*t10
      t193 = t22*t191+t56-t57
      vsigmabb(i) = -0.649176D-2*t16*t21*t193-t186
      t198 = t85*rhob*t88
      t200 = t67*t92
      t203 = t99*t8*t80
      t206 = t103*t98*t80
      t209 = t16*t108*t80
      t212 = t132*t10
      t215 = t45*t92
      t222 = t16*t21*(rhob*t143+t22*(-0.1111111111111111D0*t212
     #*sigmaaa+0.1111111111111111D0*t215*sigmaaa))
      t237 = 0.1110812266666667D0*t16/t19/t97*t63
      t239 = t88*t8
      t243 = 1/t19/t17*t85
      t245 = 1/t18
      t247 = 1/t84/t7
      t248 = t245*t247
      t271 = rhob*t245
      t285 = 0.649176D-2*t16*t21*(t22*((-0.4378024691358025D-1*t88
     #-0.6032098765432099D-1*t239+0.3157803703703704D-1*t243
     #-0.3673578308641975D-2*t248)*sigma-1.D0*(-0.6254320987654321D-2
     #*t88-0.8617283950617284D-2*t239+0.4511148148148148D-2*t243
     #-0.5247969012345679D-3*t248)*t41-0.1111111111111111D0*
     #(0.1125777777777778D0*t88+0.1551111111111111D0*t239
     #-0.8120066666666667D-1*t243+0.9446344222222222D-2*t248)*t49
     #-0.2222222222222222D0*t132*t140-0.1111111111111111D0*t45*(2.D0
     #*rhoa*t245*sigmaaa+2.D0*t271*sigmabb))-0.1333333333333333D1
     #*sigma+0.1333333333333333D1*sigmabb+0.1333333333333333D1*sigmaaa)
      t286 = t96*t17
      t287 = 1/t286
      t290 = 0.6545136693333333D-2*t103*t287*t63
      t293 = 0.4750381445333333D-2*t287*t15*t100
      t298 = 0.7628364444444444D-1*t86*rhob/t4/t18
      t300 = rhob*t21
      t302 = 0.5324598382222222D-2*t247*rhoa*t300
      t305 = 0.1096241872D-2*t99*t8*t151
      t307 = 1/t4/t286
      t308 = t307*t15
      t310 = 0.4627967769626667D-4*t308*t100
      t313 = 0.1275294711093333D-3*t308*t85*t63
      t316 = 0.151041616D-2*t103*t98*t151
      t320 = 0.1757117466133333D-3*t15*t247*t307*t63
      t323 = 0.4760624D-1*t16*t108*t151
      t325 = 0.39344D0*t9*t271
      v2rhoa2(i) = -0.4577018666666667D-1*t198+0.39344D0*t200
     #-0.1096241872D-2*t203-0.151041616D-2*t206+0.4760624D-1*t209
     #-0.1298352D-1*t222-0.649176D-2*t16*t21*(2.D0*rhob*t76
     #+0.1620551065722879D3*t71*rhob-2.D0*sigmabb)-t237-t285+t290+t293
     #+t298-t302-t305-t310-t313-t316-t320+t323-t325
      t335 = t86*t88
      t337 = t9*t92
      t340 = t99*t8*t166
      t343 = t103*t98*t166
      t346 = t16*t108*t166
      t357 = t16*t21*(rhoa*t143+t22*(-0.1111111111111111D0*t212
     #*sigmabb+0.1111111111111111D0*t215*sigmabb))
      v2rhob2(i) = -t237-t285+t290+t293+t298-t302-t305-t310-t313-t316
     #-t320+t323-0.649176D-2*t16*t21*(2.D0*rhoa*t162
     #+0.1620551065722879D3*rhoa*t158-2.D0*sigmaaa)
     #-0.4577018666666667D-1*t335+0.39344D0*t337-0.1096241872D-2*t340
     #-0.151041616D-2*t343+0.4760624D-1*t346-0.1298352D-1*t357-t325
      t377 = 0.19672D0*t337-0.2288509333333333D-1*t335-0.548120936D-3
     #*t340-0.75520808D-3*t343+0.2380312D-1*t346-0.19672D0*t8*t10
     #-0.649176D-2*t16*t21*(t27+t32+t37-t43-t51+rhob*t162+rhoa*t76)
     #-0.649176D-2*t357-0.548120936D-3*t203-0.75520808D-3*t206
     #+0.2380312D-1*t209-0.649176D-2*t222-t237
      t380 = -t285+t290+t293+t298-t302-t305-t310-t313-t316-t320+t323
     #-t325-0.2288509333333333D-1*t198+0.19672D0*t200
      v2rhoab(i) = t377+t380
      t383 = 0.1111111111111111D0*t22*t73
      t390 = 0.548120936D-3*t99*t8*t177
      t393 = 0.75520808D-3*t103*t98*t177
      t396 = 0.2380312D-1*t16*t108*t177
      t397 = 0.4690740740740741D-2*t113
      t398 = 0.6462962962962963D-2*t115
      t399 = 0.2255574074074074D-2*t119
      t407 = 0.1333333333333333D1*rho
      t411 = 0.649176D-2*t16*t21*(t22*(-t397-t398+t399
     #-0.1111111111111111D0*t132*rhoa*t10+0.1111111111111111D0*t172
     #*t92)+t407)
      t413 = t16*t300*t36
      t414 = 0.649176D-2*t413
      t416 = t99*t8*t183
      t417 = 0.548120936D-3*t416
      t419 = t103*t98*t183
      t420 = 0.75520808D-3*t419
      t422 = t16*t108*t183
      t423 = 0.2380312D-1*t422
      t428 = t16*t21*(t22*t121-0.1333333333333333D1*rho)
      t429 = 0.649176D-2*t428
      v2rhoasigmaaa(i) = -0.649176D-2*t16*t21*(rhob*t175-t383)-t390
     #-t393+t396-t411-t414-t417-t420+t423-t429
      t431 = 0.1096241872D-2*t416
      t432 = 0.151041616D-2*t419
      t433 = 0.4760624D-1*t422
      t434 = 0.1298352D-1*t428
      v2rhoasigmaab(i) = -0.1298352D-1*t413-t431-t432+t433-t434
      t443 = 0.548120936D-3*t99*t8*t193
      t446 = 0.75520808D-3*t103*t98*t193
      t449 = 0.2380312D-1*t16*t108*t193
      t460 = 0.649176D-2*t16*t21*(t22*(-t397-t398+t399
     #-0.1111111111111111D0*t132*rhob*t10+0.1111111111111111D0*t188
     #*t92)+t407)
      v2rhoasigmabb(i) = -0.649176D-2*t16*t21*(rhob*t191-2.D0*rhoa)
     #-t443-t446+t449-t460-t414-t417-t420+t423-t429
      t469 = t16*t21*rhoa*t36
      t470 = 0.649176D-2*t469
      v2rhobsigmaaa(i) = -0.649176D-2*t16*t21*(rhoa*t175-2.D0*rhob)
     #-t390-t393+t396-t411-t470-t417-t420+t423-t429
      v2rhobsigmaab(i) = -0.1298352D-1*t469-t431-t432+t433-t434
      v2rhobsigmabb(i) = -0.649176D-2*t16*t21*(rhoa*t191-t383)-t443
     #-t446+t449-t460-t470-t417-t420+t423-t429
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
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
      
      
      subroutine rks_c_lyp
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     C. Lee, W. Yang, and R.G. Parr
c     Development of the Colle-Salvetti correlation-energy formula into
c     a functional of the electron density
c     Phys. Rev. B37 (1988) 785-789
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
      t3 = 1/t2
      t6 = 1/(1.D0+0.349D0*t3)
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t12 = rho**2
      t14 = t2**2
      t20 = t3*t6
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t10*t6/t14/t12/rho*(0.25D0
     #*t12*(0.1148493600075277D2*t14*t12+(0.2611111111111111D1
     #-0.9850555555555556D-1*t3-0.1357222222222222D0*t20)*sigma-0.5D0*
     #(0.25D1-0.1407222222222222D-1*t3-0.1938888888888889D-1*t20)
     #*sigma-0.2777777777777778D-1*(t9+0.349D0*t20-11.D0)*sigma)
     #-0.4583333333333333D0*t12*sigma)
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
      t3 = 1/t2
      t5 = 1.D0+0.349D0*t3
      t6 = 1/t5
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t11 = t10*t6
      t12 = rho**2
      t14 = t2**2
      t16 = 1/t14/t12/rho
      t20 = t3*t6
      t22 = 0.2611111111111111D1-0.9850555555555556D-1*t3
     #-0.1357222222222222D0*t20
      t30 = t9+0.349D0*t20-11.D0
      t33 = 0.1148493600075277D2*t14*t12+t22*sigma-0.5D0*(0.25D1
     #-0.1407222222222222D-1*t3-0.1938888888888889D-1*t20)*sigma
     #-0.2777777777777778D-1*t30*sigma
      t38 = 0.25D0*t12*t33-0.4583333333333333D0*t12*sigma
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t11*t16*t38
      t45 = t14*rho
      t49 = t30/rho*sigma
      t54 = rho*sigma
      t60 = t5**2
      t61 = 1/t60
      t64 = t12**2
      t66 = 1/t64/rho
      t81 = 1/t2/rho
      t83 = t81*t6
      t85 = 1/t45
      t86 = t85*t61
      vrhoa(i) = -0.4918D-1*t6-0.649176D-2*t11*t16*(0.5D0*rho*t33
     #+0.25D0*t12*(0.3062649600200738D2*t45-0.2777777777777778D-1*t49)
     #-0.25D0*t54)-0.5721273333333333D-2*t61*t3-0.548120936D-3*t66*t10
     #*t6*t38-0.75520808D-3*t10*t61*t66*t38+0.2380312D-1*t11/t14/t64
     #*t38-0.649176D-2*t11*t16*(0.25D0*t12*((0.3283518518518519D-1*t81
     #+0.4524074074074074D-1*t83-0.1578901851851852D-1*t86)*sigma
     #-0.5D0*(0.4690740740740741D-2*t81+0.6462962962962963D-2*t83
     #-0.2255574074074074D-2*t86)*sigma-0.2777777777777778D-1*(
     #-0.8443333333333333D-1*t81-0.1163333333333333D0*t83
     #+0.4060033333333333D-1*t86)*sigma+0.2777777777777778D-1*t49)
     #-0.6666666666666667D0*t54)
      vsigmaaa(i) = 0.7213066666666667D-3*t11*t85-0.2596704D-1*t11*t16
     #*(0.25D0*t12*t22-0.6666666666666667D0*t12)
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
      t3 = 1/t2
      t5 = 1.D0+0.349D0*t3
      t6 = 1/t5
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t11 = t10*t6
      t12 = rho**2
      t13 = t12*rho
      t14 = t2**2
      t16 = 1/t14/t13
      t17 = t14*t12
      t18 = 0.1148493600075277D2*t17
      t20 = t3*t6
      t22 = 0.2611111111111111D1-0.9850555555555556D-1*t3
     #-0.1357222222222222D0*t20
      t23 = t22*sigma
      t28 = 0.5D0*(0.25D1-0.1407222222222222D-1*t3
     #-0.1938888888888889D-1*t20)*sigma
      t30 = t9+0.349D0*t20-11.D0
      t32 = 0.2777777777777778D-1*t30*sigma
      t33 = t18+t23-t28-t32
      t38 = 0.25D0*t12*t33-0.4583333333333333D0*t12*sigma
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t11*t16*t38
      t45 = t14*rho
      t47 = 1/rho
      t48 = t30*t47
      t49 = t48*sigma
      t51 = 0.3062649600200738D2*t45-0.2777777777777778D-1*t49
      t54 = rho*sigma
      t56 = 0.5D0*rho*t33+0.25D0*t12*t51-0.25D0*t54
      t60 = t5**2
      t61 = 1/t60
      t64 = t12**2
      t65 = t64*rho
      t66 = 1/t65
      t67 = t66*t10
      t68 = t6*t38
      t71 = t10*t61
      t76 = 1/t14/t64
      t81 = 1/t2/rho
      t83 = t81*t6
      t85 = 1/t45
      t86 = t85*t61
      t88 = 0.3283518518518519D-1*t81+0.4524074074074074D-1*t83
     #-0.1578901851851852D-1*t86
      t99 = -0.8443333333333333D-1*t81-0.1163333333333333D0*t83
     #+0.4060033333333333D-1*t86
      t103 = t88*sigma-0.5D0*(0.4690740740740741D-2*t81
     #+0.6462962962962963D-2*t83-0.2255574074074074D-2*t86)*sigma
     #-0.2777777777777778D-1*t99*sigma+0.2777777777777778D-1*t49
      t107 = 0.25D0*t12*t103-0.6666666666666667D0*t54
      vrhoa(i) = -0.4918D-1*t6-0.649176D-2*t11*t16*t56
     #-0.5721273333333333D-2*t61*t3-0.548120936D-3*t67*t68
     #-0.75520808D-3*t71*t66*t38+0.2380312D-1*t11*t76*t38-0.649176D-2
     #*t11*t16*t107
      t116 = 0.25D0*t12*t22-0.6666666666666667D0*t12
      vsigmaaa(i) = 0.7213066666666667D-3*t11*t85-0.2596704D-1*t11*t16
     #*t116
      t121 = 1/t60/t5
      t126 = t64*t12
      t127 = 1/t126
      t137 = t99*t47*sigma
      t141 = t30/t12*sigma
      t150 = rho*t51
      t165 = 1/t2/t126
      t166 = t165*t10
      t188 = 1/t2/t12
      t190 = t188*t6
      t192 = 1/t17
      t193 = t192*t61
      t195 = 1/t13
      t196 = t195*t121
      s1 = -0.2662299191111111D-2*t121*t85-0.7628364444444444D-2*t61
     #*t81+0.9500762890666667D-2*t127*t10*t68+0.9521248D-1*t11*t76*t56
     #-0.2596704D-1*t11*t16*(0.5D0*rho*t103+0.25D0*t12*(
     #-0.2777777777777778D-1*t137+0.2777777777777778D-1*t141))
     #-0.649176D-2*t11*t16*(t18+t23-t28-t32+t150)+0.9521248D-1*t11*t76
     #*t107-0.2192483744D-2*t67*t6*t56-0.302083232D-2*t71*t66*t56
      s2 = s1-0.2550589422186667D-3*t166*t61*t38-0.302083232D-2*t71
     #*t66*t107-0.3514234932266667D-3*t10*t121*t165*t38
     #-0.2192483744D-2*t67*t6*t107
      v2rhoa2(i) = s2-0.9255935539253333D-4*t166*t68
     #-0.2221624533333333D0*t11/t14/t65*t38-0.1298352D-1*t11*t16*
     #(0.25D0*t12*((-0.4378024691358025D-1*t188-0.6032098765432099D-1
     #*t190+0.3157803703703704D-1*t193-0.3673578308641975D-2*t196)
     #*sigma-0.5D0*(-0.6254320987654321D-2*t188-0.8617283950617284D-2
     #*t190+0.4511148148148148D-2*t193-0.5247969012345679D-3*t196)
     #*sigma-0.2777777777777778D-1*(0.1125777777777778D0*t188
     #+0.1551111111111111D0*t190-0.8120066666666667D-1*t193
     #+0.9446344222222222D-2*t196)*sigma+0.5555555555555556D-1*t137
     #-0.5555555555555556D-1*t141)-0.6666666666666667D0*sigma)
     #-0.649176D-2*t11*t16*(0.1D1*t150+0.2552208000167282D2*t17-0.5D0
     #*sigma)+0.1309027338666667D-1*t71*t127*t38
      v2rhoasigmaaa(i) = -0.649176D-2*t11*t16*(-0.9444444444444444D0
     #*rho-0.2777777777777778D-1*rho*t30)+0.6090232622222222D-4*t195
     #*t10*t6+0.8391200888888889D-4*t71*t195+0.9978075555555556D-2*t11
     #*t192-0.1298352D-1*t11*t16*(0.25D0*t12*(-0.3D-22*t83+0.1D-22*t86
     #+0.5555555555555556D-1*t48)+0.1333333333333333D1*rho)
     #-0.1298352D-1*t11*t192*t22-0.2192483744D-2*t67*t6*t116
     #-0.302083232D-2*t71*t66*t116+0.9521248D-1*t11*t76*t116
     #-0.2596704D-1*t11*t16*(0.25D0*t12*t88-0.1333333333333333D1*rho)
      v2sigmaaa2(i) = 0.D0
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

c:C_LYPsubrend
