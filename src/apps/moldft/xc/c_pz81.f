c:C_PZ81subrstart
c    Generated: Tue Mar  9 13:25:27 GMT 2004
      subroutine uks_c_pz81
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, and A. Zunger
c     Self-interaction correction to density-functional approximations 
c     for many-electron systems
c     Phys. Rev. B23 (1981) 5048-5079
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
      logical t4
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t5 = t1**(1.D0/6.D0)
      t11 = dlog(t3)
      t17 = piecewise(1.D0 .le. t3,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t5+0.1619735131738333D0*t2),0.1555D-1*t11
     #-0.269D-1+0.43424534362958D-3*t2*t11-0.297768235631712D-2*t2)
      zk(i) = rhob*t17
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t5 = t1**(1.D0/6.D0)
      t11 = dlog(t3)
      t17 = piecewise(1.D0 .le. t3,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t5+0.1619735131738333D0*t2),0.1555D-1*t11
     #-0.269D-1+0.43424534362958D-3*t2*t11-0.297768235631712D-2*t2)
      zk(i) = rhoa*t17
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t5 = t1**(1.D0/6.D0)
      t10 = 0.1423D0/(1.D0+0.8292885914166397D0*t5+0.20682485366586D0
     #*t2)
      t19 = (rhoa-1.D0*rhob)*t1
      t20 = 1.D0+t19
      t21 = t20**(1.D0/3.D0)
      t24 = 1.D0-1.D0*t19
      t25 = t24**(1.D0/3.D0)
      t27 = t21*t20+t25*t24-2.D0
      t31 = dlog(t3)
      t33 = t2*t31
      t43 = piecewise(1.D0 .le. t3,-t10+0.1923661050931536D1*(
     #-0.843D-1/(1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2)
     #+t10)*t27,0.311D-1*t31-0.48D-1+0.12407009817988D-2*t33
     #-0.719606569443304D-2*t2+0.1923661050931536D1*(-0.1555D-1*t31
     #+0.211D-1-0.80645563816922D-3*t33+0.421838333811592D-2*t2)*t27)
      zk(i) = rho*t43
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
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2
      t11 = dlog(t3)
      t17 = piecewise(t4,-0.843D-1/t8,0.1555D-1*t11-0.269D-1
     #+0.43424534362958D-3*t2*t11-0.297768235631712D-2*t2)
      zk(i) = rhob*t17
      vrhoa(i) = 0.D0
      t18 = t8**2
      t20 = t5**2
      t21 = t20**2
      t24 = rhob**2
      t25 = 1/t24
      t28 = t2**2
      t29 = 1/t28
      t30 = t29*t25
      t43 = piecewise(t4,0.843D-1/t18*(-0.1835293601259385D0/t21/t5
     #*t25-0.5399117105794445D-1*t30),-0.5183333333333333D-2*t1
     #-0.1447484478765267D-3*t29*t11*t25-0.1447484478765267D-3*t2*t1
     #+0.99256078543904D-3*t30)
      vrhob(i) = t17+rhob*t43
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2
      t11 = dlog(t3)
      t17 = piecewise(t4,-0.843D-1/t8,0.1555D-1*t11-0.269D-1
     #+0.43424534362958D-3*t2*t11-0.297768235631712D-2*t2)
      zk(i) = rhoa*t17
      t18 = t8**2
      t20 = t5**2
      t21 = t20**2
      t24 = rhoa**2
      t25 = 1/t24
      t28 = t2**2
      t29 = 1/t28
      t30 = t29*t25
      t43 = piecewise(t4,0.843D-1/t18*(-0.1835293601259385D0/t21/t5
     #*t25-0.5399117105794445D-1*t30),-0.5183333333333333D-2*t1
     #-0.1447484478765267D-3*t29*t11*t25-0.1447484478765267D-3*t2*t1
     #+0.99256078543904D-3*t30)
      vrhoa(i) = t17+rhoa*t43
      vrhob(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.8292885914166397D0*t5+0.20682485366586D0*t2
      t10 = 0.1423D0/t8
      t13 = 1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2
      t16 = -0.843D-1/t13+t10
      t18 = rhoa-1.D0*rhob
      t19 = t18*t1
      t20 = 1.D0+t19
      t21 = t20**(1.D0/3.D0)
      t24 = 1.D0-1.D0*t19
      t25 = t24**(1.D0/3.D0)
      t27 = t21*t20+t25*t24-2.D0
      t31 = dlog(t3)
      t33 = t2*t31
      t39 = -0.1555D-1*t31+0.211D-1-0.80645563816922D-3*t33
     #+0.421838333811592D-2*t2
      t43 = piecewise(t4,-t10+0.1923661050931536D1*t16*t27,0.311D-1
     #*t31-0.48D-1+0.12407009817988D-2*t33-0.719606569443304D-2*t2
     #+0.1923661050931536D1*t39*t27)
      zk(i) = rho*t43
      t48 = 0.1333333333333333D1*t21*t1-0.1333333333333333D1*t25*t1
      t53 = piecewise(t4,0.1923661050931536D1*t16
     #*t48,0.1923661050931536D1*t39*t48)
      t55 = t8**2
      t57 = t5**2
      t58 = t57**2
      t61 = rho**2
      t62 = 1/t61
      t63 = 1/t58/t5*t62
      t65 = t2**2
      t66 = 1/t65
      t67 = t66*t62
      t71 = 0.1423D0/t55*(-0.1382147652361066D0*t63
     #-0.6894161788861999D-1*t67)
      t72 = t13**2
      t88 = -0.1333333333333333D1*t21*t18*t62+0.1333333333333333D1*t25
     #*t18*t62
      t94 = t66*t31*t62
      t96 = t2*t1
      t109 = piecewise(t4,t71+0.1923661050931536D1*(0.843D-1/t72*(
     #-0.1835293601259385D0*t63-0.5399117105794445D-1*t67)-t71)*t27
     #+0.1923661050931536D1*t16*t88,-0.1036666666666667D-1*t1
     #-0.4135669939329333D-3*t94-0.4135669939329333D-3*t96
     #+0.2398688564811013D-2*t67+0.1923661050931536D1*
     #(0.5183333333333333D-2*t1+0.2688185460564067D-3*t94
     #+0.2688185460564067D-3*t96-0.1406127779371973D-2*t67)*t27
     #+0.1923661050931536D1*t39*t88)
      t110 = rho*t109
      vrhoa(i) = rho*t53+t43+t110
      t111 = -t48
      t116 = piecewise(t4,0.1923661050931536D1*t16
     #*t111,0.1923661050931536D1*t39*t111)
      vrhob(i) = rho*t116+t43+t110
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
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2
      t11 = dlog(t3)
      t17 = piecewise(t4,-0.843D-1/t8,0.1555D-1*t11-0.269D-1
     #+0.43424534362958D-3*t2*t11-0.297768235631712D-2*t2)
      zk(i) = rhob*t17
      vrhoa(i) = 0.D0
      t18 = t8**2
      t19 = 1/t18
      t20 = t5**2
      t21 = t20**2
      t22 = t21*t5
      t23 = 1/t22
      t24 = rhob**2
      t25 = 1/t24
      t28 = t2**2
      t29 = 1/t28
      t30 = t29*t25
      t32 = -0.1835293601259385D0*t23*t25-0.5399117105794445D-1*t30
      t36 = t29*t11
      t43 = piecewise(t4,0.843D-1*t19*t32,-0.5183333333333333D-2*t1
     #-0.1447484478765267D-3*t36*t25-0.1447484478765267D-3*t2*t1
     #+0.99256078543904D-3*t30)
      vrhob(i) = t17+rhob*t43
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t48 = t32**2
      t53 = t24**2
      t54 = 1/t53
      t58 = 1/t24/rhob
      t62 = 1/t28/t1
      t63 = t62*t54
      t65 = t29*t58
      t82 = piecewise(t4,-0.1686D0/t18/t8*t48+0.843D-1*t19*(
     #-0.1529411334382821D0/t22/t1*t54+0.367058720251877D0*t23*t58
     #-0.3599411403862963D-1*t63+0.1079823421158889D0*t65
     #),0.5183333333333333D-2*t25-0.9649896525101778D-4*t62*t11*t54
     #-0.1888622605627062D-2*t65+0.2894968957530533D-3*t36*t58
     #+0.1447484478765267D-3*t2*t25+0.6617071902926934D-3*t63)
      v2rhob2(i) = 2.D0*t43+rhob*t82
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2
      t11 = dlog(t3)
      t17 = piecewise(t4,-0.843D-1/t8,0.1555D-1*t11-0.269D-1
     #+0.43424534362958D-3*t2*t11-0.297768235631712D-2*t2)
      zk(i) = rhoa*t17
      t18 = t8**2
      t19 = 1/t18
      t20 = t5**2
      t21 = t20**2
      t22 = t21*t5
      t23 = 1/t22
      t24 = rhoa**2
      t25 = 1/t24
      t28 = t2**2
      t29 = 1/t28
      t30 = t29*t25
      t32 = -0.1835293601259385D0*t23*t25-0.5399117105794445D-1*t30
      t36 = t29*t11
      t43 = piecewise(t4,0.843D-1*t19*t32,-0.5183333333333333D-2*t1
     #-0.1447484478765267D-3*t36*t25-0.1447484478765267D-3*t2*t1
     #+0.99256078543904D-3*t30)
      vrhoa(i) = t17+rhoa*t43
      vrhob(i) = 0.D0
      t48 = t32**2
      t53 = t24**2
      t54 = 1/t53
      t58 = 1/t24/rhoa
      t62 = 1/t28/t1
      t63 = t62*t54
      t65 = t29*t58
      t82 = piecewise(t4,-0.1686D0/t18/t8*t48+0.843D-1*t19*(
     #-0.1529411334382821D0/t22/t1*t54+0.367058720251877D0*t23*t58
     #-0.3599411403862963D-1*t63+0.1079823421158889D0*t65
     #),0.5183333333333333D-2*t25-0.9649896525101778D-4*t62*t11*t54
     #-0.1888622605627062D-2*t65+0.2894968957530533D-3*t36*t58
     #+0.1447484478765267D-3*t2*t25+0.6617071902926934D-3*t63)
      v2rhoa2(i) = 2.D0*t43+rhoa*t82
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.8292885914166397D0*t5+0.20682485366586D0*t2
      t10 = 0.1423D0/t8
      t13 = 1.D0+0.1101176160755631D1*t5+0.1619735131738333D0*t2
      t16 = -0.843D-1/t13+t10
      t18 = rhoa-1.D0*rhob
      t19 = t18*t1
      t20 = 1.D0+t19
      t21 = t20**(1.D0/3.D0)
      t24 = 1.D0-1.D0*t19
      t25 = t24**(1.D0/3.D0)
      t27 = t21*t20+t25*t24-2.D0
      t31 = dlog(t3)
      t33 = t2*t31
      t39 = -0.1555D-1*t31+0.211D-1-0.80645563816922D-3*t33
     #+0.421838333811592D-2*t2
      t43 = piecewise(t4,-t10+0.1923661050931536D1*t16*t27,0.311D-1
     #*t31-0.48D-1+0.12407009817988D-2*t33-0.719606569443304D-2*t2
     #+0.1923661050931536D1*t39*t27)
      zk(i) = rho*t43
      t48 = 0.1333333333333333D1*t21*t1-0.1333333333333333D1*t25*t1
      t53 = piecewise(t4,0.1923661050931536D1*t16
     #*t48,0.1923661050931536D1*t39*t48)
      t55 = t8**2
      t56 = 1/t55
      t57 = t5**2
      t58 = t57**2
      t59 = t58*t5
      t60 = 1/t59
      t61 = rho**2
      t62 = 1/t61
      t63 = t60*t62
      t65 = t2**2
      t66 = 1/t65
      t67 = t66*t62
      t69 = -0.1382147652361066D0*t63-0.6894161788861999D-1*t67
      t71 = 0.1423D0*t56*t69
      t72 = t13**2
      t73 = 1/t72
      t76 = -0.1835293601259385D0*t63-0.5399117105794445D-1*t67
      t79 = 0.843D-1*t73*t76-t71
      t82 = t21*t18
      t85 = t25*t18
      t88 = -0.1333333333333333D1*t82*t62+0.1333333333333333D1*t85*t62
      t93 = t66*t31
      t94 = t93*t62
      t96 = t2*t1
      t103 = 0.5183333333333333D-2*t1+0.2688185460564067D-3*t94
     #+0.2688185460564067D-3*t96-0.1406127779371973D-2*t67
      t109 = piecewise(t4,t71+0.1923661050931536D1*t79*t27
     #+0.1923661050931536D1*t16*t88,-0.1036666666666667D-1*t1
     #-0.4135669939329333D-3*t94-0.4135669939329333D-3*t96
     #+0.2398688564811013D-2*t67+0.1923661050931536D1*t103*t27
     #+0.1923661050931536D1*t39*t88)
      t110 = rho*t109
      vrhoa(i) = rho*t53+t43+t110
      t111 = -t48
      t116 = piecewise(t4,0.1923661050931536D1*t16
     #*t111,0.1923661050931536D1*t39*t111)
      vrhob(i) = rho*t116+t43+t110
      t118 = t21**2
      t119 = 1/t118
      t122 = t25**2
      t123 = 1/t122
      t126 = 0.4444444444444444D0*t119*t62+0.4444444444444444D0*t123*t62
      t131 = piecewise(t4,0.1923661050931536D1*t16
     #*t126,0.1923661050931536D1*t39*t126)
      t132 = rho*t131
      t138 = 1/t61/rho
      t148 = -0.4444444444444444D0*t119*t18*t138-0.1333333333333333D1
     #*t21*t62-0.4444444444444444D0*t123*t18*t138+0.1333333333333333D1
     #*t25*t62
      t157 = piecewise(t4,0.1923661050931536D1*t79*t48
     #+0.1923661050931536D1*t16*t148,0.1923661050931536D1*t103*t48
     #+0.1923661050931536D1*t39*t148)
      t158 = rho*t157
      t160 = 2.D0*t109
      t163 = t69**2
      t165 = 0.2846D0/t55/t8*t163
      t168 = t61**2
      t169 = 1/t168
      t170 = 1/t59/t1*t169
      t172 = t60*t138
      t175 = 1/t65/t1
      t176 = t175*t169
      t178 = t66*t138
      t182 = 0.1423D0*t56*(-0.1151789710300888D0*t170
     #+0.2764295304722132D0*t172-0.4596107859241333D-1*t176
     #+0.13788323577724D0*t178)
      t185 = t76**2
      t200 = t18**2
      t211 = 0.4444444444444444D0*t119*t200*t169+0.2666666666666667D1
     #*t82*t138+0.4444444444444444D0*t123*t200*t169
     #-0.2666666666666667D1*t85*t138
      t217 = t175*t31*t169
      t220 = t93*t138
      t222 = t2*t62
      t239 = piecewise(t4,-t165+t182+0.1923661050931536D1*(-0.1686D0
     #/t72/t13*t185+0.843D-1*t73*(-0.1529411334382821D0*t170
     #+0.367058720251877D0*t172-0.3599411403862963D-1*t176
     #+0.1079823421158889D0*t178)+t165-t182)*t27+0.3847322101863073D1
     #*t79*t88+0.1923661050931536D1*t16*t211,0.1036666666666667D-1*t62
     #-0.2757113292886222D-3*t217-0.4521665800333405D-2*t178
     #+0.8271339878658667D-3*t220+0.4135669939329333D-3*t222
     #+0.1599125709874009D-2*t176+0.1923661050931536D1*(
     #-0.5183333333333333D-2*t62+0.1792123640376044D-3*t217
     #+0.2633043194706342D-2*t178-0.5376370921128133D-3*t220
     #-0.2688185460564067D-3*t222-0.9374185195813156D-3*t176)*t27
     #+0.3847322101863073D1*t103*t88+0.1923661050931536D1*t39*t211)
      t240 = rho*t239
      v2rhoa2(i) = t132+2.D0*t53+2.D0*t158+t160+t240
      t244 = -t148
      t253 = piecewise(t4,0.1923661050931536D1*t79*t111
     #+0.1923661050931536D1*t16*t244,0.1923661050931536D1*t103*t111
     #+0.1923661050931536D1*t39*t244)
      t254 = rho*t253
      v2rhob2(i) = t132+2.D0*t116+2.D0*t254+t160+t240
      t256 = -t126
      t261 = piecewise(t4,0.1923661050931536D1*t16
     #*t256,0.1923661050931536D1*t39*t256)
      v2rhoab(i) = rho*t261+t116+t254+t53+t158+t160+t240
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
      
      
      subroutine rks_c_pz81
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, and A. Zunger
c     Self-interaction correction to density-functional approximations 
c     for many-electron systems
c     Phys. Rev. B23 (1981) 5048-5079
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
      logical t4
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t5 = t1**(1.D0/6.D0)
      t11 = dlog(t3)
      t17 = piecewise(1.D0 .le. t3,-0.1423D0/(1.D0
     #+0.8292885914166397D0*t5+0.20682485366586D0*t2),0.311D-1*t11
     #-0.48D-1+0.12407009817988D-2*t2*t11-0.719606569443304D-2*t2)
      zk(i) = rho*t17
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.8292885914166397D0*t5+0.20682485366586D0*t2
      t11 = dlog(t3)
      t17 = piecewise(t4,-0.1423D0/t8,0.311D-1*t11-0.48D-1
     #+0.12407009817988D-2*t2*t11-0.719606569443304D-2*t2)
      zk(i) = rho*t17
      t18 = piecewise(t4,0.D0,0.D0)
      t20 = t8**2
      t22 = t5**2
      t23 = t22**2
      t26 = rho**2
      t27 = 1/t26
      t30 = t2**2
      t31 = 1/t30
      t32 = t31*t27
      t45 = piecewise(t4,0.1423D0/t20*(-0.1382147652361066D0/t23/t5
     #*t27-0.6894161788861999D-1*t32),-0.1036666666666667D-1*t1
     #-0.4135669939329333D-3*t31*t11*t27-0.4135669939329333D-3*t2*t1
     #+0.2398688564811013D-2*t32)
      vrhoa(i) = rho*t18+t17+rho*t45
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.8292885914166397D0*t5+0.20682485366586D0*t2
      t10 = 0.1423D0/t8
      t11 = dlog(t3)
      t13 = t2*t11
      t17 = piecewise(t4,-t10,0.311D-1*t11-0.48D-1+0.12407009817988D-2
     #*t13-0.719606569443304D-2*t2)
      zk(i) = rho*t17
      t18 = piecewise(t4,0.D0,0.D0)
      t19 = rho*t18
      t20 = t8**2
      t21 = 1/t20
      t22 = t5**2
      t23 = t22**2
      t24 = t23*t5
      t25 = 1/t24
      t26 = rho**2
      t27 = 1/t26
      t30 = t2**2
      t31 = 1/t30
      t32 = t31*t27
      t34 = -0.1382147652361066D0*t25*t27-0.6894161788861999D-1*t32
      t38 = t31*t11
      t45 = piecewise(t4,0.1423D0*t21*t34,-0.1036666666666667D-1*t1
     #-0.4135669939329333D-3*t38*t27-0.4135669939329333D-3*t2*t1
     #+0.2398688564811013D-2*t32)
      vrhoa(i) = t19+t17+rho*t45
      t53 = (-0.843D-1/(1.D0+0.1101176160755631D1*t5
     #+0.1619735131738333D0*t2)+t10)*t27
      t59 = (-0.1555D-1*t11+0.211D-1-0.80645563816922D-3*t13
     #+0.421838333811592D-2*t2)*t27
      t61 = piecewise(t4,0.1709920934161366D1*t53,0.1709920934161366D1
     #*t59)
      t68 = t34**2
      t73 = t26**2
      t74 = 1/t73
      t78 = 1/t26/rho
      t82 = 1/t30/t1
      t83 = t82*t74
      t85 = t31*t78
      t102 = piecewise(t4,-0.2846D0/t20/t8*t68+0.1423D0*t21*(
     #-0.1151789710300888D0/t24/t1*t74+0.2764295304722132D0*t25*t78
     #-0.4596107859241333D-1*t83+0.13788323577724D0*t85
     #),0.1036666666666667D-1*t27-0.2757113292886222D-3*t82*t11*t74
     #-0.4521665800333405D-2*t85+0.8271339878658667D-3*t38*t78
     #+0.4135669939329333D-3*t2*t27+0.1599125709874009D-2*t83)
      t107 = piecewise(t4,-0.1709920934161366D1*t53,
     #-0.1709920934161366D1*t59)
      v2rhoa2(i) = rho*t61+4.D0*t18+4.D0*t19+4.D0*t45+2.D0*rho*t102
     #+rho*t107
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
c:C_PZ81subrend
