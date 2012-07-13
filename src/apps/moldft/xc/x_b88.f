c:X_B88subrstart
c    Generated: Wed Jan 29 09:08:35 GMT 2003
      subroutine uks_x_b88
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     A.D. Becke
c     Density-functional exchange-energy approximation with correct 
c     asymptotic behaviour
c     Phys. Rev. A38 (1988) 3098-3100
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
      t3 = t2*rhob
      t5 = 1/t3
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t5*sigmabb/(1.D0+0.252D-1
     #*t8*t9)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t5*sigmaaa/(1.D0+0.252D-1
     #*t8*t9)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      zk(i) = -0.9305257363491D0*t5-0.42D-2*t7*sigmaaa/(1.D0+0.252D-1
     #*t10*t11)-0.9305257363491D0*t19-0.42D-2*t21*sigmabb/(1.D0
     #+0.252D-1*t24*t25)
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
      t5 = 1/t3
      t6 = t5*sigmabb
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      vrhoa(i) = 0.D0
      t17 = rhob**2
      t19 = 1/t2/t17
      t23 = t12**2
      t24 = 1/t23
      t29 = t2**2
      t34 = 1/t29/t17
      t37 = dsqrt(1.D0+sigmabb*t34)
      t38 = 1/t37
      vrhob(i) = -0.12407009817988D1*t2+0.56D-2*t19*sigmabb*t13
     #+0.42D-2*t6*t24*(-0.336D-1*t7*t19*t9-0.336D-1*sigmabb/t29/t17
     #/rhob*t38)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*(0.126D-1/t7*t5*t9
     #+0.126D-1*t34*t38)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t6 = t5*sigmaaa
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      t17 = rhoa**2
      t19 = 1/t2/t17
      t23 = t12**2
      t24 = 1/t23
      t29 = t2**2
      t34 = 1/t29/t17
      t37 = dsqrt(1.D0+sigmaaa*t34)
      t38 = 1/t37
      vrhoa(i) = -0.12407009817988D1*t2+0.56D-2*t19*sigmaaa*t13
     #+0.42D-2*t6*t24*(-0.336D-1*t7*t19*t9-0.336D-1*sigmaaa/t29/t17
     #/rhoa*t38)
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*(0.126D-1/t7*t5*t9
     #+0.126D-1*t34*t38)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t8 = t7*sigmaaa
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t14 = 1.D0+0.252D-1*t10*t11
      t15 = 1/t14
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t22 = t21*sigmabb
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t28 = 1.D0+0.252D-1*t24*t25
      t29 = 1/t28
      zk(i) = -0.9305257363491D0*t5-0.42D-2*t8*t15-0.9305257363491D0
     #*t19-0.42D-2*t22*t29
      t33 = rhoa**2
      t35 = 1/t4/t33
      t39 = t14**2
      t40 = 1/t39
      t45 = t4**2
      t50 = 1/t45/t33
      t53 = dsqrt(1.D0+sigmaaa*t50)
      t54 = 1/t53
      vrhoa(i) = -0.12407009817988D1*t4+0.56D-2*t35*sigmaaa*t15
     #+0.42D-2*t8*t40*(-0.336D-1*t9*t35*t11-0.336D-1*sigmaaa/t45/t33
     #/rhoa*t54)
      t62 = rhob**2
      t64 = 1/t18/t62
      t68 = t28**2
      t69 = 1/t68
      t74 = t18**2
      t79 = 1/t74/t62
      t82 = dsqrt(1.D0+sigmabb*t79)
      t83 = 1/t82
      vrhob(i) = -0.12407009817988D1*t18+0.56D-2*t64*sigmabb*t29
     #+0.42D-2*t22*t69*(-0.336D-1*t23*t64*t25-0.336D-1*sigmabb/t74/t62
     #/rhob*t83)
      vsigmaaa(i) = -0.42D-2*t7*t15+0.42D-2*t8*t40*(0.126D-1/t9*t7*t11
     #+0.126D-1*t50*t54)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.42D-2*t21*t29+0.42D-2*t22*t69*(0.126D-1/t23*t21
     #*t25+0.126D-1*t79*t83)
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
      t5 = 1/t3
      t6 = t5*sigmabb
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      vrhoa(i) = 0.D0
      t17 = rhob**2
      t19 = 1/t2/t17
      t20 = t19*sigmabb
      t23 = t12**2
      t24 = 1/t23
      t28 = t17*rhob
      t29 = t2**2
      t34 = 1/t29/t17
      t36 = 1.D0+sigmabb*t34
      t37 = dsqrt(t36)
      t38 = 1/t37
      t41 = -0.336D-1*t7*t19*t9-0.336D-1*sigmabb/t29/t28*t38
      t42 = t24*t41
      vrhob(i) = -0.12407009817988D1*t2+0.56D-2*t20*t13+0.42D-2*t6*t42
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t53 = 0.126D-1/t7*t5*t9+0.126D-1*t34*t38
      vsigmabb(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*t53
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t60 = 1/t2/t28
      t67 = 1/t23/t12
      t68 = t41**2
      t75 = t17**2
      t81 = sigmabb**2
      t87 = 1/t37/t36
      v2rhob2(i) = -0.4135669939329333D0/t29-0.1306666666666667D-1*t60
     #*sigmabb*t13-0.112D-1*t20*t42-0.84D-2*t6*t67*t68+0.42D-2*t6*t24*
     #(0.784D-1*t7*t60*t9+0.168D0*sigmabb/t29/t75*t38-0.448D-1*t81/t2
     #/t75/t28*t87)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t97 = t53**2
      v2sigmabb2(i) = 0.84D-2*t5*t24*t53-0.84D-2*t6*t67*t97+0.42D-2*t6
     #*t24*(-0.63D-2/t7/sigmabb*t5*t9+0.63D-2/sigmabb*t34*t38-0.63D-2
     #/t2/t75/rhob*t87)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t6 = t5*sigmaaa
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      t17 = rhoa**2
      t19 = 1/t2/t17
      t20 = t19*sigmaaa
      t23 = t12**2
      t24 = 1/t23
      t28 = t17*rhoa
      t29 = t2**2
      t34 = 1/t29/t17
      t36 = 1.D0+sigmaaa*t34
      t37 = dsqrt(t36)
      t38 = 1/t37
      t41 = -0.336D-1*t7*t19*t9-0.336D-1*sigmaaa/t29/t28*t38
      t42 = t24*t41
      vrhoa(i) = -0.12407009817988D1*t2+0.56D-2*t20*t13+0.42D-2*t6*t42
      vrhob(i) = 0.D0
      t53 = 0.126D-1/t7*t5*t9+0.126D-1*t34*t38
      vsigmaaa(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*t53
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t60 = 1/t2/t28
      t67 = 1/t23/t12
      t68 = t41**2
      t75 = t17**2
      t81 = sigmaaa**2
      t87 = 1/t37/t36
      v2rhoa2(i) = -0.4135669939329333D0/t29-0.1306666666666667D-1*t60
     #*sigmaaa*t13-0.112D-1*t20*t42-0.84D-2*t6*t67*t68+0.42D-2*t6*t24*
     #(0.784D-1*t7*t60*t9+0.168D0*sigmaaa/t29/t75*t38-0.448D-1*t81/t2
     #/t75/t28*t87)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t97 = t53**2
      v2sigmaaa2(i) = 0.84D-2*t5*t24*t53-0.84D-2*t6*t67*t97+0.42D-2*t6
     #*t24*(-0.63D-2/t7/sigmaaa*t5*t9+0.63D-2/sigmaaa*t34*t38-0.63D-2
     #/t2/t75/rhoa*t87)
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
      t7 = 1/t5
      t8 = t7*sigmaaa
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t14 = 1.D0+0.252D-1*t10*t11
      t15 = 1/t14
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t22 = t21*sigmabb
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t28 = 1.D0+0.252D-1*t24*t25
      t29 = 1/t28
      zk(i) = -0.9305257363491D0*t5-0.42D-2*t8*t15-0.9305257363491D0
     #*t19-0.42D-2*t22*t29
      t33 = rhoa**2
      t35 = 1/t4/t33
      t36 = t35*sigmaaa
      t39 = t14**2
      t40 = 1/t39
      t44 = t33*rhoa
      t45 = t4**2
      t47 = 1/t45/t44
      t50 = 1/t45/t33
      t52 = 1.D0+sigmaaa*t50
      t53 = dsqrt(t52)
      t54 = 1/t53
      t57 = -0.336D-1*t9*t35*t11-0.336D-1*sigmaaa*t47*t54
      t58 = t40*t57
      vrhoa(i) = -0.12407009817988D1*t4+0.56D-2*t36*t15+0.42D-2*t8*t58
      t62 = rhob**2
      t64 = 1/t18/t62
      t65 = t64*sigmabb
      t68 = t28**2
      t69 = 1/t68
      t73 = t62*rhob
      t74 = t18**2
      t76 = 1/t74/t73
      t79 = 1/t74/t62
      t81 = 1.D0+sigmabb*t79
      t82 = dsqrt(t81)
      t83 = 1/t82
      t86 = -0.336D-1*t23*t64*t25-0.336D-1*sigmabb*t76*t83
      t87 = t69*t86
      vrhob(i) = -0.12407009817988D1*t18+0.56D-2*t65*t29+0.42D-2*t22*t87
      t92 = 1/t9
      t98 = 0.126D-1*t92*t7*t11+0.126D-1*t50*t54
      t99 = t40*t98
      vsigmaaa(i) = -0.42D-2*t7*t15+0.42D-2*t8*t99
      vsigmaab(i) = 0.D0
      t104 = 1/t23
      t110 = 0.126D-1*t104*t21*t25+0.126D-1*t79*t83
      t111 = t69*t110
      vsigmabb(i) = -0.42D-2*t21*t29+0.42D-2*t22*t111
      t117 = 1/t4/t44
      t124 = 1/t39/t14
      t125 = t57**2
      t132 = t33**2
      t138 = sigmaaa**2
      t144 = 1/t53/t52
      v2rhoa2(i) = -0.4135669939329333D0/t45-0.1306666666666667D-1
     #*t117*sigmaaa*t15-0.112D-1*t36*t58-0.84D-2*t8*t124*t125+0.42D-2
     #*t8*t40*(0.784D-1*t9*t117*t11+0.168D0*sigmaaa/t45/t132*t54
     #-0.448D-1*t138/t4/t132/t44*t144)
      t154 = 1/t18/t73
      t161 = 1/t68/t28
      t162 = t86**2
      t169 = t62**2
      t175 = sigmabb**2
      t181 = 1/t82/t81
      v2rhob2(i) = -0.4135669939329333D0/t74-0.1306666666666667D-1
     #*t154*sigmabb*t29-0.112D-1*t65*t87-0.84D-2*t22*t161*t162+0.42D-2
     #*t22*t69*(0.784D-1*t23*t154*t25+0.168D0*sigmabb/t74/t169*t83
     #-0.448D-1*t175/t18/t169/t73*t181)
      v2rhoab(i) = 0.D0
      t192 = t7*t40
      v2rhoasigmaaa(i) = 0.56D-2*t35*t15-0.56D-2*t36*t99+0.42D-2*t192
     #*t57-0.84D-2*t8*t124*t57*t98+0.42D-2*t8*t40*(-0.168D-1*t92*t35
     #*t11-0.504D-1*t47*t54+0.168D-1*sigmaaa/t4/t132/t33*t144)
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      t218 = t21*t69
      v2rhobsigmabb(i) = 0.56D-2*t64*t29-0.56D-2*t65*t111+0.42D-2*t218
     #*t86-0.84D-2*t22*t161*t86*t110+0.42D-2*t22*t69*(-0.168D-1*t104
     #*t64*t25-0.504D-1*t76*t83+0.168D-1*sigmabb/t18/t169/t62*t181)
      t242 = t98**2
      v2sigmaaa2(i) = 0.84D-2*t192*t98-0.84D-2*t8*t124*t242+0.42D-2*t8
     #*t40*(-0.63D-2/t9/sigmaaa*t7*t11+0.63D-2/sigmaaa*t50*t54-0.63D-2
     #/t4/t132/rhoa*t144)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t266 = t110**2
      v2sigmabb2(i) = 0.84D-2*t218*t110-0.84D-2*t22*t161*t266+0.42D-2
     #*t22*t69*(-0.63D-2/t23/sigmabb*t21*t25+0.63D-2/sigmabb*t79*t83
     #-0.63D-2/t18/t169/rhob*t181)
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
      
      
      subroutine rks_x_b88
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     A.D. Becke
c     Density-functional exchange-energy approximation with correct 
c     asymptotic behaviour
c     Phys. Rev. A38 (1988) 3098-3100
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
      t3 = t2*rho
      t5 = 1/t3
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t5*sigma/
     #(1.D0+0.317500104573508D-1*t8*t10)
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
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t6*t14
      t18 = rho**2
      t20 = 1/t2/t18
      t24 = t13**2
      t25 = 1/t24
      t30 = t2**2
      t35 = 1/t30/t18
      t39 = dsqrt(1.D0+0.1587401051968199D1*sigma*t35)
      t40 = 1/t39
      vrhoa(i) = -0.9847450218426965D0*t2+0.705555787941129D-2*t20
     #*sigma*t14+0.2645834204779234D-2*t6*t25*(-0.8466669455293548D-1
     #*t7*t20*t10-0.106673350692263D0*sigma/t30/t18/rho*t40)
      vsigmaaa(i) = -0.2116667363823387D-1*t5*t14
     #+0.5291668409558467D-2*t6*t25*(0.6350002091470161D-1/t7*t5*t10
     #+0.8000501301919725D-1*t35*t40)
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
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t6*t14
      t18 = rho**2
      t20 = 1/t2/t18
      t21 = t20*sigma
      t24 = t13**2
      t25 = 1/t24
      t29 = t18*rho
      t30 = t2**2
      t32 = 1/t30/t29
      t35 = 1/t30/t18
      t38 = 1.D0+0.1587401051968199D1*sigma*t35
      t39 = dsqrt(t38)
      t40 = 1/t39
      t43 = -0.8466669455293548D-1*t7*t20*t10-0.106673350692263D0
     #*sigma*t32*t40
      t44 = t25*t43
      vrhoa(i) = -0.9847450218426965D0*t2+0.705555787941129D-2*t21*t14
     #+0.2645834204779234D-2*t6*t44
      t49 = 1/t7
      t55 = 0.6350002091470161D-1*t49*t5*t10+0.8000501301919725D-1*t35
     #*t40
      t56 = t25*t55
      vsigmaaa(i) = -0.2116667363823387D-1*t5*t14
     #+0.5291668409558467D-2*t6*t56
      t62 = 1/t2/t29
      t69 = 1/t24/t13
      t70 = t43**2
      t77 = t18**2
      t83 = sigma**2
      t89 = 1/t39/t38
      v2rhoa2(i) = -0.6564966812284644D0/t30-0.3292593677058602D-1*t62
     #*sigma*t14-0.1411111575882258D-1*t21*t44-0.5291668409558467D-2
     #*t6*t69*t70+0.2645834204779234D-2*t6*t25*(0.3951112412470322D0
     #*t7*t62*t10+0.106673350692263D1*sigma/t30/t77*t40
     #-0.4515557042823225D0*t83/t2/t77/t29*t89)
      t100 = t5*t25
      v2rhoasigmaaa(i) = 0.2822223151764516D-1*t20*t14
     #-0.705555787941129D-2*t21*t56+0.1058333681911693D-1*t100*t43
     #-0.5291668409558467D-2*t6*t69*t43*t55+0.2645834204779234D-2*t6
     #*t25*(-0.169333389105871D0*t49*t20*t10-0.640040104153578D0*t32
     #*t40+0.3386667782117419D0*sigma/t2/t77/t18*t89)
      t124 = t55**2
      v2sigmaaa2(i) = 0.4233334727646774D-1*t100*t55
     #-0.1058333681911693D-1*t6*t69*t124+0.5291668409558467D-2*t6*t25*
     #(-0.1270000418294032D0/t7/sigma*t5*t10+0.1600100260383945D0
     #/sigma*t35*t40-0.2540000836588064D0/t2/t77/rho*t89)
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
c:X_B88subrend
