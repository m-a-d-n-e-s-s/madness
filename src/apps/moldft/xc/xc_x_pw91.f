c:X_PW91subrstart

c    Generated: Sun Oct 24 15:07:52 BST 2004

      subroutine uks_x_pw91
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
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
      t4 = dsqrt(sigmabb)
      t6 = t4/t3
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.1645530784602056D1*sigmabb*t14)
      t25 = sigmabb**2
      t26 = t11**2
      zk(i) = -0.9305257363491D0*t3*(1.D0+t10+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t17)*sigmabb*t14)/(1.D0+t10
     #+0.1083108625229223D-5*t25/t2/t26/rhob)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t6 = t4/t3
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.1645530784602056D1*sigmaaa*t14)
      t25 = sigmaaa**2
      t26 = t11**2
      zk(i) = -0.9305257363491D0*t3*(1.D0+t10+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t17)*sigmaaa*t14)/(1.D0+t10
     #+0.1083108625229223D-5*t25/t2/t26/rhoa)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t8 = t6/t5
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t19 = dexp(-0.1645530784602056D1*sigmaaa*t16)
      t27 = sigmaaa**2
      t28 = t13**2
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t42 = t40/t39
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t53 = dexp(-0.1645530784602056D1*sigmabb*t50)
      t61 = sigmabb**2
      t62 = t47**2
      zk(i) = -0.9305257363491D0*t5*(1.D0+t12+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t19)*sigmaaa*t16)/(1.D0+t12
     #+0.1083108625229223D-5*t27/t4/t28/rhoa)-0.9305257363491D0*t39*
     #(1.D0+t46+0.1645530784602056D-1*(0.2743D0-0.1508D0*t53)*sigmabb
     #*t50)/(1.D0+t46+0.1083108625229223D-5*t61/t38/t62/rhob)
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
      t4 = dsqrt(sigmabb)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmabb*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmabb
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmabb**2
      t26 = t11**2
      t29 = 1/t2/t26/rhob
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      vrhoa(i) = 0.D0
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rhob
      t50 = dsqrt(1.D0+0.10000117555961D1*t15)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmabb*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhob(i) = -0.12407009817988D1*t2*t23*t33-0.9305257363491D0*t3*(
     #-t43-t53-0.1088885204563779D-1*t57*t17-0.4388082092272149D-1*t20
     #*t46)*t33+0.9305257363491D0*t24*t67*(-t43-t53
     #-0.5776579334555855D-5*t57)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      vsigmabb(i) = -0.9305257363491D0*t3*(t76+t78
     #+0.408331951711417D-2*t29*t17*sigmabb+0.1645530784602056D-1*t19
     #*t14)*t33+0.9305257363491D0*t24*t67*(t76+t78
     #+0.2166217250458446D-5*sigmabb*t29)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmaaa*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmaaa
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmaaa**2
      t26 = t11**2
      t29 = 1/t2/t26/rhoa
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rhoa
      t50 = dsqrt(1.D0+0.10000117555961D1*t15)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmaaa*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhoa(i) = -0.12407009817988D1*t2*t23*t33-0.9305257363491D0*t3*(
     #-t43-t53-0.1088885204563779D-1*t57*t17-0.4388082092272149D-1*t20
     #*t46)*t33+0.9305257363491D0*t24*t67*(-t43-t53
     #-0.5776579334555855D-5*t57)
      vrhob(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      vsigmaaa(i) = -0.9305257363491D0*t3*(t76+t78
     #+0.408331951711417D-2*t29*t17*sigmaaa+0.1645530784602056D-1*t19
     #*t14)*t33+0.9305257363491D0*t24*t67*(t76+t78
     #+0.2166217250458446D-5*sigmaaa*t29)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t7 = 1/t5
      t8 = t6*t7
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t17 = sigmaaa*t16
      t19 = dexp(-0.1645530784602056D1*t17)
      t21 = 0.2743D0-0.1508D0*t19
      t22 = t21*sigmaaa
      t25 = 1.D0+t12+0.1645530784602056D-1*t22*t16
      t26 = t5*t25
      t27 = sigmaaa**2
      t28 = t13**2
      t31 = 1/t4/t28/rhoa
      t34 = 1.D0+t12+0.1083108625229223D-5*t27*t31
      t35 = 1/t34
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t41 = 1/t39
      t42 = t40*t41
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t51 = sigmabb*t50
      t53 = dexp(-0.1645530784602056D1*t51)
      t55 = 0.2743D0-0.1508D0*t53
      t56 = t55*sigmabb
      t59 = 1.D0+t46+0.1645530784602056D-1*t56*t50
      t60 = t39*t59
      t61 = sigmabb**2
      t62 = t47**2
      t65 = 1/t38/t62/rhob
      t68 = 1.D0+t46+0.1083108625229223D-5*t61*t65
      t69 = 1/t68
      zk(i) = -0.9305257363491D0*t26*t35-0.9305257363491D0*t60*t69
      t79 = 0.3360034800657352D-1*t6/t4/t13*t10
      t82 = 1/t14/t13/rhoa
      t86 = dsqrt(1.D0+0.10000117555961D1*t17)
      t87 = 1/t86
      t89 = 0.3360054550205309D-1*sigmaaa*t82*t87
      t93 = t27/t4/t28/t13
      t102 = t34**2
      t103 = 1/t102
      vrhoa(i) = -0.12407009817988D1*t4*t25*t35-0.9305257363491D0*t5*(
     #-t79-t89-0.1088885204563779D-1*t93*t19-0.4388082092272149D-1*t22
     #*t82)*t35+0.9305257363491D0*t26*t103*(-t79-t89
     #-0.5776579334555855D-5*t93)
      t116 = 0.3360034800657352D-1*t40/t38/t47*t44
      t119 = 1/t48/t47/rhob
      t123 = dsqrt(1.D0+0.10000117555961D1*t51)
      t124 = 1/t123
      t126 = 0.3360054550205309D-1*sigmabb*t119*t124
      t130 = t61/t38/t62/t47
      t139 = t68**2
      t140 = 1/t139
      vrhob(i) = -0.12407009817988D1*t38*t59*t69-0.9305257363491D0*t39
     #*(-t116-t126-0.1088885204563779D-1*t130*t53
     #-0.4388082092272149D-1*t56*t119)*t69+0.9305257363491D0*t60*t140*
     #(-t116-t126-0.5776579334555855D-5*t130)
      t149 = 0.1260013050246507D-1/t6*t7*t10
      t151 = 0.1260020456326991D-1*t16*t87
      vsigmaaa(i) = -0.9305257363491D0*t5*(t149+t151
     #+0.408331951711417D-2*t31*t19*sigmaaa+0.1645530784602056D-1*t21
     #*t16)*t35+0.9305257363491D0*t26*t103*(t149+t151
     #+0.2166217250458446D-5*sigmaaa*t31)
      vsigmaab(i) = 0.D0
      t170 = 0.1260013050246507D-1/t40*t41*t44
      t172 = 0.1260020456326991D-1*t50*t124
      vsigmabb(i) = -0.9305257363491D0*t39*(t170+t172
     #+0.408331951711417D-2*t65*t53*sigmabb+0.1645530784602056D-1*t55
     #*t50)*t69+0.9305257363491D0*t60*t140*(t170+t172
     #+0.2166217250458446D-5*sigmabb*t65)
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
      t4 = dsqrt(sigmabb)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmabb*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmabb
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmabb**2
      t26 = t11**2
      t29 = 1/t2/t26/rhob
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      vrhoa(i) = 0.D0
      t36 = t2*t23
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t44 = t11*rhob
      t46 = 1/t12/t44
      t49 = 1.D0+0.10000117555961D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmabb*t46*t51
      t57 = t25/t2/t26/t11
      t62 = -t43-t53-0.1088885204563779D-1*t57*t17
     #-0.4388082092272149D-1*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.5776579334555855D-5*t57
      t70 = t67*t69
      vrhob(i) = -0.12407009817988D1*t36*t33-0.9305257363491D0*t63*t33
     #+0.9305257363491D0*t24*t70
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      t79 = t29*t17
      t85 = t3*(t76+t78+0.408331951711417D-2*t79*sigmabb
     #+0.1645530784602056D-1*t19*t14)
      t90 = t76+t78+0.2166217250458446D-5*sigmabb*t29
      t91 = t67*t90
      vsigmabb(i) = -0.9305257363491D0*t85*t33+0.9305257363491D0*t24*t91
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t107 = 0.7840081201533821D-1*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1680027275102654D0*sigmabb*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4480125399532632D-1*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhob2(i) = -0.4135669939329333D0/t12*t23*t33
     #-0.24814019635976D1*t2*t62*t33+0.24814019635976D1*t36*t70
     #-0.9305257363491D0*t3*(t107+t112-t120+0.9799966841074009D-1*t116
     #*t17-0.4778117666686413D-1*t25*sigmabb/t124/t11*t17
     #+0.1608963433833121D0*t20*t109)*t33+0.18610514726982D1*t63*t70
     #-0.18610514726982D1*t24*t139*t140+0.9305257363491D0*t24*t67*
     #(t107+t112-t120+0.3658500245218708D-4*t116)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t153 = 0.6300065251232535D-2/t4/sigmabb*t5*t8
      t157 = 0.6300102281634954D-2/sigmabb*t14*t51
      t159 = 0.6300176343092764D-2*t29*t118
      t171 = t90**2
      v2sigmabb2(i) = -0.9305257363491D0*t3*(-t153+t157-t159
     #-0.6719227968777768D-2/t124*t17*sigmabb+0.8166639034228341D-2
     #*t79)*t33+0.18610514726982D1*t85*t91-0.18610514726982D1*t24*t139
     #*t171+0.9305257363491D0*t24*t67*(-t153+t157-t159
     #+0.2166217250458446D-5*t29)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmaaa*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmaaa
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmaaa**2
      t26 = t11**2
      t29 = 1/t2/t26/rhoa
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      t36 = t2*t23
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t44 = t11*rhoa
      t46 = 1/t12/t44
      t49 = 1.D0+0.10000117555961D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmaaa*t46*t51
      t57 = t25/t2/t26/t11
      t62 = -t43-t53-0.1088885204563779D-1*t57*t17
     #-0.4388082092272149D-1*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.5776579334555855D-5*t57
      t70 = t67*t69
      vrhoa(i) = -0.12407009817988D1*t36*t33-0.9305257363491D0*t63*t33
     #+0.9305257363491D0*t24*t70
      vrhob(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      t79 = t29*t17
      t85 = t3*(t76+t78+0.408331951711417D-2*t79*sigmaaa
     #+0.1645530784602056D-1*t19*t14)
      t90 = t76+t78+0.2166217250458446D-5*sigmaaa*t29
      t91 = t67*t90
      vsigmaaa(i) = -0.9305257363491D0*t85*t33+0.9305257363491D0*t24*t91
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t107 = 0.7840081201533821D-1*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1680027275102654D0*sigmaaa*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4480125399532632D-1*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhoa2(i) = -0.4135669939329333D0/t12*t23*t33
     #-0.24814019635976D1*t2*t62*t33+0.24814019635976D1*t36*t70
     #-0.9305257363491D0*t3*(t107+t112-t120+0.9799966841074009D-1*t116
     #*t17-0.4778117666686413D-1*t25*sigmaaa/t124/t11*t17
     #+0.1608963433833121D0*t20*t109)*t33+0.18610514726982D1*t63*t70
     #-0.18610514726982D1*t24*t139*t140+0.9305257363491D0*t24*t67*
     #(t107+t112-t120+0.3658500245218708D-4*t116)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t153 = 0.6300065251232535D-2/t4/sigmaaa*t5*t8
      t157 = 0.6300102281634954D-2/sigmaaa*t14*t51
      t159 = 0.6300176343092764D-2*t29*t118
      t171 = t90**2
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(-t153+t157-t159
     #-0.6719227968777768D-2/t124*t17*sigmaaa+0.8166639034228341D-2
     #*t79)*t33+0.18610514726982D1*t85*t91-0.18610514726982D1*t24*t139
     #*t171+0.9305257363491D0*t24*t67*(-t153+t157-t159
     #+0.2166217250458446D-5*t29)
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
      t6 = dsqrt(sigmaaa)
      t7 = 1/t5
      t8 = t6*t7
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t17 = sigmaaa*t16
      t19 = dexp(-0.1645530784602056D1*t17)
      t21 = 0.2743D0-0.1508D0*t19
      t22 = t21*sigmaaa
      t25 = 1.D0+t12+0.1645530784602056D-1*t22*t16
      t26 = t5*t25
      t27 = sigmaaa**2
      t28 = t13**2
      t31 = 1/t4/t28/rhoa
      t34 = 1.D0+t12+0.1083108625229223D-5*t27*t31
      t35 = 1/t34
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t41 = 1/t39
      t42 = t40*t41
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t51 = sigmabb*t50
      t53 = dexp(-0.1645530784602056D1*t51)
      t55 = 0.2743D0-0.1508D0*t53
      t56 = t55*sigmabb
      t59 = 1.D0+t46+0.1645530784602056D-1*t56*t50
      t60 = t39*t59
      t61 = sigmabb**2
      t62 = t47**2
      t65 = 1/t38/t62/rhob
      t68 = 1.D0+t46+0.1083108625229223D-5*t61*t65
      t69 = 1/t68
      zk(i) = -0.9305257363491D0*t26*t35-0.9305257363491D0*t60*t69
      t72 = t4*t25
      t76 = 1/t4/t13
      t79 = 0.3360034800657352D-1*t6*t76*t10
      t80 = t13*rhoa
      t82 = 1/t14/t80
      t85 = 1.D0+0.10000117555961D1*t17
      t86 = dsqrt(t85)
      t87 = 1/t86
      t89 = 0.3360054550205309D-1*sigmaaa*t82*t87
      t92 = 1/t4/t28/t13
      t93 = t27*t92
      t98 = -t79-t89-0.1088885204563779D-1*t93*t19
     #-0.4388082092272149D-1*t22*t82
      t99 = t5*t98
      t102 = t34**2
      t103 = 1/t102
      t105 = -t79-t89-0.5776579334555855D-5*t93
      t106 = t103*t105
      vrhoa(i) = -0.12407009817988D1*t72*t35-0.9305257363491D0*t99*t35
     #+0.9305257363491D0*t26*t106
      t109 = t38*t59
      t113 = 1/t38/t47
      t116 = 0.3360034800657352D-1*t40*t113*t44
      t117 = t47*rhob
      t119 = 1/t48/t117
      t122 = 1.D0+0.10000117555961D1*t51
      t123 = dsqrt(t122)
      t124 = 1/t123
      t126 = 0.3360054550205309D-1*sigmabb*t119*t124
      t129 = 1/t38/t62/t47
      t130 = t61*t129
      t135 = -t116-t126-0.1088885204563779D-1*t130*t53
     #-0.4388082092272149D-1*t56*t119
      t136 = t39*t135
      t139 = t68**2
      t140 = 1/t139
      t142 = -t116-t126-0.5776579334555855D-5*t130
      t143 = t140*t142
      vrhob(i) = -0.12407009817988D1*t109*t69-0.9305257363491D0*t136
     #*t69+0.9305257363491D0*t60*t143
      t146 = 1/t6
      t149 = 0.1260013050246507D-1*t146*t7*t10
      t151 = 0.1260020456326991D-1*t16*t87
      t152 = t31*t19
      t157 = t149+t151+0.408331951711417D-2*t152*sigmaaa
     #+0.1645530784602056D-1*t21*t16
      t158 = t5*t157
      t163 = t149+t151+0.2166217250458446D-5*sigmaaa*t31
      t164 = t103*t163
      vsigmaaa(i) = -0.9305257363491D0*t158*t35+0.9305257363491D0*t26
     #*t164
      vsigmaab(i) = 0.D0
      t167 = 1/t40
      t170 = 0.1260013050246507D-1*t167*t41*t44
      t172 = 0.1260020456326991D-1*t50*t124
      t173 = t65*t53
      t178 = t170+t172+0.408331951711417D-2*t173*sigmabb
     #+0.1645530784602056D-1*t55*t50
      t179 = t39*t178
      t184 = t170+t172+0.2166217250458446D-5*sigmabb*t65
      t185 = t140*t184
      vsigmabb(i) = -0.9305257363491D0*t179*t69+0.9305257363491D0*t60
     #*t185
      t201 = 0.7840081201533821D-1*t6/t4/t80*t10
      t203 = 1/t14/t28
      t206 = 0.1680027275102654D0*sigmaaa*t203*t87
      t210 = t27/t4/t28/t80
      t212 = 1/t86/t85
      t214 = 0.4480125399532632D-1*t210*t212
      t218 = t28**2
      t233 = 1/t102/t34
      t234 = t105**2
      v2rhoa2(i) = -0.4135669939329333D0/t14*t25*t35
     #-0.24814019635976D1*t4*t98*t35+0.24814019635976D1*t72*t106
     #-0.9305257363491D0*t5*(t201+t206-t214+0.9799966841074009D-1*t210
     #*t19-0.4778117666686413D-1*t27*sigmaaa/t218/t13*t19
     #+0.1608963433833121D0*t22*t203)*t35+0.18610514726982D1*t99*t106
     #-0.18610514726982D1*t26*t233*t234+0.9305257363491D0*t26*t103*
     #(t201+t206-t214+0.3658500245218708D-4*t210)
      t256 = 0.7840081201533821D-1*t40/t38/t117*t44
      t258 = 1/t48/t62
      t261 = 0.1680027275102654D0*sigmabb*t258*t124
      t265 = t61/t38/t62/t117
      t267 = 1/t123/t122
      t269 = 0.4480125399532632D-1*t265*t267
      t273 = t62**2
      t288 = 1/t139/t68
      t289 = t142**2
      v2rhob2(i) = -0.4135669939329333D0/t48*t59*t69
     #-0.24814019635976D1*t38*t135*t69+0.24814019635976D1*t109*t143
     #-0.9305257363491D0*t39*(t256+t261-t269+0.9799966841074009D-1
     #*t265*t53-0.4778117666686413D-1*t61*sigmabb/t273/t47*t53
     #+0.1608963433833121D0*t56*t258)*t69+0.18610514726982D1*t136*t143
     #-0.18610514726982D1*t60*t288*t289+0.9305257363491D0*t60*t140*
     #(t256+t261-t269+0.3658500245218708D-4*t265)
      v2rhoab(i) = 0.D0
      t305 = 0.1680017400328676D-1*t146*t76*t10
      t307 = 0.5040081825307963D-1*t82*t87
      t308 = sigmaaa*t92
      t310 = 0.1680047024824737D-1*t308*t212
      v2rhoasigmaaa(i) = -0.12407009817988D1*t4*t157*t35
     #+0.12407009817988D1*t72*t164-0.9305257363491D0*t5*(-t305-t307
     #+t310-0.3266655613691336D-1*t92*t19*sigmaaa
     #+0.1791794125007405D-1*t27/t218/rhoa*t19-0.4388082092272149D-1
     #*t21*t82)*t35+0.9305257363491D0*t99*t164+0.9305257363491D0*t158
     #*t106-0.18610514726982D1*t26*t233*t105*t163+0.9305257363491D0
     #*t26*t103*(-t305-t307+t310-0.1155315866911171D-4*t308)
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      t345 = 0.1680017400328676D-1*t167*t113*t44
      t347 = 0.5040081825307963D-1*t119*t124
      t348 = sigmabb*t129
      t350 = 0.1680047024824737D-1*t348*t267
      v2rhobsigmabb(i) = -0.12407009817988D1*t38*t178*t69
     #+0.12407009817988D1*t109*t185-0.9305257363491D0*t39*(-t345-t347
     #+t350-0.3266655613691336D-1*t129*t53*sigmabb
     #+0.1791794125007405D-1*t61/t273/rhob*t53-0.4388082092272149D-1
     #*t55*t119)*t69+0.9305257363491D0*t136*t185+0.9305257363491D0
     #*t179*t143-0.18610514726982D1*t60*t288*t142*t184
     #+0.9305257363491D0*t60*t140*(-t345-t347+t350
     #-0.1155315866911171D-4*t348)
      t382 = 0.6300065251232535D-2/t6/sigmaaa*t7*t10
      t386 = 0.6300102281634954D-2/sigmaaa*t16*t87
      t388 = 0.6300176343092764D-2*t31*t212
      t400 = t163**2
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(-t382+t386-t388
     #-0.6719227968777768D-2/t218*t19*sigmaaa+0.8166639034228341D-2
     #*t152)*t35+0.18610514726982D1*t158*t164-0.18610514726982D1*t26
     #*t233*t400+0.9305257363491D0*t26*t103*(-t382+t386-t388
     #+0.2166217250458446D-5*t31)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t413 = 0.6300065251232535D-2/t40/sigmabb*t41*t44
      t417 = 0.6300102281634954D-2/sigmabb*t50*t124
      t419 = 0.6300176343092764D-2*t65*t267
      t431 = t184**2
      v2sigmabb2(i) = -0.9305257363491D0*t39*(-t413+t417-t419
     #-0.6719227968777768D-2/t273*t53*sigmabb+0.8166639034228341D-2
     #*t173)*t69+0.18610514726982D1*t179*t185-0.18610514726982D1*t60
     #*t288*t431+0.9305257363491D0*t60*t140*(-t413+t417-t419
     #+0.2166217250458446D-5*t65)
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
      
      
      subroutine rks_x_pw91
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
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
      t4 = dsqrt(sigma)
      t6 = t4/t3
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.261211729852336D1*sigma*t14)
      t25 = sigma**2
      t26 = t11**2
      zk(i) = -0.7385587663820224D0*t3*(1.D0+t10+0.261211729852336D-1*
     #(0.2743D0-0.1508D0*t17)*sigma*t14)/(1.D0+t10
     #+0.272926271249799D-5*t25/t2/t26/rho)
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
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigma*t14
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t29 = 1/t2/t26/rho
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.7385587663820224D0*t24*t33
      t43 = 0.8466757147455043D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rho
      t50 = dsqrt(1.D0+0.1587419712813815D1*t15)
      t51 = 1/t50
      t53 = 0.1066750825533289D0*sigma*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhoa(i) = -0.9847450218426965D0*t2*t23*t33-0.3692793831910112D0
     #*t3*(-t43-t53-0.5487637560595959D-1*t57*t17-0.1393129225879125D0
     #*t20*t46)*t33+0.3692793831910112D0*t24*t67*(-t43-t53
     #-0.2911213559997856D-4*t57)
      t76 = 0.6350067860591282D-1/t4*t5*t8
      t78 = 0.8000631191499664D-1*t14*t51
      vsigmaaa(i) = -0.7385587663820224D0*t3*(t76+t78
     #+0.411572817044697D-1*t29*t17*sigma+0.1044846919409344D0*t19*t14
     #)*t33+0.7385587663820224D0*t24*t67*(t76+t78
     #+0.2183410169998392D-4*sigma*t29)
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
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigma*t14
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t29 = 1/t2/t26/rho
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.7385587663820224D0*t24*t33
      t36 = t2*t23
      t40 = 1/t2/t11
      t43 = 0.8466757147455043D-1*t4*t40*t8
      t44 = t11*rho
      t46 = 1/t12/t44
      t49 = 1.D0+0.1587419712813815D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.1066750825533289D0*sigma*t46*t51
      t56 = 1/t2/t26/t11
      t57 = t25*t56
      t62 = -t43-t53-0.5487637560595959D-1*t57*t17
     #-0.1393129225879125D0*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.2911213559997856D-4*t57
      t70 = t67*t69
      vrhoa(i) = -0.9847450218426965D0*t36*t33-0.3692793831910112D0
     #*t63*t33+0.3692793831910112D0*t24*t70
      t73 = 1/t4
      t76 = 0.6350067860591282D-1*t73*t5*t8
      t78 = 0.8000631191499664D-1*t14*t51
      t79 = t29*t17
      t84 = t76+t78+0.411572817044697D-1*t79*sigma
     #+0.1044846919409344D0*t19*t14
      t85 = t3*t84
      t90 = t76+t78+0.2183410169998392D-4*sigma*t29
      t91 = t67*t90
      vsigmaaa(i) = -0.7385587663820224D0*t85*t33+0.7385587663820224D0
     #*t24*t91
      t107 = 0.395115333547902D0*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1066750825533289D1*sigma*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4515683437631874D0*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhoa2(i) = -0.6564966812284644D0/t12*t23*t33
     #-0.1969490043685393D1*t2*t62*t33+0.1969490043685393D1*t36*t70
     #-0.3692793831910112D0*t3*(t107+t112-t120+0.9877747609072727D0
     #*t116*t17-0.764498826669826D0*t25*sigma/t124/t11*t17
     #+0.1021628098978025D1*t20*t109)*t33+0.7385587663820224D0*t63*t70
     #-0.7385587663820224D0*t24*t139*t140+0.3692793831910112D0*t24*t67
     #*(t107+t112-t120+0.3687537175997285D-3*t116)
      t156 = 0.1693351429491009D0*t73*t40*t8
      t158 = 0.6400504953199731D0*t46*t51
      t159 = sigma*t56
      t161 = 0.3386762578223905D0*t159*t118
      v2rhoasigmaaa(i) = -0.9847450218426965D0*t2*t84*t33
     #+0.9847450218426965D0*t36*t91-0.3692793831910112D0*t3*(-t156
     #-t158+t161-0.6585165072715151D0*t56*t17*sigma
     #+0.5733741200023695D0*t25/t124/rho*t17-0.5572516903516501D0*t19
     #*t46)*t33+0.3692793831910112D0*t63*t91+0.3692793831910112D0*t85
     #*t70-0.7385587663820224D0*t24*t139*t69*t90+0.3692793831910112D0
     #*t24*t67*(-t156-t158+t161-0.2328970847998285D-3*t159)
      t193 = 0.1270013572118256D0/t4/sigma*t5*t8
      t197 = 0.1600126238299933D0/sigma*t14*t51
      t199 = 0.2540071933667929D0*t29*t118
      t211 = t90**2
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(-t193+t197-t199
     #-0.4300305900017772D0/t124*t17*sigma+0.3292582536357576D0*t79)
     #*t33+0.1477117532764045D1*t85*t91-0.1477117532764045D1*t24*t139
     #*t211+0.7385587663820224D0*t24*t67*(-t193+t197-t199
     #+0.8733640679993569D-4*t29)
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

c:X_PW91subrend
