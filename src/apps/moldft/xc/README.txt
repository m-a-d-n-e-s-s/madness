     CITATION:

     Functionals were obtained from the Density Functional Repository 
     as developed and distributed by the Quantum Chemistry Group, 
     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
     Paul Sherwood for further information.

     COPYRIGHT:

     Users may incorporate the source code into software packages and
     redistribute the source code provided the source code is not
     changed in anyway and is properly cited in any documentation or
     publication related to its use.

     ACKNOWLEDGEMENT:

     The source code was generated using Maple 8 through a modified
     version of the dfauto script published in:

        R. Strange, F.R. Manby, P.J. Knowles
        Automatic code generation in density functional theory
        Comp. Phys. Comm. 136 (2001) 310-318.

     FTP:

         ftp://ftp.dl.ac.uk/qcg/dft_library

----------------------------------------------------------------------

      subroutine uks_c_DFTfun
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)

      subroutine rks_c_DFTfun
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)

      DFTfun:

      
      x_b3
      x_b88
      x_ft97b
      x_lda
      x_pbe
      x_pw91

      c_ft97
      c_lyp
      c_p86
      c_pbe
      c_pw91
      c_pw92
      c_pz81
      c_vwn5
      c_vwn5rpa

      xc_b3lyp
      xc_b97
      xc_b97_1
      xc_b97_2
      xc_edf1
      xc_ft97
      xc_hcth
      xc_hcth147
      xc_hcth407
      xc_pbe
      xc_pw91

