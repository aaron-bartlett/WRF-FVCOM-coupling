!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!==============================================================================!
!  VERSION 4.1
!==============================================================================!

PROGRAM FVCOM

! CODING STATUS: Beta release
! Remaining Work:
! Add variable, 1/ART, 1/ART2... etc
! Fix bottom boundary condition for scalars
! Add sediment 
! Add assimilation modes
! Test GOTM Interface
!

  !================================================================================!
  !                                                                                !
  !                           USG-FVCOM                                            !
  !    The Unstructured Grid Finite Volume Coastal Ocean Model                     !
  !                                                                                !
  !    The USG-FVCOM (publicly called FVCOM) was developed by Drs. Changsheng Chen !
  !  and Hedong Liu at the Marine Ecosystem Dynamics Modeling Laboratory at the    !
  !  School of Marine Science and Technology (SMAST), University of Massachusetts- !
  !  Dartmouth (UMASSD) and Dr. Robert C. Beardsley at the Department of Physical  !
  !  Oceanography, Woods Hole Oceanographic Institution (WHOI). This code was      !
  !  rewritten in FORTRAN 90/2K, modularized, rendered somewhat understandable,    !
  !  and parallelized by Geoff Cowles at SMAST/UMASSD. FVCOM is being upgraded by  !
  !  the UMASSD-WHOI joint team led by Dr. Chen and Dr. Beardsley.                 ! 
  !                                                                                !
  !    The Development was initially supported by the Georgia Sea Grant College    !
  !  Program for the study of the complex dynamic system in Georgia estuaries. The !
  !  code improvement has been supported by Dr. Chen's research grants received    !
  !  from NSF and NOAA Coastal Ocean research grants and SMAST fishery research    !
  !  grants.                                                                       !
  !                                                                                !
  !    FVCOM is a prognostic, unstructured grid, Finite-Volume free-surface three- !
  !  dimensional hydrostatic/non-hydrostatic primitive equations, Coastal Ocean    !
  !  circulation Model (Chen et al., 2003; 2006, 2007). The model computes the     !
  !  momentum, continuity, temperature salinity, and density equations and is      !
  !  closed physically and mathematically using user-specified turbulent parameters!
  !  or turbulent closure submodel. The irregular bottom slope is represented using!
  !  a generalized terrain-following coordinate transformation with spatially      !
  !  variable vertical distribution (Pietrzak et al., 2002). The horizontal grids  !
  !  comprise unstructured triangular cells. The spatial fluxes of momentum are    !
  !  discretized using a second-order accurate finite-volume method (Kobayashi et  !
  !  al. 1999) in conjunction with a vertical velocity adjustment enforce exact    !
  !  volume conservation in individual control volumes. The flux of scalars (e. g. !
  !  temperature, salinity) is calculated by the second-order accurate upwind      !
  !  scheme. The finite-volume method (FVM) used in this model combines the        !
  !  advantages of the finite-element method (FEM) for geometric flexibility and   !
  !  the finite-difference method (FDM) for simple discrete computation. Current,  !
  !  temperature, and salinity in the model are computed in the integral form of   !
  !  the equations, which provides a better representation of the conservation laws!
  !  for mass, momentum, and heat in the coastal region with complex geometry.     !
  !                                                                                !
  !    FVCOM was originally developed for the coastal and estuarine applications.  !
  !  With implementation of spherical coordinates, this model is capable of using  !
  !  for the global ocean.                                                         !
  !                                                                                !
  !    All users should read this agreement carefully.  A user, who receives any   !
  !  version of the source code of FVCOM, must accept all the terms and conditions !
  !  of this agreement and also agree that this agreement is like any written      !
  !  negotiated agreement signed by you. You may be required to have another       !
  !  written agreement directly with Dr. Changsheng Chen at SMAST/UMASSD that      !
  !  supplements all or portions of this agreement. Dr. Changsheng Chen, leader of !
  !  the FVCOM development team, owns all intellectual property rights to the      !
  !  software. The University of Massachusetts-Dartmouth owns the copyright of the !
  !  software. All copyrights are reserved. Unauthorized reproduction and re-      !
  !  distribution of this program are expressly prohibited. This program is only   !
  !  permitted for use in non-commercial academic research and education.          !
  !  Commercial use must be approved by Dr. Chen. Registration is required for all !
  !  new users. Users should realize that this model software is a research product!
  !  without any warranty. Users must take full responsibility for any mistakes    !
  !  that might be caused by any inappropriate modification of the source code.    !
  !  Modification is not encouraged for users who do not have a deep understanding !
  !  of the finite-volume numerical methods used in FVCOM. Contributions made to   !
  !  correcting and modifying the programs will be credited, but will not affect   !
  !  copyrights. No duplicate configurations of FVCOM are allowed in the same      !
  !  geographical region, particularly in the regions where FVCOM has been already !
  !  been applied. Users who want to use FVCOM in a region that the SMAST/UMASS    !
  !  Marine Ecosystem Dynamics Modeling (MEDM) group (led by Dr. Chen) is working  !
  !  on must request permission from Dr. Chen. No competition is allowed in the    !
  !  same region using FVCOM, especially with Dr. Chen's group. FVCOM has been     !
  !  validated for many standard model test cases.  Users are welcome to do any    !
  !  further model validation experiments. These experiments shall be carried out  !
  !  in collaboration with the SMAST/UMASSD model development team. To avoid or    !
  !  reduce deriving any incorrect conclusions due to an inappropriate use of      !
  !  FVCOM, users are required to contact the scientific leaders of the FVCOM      !
  !  development team (Dr. Chen at SMAST/UMASSD and Dr. Beardsley at WHOI) before  !
  !  any formal publications are prepared for model validation.                    !
  !                                                                                !
  !    For public use, all users should name this model as "FVCOM". In any         !
  !  publications with the use of FVCOM, acknowledgement must be included. The     !
  !  rationale behind this FVCOM distribution policy is straightforward.  New      !
  !  researchers and educators who want to use FVCOM and agree to the above        !
  !  requirements get free access to the latest version of FVCOM and the collective!
  !  wisdom and experience of the FVCOM development team and existing users.       !
  !  Problems arising in new FVCOM applications, both related to conceptual as well!
  !  as numerical and coding issues, can be shared with the development team and   !
  !  other users who can work together on physics and code improvements that over  !
  !  time will lead to a better FVCOM.                                             !
  !                                                                                !
  !    FVCOM has been developed to date with state and federal funding with the    !
  !  idea that FVCOM will become a community model that new users start to use the !
  !  model and its scientific usefulness and numerical accuracy and efficiency     !
  !  continue to improve. The FVCOM distribution policy is designed to encourage   !
  !  this transition while maintaining a central core group responsible for overall!
  !  FVCOM development direction, implementing official code improvements, and     !
  !  maintaining well tested and documented updated code versions.                 !
  !                                                                                !
  !                                                                                !
  !  External forces used to drive this model:                                     !
  !                                                                                !
  !  1) Tidal amplitudes and phases at open boundaries (initial designs include 8  !
  !         tidal constituents, more can be added as needed) plus astronomical     !
  !         forcing implemented via gradients of tidal-generating potential for    !
  !         each tidal constituent and various corrections due to the Earth tidal  !
  !         and ocean loading.                                                     !
  !  2) Wind stress [3 ways: a) uniform wind speed and direction, b) spatially     !
  !         distributed wind velocity field, and c) the weather model output wind  !
  !         fields (NCEP, MM5, WRF, ECMWF).                                        !
  !  3) Air pressure gradients [2 ways: weather model output and Hurricane or      !
  !         typhoon model predicted.                                               !
  !  4) Surface heat flux [3 ways: a) uniform heat flux, b) spatially distributed  !
  !         heat flux, and c) the weather model-output heat flux fields. All the   !
  !         surface net heat flux and short-wave radiation are needed in the input !
  !         file.                                                                  !
  !  5) Precipitation via evaporation.                                             !
  !  6) River discharges: specify the location and discharge volume, temperature,  !
  !         and salinity.                                                          !
  !  7) Groundwater input: currently diffused bottom flux only.                    !
  !                                                                                !
  !  Note: FVCOM can be run directly by coupling with the weather model (WRF or    !
  !        MM5). This is used in the US Northeast Coastal Ocean Forecast System    !
  !        (NECOFS) (http://fvcom.smast.umassd.edu). Modes use for air-sea coupling!
  !        is not included in this release.                                        !
  !                                                                                !
  !  Hydrodynamics initial conditions:                                             !
  !                                                                                !
  !  The model can be prognostically run for both barotropic and baroclinic cases. !
  !                                                                                !
  !  Tidal forcing can be added into the system with zero velocity field at initial!
  !  or specified the 3-D tidal initial velocity field using the model-predicted   !
  !  harmonic tidal currents.                                                      !
  !                                                                                !
  !  Initial fields of temperature and salinity needed to be specified by using    !
  !  either limatological field, real-time observed field or idealized functions.  !
  !  The model has included Gregorian time for the time simulation for tidal       !
  !  currents.                                                                     !
  !                                                                                !
  !  For the purpose of interdisciplinary studies, biological, chemical, and       !
  !  sediment suspension models are available for FVCOM.  These submodels are      !
  !  directly driven by the FVCOM physical model. These submodels are called       !
  !  offline FVCOM modes.  A description of these submodels follows.               !
  !                                                                                !
  !  Generalized Biological Model (GBM)-a software platform that allows users to   !
  !        build his own biological model in FVCOM (developed by a team including  !
  !        Drs. Chen, Tian and Qi)                                                 !
  !                                                                                !
  !  NPZ model-built using GBM for 3 component nutrient-phytoplankton-zooplankton. !
  !                                                                                ! 
  !  NPZD model-an 4 and 8 component nutrient-phytoplankton-zooplankton-detritus   !
  !        model constructed by GBM                                                !
  !                                                                                !
  !  NPZDB-model-a 9 phosphorus-controlled component nutrient-phytoplankton-       !
  !        zooplankton-detritus-bacteria model constructed by GBM                  !
  !                                                                                !
  !  FVCOM-WQM: The water quality model modified from US-EPA WASP with inclusion of!
  !        the benthic process. This code was converted from the structured grid   !
  !        WQM (developed by Zheng and Chen) to the unstructured grid version.     !
  !                                                                                !
  !  UG-RCA: The unstructured grid Row/column version of the water quality model   !
  !        ESOP (RCA). UG-RCA developed by the FVCOM development team (Drs. Chen,  !
  !        J. Qi and R. Tian) on converting HydroQual's structured grid RCA to     !
  !        unstructured grid finite-volume version.                                !
  !                                                                                !
  !  UG-CE-QUAL-ICM: A unstructured grid version of the army corp water quality    !
  !        model. The code was written by Drs. Qi and Chen at UMASSD and modified  !
  !        and tested by Drs. Kim, Labiosa, Khanhaonkar and Yang at PNL.           !
  !                                                                                !
  !  FVCOM-SED: The 3-D sediment module (developed by Cowles based on the U.S.G.S. !
  !        national community sediment transport model and modified by the FVCOM   !
  !        development team staff-Ge, Wu, Chen and FVCOM users to include cohesive !
  !        process). This code can run online coupled with FVCOM or offline driven !
  !        by the FVCOM output                                                     ! 
  !                                                                                ! 
  !  FVCOM-SWAVE: The unstructured-grid (UG) version of the Simulating Wave        !
  !        Nearshore (SWAN) surface wave model. The code was developed by Drs. Qi  !
  !        and Chen at SMAST/UMASSD (Qi et al. 2008)                               !
  !                                                                                !
  !  Lagrangian particle tracking:                                                 !
  !                                                                                !
  !  FVCOM-LAG: A bilinear interpolation scheme is used to determine the particle  !
  !        velocity for the Lagrangian particle tracking. A random walk process    !
  !        also could be included with a specified function related to horizontal  !
  !        and vertical diffusion coefficients                                     !
  !                                                                                !
  !  Key references:                                                               !
  !                                                                                !
  !  Chen, C., H. Liu, and R. C. Beardsley, 2003. An unstructured grid, finite-    !
  !      volume, three-dimensional, primitive equations ocean model: application to!
  !      coastal ocean and estuaries, Journal of Atmospheric and Oceanic           !
  !      Technology,  20, 159-186.                                                 !
  !                                                                                !
  !  Chen, C, R. C. Beardsley and G. Cowles, 2006. An unstructured grid, finite-   !
  !      volume coastal ocean model (FVCOM) system. Special Issue entitled Advance !
  !      in Computational Oceanography, Oceanography, 19(1), 78-89.                !
  !                                                                                !
  !  Chen, C. H. Huang, R. C. Beardsley, H. Liu, Q. Xu and G. Cowles, 2007. A      !
  !      finite-volume numerical approach for coastal ocean circulation studies:   !
  !      comparisons with finite difference models. J. Geophys. Res. 112, C03018,  !
  !      doi:10.1029/2006JC003485.                                                 !
  !                                                                                !
  !  Chen, C., P. Malanotte-Rizzoli, J. Wei, R. C. Beardsley, Z. Lai, P. Xue, S.   !
  !      Lyu, Q. Xu, J. Qi, and G. W. Cowles, 2009. Application and comparison of  !
  !      Kalman filters for coastal ocean problems: An experiment with FVCOM, J.   !
  !      Geophys. Res., 114, C05011, doi:10.1029/2007JC004548.                     !
  !                                                                                !
  !  Cowles, G., 2008. Parallelization of the FVCOM Coastal Ocean Model,           !
  !      International Journal of High Performance Computing Applications, 22(2),  !
  !      177-193.                                                                  !
  !                                                                                !
  !  Huang, H. C. Chen, G. Cowles, C. D. Winant, R. C. Beardsley, K. Hedstrom, D.  !
  !      B. Haidvogel, 2008. FVCOM validation experiments: comparisons with ROMS   !
  !      for three idealized test problems. J. Geophys. Res, 113, C07042, doi:     !
  !      10.1029/2007JC004557.                                                     !
  !                                                                                !
  !  Lai, Z, C. Chen, G. Cowles and R. C. Beardsley, 2009. A non-hydrostatic       !
  !      version of FVCOM, part I: validation experiments. J. Geophys. Res.,       !
  !      in revision.                                                              !
  !                                                                                !
  !  Lai, Z., C. Chen, G. Cowles and R. C. Beardsley, 2009. A non-hydrostatic      !
  !      version of FVCOM, part II: mechanistic study of tidally generated         !
  !      nonlinear internal waves in Massachusetts Bay. J. Geophys. Res.,          !
  !      submitted.                                                                !
  !                                                                                !
  !  Qi, J. C. Chen, R. C. Beardsley, W. Perrie G. Cowles and Z. Lai, 2009. An     !
  !      unstructured-grid finite-volume surface wave model (FVCOM-SWAVE):         !
  !      implementation, validations and applications. Ocean Modelling, 28,        !
  !      153-166. doi:10.1016/j.ocemod.2009.01.007.                                !
  !                                                                                !
  !  Please direct criticisms and suggestions to                                   !
  !                                                                                !
  !               Changsheng Chen                                                  !
  !               School for Marine Science and Technology                         ! 
  !               University of Massachusetts-Dartmouth                            !
  !               New Bedford, MA 02742                                            !
  !               Phone: 508-910-6388, Fax: 508-910-6371                           !
  !               E-mail: c1chen@umassd.edu, cchen@whoi.edu                        !
  !               Web: http://fvcom.smast.umassd.edu                               !
  !                                                                                !
  !  What are new for version 3.1?                                                 !
  !                                                                                ! 
  !  1) Fully sea-ice coupling. The sea ice model was developed by converting the  !
  !     structured grid CICE into the unstructured grid version called UG-CICE).   !
  !     This is one component of G. Gao's Ph.D. thesis work supervised by C. Chen, !
  !     R. C. Beardsley and A. Proshutinsky. The code was written by G. Gao with   !
  !     assistance of J. Qi and C. Chen and tested by G. Gao and C. Chen.          !
  !                                                                                !
  !  2) Non-hydrostatic version. The non-hydrostatic dynamics is added into FVCOM  !
  !     as one component of Z. Lai's Ph.D. thesis work supervised by C. Chen, G.   !
  !     Cowles and R. C. Beardsley.  The code was written by Z. Lai with assistance!
  !     of C. Chen and G. Cowles, and tested by Z. Lai and C. Chen.                !
  !                                                                                ! 
  !  3) Fully current-wave interaction. An unstructured grid surface wave model    !
  !     (called FVCOM-SWAVE) was developed by Qi et al. (2008). This wave model    !
  !     is modified from SWAN by converting it to unstructured grid version using  !
  !     the second-order accurate FVCOM algorithms. The implementing FVCOM-SWAVE   !
  !     into FVCOM was initialized by A. Wu at SMAST/UMASSD as one component of his!
  !     Ph.D. thesis work. The work is supervised by C. Chen at SMAST/UMASSD. The  !
  !     code is re-organized, modified and re-debugged to fit more flexible and    !
  !     general setup of FVCOM by J. Qi and C. Chen.                               !
  !                                                                                !
  !  4) Fully current-wave-sediment interaction. Coupling of wave-sediment was     !
  !     initialized by A. Wu and J. Ge at SMAST/UMASSD as a visiting student       !
  !     project. The code was re-organized, modified and re-debugged by J. Qi and  !
  !     C. Chen when it is implemented into version 3.1.                           !
  !                                                                                !
  !  5) NetCDF Input and output. Version 3.1 uses the NetCDF input and output      !
  !     implemented into FVCOM by David Stuebe at SMAST/UMASSD. The code is also   !
  !     modified to improve the efficiency of inter-node data exchange and model   !
  !     output writing.                                                            !
  !                                                                                !
  !  6) Semi-implicit solver. In addition to the mode-split version, a semi-       !
  !     implicit solver was implemented into FVCOM by Z. Lai and C. Chen. Using    !
  !     this option requires the installation of PETCs.  With this implementation, !
  !     FVCOM can be run by choosing either mode-split scheme or semi-implicit     !
  !     scheme. The semi-implicit version of FVCOM can significantly increase the  !
  !     computational efficiency for regional and global scale application.        !
  !                                                                                !
  !  7) Nesting software module.  This module is built in FVCOM to allow multi-sub-!
  !     regional FVCOM domains to be run simultaneously through nesting approach.  !
  !     The module was originally written by P. Xue at SMAST/UMASSD when he and    !
  !     Chen applied FVCOM Kalman Filter to conduct The Observing System Simulation!
  !     Experiments (OSSEs) in Nantucket Sound. This module was modified by Dave   !
  !     Stuebe at SMAST/UMASSD and implemented into the version 3.0 of FVCOM, and  !
  !     graded by J. Qi and C. Chen at SMAST/UMASSD for current-wave interaction   !
  !     modules in in version 3.1.                                                 !
  !                                                                                !
  !  8) Dike-Groyne module. This module is implemented into FVCOM to deal with the !
  !     vertical wall under or above the sea surface. The algorithm was derived by ! 
  !     C. Chen at SMAST/UMASSD, and the parallelized code was written and tested  !
  !     by J. Qi, J. Ge and C. Chen at SMAST/UMASSD.                               !
  !                                                                                !
  !  9) SST/SSH assimilation module.  A SSH assimilation module is added to include!
  !     the satellite-altimeter data into FVCOM. The algorithm is the same as that !
  !     used for SST assimilation. The SST/SSH is merged together by Q. Xu, Z. Lai !
  !     and C. Chen at SMAST/UMASSD and implemented into version 3.1. C. Chen      !
  !     improved the SST algorithm to avoid the unreal thin layer caused by SST    !
  !     nudging.                                                                   !
  !                                                                                ! 
  !  10) Kalman Filter Assimilation modules. The Kalman Filter package is upgraded !
  !     by adding a Signular Evolutive Interpolated Kalman filter (SEIK). The      !
  !     source code of SEIK code is provided by Dr. Lars Nerger at AWI, Germany.   !
  !     In addition to Reduced Rank Kalman Filter (RRKF), Ensemble Kalman Filter   !
  !     (EnKF), we also reconstruct Ensemble Transform Kalman Filter (EnTKF) for   !
  !     the OSSEs application. The Kalman Filter assimilation modules were         !
  !     originally developed by a team effort led by C. Chen at SMAST/UMASSD and P.!
  !     Rizzoli/MIT. The team members include P. Xue, Z. Lai, Q. Xu and J. Qi at   !
  !     SMAST/UMASSD, J. Wei and S. Lyu at MIT and R. C. Beardsley at WHOI. The    !
  !     package is upgraded and implemented into version 3.1 by P. Xue, J.Qi and   !
  !     C.Chen. Dr. Lars Nerger owns the copyright of SEIK method used in FVCOM.   !                                   !
  !                                                                                !
  !  11) Optimal Interpolation Assimilation Module. This module is implemented into!
  !     FVCOM to integrate vertical profiles of hydrographic/current data into the !
  !     assimilation. The code was written by Q. Xu, J. Qi and C. Chen at          !
  !     SMAST/UMASSD.                                                              !
  !                                                                                ! 
  !  12) The pre-process programs to convert the input from the old version of     !
  !     FVCOM to version 3.0 or up.                                                !
  !                                                                                !
  !  13) Bugs in the wet/dry treatments on heat flux at the wet/dry boundary are   !
  !     corrected. The bug was detected by T. Kim and fixed by C. Chen and J. Qi at!
  !     SMAST/UMASSD. The second-order limiter for advection terms in the wet/dry  !
  !     treatment was originally coded in Fortran 77 version of FVCOM. This limiter!
  !     is added back to updated version 2.6 or up of FVCOM by C. Chen and J. Qi at!
  !     SMAST/UMASSD.                                                              !
  !                                                                                !
  !     Bugs in selecting General Turbulence Module is corrected by J. Qi and C.   !
  !  Chen at SMAST/UMASSD.                                                         !
  !                                                                                !
  !     Non-general definitions and commands in arrays, variables and NetCDF       !
  !  modules are changed by G. Cowles and Z. Lai to allow to use the version 3.1 of!
  !  FVCOM on IBM super computer. These problems only occur in version 3.0 or up.  !
  !  The bugs in PETCs were detected by Z. Lai, C. Chen and J. Qi at SMAST/UMASSD  !
  !  and corrected by the PETCs development team.                                  ! 
  !                                                                                ! 
  !    Bugs in the Lagrangian-tracking in the generalized terrain-following        !
  !  coordinates are corrected. Many minor bugs are also corrected.                !
  !                                                                                !
  !  Enjoy!      C. Chen at SMAST/UMASSD at Dec. 8 2009.                           !
  !================================================================================!

  !==============================================================================!
  !  INCLUDE MODULES                                                             !
  !==============================================================================!

  USE MOD_UTILS
  USE CONTROL

  USE MOD_PAR  

  USE MOD_STARTUP

  USE MOD_TIME
  USE MOD_CLOCK

  USE MOD_INPUT
  USE MOD_NCDIO
  USE MOD_NCLL

  USE MOD_SETUP
  USE MOD_SET_TIME

  USE MOD_FORCE
  USE MOD_OBCS

# if defined (DATA_ASSIM)
  USE MOD_ASSIM
# endif  
  USE MOD_NESTING

  USE MOD_REPORT
# if defined (LAG_PARTICLE)
  USE MOD_LAG
# endif
  USE PROBES
  USE MOD_BOUNDSCHK !bounds checking
# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT
# endif
# if defined (BALANCE_2D)
  USE MOD_BALANCE_2D
# endif  
# if defined (ONE_D_MODEL)
  USE MOD_ONEDTIDE
# endif
# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif    
# if defined (BioGen)
  USE MOD_BIO_3D
# endif    
# if defined (DYE_RELEASE)
  USE MOD_DYE
# endif
# if defined (RRKF)
   USE RRKVAL
   USE MOD_RRK
   USE MOD_RRKA
   USE MOD_RRKF_OBS
# endif
# if defined (ENKF)  
   USE ENKFVAL 
   USE MOD_ENKF
   USE MOD_ENKFA
   USE MOD_ENKF_OBS
# endif


# if defined (GOTM)
  USE MOD_GOTM 
# endif  
  
# if defined (SEDIMENT)
# if defined (ORIG_SED)
  USE MOD_SED 
# elif defined (CSTMS_SED)
  USE MOD_SED_CSTMS 
# endif
# endif

# if defined (MEAN_FLOW)
  USE MOD_MEANFLOW
  USE MOD_OBCS2
  USE MOD_OBCS3
# endif

# if defined (ICE)
  USE MOD_ICE
# endif

# if defined (NH)
  USE NON_HYDRO
# endif

# if defined (SEMI_IMPLICIT) || defined (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
  USE MOD_PETSc, ONLY : PETSc_SET, PETSc_CLEANUP
# endif

! for periodic lateral boundary conditions
# if defined (PLBC)
  USE MOD_PERIODIC_LBC
# endif

# if defined (WAVE_CURRENT_INTERACTION)
  USE TIMECOMM
  USE SWCOMM3
  USE VARS_WAVE
  USE MOD_WAVE_CURRENT_INTERACTION
# endif 


# if defined (WAVE_CURRENT_INTERACTION) && (SEDIMENT)
  USE MOD_BBL
# endif

# if defined (THIN_DAM)
  USE MOD_DAM
# endif

# if defined (PWP)
  USE MOD_PWP
# endif

  USE MOD_STATION_TIMESERIES 
  USE MOD_SPARSE_TIMESERIES

!----ch201610_oasis
# if defined (oasis_coupler)
  use mod_var_cpl
  use ALL_VARS
  use mod_oasis
# endif
!----end

  !------------------------------------------------------------------------------|
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string

# if defined (RRKF)
  CHARACTER(LEN=120):: RRKFILES
  CHARACTER(LEN=8)  :: RRKINP1
  CHARACTER(LEN=4)  :: RRKINP2
  CHARACTER(LEN=100)::  MKDIR
# endif
  INTEGER :: IERR

  type(watch) Timer
  
  type(TIME) ::GET_BEGIN
  integer status
!--ch2018
#  if defined (oasis_coupler)
  integer :: tmpwid,iwd
  integer :: count_sst,count_aice
  !character(len=*),parameter :: cpl_outfmt="(a,i10,i10,i5,i5)"
  character(len=*),parameter :: cpl_outfmt="(a,f12.2,i10,i5,i5)"
# endif
!--end
  !==============================================================================!
  ! INITIALIZE ALL CONTROL VARIABLES
  !==============================================================================!
  CALL INITIALIZE_CONTROL("FVCOM")

#  if defined (MULTIPROCESSOR)
  ! INTIALIZE MPI CONTROL VARIABLES
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = il_commlocal ! FOR NOW MAKE THEM EQUAL
#  endif

  !==============================================================================!
  !   INITIALIZE A STOP WATCH TIMER FOR TESTING SUBROUTINE EFFICENCY             !
  !==============================================================================!
  CALL WATCH_INIT(TIMER)

  !==============================================================================!
  !   IMPORT CASENAME AND COMMAND LINE ARGUMENTS AND START LOG FILE              !
  !==============================================================================!
  CALL COMMANDLINEIO(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)       
  if(DBG_SET(dbg_log)) Call WRITE_BANNER(PAR,NPROCS,MYID)

  !==============================================================================!
  ! SET DEFAULT VALUES AND READ NAME LISTS                                            
  !==============================================================================!

  CALL NAMELIST

  !==============================================================================!
  !   SET MODEL CONTROL PARAMTERS BASED ON NAME LIST HERE                        !
  !==============================================================================!
  CALL CNTRL_PRMTRS

  !==============================================================================!
  !   SET THE STARTUP TYPE TO BE USED!                                           !
  !==============================================================================!
  CALL SET_STARTUP_TYPE ! see: startup_type.F

  !==============================================================================!
  !   OPEN ALL FILES NEEDED BASED ON THE RUN PARAMETERS                          !
  !==============================================================================!
  CALL OPEN_ALL

  !==============================================================================!
  !   SET MODEL TIME BASED ON THE NAMELIST TIME STRINGS OR RESTART FILE          !
  !==============================================================================!
  CALL SETUP_TIME

  !==============================================================================!
  !   LOAD GRID CONNECTIVITY AND OBC LIST FOR METIS DECOMPOSITION                !
  !==============================================================================!
  CALL LOAD_GRID

  !==============================================================================!
  !   SETUP THE DOMAIN FOR PARALLEL OR SERIAL RUNNING                            !
  !==============================================================================!
  CALL SETUP_DOMAIN

  !==============================================================================!
  !   ALLOCATE ALL DOMAIN SIZE VARIABLES HERE                                    !
  !==============================================================================!
  CALL ALLOCATE_ALL
# if defined (DYE_RELEASE)
  CALL ALLOC_VARS_DYE
# endif
# if defined (BALANCE_2D)
  CALL ALLOC_BALANCE_VARS
# endif    
  !==============================================================================!
  !   LOAD/SETUP PHYSICAL QUANTITIES (CORIOLIS, GRAVITY, SPONGE LAYER, XY/LATLON)!
  !==============================================================================!
  CALL COORDS_N_CONST

  !==============================================================================!
  ! CALCULATE GRID METRICS - NEIGHBORS, GRADIENTS, CELL AREA, INTERP COEFF'S     !
  !==============================================================================!
  CALL GRID_METRICS

  !==============================================================================!
  ! SETUP THE SEDIMENT MODEL (MUST COME BEFORE SETUP_FORCING)                    ! 
  !==============================================================================!
#  if defined (SEDIMENT)
   IF(SEDIMENT_MODEL) &
  &  CALL SETUP_SED(SEDIMENT_MODEL_FILE,IOBCN,IOBCN_GL,N_SED,N_SED_MAX,SED_NAMES)
#  endif


!JQI  !==============================================================================!
!JQI  !  SETUP THE MODEL FORCING                                                     !
!JQI  !==============================================================================!
!JQI  CALL SETUP_FORCING

  !==============================================================================!
  !  GET THE PARAMETERS OF BIOLOGICAL MODEL                                      !
  !==============================================================================!
#  if defined (BioGen)
   IF(BIOLOGICAL_MODEL)THEN
     KBV=KB
!     CALL GET_PARAMETER("./"//trim(input_dir)//"/"//trim(BIOLOGICAL_MODEL_FILE))
     CALL GET_PARAMETER_NEW("./"//trim(input_dir)//"/"//trim(BIOLOGICAL_MODEL_FILE))
     CALL BIO_PARAMETER_REPORT
   END IF
!   RIVER_TS_SETTING = 'NONE' !No bio river yet
#  endif

  !==============================================================================!
  !  SETUP THE MODEL FORCING                                                     !
  !==============================================================================!
  CALL SETUP_FORCING

  !==============================================================================!
  !  SETUP OTHER TOOLS, MODELS AND DATA ASSIMILATION                             !
  !==============================================================================!

# if defined (ICE)
  IF(ICE_MODEL) CALL ICE_INIT_0
# endif

  !  SETUP PETSc FOR SEMI_IMPLICIT AND NON-HYDROSTATIC MODULE
# if defined (SEMI_IMPLICIT) || (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
  CALL PETSc_SET
# endif

#  if defined (WAVE_CURRENT_INTERACTION) 
#  if defined (WAVE_OFFLINE)
   CALL WAVE_CURRENT_SETUP
#  else
   CALL SWMAIN_SETUP
   CALL WAVE_CURRENT_SETUP
#  endif
#  endif  

#  if defined (PLBC)
   CALL FIND_NODE_CELL
#  endif

# if defined (WAVE_CURRENT_INTERACTION) && (SEDIMENT)
   CALL INIT_BBL
# endif


  ! SETUP DATA ASSIMILATION MODE
# if defined (DATA_ASSIM) 
  CALL SETUP_DATA_ASSIMILATION
# if defined (PWP)
  CALL SETUP_PWP
# endif
# endif
#  if defined (ENKF)
   if(msr)  print *, 'before enkf_set_time'
     CALL ENKF_SET_TIME
     CALL SET_ASSIM_ENKF_EVE
        IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
   if(msr)  print *, 'finish enkf_set_time'
     !check point done
#  endif
! New Open Boundary Condition ----2
#  if defined (MEAN_FLOW)
     CALL FIND_OBSIDE
     IF(MSR)WRITE(IPT,*)'FIND_OBSIDE COMPLETED....'
     CALL ALLOC_OBC3_DATA
     IF(MSR)WRITE(IPT,*)'OBC3_DATA COMPLETED....'
     CALL SETUP_OBC3
     IF(MSR)WRITE(IPT,*)'SETUP_OBC3 COMPLETED....'
     CALL ALLOC_OBC2_DATA
     IF(MSR)WRITE(IPT,*)'OBC2_DATA COMPLETED....'
     MF_RST_STCNT = 0
#    endif
#  if defined (RRKF)

     CALL RRK_SET_TIME
     print *, 'finish rrkf_set_time'
!=============================================================
! THIS IS NOT STANDARD FORTRAN - REMOVE BEFORE PUBLIC RELEASE!
     IF(MSR) THEN 
        MKDIR = "mkdir -p "//TRIM(OUTPUT_DIR)//"/rrktemp"
        CALL SYSTEM( TRIM(MKDIR) )         
     ENDIF
# if defined(MULTIPROCESSOR)     
     CALL MPI_BARRIER(MPI_FVCOM_GROUP,ierr)
# endif
!=============================================================
! IF (RRK_ON) then
    print *, 'are you do here?' ! should be marked
     CALL RRK_SET_STARTUP
    print *, 'finish rrk_set_startup'
     CALL READ_SSH
    print *, 'finish read_ssh' 
     IF(WETTING_DRYING_ON) CALL READ_WETDRY
    print *, 'finish read_read_wetdry' 
     CALL READ_UV
    print *, 'finish read_uv' 
     CALL READ_TURB
    print *, 'finish read_turb' 
     CALL READ_TS
    print *, 'finish read_TS'
!  END IF
# else

  !==============================================================================!
  !  SET THE INITIAL CONDITIONS FOR THE MODEL                                    !
  !==============================================================================!
  CALL STARTUP
# endif

!qxu{ test function READ_DATETIME
!  GET_BEGIN = READ_DATETIME('2009-01-01T00:00:00.0','ymd','UTC',status)
!  CALL PRINT_TIME(GET_BEGIN,IPT,'2009-01-01T00:00:00.0')
!qxu}

  !==============================================================================!
  !  CALL ARCHIVE TO SETUP THE OUTPUT AND DUMP CONSTANT VALUES                   !
  !==============================================================================!
!-------------- ch201610 added for coupler-----------------!
#  if defined (oasis_coupler)
  IF(USE_OASIS_COUPLER) THEN
  !---
     print *,'coupler -- fvcom'

      call ini_coupler('ocean1')
      !------ ch201611 added couping ----------------!
!      print *,'coupling period from oasis',namflddti
!      print *,'coupling lag from oasis',namfldlag
!      coupler_period=namflddti(1) ! coupling period (regcm to fvcom) 
!      coupler_period_fvcom=namflddti(nnamcpl) !coupling period(fvcom to regcm)
!      coupler_lag_fvcom=namfldlag(nnamcpl) !lag for last field (fvcom)

   IF(STARTUP_TYPE==STARTUP_TYPE_COLDSTART) THEN
       itap_sec=0 !for initial
       itap_sec_save=itap_sec !for initial
       itap_sec_float=0.0 !for initial

!-- REORDERED (was after the fromcoupler/RECV block below): FVCOM's initial
!-- fromcoupler calls block waiting on data from WRF/waves, but WRF's own
!-- first coupling action is itself a blocking RECV of FVCOM's SST. With
!-- RECV-before-SEND on both sides, neither can ever unblock the other --
!-- a permanent mutual deadlock at t=0. Sending FVCOM's initial state FIRST
!-- (using STARTUP's already-initialized T/T1/AICE/U/V/EL1, which do not
!-- depend on anything received below) breaks the deadlock: WRF's blocking
!-- SST receive can complete immediately once this send happens.
        sst_on_cell(1:N)=T(1:N,1)
        sst_on_node(1:M)=T1(1:M,1)
# if defined (ICE)
       IF(ICE_MODEL) THEN
        aice_on_node(1:M)=AICE(1:M,1)
       END IF
# endif

# if defined (MULTIPROCESSOR)
     IF (PAR) THEN
        CALL AEXCHANGE(NC,MYID,NPROCS,sst_on_node) !--
        CALL AEXCHANGE(NC,MYID,NPROCS,aice_on_node)
        CALL AEXCHANGE(EC,MYID,NPROCS,sst_on_cell)
        CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
        CALL AEXCHANGE(EC,MYID,NPROCS,EL1)
     END IF
# endif
       CALL N2E2D(aice_on_node,aice_on_cell)

             sst_sum=sst_on_node !ch20170102
             aice_sum=aice_on_node
             sst_avg=sst_sum
             aice_avg=aice_sum

        do iwd=1,oasis_nnn2
           field_fvcom_send_sst(iwd,1)=sst_on_node(local_nn(iwd))+273.15 !C to K
           field_fvcom_send_aice(iwd,1)=aice_on_node(local_nn(iwd)) !
        end do
           field_fvcom_send_u(:,1)=u(:,1)
           field_fvcom_send_v(:,1)=v(:,1)
           field_fvcom_send_zeta(:,1)=el1(:)

      if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_ocn2atm-dt/lag",itap_sec_float,0,cpl_dt_ocn2atm,cpl_lag_ocn2atm
        call tocoupler(itap_sec,id_fvcom_send_sst,oasis_nnn2,1,field_fvcom_send_sst)
        call tocoupler(itap_sec,id_fvcom_send_aice,oasis_nnn2,1,field_fvcom_send_aice)
      if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_ocn2wav-dt/lag",itap_sec_float,0,cpl_dt_ocn2wav,cpl_lag_ocn2wav
        call tocoupler(itap_sec,id_fvcom_send_u,oasis_nnn,1,field_fvcom_send_u)
        call tocoupler(itap_sec,id_fvcom_send_v,oasis_nnn,1,field_fvcom_send_v)
        call tocoupler(itap_sec,id_fvcom_send_zeta,oasis_nnn,1,field_fvcom_send_zeta)

     if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_atm2ocn-dt/lag",itap_sec_float,0,cpl_dt_atm2ocn,cpl_lag_atm2ocn
       !--meteorological forcing
       call fromcoupler(itap_sec,id_fvcom_recv_evap,oasis_nnn,1,field_fvcom_recv_evap)
       call fromcoupler(itap_sec,id_fvcom_recv_precip,oasis_nnn,1,field_fvcom_recv_precip)
       call fromcoupler(itap_sec,id_fvcom_recv_latent,oasis_nnn,1,field_fvcom_recv_latent)
       call fromcoupler(itap_sec,id_fvcom_recv_sensible,oasis_nnn,1,field_fvcom_recv_sensible)
       call fromcoupler(itap_sec,id_fvcom_recv_shortwave,oasis_nnn,1,field_fvcom_recv_shortwave)
       call fromcoupler(itap_sec,id_fvcom_recv_longwave,oasis_nnn,1,field_fvcom_recv_longwave)
       call fromcoupler(itap_sec,id_fvcom_recv_netheat,oasis_nnn,1,field_fvcom_recv_netheat)
       call fromcoupler(itap_sec,id_fvcom_recv_spq,oasis_nnn,1,field_fvcom_recv_spq)
       call fromcoupler(itap_sec,id_fvcom_recv_rhumidity,oasis_nnn,1,field_fvcom_recv_rhumidity)
       call fromcoupler(itap_sec,id_fvcom_recv_u10,oasis_nnn,1,field_fvcom_recv_u10)
       call fromcoupler(itap_sec,id_fvcom_recv_v10,oasis_nnn,1,field_fvcom_recv_v10)
       call fromcoupler(itap_sec,id_fvcom_recv_airt,oasis_nnn,1,field_fvcom_recv_airt)
       call fromcoupler(itap_sec,id_fvcom_recv_airp,oasis_nnn,1,field_fvcom_recv_airp)
       call fromcoupler(itap_sec,id_fvcom_recv_cloud,oasis_nnn,1,field_fvcom_recv_cloud)
       call fromcoupler(itap_sec,id_fvcom_recv_dew,oasis_nnn,1,field_fvcom_recv_dew)
     if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_wav2ocn-dt/lag",itap_sec_float,0,cpl_dt_wav2ocn,cpl_lag_wav2ocn
       !--wave radiations
       call fromcoupler(itap_sec,id_fvcom_recv_fx,oasis_nnn,1,field_fvcom_recv_fx)
       call fromcoupler(itap_sec,id_fvcom_recv_fy,oasis_nnn,1,field_fvcom_recv_fy)
       !--wave parameters
       call fromcoupler(itap_sec,id_fvcom_recv_hs,oasis_nnn2,1,field_fvcom_recv_hs)
       call fromcoupler(itap_sec,id_fvcom_recv_wlen,oasis_nnn2,1,field_fvcom_recv_wlen)
       call fromcoupler(itap_sec,id_fvcom_recv_dir,oasis_nnn2,1,field_fvcom_recv_dir)
       call fromcoupler(itap_sec,id_fvcom_recv_rtp,oasis_nnn2,1,field_fvcom_recv_rtp)
       call fromcoupler(itap_sec,id_fvcom_recv_tmbot,oasis_nnn2,1,field_fvcom_recv_tmbot)
       call fromcoupler(itap_sec,id_fvcom_recv_ubot,oasis_nnn2,1,field_fvcom_recv_ubot)

   END IF !(STARTUP_TYPE==STARTUP_TYPE_COLDSTART)

  END IF !end if USE_OASIS_COUPLER
# endif
!================ end ch201610 added for coupler=====================!

  CALL BCOND_GCN(8,0)
  
  CALL ARCHIVE
  ! ORDER MATTERS - ARCHIVE_NEST MUST GO AFTER ARCHIVE DURING SETUP!
  CALL ARCHIVE_NEST
# if defined(ENKF)
   enkf_out_Time=NC_DAT%FTIME%NEXT_IO
   DO ENKF_memberCONTR =1, ENKF_NENS
   write(enkf_num,'(i2.2)')  ENKF_memberCONTR
    IF(NC_ON)THEN
             call    ENKF_SETUP_NC(NC_DAT)
             IF (ENKF_memberCONTR == ENKF_NENS) ARCHIVE_TIMES=ARCHIVE_TIMES+1

    END IF
   END DO
# endif
# if defined (WAVE_CURRENT_INTERACTION)
  CALL ARCHIVE_NEST_WAVE
# endif

  CALL SET_PROBES(PROBES_ON,PROBES_NUMBER,PROBES_FILE)
  IF(OUT_STATION_TIMESERIES_ON)CALL READ_STATION_FILE
# if defined (WAVE_CURRENT_INTERACTION)
  IF(OUT_WAVE_SPARSE_TIMESERIES_ON)CALL SPARSE_STATION
# endif

  ! Setup Bounds checking (shutdown if variables exceed threshold)
  CALL SETUP_BOUNDSCHK !bounds checking

# if defined (LAG_PARTICLE)
  CALL SET_LAG
# endif

# if defined (GOTM)
  CALL INIT_GOTM
# endif
  
  IF(OUT_STATION_TIMESERIES_ON)THEN
    CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL),INTERVAL_TIME_SERIES)
    CALL OUT_STATION_TIMESERIES
    TIME_SERIES = STARTTIME + INTERVAL_TIME_SERIES
  END IF
# if defined (WAVE_CURRENT_INTERACTION)
  IF(OUT_WAVE_SPARSE_TIMESERIES_ON)THEN
    CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL_SPARSE),WAVE_INTERVAL_TIME_SERIES)
    CALL OUT_WAVE_SPARSE_TIMESERIES
    WAVE_TIME_SERIES = STARTTIME + WAVE_INTERVAL_TIME_SERIES
  END IF
# endif
  !==============================================================================!
  !  SELECT THE RUN MODE AND EXECUTE THE MAIN LOOP
  !==============================================================================!
  SELECT CASE(FVCOM_RUN_MODE)
     ! RUN MODE SET IN mod_assim.F(set_assim_param)

     ! =============================================================================!
     ! == PURE SIMULATION MODE - Instantanious data assimilation only ==============!
  CASE(FVCOM_PURE_SIM)
     ! =============================================================================!


     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
        write(IPT,*) "!============== STARTING MAIN LOOP AT:==============="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "!===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
     end if

     !////////////////////////// MAIN LOOP //////////////////////////////////////////
#  if defined (oasis_coupler)
   tmpwid=0 !-ch2016
#  endif
     DO IINT=ISTART,IEND

        IntTime=IntTime + IMDTI

!-------------- ch201610 added for coupler-----------------!
#  if defined (oasis_coupler)
  IF(USE_OASIS_COUPLER) THEN
   IF(STARTUP_TYPE==STARTUP_TYPE_COLDSTART) THEN
        itap_sec=seconds(IMDTI)*(IINT)
   ELSEIF(STARTUP_TYPE==STARTUP_TYPE_HOTSTART) THEN
        !itap_sec=seconds(IMDTI)*(IINT-1)
        itap_sec=seconds(IMDTI)*(IINT-ISTART) !ch202001 solve coupled restart issue
        itap_sec_float=real(seconds(IMDTI)*(IINT-ISTART)) !ch202012 solve the frication dt issue (dt<1)
   END IF
        itap_sec_save=itap_sec !
!-- REORDERED (was after the fromcoupler/RECV block below): same deadlock
!-- class as the initial t=0 exchange -- FVCOM's RECV-before-SEND ordering
!-- here means both models can block waiting on each other's data at the
!-- first *real* coupling boundary. sst_on_cell/node, aice_on_node, u, v,
!-- and el1 below are FVCOM's own prognostic state as of the END of the
!-- PREVIOUS timestep (INTERNAL_STEP, which updates them, has not run yet
!-- this iteration) -- they do not depend on anything in the RECV block
!-- that follows, so sending first is safe and introduces no new lag beyond
!-- what was already implicit in this block's position (before INTERNAL_STEP).
!----collect and send data to cpl
  IF (OASIS_ocn2atm) THEN
   IF (cpl_dt_ocn2atm>0) then
        sst_on_cell(1:N)=T(1:N,1)
        sst_on_node(1:M)=T1(1:M,1)
# if defined (ICE)
      if (ice_model) then
        aice_on_node(1:M)=AICE(1:M,1)
      end if
# endif
# if defined (MULTIPROCESSOR)
     IF (PAR) THEN
        CALL AEXCHANGE(EC,MYID,NPROCS,sst_on_cell)
        CALL AEXCHANGE(NC,MYID,NPROCS,sst_on_node)
        CALL AEXCHANGE(NC,MYID,NPROCS,aice_on_node)
     END IF
# endif
       CALL N2E2D(aice_on_node,aice_on_cell)

      !--calc avg sst/aice during period (itself if coupled on every step)
      tmpwid=tmpwid+1 !count # of each period
      sst_sum=sst_sum+sst_on_node
      aice_sum=aice_sum+aice_on_node
      !
      if(mod(itap_sec+cpl_lag_ocn2atm,cpl_dt_ocn2atm)==0) then
         if(MSR) print *,"taking average value of sst,aice during coupling period"

         sst_avg=sst_sum/cpl_dt_ocn2atm*seconds(IMDTI)
         sst_sum=0. !reset to 0 for next period

         aice_avg=aice_sum/cpl_dt_ocn2atm*seconds(IMDTI)
         aice_sum=0. !reset to 0 for next period

         if(MSR) print *,"# in each exchange period = ",tmpwid
         tmpwid=0 !reset count number
      end if

      !--check if sst/aice is set to be controlled (for ice)
      if(mod(itap_sec+cpl_lag_ocn2atm,cpl_dt_ocn2atm)==0) then
         count_sst=0
         count_aice=0
         !do iwd=1,N
         do iwd=1,M
          IF(OASIS_ICE_CTRL_SST) THEN
           if (sst_avg(iwd)<=OASIS_CTRL_VAL_SST) count_sst=count_sst+1
          END IF
          !
          IF(OASIS_ICE_CTRL_AICE) THEN
           if (aice_avg(iwd)>=OASIS_CTRL_VAL_AICE) count_aice=count_aice+1
          END IF
         end do
         !
         IF(OASIS_ICE_CTRL_SST.or.OASIS_ICE_CTRL_AICE) THEN
          if (count_sst>0.or.count_aice>0) &
            print *,"aice_control#,sst_control#,myid",count_aice,count_sst,myid
         END IF
      end if

      do iwd=1,oasis_nnn2
         field_fvcom_send_sst(iwd,1)=sst_avg(local_nn(iwd))+273.15 ! C to K
         field_fvcom_send_aice(iwd,1)=aice_avg(local_nn(iwd)) !
      end do

      !--ch201611 added aice/sst control--!
      IF(OASIS_ICE_CTRL_AICE) THEN
       if (MSR.and.mod(itap_sec+cpl_lag_ocn2atm,cpl_dt_ocn2atm)==0) &
            print *,"ice-control(aice) in fvcom,oasis_ctrl_val_aice",&
                        OASIS_CTRL_VAL_AICE,OASIS_freezing_temp
       do iwd=1,oasis_nnn2
         if (aice_avg(local_nn(iwd))>=OASIS_CTRL_VAL_AICE) then
            field_fvcom_send_sst(iwd,1)=OASIS_freezing_temp+273.15
         end if
       end do
      ELSE
       if (MSR.and.mod(itap_sec+cpl_lag_ocn2atm,cpl_dt_ocn2atm)==0) &
            print *,"*NO* ice-control(aice) in fvcom"
      END IF
      !
      IF(OASIS_ICE_CTRL_SST) THEN
       if (MSR.and.mod(itap_sec+cpl_lag_ocn2atm,cpl_dt_ocn2atm)==0) &
            print *,"ice-control(sst) in fvcom,oasis_ctrl_val_sst",&
                        OASIS_CTRL_VAL_SST,OASIS_freezing_temp
       do iwd=1,oasis_nnn2
         if (sst_avg(local_nn(iwd))<=OASIS_CTRL_VAL_SST) then
            field_fvcom_send_sst(iwd,1)=OASIS_freezing_temp+273.15
         end if
       end do
      ELSE
       if (MSR.and.mod(itap_sec+cpl_lag_ocn2atm,cpl_dt_ocn2atm)==0) &
            print *,"*NO* ice-control(sst) in fvcom"
      END IF

      !---
      if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_ocn2atm-dt/lag",&
                                itap_sec_float,IINT,cpl_dt_ocn2atm,cpl_lag_ocn2atm
      !if (mod(itap_sec_float,1.0)==0) then !-to aviod fractional dt (i.e., dt<1) issue
        call tocoupler(itap_sec,id_fvcom_send_sst,oasis_nnn2,1,field_fvcom_send_sst)
        call tocoupler(itap_sec,id_fvcom_send_aice,oasis_nnn2,1,field_fvcom_send_aice)
      !end if
   END IF !(cpl_dt_ocn2atm>0)
  END IF !(OASIS_ocn2atm)

  IF (OASIS_ocn2wav) THEN !-send zeta
   IF (cpl_dt_ocn2wav>0) then
# if defined (MULTIPROCESSOR)
     IF (PAR) THEN
        CALL AEXCHANGE(EC,MYID,NPROCS,el1)
     END IF
# endif
        field_fvcom_send_zeta(:,1)=el1(:)

      if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_ocn2wav dt,lag",&
                                itap_sec_float,IINT,cpl_dt_ocn2wav,cpl_lag_ocn2wav
      !if (mod(itap_sec_float,1.0)==0) then !-to aviod fractional dt (i.e., dt<1) issue
        call tocoupler(itap_sec,id_fvcom_send_zeta,oasis_nnn,1,field_fvcom_send_zeta)
      !end if
   END IF !(cpl_dt_ocn2wav>0)
  END IF !(OASIS_ocn2wav)

  IF (OASIS_ocn2wav.or.OASIS_ocn2atm) THEN !-send u,v
   IF (cpl_dt_ocn2wav>0.or.cpl_dt_ocn2atm>0) then
# if defined (MULTIPROCESSOR)
     IF (PAR) THEN
        CALL AEXCHANGE(EC,MYID,NPROCS,u,v)
     END IF
# endif
        field_fvcom_send_u(:,1)=u(:,1)
        field_fvcom_send_v(:,1)=v(:,1)

      if(MSR.and.cpl_dt_ocn2wav>0) &
              print cpl_outfmt,"itap_sec,IINT,cpl_ocn2wav dt,lag",&
                                itap_sec_float,IINT,cpl_dt_ocn2wav,cpl_lag_ocn2wav
      if(MSR.and.cpl_dt_ocn2atm>0) &
              print cpl_outfmt,"itap_sec_float,IINT,cpl_ocn2atm dt,lag",&
                                itap_sec_float,IINT,cpl_dt_ocn2atm,cpl_lag_ocn2atm
      !if (mod(itap_sec_float,1.0)==0) then !-to aviod fractional dt (i.e., dt<1) issue
        call tocoupler(itap_sec,id_fvcom_send_u,oasis_nnn,1,field_fvcom_send_u)
        call tocoupler(itap_sec,id_fvcom_send_v,oasis_nnn,1,field_fvcom_send_v)
      !end if
   END IF !(ocn2wav or ocn2atm)
  END IF !(OASIS_ocn2wav/ocn2atm)

     if(cpl_dt_atm2ocn>0) then
         !if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_atm2ocn-dt/lag",&
         !                  itap_sec,IINT,cpl_dt_atm2ocn,cpl_lag_atm2ocn
         if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_atm2ocn-dt/lag",&
                           itap_sec_float,IINT,cpl_dt_atm2ocn,cpl_lag_atm2ocn
      if (IINT<IEND) then !-ch201808
      !if (mod(itap_sec_float,1.0)==0) then !-to aviod fractional dt (i.e., dt<1) issue
       call fromcoupler(itap_sec,id_fvcom_recv_evap,oasis_nnn,1,field_fvcom_recv_evap)
       call fromcoupler(itap_sec,id_fvcom_recv_precip,oasis_nnn,1,field_fvcom_recv_precip)
       call fromcoupler(itap_sec,id_fvcom_recv_latent,oasis_nnn,1,field_fvcom_recv_latent)
       call fromcoupler(itap_sec,id_fvcom_recv_sensible,oasis_nnn,1,field_fvcom_recv_sensible)
       call fromcoupler(itap_sec,id_fvcom_recv_shortwave,oasis_nnn,1,field_fvcom_recv_shortwave)
       call fromcoupler(itap_sec,id_fvcom_recv_longwave,oasis_nnn,1,field_fvcom_recv_longwave)
       call fromcoupler(itap_sec,id_fvcom_recv_netheat,oasis_nnn,1,field_fvcom_recv_netheat)
       call fromcoupler(itap_sec,id_fvcom_recv_spq,oasis_nnn,1,field_fvcom_recv_spq)
       call fromcoupler(itap_sec,id_fvcom_recv_rhumidity,oasis_nnn,1,field_fvcom_recv_rhumidity)
       call fromcoupler(itap_sec,id_fvcom_recv_u10,oasis_nnn,1,field_fvcom_recv_u10)
       call fromcoupler(itap_sec,id_fvcom_recv_v10,oasis_nnn,1,field_fvcom_recv_v10)
       call fromcoupler(itap_sec,id_fvcom_recv_airt,oasis_nnn,1,field_fvcom_recv_airt)
       call fromcoupler(itap_sec,id_fvcom_recv_airp,oasis_nnn,1,field_fvcom_recv_airp)
       call fromcoupler(itap_sec,id_fvcom_recv_cloud,oasis_nnn,1,field_fvcom_recv_cloud)
       call fromcoupler(itap_sec,id_fvcom_recv_dew,oasis_nnn,1,field_fvcom_recv_dew)
      !end if
      end if
     end if
     !
     if(cpl_dt_wav2ocn>0) then
         if(MSR) print cpl_outfmt,"itap_sec,IINT,cpl_wav2ocn-dt/lag",&
                           itap_sec_float,IINT,cpl_dt_wav2ocn,cpl_lag_wav2ocn
      if (IINT<IEND) then !-ch201808
      !if (mod(itap_sec_float,1.0)==0) then !-to aviod fractional dt (i.e., dt<1) issue
         if (MSR) print *, "check integer itap_sec for fromcpl subroutine"
       !--wave radiation stress gradient
       call fromcoupler(itap_sec,id_fvcom_recv_fx,oasis_nnn,1,field_fvcom_recv_fx)
       call fromcoupler(itap_sec,id_fvcom_recv_fy,oasis_nnn,1,field_fvcom_recv_fy)
       !--wave parameters
       call fromcoupler(itap_sec,id_fvcom_recv_hs,oasis_nnn2,1,field_fvcom_recv_hs)
       call fromcoupler(itap_sec,id_fvcom_recv_wlen,oasis_nnn2,1,field_fvcom_recv_wlen)
       call fromcoupler(itap_sec,id_fvcom_recv_dir,oasis_nnn2,1,field_fvcom_recv_dir)
       call fromcoupler(itap_sec,id_fvcom_recv_rtp,oasis_nnn2,1,field_fvcom_recv_rtp)
       call fromcoupler(itap_sec,id_fvcom_recv_tmbot,oasis_nnn2,1,field_fvcom_recv_tmbot)
       call fromcoupler(itap_sec,id_fvcom_recv_ubot,oasis_nnn2,1,field_fvcom_recv_ubot)
      !end if
      end if
     end if
     !

       !---check value (units)--!
        do iwd=1,N
           if(field_fvcom_recv_airt(iwd,1)>100.0) then
              print *,"Tair>100 degC",field_fvcom_recv_airt(iwd,1)
              field_fvcom_recv_airt(iwd,1)=field_fvcom_recv_airt(iwd,1)-273.15
           end if

           if(field_fvcom_recv_airp(iwd,1)<2000.0) then
              print *,"Ps<2000 Pa",field_fvcom_recv_airp(iwd,1)
              field_fvcom_recv_airp(iwd,1)=field_fvcom_recv_airp(iwd,1)*100.0
           end if
        end do !iwd=1,N
       !---end check---------------!

  IF (OASIS_atm2ocn) THEN !-ch2018 do only if active atm2ocn
   IF (cpl_dt_atm2ocn<0) then
      CALL FATAL_ERROR("OASIS_atm2ocn=T,but no cpl field,please check namelist")
   ELSE
   !!---- save received data into cell-based variables
      do iwd=1,N
# if defined (ICE) || (HEATING_CALCULATED) || (HEATING_CALCULATED_GL) || (HEATING_SOLAR)
!--this is new for fvcom4.1 due to heating_solar option
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then !save IO_previous every dt
               T_AIR_cell0(iwd)=T_AIR_cell2(iwd)
            end if
               T_AIR_cell2(iwd)=field_fvcom_recv_airt(iwd,1)

            !- make IO_next and IO_previous equal for first 2 coupling steps
            if(itap_sec<=cpl_dt_atm2ocn) then
               T_AIR_cell0(iwd)=T_AIR_cell2(iwd)
            end if
# endif
# if defined (ICE) || (HEATING_CALCULATED) || (HEATING_CALCULATED_GL)
      !!-- ch201611
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then
               !T_AIR_cell0(iwd)=T_AIR_cell2(iwd) !-moved to single section above
               RH_AIR_cell0(iwd)=RH_AIR_cell2(iwd)
               DLW_AIR_cell0(iwd)=DLW_AIR_cell2(iwd)
               DSW_AIR_cell0(iwd)=DSW_AIR_cell2(iwd)
               HEAT_NET_cell0(iwd)=HEAT_NET_cell2(iwd)
            end if
               !T_AIR_cell2(iwd)=field_fvcom_recv_airt(iwd,1)
               RH_AIR_cell2(iwd)=field_fvcom_recv_rhumidity(iwd,1)
               DLW_AIR_cell2(iwd)=field_fvcom_recv_longwave(iwd,1)
               DSW_AIR_cell2(iwd)=field_fvcom_recv_shortwave(iwd,1)
               HEAT_NET_cell2(iwd)=field_fvcom_recv_netheat(iwd,1)

            if(itap_sec<=cpl_dt_atm2ocn) then
               !T_AIR_cell0(iwd)=T_AIR_cell2(iwd)
               RH_AIR_cell0(iwd)=RH_AIR_cell2(iwd)
               DLW_AIR_cell0(iwd)=DLW_AIR_cell2(iwd)
               DSW_AIR_cell0(iwd)=DSW_AIR_cell2(iwd)
               HEAT_NET_cell0(iwd)=HEAT_NET_cell2(iwd)
            end if
# endif
# if defined (HEATING_SOLAR) || defined (ICE)
!-ch2018
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then
               CLOUD_cell0(iwd)=CLOUD_cell2(iwd)
            end if
               CLOUD_cell2(iwd)=field_fvcom_recv_cloud(iwd,1)

            if(itap_sec<=cpl_dt_atm2ocn) then
               CLOUD_cell0(iwd)=CLOUD_cell2(iwd)
            end if
# endif
# if defined (HEATING_SOLAR)
!-ch2019
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then
               TDEW_cell0(iwd)=TDEW_cell2(iwd) !-dew point temp
            end if
               TDEW_cell2(iwd)=field_fvcom_recv_dew(iwd,1)

            if(itap_sec<=cpl_dt_atm2ocn) then
               TDEW_cell0(iwd)=TDEW_cell2(iwd)
            end if
# endif
# if defined (ICE)
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then
               QA_AIR_cell0(iwd)=QA_AIR_cell2(iwd)
            end if
               QA_AIR_cell2(iwd)=field_fvcom_recv_spq(iwd,1)

            if(itap_sec<=cpl_dt_atm2ocn) then
               QA_AIR_cell0(iwd)=QA_AIR_cell2(iwd)
            end if
# endif
            !-u10,v10,evap,precip
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then
               UUWIND_oasis0(iwd)=UUWIND_oasis2(iwd)
               VVWIND_oasis0(iwd)=VVWIND_oasis2(iwd)
               if (PRECIPITATION_ON) then
                  QPREC_cell0(iwd)=QPREC_cell2(iwd)
                  QEVAP_cell0(iwd)=QEVAP_cell2(iwd)
                  QPREC2_cell0(iwd)=QPREC2_cell2(iwd)
                  QEVAP2_cell0(iwd)=QEVAP2_cell2(iwd)
               end if
            end if
               UUWIND_oasis2(iwd)=field_fvcom_recv_u10(iwd,1)
               VVWIND_oasis2(iwd)=field_fvcom_recv_v10(iwd,1)
               if (PRECIPITATION_ON) then
                  QPREC_cell2(iwd)=field_fvcom_recv_precip(iwd,1)
                  QEVAP_cell2(iwd)=-field_fvcom_recv_evap(iwd,1) !-evap is negative in fvcom
                  QPREC2_cell2(iwd)=field_fvcom_recv_precip(iwd,1)
                  QEVAP2_cell2(iwd)=-field_fvcom_recv_evap(iwd,1) !-evap is negative in fvcom
               end if

            if(itap_sec<=cpl_dt_atm2ocn) then
               UUWIND_oasis0(iwd)=UUWIND_oasis2(iwd)
               VVWIND_oasis0(iwd)=VVWIND_oasis2(iwd)
               if (PRECIPITATION_ON) then
                  QPREC_cell0(iwd)=QPREC_cell2(iwd)
                  QEVAP_cell0(iwd)=QEVAP_cell2(iwd)
                  QPREC2_cell0(iwd)=QPREC2_cell2(iwd)
                  QEVAP2_cell0(iwd)=QEVAP2_cell2(iwd)
               end if
            end if

# if defined (AIR_PRESSURE) || (HEATING_CALCULATED) || (HEATING_CALCULATED_GL)
            if(mod(itap_sec,cpl_dt_atm2ocn)==0) then
               PA_AIR_cell0(iwd)=PA_AIR_cell2(iwd)
            end if
               PA_AIR_cell2(iwd)=field_fvcom_recv_airp(iwd,1)

            if(itap_sec<=cpl_dt_atm2ocn) then
               PA_AIR_cell0(iwd)=PA_AIR_cell2(iwd)
            end if
# endif
      end do
      !--------- end save received data into cell-based variables

   !-------ch201611 added for partition boundary issue---------!
# if defined (MULTIPROCESSOR)
     IF (PAR) THEN
        CALL AEXCHANGE(EC,MYID,NPROCS,UUWIND_oasis2,VVWIND_oasis2)
        CALL AEXCHANGE(EC,MYID,NPROCS,UUWIND_oasis0,VVWIND_oasis0)
       if (PRECIPITATION_ON) then
        CALL AEXCHANGE(EC,MYID,NPROCS,QPREC_cell2,QEVAP_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,QPREC2_cell2,QEVAP2_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,QPREC_cell0,QEVAP_cell0)
        CALL AEXCHANGE(EC,MYID,NPROCS,QPREC2_cell0,QEVAP2_cell0)
       end if
!
# if defined (ICE) || (HEATING_CALCULATED) || (HEATING_CALCULATED_GL) || (HEATING_SOLAR)
        CALL AEXCHANGE(EC,MYID,NPROCS,T_AIR_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,T_AIR_cell0)
# endif
!
# if defined (ICE) || (HEATING_CALCULATED) || (HEATING_CALCULATED_GL)
        !CALL AEXCHANGE(EC,MYID,NPROCS,T_AIR_cell2,RH_AIR_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,RH_AIR_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,DLW_AIR_cell2,DSW_AIR_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,HEAT_NET_cell2)
!
        !CALL AEXCHANGE(EC,MYID,NPROCS,T_AIR_cell0,RH_AIR_cell0)
        CALL AEXCHANGE(EC,MYID,NPROCS,RH_AIR_cell0)
        CALL AEXCHANGE(EC,MYID,NPROCS,DLW_AIR_cell0,DSW_AIR_cell0)
        CALL AEXCHANGE(EC,MYID,NPROCS,HEAT_NET_cell0)
# endif
!
# if defined (HEATING_SOLAR) || defined (ICE)
!-ch2018
        CALL AEXCHANGE(EC,MYID,NPROCS,CLOUD_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,CLOUD_cell0)
# endif
!
# if defined (HEATING_SOLAR)
!-ch2019
        CALL AEXCHANGE(EC,MYID,NPROCS,TDEW_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,TDEW_cell0)
# endif
!
# if defined (ICE)
        CALL AEXCHANGE(EC,MYID,NPROCS,QA_AIR_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,QA_AIR_cell0)
# endif
!
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED) || (HEATING_CALCULATED_GL)
        CALL AEXCHANGE(EC,MYID,NPROCS,PA_AIR_cell2)
        CALL AEXCHANGE(EC,MYID,NPROCS,PA_AIR_cell0)
# endif
     END IF
# endif
   !-------end ch201611 added for partition boundary issue---------!
   END IF !--cpl_dt_atm2ocn<0
  END IF !(OASIS_atm2ocn)

  END IF !end if USE_OASIS_COUPLER
# endif
!================ end ch201610 added for coupler=====================!

        CALL INTERNAL_STEP

        !==============================================================================!
        !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
        !==============================================================================!
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        
        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

        !==============================================================================!
        !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
        !==============================================================================!

!for Eulerian velocity output
# if defined (WAVE_CURRENT_INTERACTION)
         !U = U-U_STOKES_3D
         !V = V-V_STOKES_3D
         !UA=UA-U_STOKES_2D
         !VA=VA-V_STOKES_2D
# endif

        CALL ARCHIVE

#  if defined (WAVE_CURRENT_INTERACTION)
!Convert Eulerian velocity back to Lagrangian velocity to maintain consistent with Mellor's equation
         !U = U+U_STOKES_3D
         !V = V+V_STOKES_3D
         !UA=UA+U_STOKES_2D
         !VA=VA+V_STOKES_2D
#  endif

        CALL DUMP_PROBE_DATA
        IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
# if defined (WAVE_CURRENT_INTERACTION)
        IF(OUT_WAVE_SPARSE_TIMESERIES_ON) CALL OUT_WAVE_SPARSE_TIMESERIES
# endif
        !==============================================================================!
        !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
        !==============================================================================!
        CALL SHUTDOWN_CHECK(D1)

        !==============================================================================!
        !  CALL BOUNDS CHECK TO SEE IF VARIABLES EXCEED USER-DEFINED THRESHOLDS 
        !==============================================================================!
        CALL BOUNDSCHK  !bounds checking

        !==============================================================================!
        !    LAGRANGIAN PARTICLE TRACKING                                              |
        !==============================================================================!
# if defined (LAG_PARTICLE)
        CALL LAG_UPDATE
# endif

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
        IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
        IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif

     END DO
     !////////////////////////// END MAIN LOOP //////////////////////////////////////

  ! ================================================================================!
  ! == THIS MAIN LOOP USES NUDGING OR OI METHODS TO ASSIMILATE =====================!
  CASE(FVCOM_NUDGE_OI_ASSIM)
  ! ================================================================================!
# if defined (DATA_ASSIM) 

     CALL ALLOC_BUFFER

     ! SET THE ASSIMILATION/SIMULATION RESET TIME
     ASSIM_RESTART_TIME = StartTime
     IF(TSGRD_ASSIM) THEN
       CALL Now_2_month_days(StartTime,Pyear,Pmonth,Pmdays)
       TSC_SAVE_TIME2  = IntTime + days2time(real(Pmdays,kind=DP))
     ENDIF

     IF(SSHGRD_ASSIM) SSH_SAVED(1:M,0) = EL(1:M)

     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
     INT_START = ISTART
     INT_COUNT = ISTART
     IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
       INT_END = ISTART + 2*(IEND-ISTART)
     ELSEIF(TS_NGASSIM .OR. CUR_NGASSIM .OR.      &
            TS_OIASSIM .OR. CUR_OIASSIM)THEN
       INT_END = IEND
     END IF

     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
        Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     IF(TSGRD_ASSIM) THEN
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "SETUP ASSIMILATION RSTFILE"
       CALL SETUP_RSTFILE_ASSIM
     ENDIF

     DO WHILE(IntTime < EndTime)
       IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
         ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY
       ELSEIF(TS_NGASSIM.OR.CUR_NGASSIM .OR.       &
              TS_OIASSIM.OR.CUR_OIASSIM)THEN
         ASSIM_RESTART_TIME = EndTime
       END IF

        IF(SST_ASSIM.OR.SSTGRD_ASSIM) THEN
          SST_SAVE_INDEX = 0
          SST_SAVE_TIME  = IntTime + sst_save_interval
        ENDIF

        IF(TSGRD_ASSIM) THEN
          TSC_SAVE_INDEX = 0
          TSC_SAVE_TIME = IntTime + TSC_SAVE_INTERVAL 

          CALL Now_2_month_days(IntTime+IMDTI,Pyear,Pmonth,Pmdays)

          IF(ASSIM_RESTART_TIME > TSC_SAVE_TIME2) THEN
            TSC_SAVE_INDEX2 = 1
            TSC_SAVE_TIME2 = TSC_SAVE_TIME2 + days2time(real(Pmdays,kind=DP))
          ELSE
            TSC_SAVE_INDEX2 = 0
          ENDIF
        ENDIF

        IF(SSHGRD_ASSIM) THEN
          SSH_SAVE_INDEX = 0
          SSH_SAVE_TIME  = IntTime + ssh_save_interval
        ENDIF

        CALL SAVE_STATE

        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Simulation  ===================="
        end if

        !==============================================================================!
        !    RUN PURE SIMULATION MODE:
        !==============================================================================!

        ASSIM_FLAG = 0
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
             ! do nothing
           ELSE
             CALL ARCHIVE

             CALL DUMP_PROBE_DATA
             IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
#            if defined (WAVE_CURRENT_INTERACTION)
             IF(OUT_WAVE_SPARSE_TIMESERIES_ON) CALL OUT_WAVE_SPARSE_TIMESERIES
#            endif


#            if defined (LAG_PARTICLE)
             CALL LAG_UPDATE
#            endif

             IF(NCNEST_ON)      CALL ARCHIVE_NEST     
#            if defined (WAVE_CURRENT_INTERACTION)
             IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#            endif
           END IF

           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

        END DO

        !==============================================================================!
        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
        !==============================================================================!
        IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
         if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
         CALL RESTORE_STATE
         CALL RHO_PMEAN

         IF(SST_ASSIM.OR.SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE
         IF(TSGRD_ASSIM) CALL TSGRD_OBSERVATION_UPDATE(IntTime+(days2time(1.0_DP))/2)
         IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

         if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="

         if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Assimilation  ===================="
         end if

        ELSEIF(TS_NGASSIM .OR. CUR_NGASSIM .OR. TS_OIASSIM .OR. CUR_OIASSIM)THEN
          IntTime = ASSIM_RESTART_TIME
        END IF
        !==============================================================================!
        !    RUN DATA ASSIMILATION MODE:
        !==============================================================================!
        ASSIM_FLAG = 1
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
           !==============================================================================!
           CALL ARCHIVE

           CALL DUMP_PROBE_DATA
           IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
#          if defined (WAVE_CURRENT_INTERACTION)
           IF(OUT_WAVE_SPARSE_TIMESERIES_ON)CALL OUT_WAVE_SPARSE_TIMESERIES
#          endif
           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

           !==============================================================================!
           !    LAGRANGIAN PARTICLE TRACKING                                              |
           !==============================================================================!
# if defined (LAG_PARTICLE)
           CALL LAG_UPDATE
# endif

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
        IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
        IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif

        END DO

        ! RUN THE NEXT DAY

        IF(SSHGRD_ASSIM) SSH_SAVED(1:M,0) = SSH_SAVED(1:M,SSHGRD_N_PER_INTERVAL)

     END DO
# endif

!  ! ================================================================================!
!  ! == THIS MAIN LOOP USES OPTIMAL INTERPOLATION  METHODS TO ASSIMILATE ============!
!  CASE(FVCOM_OI_ASSIM)
!  ! ================================================================================!
!
!     CALL SET_OIASSIM_INTERVALS
!     CALL ALLOC_BUFFER_OI
!
!     ! SET THE ASSIMILATION/SIMULATION RESET TIME
!     ASSIM_RESTART_TIME = StartTime
!
!     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
!     INT_START = ISTART
!     INT_COUNT = ISTART
!     IF(SST_OIASSIM.OR.SSTGRD_OIASSIM)THEN
!       INT_END = ISTART + 2*(IEND-ISTART)
!     ELSEIF(TS_OIASSIM.OR.CUR_OIASSIM)THEN
!       INT_END = IEND
!     END IF
!     !==============================================================================!
!     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
!     !==============================================================================!
!     if(DBG_SET(dbg_log)) THEN
!        write(IPT,*) "===================================================="
!        write(IPT,*) "===================================================="
!        write(IPT,*) "======STARTING MAIN LOOP OIASSIMILATION MODE:======="
!        Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!        write(IPT,*) "===================================================="
!     end if
!
!     CALL REPORT('INITIAL CONDITIONS')
!     
!     if(DBG_SET(dbg_log)) THEN
!        write(IPT,*) "===================================================="
!        write(IPT,*) "===================================================="
!        write(IPT,*) "===================================================="
!     end if

!     DO WHILE(IntTime < EndTime)
!        IF(SST_OIASSIM.OR.SSTGRD_OIASSIM)THEN
!          ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY
!        ELSEIF(TS_OIASSIM.OR.CUR_OIASSIM)THEN
!          ASSIM_RESTART_TIME = EndTime
!        END IF
!        SST_SAVE_INDEX = 0
!        SST_SAVE_TIME  = IntTime + sst_save_interval
!
!        CALL OI_SAVE_STATE
!
!
!        if(DBG_SET(dbg_log)) THEN
!           Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!           write(IPT,*) "======= Start 1 Day Simulation  ===================="
!        end if
!        !==============================================================================!
!        !    RUN PURE SIMULATION MODE:
!        !==============================================================================!
!
!        OIASSIM_FLAG = 0
!        DO WHILE(IntTime < ASSIM_RESTART_TIME)
!
!           IINT = IINT + 1
!           IntTime=IntTime + IMDTI
!
!           CALL INTERNAL_STEP
!
!           !==============================================================================!
!           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!           !==============================================================================!
!           if(DBG_SET(dbg_log)) &
!                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!
!           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')
!
!           !==============================================================================!
!           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
!           !==============================================================================!
!           IF(TS_OIASSIM.OR.CUR_OIASSIM)THEN
!
!             CALL ARCHIVE
!
!             CALL DUMP_PROBE_DATA
!
!           END IF
!
!           CALL SHUTDOWN_CHECK(D1)
!
!        END DO
!
!        IF(SST_OIASSIM.OR.SSTGRD_OIASSIM)THEN
!        !==============================================================================!
!        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
!        !==============================================================================!
!          CALL OI_RESTORE_STATE
!
!          if(DBG_SET(dbg_log)) THEN
!             Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
!             write(IPT,*) "======= Start 1 Day Assimilation  ===================="
!          end if
!
!        !==============================================================================!
!        !    RUN DATA ASSIMILATION MODE:
!        !==============================================================================!
!          OIASSIM_FLAG = 1
!          DO WHILE(IntTime < ASSIM_RESTART_TIME)
!
!             IINT = IINT + 1
!             IntTime=IntTime + IMDTI
!
!             CALL INTERNAL_STEP
!
!           !==============================================================================!
!           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!           !==============================================================================!
!             if(DBG_SET(dbg_log)) &
!                  & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!
!             IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')
!
!           !==============================================================================!
!           !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
!           !==============================================================================!
!             CALL ARCHIVE
!
!             CALL DUMP_PROBE_DATA
!           !==============================================================================!
!           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
!           !==============================================================================!
!             CALL SHUTDOWN_CHECK(D1)
!
!
!           !==============================================================================!
!           !    LAGRANGIAN PARTICLE TRACKING                                              |
!           !==============================================================================!
!# if defined (LAG_PARTICLE)
!             CALL LAG_UPDATE
!# endif
!
!          END DO
!        END IF
!        ! RUN THE NEXT DAY
!
!     END DO
    
  CASE(FVCOM_RRKF_WITHOUT_SSA)
# if defined (RRKF)     

    IF(MSR) WRITE(IPT,*) 
    IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'

!    IF(MSR) CALL RRK_REF
     CALL RRK_REF
     print *, 'finish rrk_ref'
     IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
    print *, 'call setup_rrk_restart'   
!    IF(LOCAL_DISK)THEN
!       NCF_RRKRE => SETUP_RESTART(LOCAL_FILE)
!    ELSE
       NCF_RRKRE => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE))
!        NCF_RRKRE=>new_file()
!        NCF_RRKRE%FNAME=TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE)
!        call NC_LOAD(NCF_RRKRE)
!        NCF_rrkre%FTIME=>NEW_FTIME()
!        print *, 'finish nc_load'
!    END IF

!    IF(LOCAL_DISK)THEN
!       NCF_RRKFCT => SETUP_RESTART(LOCAL_FILE1)
!    ELSE
       NCF_RRKFCT => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE1))
!    END IF
    print *, 'call rrk_rrk' 
    CALL RRK_RRK(1)    !  perturbation
!        NCF_RRKRE=>new_file()
!        NCF_RRKRE%FNAME=TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE)
!        call NC_LOAD(NCF_RRKRE)
!        NCF_rrkre%FTIME=>NEW_FTIME()
        print *, 'finish nc_load'
    print *, 'finish rrk_rrki'
    print *, 'before nc open'
    CALL NC_OPEN(NCF_RRKRE)
!    CALL NC_LOAD(NCF_RRKRE)
    print *, 'finish nc_open'
    DO RRK_EOFCONTR =1, RRK_NEOF
           print *, 'before readrestart' 
              IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
              print *, NCF_RRKRE%FNAME
         CALL READRESTART(NCF_RRKRE,RRK_EOFCONTR)
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
           print *, 'finish readrestart'
         IntTime = REF_START_TIME
         ExtTime = REF_START_TIME

         iint = rkint_start         
  
!==============================================================================!
!  PREPARE TO START FVCOM'S MAIN LOOP                                          !
!==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "============   ", 'RRK_NEOF =', RRK_EOFCONTR,"=============="
        write(IPT,*) "===================================================="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     

     !////////////////////////// MAIN LOOP //////////////////////////////////////////
     DO WHILE(IntTime < (REF_START_TIME + RRK_INTERVAL) )

        IntTime=IntTime + IMDTI
        iint = iint +1
        
        CALL INTERNAL_STEP

        !==============================================================================!
        !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
        !==============================================================================!
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,rkInt_start,rkInt_END,IntTime)
        
!        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS: Perturbation Mode')

        !==============================================================================!
        !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
        !==============================================================================!
        CALL SHUTDOWN_CHECK(D1)

     END DO
     !////////////////////////// END MAIN LOOP /////////////////////
     !/////////////////


     CALL WRITERESTART(NCF_RRKFCT,RRK_EOFCONTR)  
     
   END DO  ! RRK_NEOF
    print *, 'before rrk_rrk(2)'  
   CALL RRK_RRK(2)   ! calculate linearized model matrix
    print *, 'before rrk_rrk(4)'
   
   IF(MSR) CALL RRK_RRK(4) ! doubling algorithm

   
!JQI4000 CONTINUE   

   IF(MSR) WRITE(IPT,*)
   IF(MSR) WRITE(IPT,*) 'Finish Kalman filter preparation! '
   IF(MSR) WRITE(IPT,*) 'Start data assimilation! '
   IF(MSR) WRITE(IPT,*)
   
   CALL NC_OPEN(NC_START)  

   
   CALL READRESTART(NC_START,REF_START_TIME)
   
   
   IntTime = RRK_START_TIME
   ExtTime = RRK_START_TIME

   iint = rkint_start 
   

!==============================================================================!
!  PREPARE TO START FVCOM'S MAIN LOOP                                          !
!==============================================================================!
   if(DBG_SET(dbg_log)) THEN
      write(IPT,*) "===================================================="
      write(IPT,*) "===================================================="
      write(IPT,*) "============== STARTING MAIN LOOP AT:==============="
      if(DBG_SET(dbg_log)) &
           & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
      write(IPT,*) "===================================================="
   end if
   
   CALL REPORT('INITIAL CONDITIONS')
   
   if(DBG_SET(dbg_log)) THEN
      write(IPT,*) "===================================================="
      write(IPT,*) "===================================================="
      write(IPT,*) "===================================================="
   end if
   
   ICYC = 1
   RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)

   !////////////////////////// MAIN LOOP //////////////////////////////////////////
   DO WHILE(IntTime < ENDTIME)
      
      IntTime=IntTime + IMDTI
      iint = iint +1
      
      CALL INTERNAL_STEP    
           
      !==============================================================================!
      !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
      !==============================================================================!
      if(DBG_SET(dbg_log)) &
           & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
      
!      IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS: Assimilation Mode')
      
      !==============================================================================!
      !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
      !==============================================================================!
!      CALL ARCHIVE
      
!      CALL DUMP_PROBE_DATA
      !==============================================================================!
      !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
      !==============================================================================!
!      CALL SHUTDOWN_CHECK(D1)
      
      !==============================================================================!
      !    LAGRANGIAN PARTICLE TRACKING                                              |
      !==============================================================================!
!      CALL LAG_UPDATE
      
      !==============================================================================!
      !    NESTING OUTPUT                                                            |
      !==============================================================================!
!      IF(NESTING_ON .and. NESTING_MODE == "large_domain") CALL ARCHIVE_NEST

      !==============================================================================!
      !    KALMAN FILTERING                                                          |
      !==============================================================================!
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
                if (msr) print *, 'finish the filter_obs_data'
      IF(IntTime>=RRK_START_TIME .AND. IntTime<=RRK_END_TIME) THEN
         IF(IntTime == RRK_CYC) THEN
            CALL PRINT_TIME(IntTime,IPT,"DOING RRK_RRKF")
#if defined (MULTIPROCESSOR)
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
#endif 
                CALL filter_obs_data
#if defined (MULTIPROCESSOR)
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
#endif
            CALL set_rrkf_assim_data
               print *, 'kf observation location is:', nloc
               if (nloc>0) then
            print *, 'rrk_rrkf'
            CALL RRK_RRKF
               end if
           call deallocate_obs_data
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
           print *, 'finish rrkf'
	    ICYC = ICYC+1
            RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)	    
         END IF
      END IF
            
      
   END DO
   !////////////////////////// END MAIN LOOP ////////////////////////////////////// 


! #####################################################################################
  
# endif
  CASE(FVCOM_RRKF_WITH_SSA)
# if defined (RRKF)
# if defined (DATA_ASSIM) 
   IF (SSTGRD_ASSIM .or. SSHGRD_ASSIM) THEN
   ELSE
   PRINT *, 'STOP, RRKF ONLY WITH SSTGRD_ASSIM OR SSHGRD_ASSIM'
   CALL PSTOP
   END IF
  IF(MSR) WRITE(IPT,*) 'Starting RRKF input data processing......'
     CALL RRK_REF
     print *, 'finish rrk_ref'
     IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
    print *, 'call setup_rrk_restart'   

       NCF_RRKRE => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE))

       NCF_RRKFCT => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE1))
    print *, 'call rrk_rrk' 
    CALL RRK_RRK(1)    !  perturbation


    CALL NC_OPEN(NCF_RRKRE)
    print *, 'finish nc_open'
    DO RRK_EOFCONTR =1, RRK_NEOF
           print *, 'before readrestart' 
              IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
              print *, NCF_RRKRE%FNAME
         CALL READRESTART(NCF_RRKRE,RRK_EOFCONTR)
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
           print *, 'finish readrestart'
         IntTime = REF_START_TIME
         ExtTime = REF_START_TIME

         iint = rkint_start         
  
!==============================================================================!
!  PREPARE TO START FVCOM'S MAIN LOOP                                          !
!==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "============   ", 'RRK_NEOF =', RRK_EOFCONTR,"=============="
        write(IPT,*) "===================================================="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if
    IF (starttime%MuSOD ==0 )  THEN
     IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
     CALL ALLOC_BUFFER

     ! SET THE ASSIMILATION/SIMULATION RESET TIME
     ASSIM_RESTART_TIME = StartTime


     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
     INT_START = ISTART
     INT_COUNT = ISTART
     INT_END   = ISTART + 2*(IEND-ISTART)

     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
        Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     end if  !(starttime%MuSOD==0)

     DO WHILE(IntTime < EndTime)
       IF (starttime%MuSOD ==0 )  THEN
          IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
	  
        ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY

        IF(SSTGRD_ASSIM) THEN
          SST_SAVE_INDEX = 0
          SST_SAVE_TIME  = IntTime + sst_save_interval
        ENDIF


        IF(SSHGRD_ASSIM) THEN
          SSH_SAVE_INDEX = 0
          SSH_SAVE_TIME  = IntTime + ssh_save_interval
        ENDIF

        CALL SAVE_STATE

        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Simulation  ===================="
        end if

        !==============================================================================!
        !    RUN PURE SIMULATION MODE:
        !==============================================================================!

        ASSIM_FLAG = 0
        
        IF(SSHGRD_ASSIM) THEN
           CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        IF(SSTGRD_ASSIM) THEN
           CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        RECALC_RHO_MEAN%mjd = Inttime%mjd
        RECALC_RHO_MEAN%musod=0
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

        END DO

        !==============================================================================!
        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
        !==============================================================================!
        if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
        CALL RESTORE_STATE
        CALL RHO_PMEAN

        IF(SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE

        IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

        if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="
        ELSE  !(starttime%MuSOD  .NOT 0 )
          IF (MSR) PRINT *, 'DO NOT NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
        END IF
        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Assimilation  ===================="
        end if

        !==============================================================================!
        !    RUN DATA ASSIMILATION MODE:
        !==============================================================================!
        ASSIM_FLAG = 1
        IF(SSHGRD_ASSIM) THEN
           CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        IF(SSTGRD_ASSIM) THEN
           CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        RECALC_RHO_MEAN%mjd = Inttime%mjd
        RECALC_RHO_MEAN%musod=0

        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')


           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

        END DO ! ASSIMILATION STAGE LOOP
     

     END DO ! ASSIMILATION TOTAL DAY CONTROL LOOP

     CALL WRITERESTART(NCF_RRKFCT,RRK_EOFCONTR)  
     
   END DO  ! RRK_NEOF

         
    print *, 'before rrk_rrk(2)'  
   CALL RRK_RRK(2)   ! calculate linearized model matrix
    print *, 'before rrk_rrk(4)'
   
   IF(MSR) CALL RRK_RRK(4) ! doubling algorithm

   
!JQI4000 CONTINUE   

   IF(MSR) WRITE(IPT,*)
   IF(MSR) WRITE(IPT,*) 'Finish Kalman filter preparation! '
   IF(MSR) WRITE(IPT,*) 'Start data assimilation! '
   IF(MSR) WRITE(IPT,*)
   
   CALL NC_OPEN(NC_START)  

   
   CALL READRESTART(NC_START,REF_START_TIME)
   
   
   IntTime = RRK_START_TIME
   ExtTime = RRK_START_TIME

   iint = rkint_start 
   

!==============================================================================!
!  PREPARE TO START FVCOM'S MAIN LOOP                                          !
!==============================================================================!
   if(DBG_SET(dbg_log)) THEN
      write(IPT,*) "===================================================="
      write(IPT,*) "===================================================="
      write(IPT,*) "============== STARTING MAIN LOOP AT:==============="
      if(DBG_SET(dbg_log)) &
           & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
      write(IPT,*) "===================================================="
   end if
   
   CALL REPORT('INITIAL CONDITIONS')
   
   if(DBG_SET(dbg_log)) THEN
      write(IPT,*) "===================================================="
      write(IPT,*) "===================================================="
      write(IPT,*) "===================================================="
   end if
   
   ICYC = 1
   RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)
     CALL ALLOC_BUFFER

     ! SET THE ASSIMILATION/SIMULATION RESET TIME
     ASSIM_RESTART_TIME = StartTime
     IF(TSGRD_ASSIM) THEN
       CALL Now_2_month_days(StartTime,Pyear,Pmonth,Pmdays)
       TSC_SAVE_TIME2  = IntTime + days2time(real(Pmdays,kind=DP))
     ENDIF

     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
     INT_START = ISTART
     INT_COUNT = ISTART
     INT_END   = ISTART + 2*(IEND-ISTART)

     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
        Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     IF(TSGRD_ASSIM) THEN
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "SETUP ASSIMILATION RSTFILE"
       CALL SETUP_RSTFILE_ASSIM
     ENDIF

     DO WHILE(IntTime < EndTime)

        ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY

        IF(SSTGRD_ASSIM) THEN
          SST_SAVE_INDEX = 0
          SST_SAVE_TIME  = IntTime + sst_save_interval
        ENDIF

        IF(TSGRD_ASSIM) THEN
          TSC_SAVE_INDEX = 0
          TSC_SAVE_TIME = IntTime + TSC_SAVE_INTERVAL 

          CALL Now_2_month_days(IntTime+IMDTI,Pyear,Pmonth,Pmdays)

          IF(ASSIM_RESTART_TIME > TSC_SAVE_TIME2) THEN
            TSC_SAVE_INDEX2 = 1
            TSC_SAVE_TIME2 = TSC_SAVE_TIME2 + days2time(real(Pmdays,kind=DP))
          ELSE
            TSC_SAVE_INDEX2 = 0
          ENDIF
        ENDIF

        IF(SSHGRD_ASSIM) THEN
          SSH_SAVE_INDEX = 0
          SSH_SAVE_TIME  = IntTime + ssh_save_interval
        ENDIF

        CALL SAVE_STATE

        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Simulation  ===================="
        end if

        !==============================================================================!
        !    RUN PURE SIMULATION MODE:
        !==============================================================================!

        ASSIM_FLAG = 0
        
        IF(SSHGRD_ASSIM) THEN
           CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        IF(SSTGRD_ASSIM) THEN
           CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        RECALC_RHO_MEAN%mjd = Inttime%mjd
        RECALC_RHO_MEAN%musod=0
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

        END DO

        !==============================================================================!
        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
        !==============================================================================!
        if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
        CALL RESTORE_STATE
        CALL RHO_PMEAN

        IF(SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE
        IF(TSGRD_ASSIM) CALL TSGRD_OBSERVATION_UPDATE(IntTime+(days2time(1.0_DP))/2)
        IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

        if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="

        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Assimilation  ===================="
        end if

        !==============================================================================!
        !    RUN DATA ASSIMILATION MODE:
        !==============================================================================!
        ASSIM_FLAG = 1
        IF(SSHGRD_ASSIM) THEN
           CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        IF(SSTGRD_ASSIM) THEN
           CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        RECALC_RHO_MEAN%mjd = Inttime%mjd
        RECALC_RHO_MEAN%musod=0

        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
           !==============================================================================!
!           CALL ARCHIVE

!           CALL DUMP_PROBE_DATA
           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

           !==============================================================================!
           !    LAGRANGIAN PARTICLE TRACKING                                              |
           !==============================================================================!
!# if defined (LAG_PARTICLE)
!           CALL LAG_UPDATE
!# endif

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
!        IF(NCNEST_ON) CALL ARCHIVE_NEST

      !==============================================================================!
      !    KALMAN FILTERING                                                          |
      !==============================================================================!
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
                if (msr) print *, 'finish the filter_obs_data'
      IF(IntTime>=RRK_START_TIME .AND. IntTime<=RRK_END_TIME) THEN
         IF(IntTime == RRK_CYC) THEN
            CALL PRINT_TIME(IntTime,IPT,"DOING RRK_RRKF")
#if defined (MULTIPROCESSOR)
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
#endif 
                CALL filter_obs_data
#if defined (MULTIPROCESSOR)
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
#endif
            CALL set_rrkf_assim_data
               print *, 'kf observation location is:', nloc
               if (nloc>0) then
            print *, 'rrk_rrkf'
            CALL RRK_RRKF
               end if
           call deallocate_obs_data
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
           print *, 'finish rrkf'
	    ICYC = ICYC+1
            RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)	    
         END IF
      END IF


        END DO

        ! RUN THE NEXT DAY

     END DO
!qxu{
# endif
# endif  
    
!============================================END RRKF_WITH_SSA================================| 





  CASE(FVCOM_ENKF_WITHOUT_SSA)
# if defined (ENKF)
  IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'
  ICYC = 1
  enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
  DO WHILE(IntTime < ENDTIME) ! biggest control runtime
   if(msr)   print *,'start a forecast cycle'

   kf_times=kf_times+1 ! 


   DO ENKF_memberCONTR =1, ENKF_NENS 
    
    CALL ENKF_RESTART
    if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "============   ", 'ENKF_NENS =', ENKF_memberCONTR,"=============="
        write(IPT,*) "===================================================="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "===================================================="
     end if
!     istart=iint+1
!     NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)
!     IEND = ISTART + NSTEPS
     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     !////////////////////////// MAIN LOOP //////////////////////////////////////////
     DO while (IntTime < EndTime)

        IntTime=IntTime + IMDTI
        IINT = IINT + 1
        CALL INTERNAL_STEP

        !==============================================================================!
        !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
        !==============================================================================!
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        
        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

        !==============================================================================!
        !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
        !==============================================================================!

!for Eulerian velocity output
# if defined (WAVE_CURRENT_INTERACTION)
         !U = U-U_STOKES_3D
         !V = V-V_STOKES_3D
         !UA=UA-U_STOKES_2D
         !VA=VA-V_STOKES_2D
# endif


IF (ENKF_memberCONTR ==1 ) THEN
           if (msr) print *, "call archive for ensemble member: ", ENKF_memberCONTR
           CALL ARCHIVE

           END IF

    IF(NC_ON)THEN
       ! bounds checking
       !IF(NC_DAT%FTIME%NEXT_IO == IntTime .or. FORCE_ARCHIVE) THEN
                 if (inttime > enkf_out_Time+0.1_SP*IMDTI) then
                 print *, 'inttime',inttime
                 print *, 'enkf_out_time',enkf_out_time
                 Call Fatal_Error("The time of the output &
               &must be greater than or equal to the current time")
                 end if
       IF(abs(enkf_out_Time -IntTime)<0.1_SP*IMDTI) THEN
             call    ENKF_SETUP_NC(NC_DAT)
             IF (ENKF_memberCONTR == ENKF_NENS) then
                ARCHIVE_TIMES=ARCHIVE_TIMES+1
                enkf_out_Time=enkf_out_Time+NC_DAT%FTIME%INTERVAL
             end if

       END IF
    END IF


#  if defined (WAVE_CURRENT_INTERACTION)
!Convert Eulerian velocity back to Lagrangian velocity to maintain consistent with Mellor's equation
         !U = U+U_STOKES_3D
         !V = V+V_STOKES_3D
         !UA=UA+U_STOKES_2D
         !VA=VA+V_STOKES_2D
#  endif

        CALL DUMP_PROBE_DATA
        IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
# if defined (WAVE_CURRENT_INTERACTION)
        IF(OUT_WAVE_SPARSE_TIMESERIES_ON)CALL OUT_WAVE_SPARSE_TIMESERIES
# endif
        !==============================================================================!
        !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
        !==============================================================================!
        CALL SHUTDOWN_CHECK(D1)

        !==============================================================================!
        !  CALL BOUNDS CHECK TO SEE IF VARIABLES EXCEED USER-DEFINED THRESHOLDS 
        !==============================================================================!
        CALL BOUNDSCHK  !bounds checking

        !==============================================================================!
        !    LAGRANGIAN PARTICLE TRACKING                                              |
        !==============================================================================!
# if defined (LAG_PARTICLE)
        CALL LAG_UPDATE
# endif

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
        IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
        IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif
          if (msr ) print *,  'inttime',inttime, 'enkf_cyc',enkf_cyc
          IF(IntTime == enkf_cyc) EXIT ! stop run when KF launch , exit from assimilation stage
     END DO
     !////////////////////////// END MAIN LOOP //////////////////////////////////////
     write(enkf_num,'(i2.2)')  ENKF_memberCONTR
     NCF_ENKFfct => SETUP_RESTART(TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//enkf_num//".nc")


     CALL ENKF_WRITERESTART(NCF_ENKFfct)
    end do ! enkf_nens
         IF(IntTime == enkf_cyc) THEN
            CALL PRINT_TIME(IntTime,IPT,"DOING ENKF")
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
              if (enkf_test) then
               call getnloc
              else 
               CALL filter_obs_data
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
                if (msr) print *, 'finish the filter_obs_data'
                CALL set_enkf_assim_data
               end if
!           IF (MSR) THEN
               print *, 'kf observation location is:', nloc
               if (nloc>0) then
               CALL ENKF_ASS_EVE
               end if
!           END IF
           IF (MSR) THEN
               if (ICYC==1 .AND. NLOC==0) THEN
                    open(435,file='idobs.dat')
                    write(435,*) -31,-711
                    close(435)

               END IF
               PRINT *, 'FINISH ENKF, ',inttime
           END IF
           if (.not. enkf_test) then
           call deallocate_obs_data
           end if 
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
           ICYC = ICYC+1
            enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
          END IF

   END DO  ! end biggest run control


   IF(MSR) WRITE(IPT,*)
   IF(MSR) WRITE(IPT,*) 'Finish Ensemble Kalman filter! '
   IF(MSR) WRITE(IPT,*)
   IF(MSR) WRITE(IPT,*)

# endif
    
!============================================END ENKF_WITHOUT_SSA================================|     

  CASE(FVCOM_ENKF_WITH_SSA)
# if defined (ENKF)
# if defined (DATA_ASSIM) 
   IF (SSTGRD_ASSIM .or. SSHGRD_ASSIM) THEN
   ELSE
   PRINT *, 'STOP, ENKF ONLY WITH SSTGRD_ASSIM OR SSHGRD_ASSIM'
   CALL PSTOP
   END IF
  IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'
  ICYC = 1
  enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
  if(msr)  print *, "inttime",inttime
  if(msr)  print *,  "endtime",endtime
  DO WHILE(IntTime < ENDTIME) ! biggest control runtime
   if(msr)   print *,'start a forecast cycle'

   kf_times=kf_times+1 ! 


   DO ENKF_memberCONTR =1, ENKF_NENS 
    
    CALL ENKF_RESTART
    if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "============   ", 'ENKF_NENS =', ENKF_memberCONTR,"=============="
        write(IPT,*) "===================================================="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if
!============================SSH AND SST ASSIMILATION CODE BELOW=============== 
     IF (starttime%MuSOD ==0 )  THEN
     IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
     CALL ALLOC_BUFFER

     ! SET THE ASSIMILATION/SIMULATION RESET TIME
     ASSIM_RESTART_TIME = StartTime


     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
     INT_START = ISTART
     INT_COUNT = ISTART
     INT_END   = ISTART + 2*(IEND-ISTART)

     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
        Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     end if  !(starttime%MuSOD==0)

     DO WHILE(IntTime < EndTime)
       IF (starttime%MuSOD ==0 )  THEN
          IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
	  
        ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY

        IF(SSTGRD_ASSIM) THEN
          SST_SAVE_INDEX = 0
          SST_SAVE_TIME  = IntTime + sst_save_interval
        ENDIF


        IF(SSHGRD_ASSIM) THEN
          SSH_SAVE_INDEX = 0
          SSH_SAVE_TIME  = IntTime + ssh_save_interval
        ENDIF

        CALL SAVE_STATE

        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Simulation  ===================="
        end if

        !==============================================================================!
        !    RUN PURE SIMULATION MODE:
        !==============================================================================!

        ASSIM_FLAG = 0
        
        IF(SSHGRD_ASSIM) THEN
           CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        IF(SSTGRD_ASSIM) THEN
           CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        RECALC_RHO_MEAN%mjd = Inttime%mjd
        RECALC_RHO_MEAN%musod=0
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

        END DO

        !==============================================================================!
        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
        !==============================================================================!
        if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
        CALL RESTORE_STATE
        CALL RHO_PMEAN

        IF(SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE

        IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

        if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="
        ELSE  !(starttime%MuSOD  .NOT 0 )
          IF (MSR) PRINT *, 'DO NOT NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
        END IF
        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Assimilation  ===================="
        end if

        !==============================================================================!
        !    RUN DATA ASSIMILATION MODE:
        !==============================================================================!
        ASSIM_FLAG = 1
        IF(SSHGRD_ASSIM) THEN
           CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        IF(SSTGRD_ASSIM) THEN
           CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
        ENDIF
        RECALC_RHO_MEAN%mjd = Inttime%mjd
        RECALC_RHO_MEAN%musod=0

        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
           !==============================================================================!
	   IF (ENKF_memberCONTR ==1 ) THEN
	   if (msr) print *, "call archive for ensemble member: ", ENKF_memberCONTR
           CALL ARCHIVE

       END IF
           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

           !==============================================================================!
           !    LAGRANGIAN PARTICLE TRACKING                                              |
           !==============================================================================!
# if defined (LAG_PARTICLE)
           CALL LAG_UPDATE
# endif

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
        IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
        IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif	
          if (msr ) print *,  'inttime',inttime, 'enkf_cyc',enkf_cyc
          IF(IntTime == enkf_cyc) EXIT ! stop run when KF launch , exit from assimilation stage
        END DO ! ASSIMILATION STAGE LOOP
     IF(IntTime == enkf_cyc) EXIT ! not run one day, exit from assimilation loop when KF should laungh
        

     END DO ! ASSIMILATION TOTAL DAY CONTROL LOOP
     write(enkf_num,'(i2.2)')  ENKF_memberCONTR
     NCF_ENKFfct => SETUP_RESTART(TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//enkf_num//".nc")


     CALL ENKF_WRITERESTART(NCF_ENKFfct)
    end do ! enkf_nens
         IF(IntTime == enkf_cyc) THEN
 
            CALL PRINT_TIME(IntTime,IPT,"DOING ENKF")
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
!               call getobsloc
               CALL filter_obs_data
                 IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
                if (msr) print *, 'finish the filter_obs_data'
                CALL set_enkf_assim_data
           IF (MSR) THEN
               print *, 'kf observation location is:', nloc
               if (nloc>0) then
               CALL ENKF_ASS_EVE
               end if
               PRINT *, 'FINISH ENKF, ',inttime
           END IF
           call deallocate_obs_data
           IF(PAR)CALL MPI_BARRIER(il_commlocal,IERR)
           ICYC = ICYC+1
            enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
          END IF

   END DO  ! end biggest run control


   IF(MSR) WRITE(IPT,*)
   IF(MSR) WRITE(IPT,*) 'Finish Ensemble Kalman filter! '
   IF(MSR) WRITE(IPT,*)
   IF(MSR) WRITE(IPT,*)

# endif
# endif
    
!============================================END ENKF_WITH_SSA================================|     
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN FVCOM_RUN_MODE :'"//TRIM(FVCOM_RUN_MODE),&
          & "Options are the following: '"//TRIM(FVCOM_PURE_SIM)//"&
          &' OR '"//TRIM(FVCOM_NUDGE_OI_ASSIM) )
  END SELECT

#  if defined (SEMI_IMPLICIT) || defined (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
   CALL PETSc_CLEANUP
#  endif


#  if defined (DATA_ASSIM)   
   IF(TSGRD_ASSIM) THEN
     TSC_SAVE_INDEX2 = 1
     CALL SAVE_STATE  !! save the information netcdf file
   ENDIF
#  endif

#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
   CALL SWMAIN_CLOSE
#  endif

  
  if(DBG_SET(dbg_log)) write(IPT,*)"TADA!"
  CALL PSHUTDOWN

END PROGRAM FVCOM

