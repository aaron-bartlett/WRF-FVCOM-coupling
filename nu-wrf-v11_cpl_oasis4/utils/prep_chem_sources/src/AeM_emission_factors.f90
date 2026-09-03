!###############################################################################!
!  CCATT-BRAMS/MCGA-CPTEC/WRF emission model    CPTEC/INPE                      !
!  Version 1.0.0: 12/nov/2010                                                   !
!  Coded by Saulo Freitas and Karla Longo                                       !
!  Contact: gmai@cptec.inpe.br - http://meioambiente.cptec.inpe.br              !
!###############################################################################!

module AeM_emission_factors
! ------------------------------------------------
! Índices para os dados da matriz emission-factors
! ------------------------------------------------
! Coluna - Bioma
! ------------------------------------------------

integer	:: &
	TropFor       = 01, & ! TropFor
	ExtratropF    = 02, & ! ExtratropF 
	Savanna       = 03, & ! Savanna
	pasture       = 04, & ! pasture/cropland       
	Biofuel       = 05, & ! Biofuel       
	CharcMak      = 06, & ! CharcMak
	CharcBurn     = 07, & ! CharcBurn
	AgResid       = 08, & ! AgResid
	Laboratory    = 09, & ! Laboratory
        Molec_Weight  = 10, & ! Molec. Weight
	
! ------------------------------------------------
! Linha - Gas
! ------------------------------------------------
	
	CO2    = 01, &        !CO2
	CO     = 02, &        !CO
	CH4    = 03, &        !CH4
	NHMC   = 04, &        !NHMC
	
	C2H2   = 05, &        !C2H2
	C2H4   = 06, &        !C2H4
	C2H6   = 07, &        !C2H6
	C3H4   = 08, &        !C3H4
	C3H6   = 09, &        !C3H6
	C3H8   = 10, &        !C3H8
	
    	butene_1     = 11, &  !1-butene
	butene_i     = 12, &  !i-butene
	butene_tr_2  = 13, &  !tr-2-butene
	butene_cis_2 = 14, &  !cis-2-butene
	butadiene    = 15, &  !butadiene
	butane_n     = 16, &  !n-butane
	butane_i     = 17, &  !i-butane
	
    	pentene_1    = 18, &  !1-pentene
	pentene_2    = 19, &  !2-pentene (cis&trans)
	pentane_n    = 20, &  !n-pentane
	Butene_2_Me  = 21, &  !2-Me-Butene
	butane_2_Me  = 22, &  !2-Me-butane
	pentadienes  = 23, &  !pentadienes
	Isoprene     = 24, &  !Isoprene
	
    	cyclopentene	= 25, & !cyclopentene
	cyclopentadiene = 26, & !cyclopentadiene
	pentene_4_me_1  = 27, & !4-me-1-pentene
	pentene_2_me_1  = 28, & !2-me-1-pentene
	
        hexene_1    = 29, &   !1-hexene
	hexadienes  = 30, &   !hexadienes
	hexane_n    = 31, &   !n-hexane
	isohexanes  = 32, &   !isohexanes
	
        heptane       = 33, & !heptane
	octenes       = 34, & !octenes
	terpenes      = 35, & !terpenes
	benzene       = 36, & !benzene
	toluene       = 37, & !toluene
	xylenes       = 38, & !xylenes
	ethylbenzene  = 39, & !ethylbenzene
	styrene       = 40, & !styrene
	PAH	      = 41, & !PAH
	
        Methanol       = 42, & !Methanol
	Ethanol        = 43, & !Ethanol
	Propanol_1     = 44, & !1-Propanol
	propanol_2     = 45, & !2-propanol
	Butanols       = 46, & !Butanols
	cyclopentanol  = 47, & !cyclopentanol
	phenol         = 48, & !phenol
	
        Formaldehyde	    = 49, &!Formaldehyde
	Acetald 	    = 50, &!Acetald
	Hydroxyacetaldehyde = 51, &!Hydroxyacetaldehyde
	Acrolein	    = 52, &!Acrolein
	
        Propanal  = 53, &     !Propanal
	Butanals  = 54, &     !Butanals
	Hexanals  = 55, &     !Hexanals
	Heptanals = 56, &     !Heptanals
	
        Acetone 	= 57, &  !Acetone
	Butanone_2	= 58, &  !2-Butanone
	Butanedione_2_3 = 59, &  !2,3-Butanedione
	Pentanones	= 60, &  !Pentanones
	Hexanones	= 61, &  !Hexanones
	Heptanones	= 62, &  !Heptanones
	Octanones	= 63, &  !Octanones
	
        Benzaldehyde	 = 64, &  !Benzaldehyde   
	Furan		 = 65, &  !Furan
	Me_Furan_2       = 66, &  !2-Me-Furan
	Me_Furan_3	 = 67, &  !3-Me-Furan
	ethylfuran_2     = 68, &  !2-ethylfuran
	dime_furan_2_4   = 69, &  !2,4-dime-furan
	Dime_furan_2_5   = 70, &  !2,5-Dime-furan
	Tetrahydrofuran  = 71, &  !Tetrahydrofuran
	dihydrofuran_2_3 = 72, &  !2,3-dihydrofuran
	benzofuran	 = 73, &  !benzofuran
	Furfural	 = 74, &  !Furfural
	
        Me_format        = 75, &  !Me-format
	Me_Acetate       = 76, &  !Me-Acetate
	Acetonitrile     = 77, &  !Acetonitrile
	Acrylonitrile    = 78, &  !Acrylonitrile
	Propionitrile    = 79, &  !Propionitrile
	
        pyrrole 	  = 80, & !pyrrole
	trimethylpyrazole = 81, & !trimethylpyrazole
	methylamine	  = 82, & !methylamine
	dimethylamine	  = 83, & !dimethylamine
	ethylamine	  = 84, & !ethylamine
	trimethylamine    = 85, & !trimethylamine
	n_pentylamine	  = 86, & !n-pentylamine
	me_1_butylamine_2 = 87, & !2-me-1-butylamine
	
        HFo	  = 88, &  !HFo
	HAc	  = 89, &  !HAc
	Propanoic = 90, &  !Propanoic acid
	H2	  = 91, &  !H2
	NOx	  = 92, &  !NOx
	NOy	  = 93, &  !NOy
	
        N2O       = 94, &  !N2O
	NH3       = 95, &  !NH3
	HCN       = 96, &  !HCN
	cyanogen  = 97, &  !cyanogen
	N2	  = 98, &  !N2
	SO2	  = 99, &  !SO2
	DMS	  =100, &  !DMS
	
        COS_	  =101, &  !COS
	CH3Cl	  =102, &  !CH3Cl
	CH3Br	  =103, &  !CH3Br
	CH3I	  =104, &  !CH3I   
        Hg	  =105, &  !Hg

	BBURN2	  =106, &  !PM2.5
		  
    BBURN3	  =107, &  !old - TPM

	TC	  =108, &  !TC
	OC	  =109, &  !OC
!	BC	  =110     !BC
! NUWRF EMK...Added dry mass burned (DM), which will be passed to WRF-Chem. 
	BC	  =110, &  !BC
        DM        =111     !DM

!---------------------------------------------------------------------------------------------------
! NUWRF EMK...Added dry mass burned (DM), which will be passed to WRF-Chem. 
!integer, parameter :: AeM_nspecies=110
integer, parameter :: AeM_nspecies=111
character(LEN=20),dimension(AeM_nspecies),parameter :: AeM_spc_name= &
! avoid tab character
!'12345678901234567890'
(/                     &
 'CO2                 '&     
,'CO                  '&     
,'CH4                 '&     
,'NHMC                '&     
,'C2H2                '&     
,'C2H4                '&     
,'C2H6                '&     
,'C3H4                '&     
,'C3H6                '&     
,'C3H8                '&     
,'1_butene            '&  
,'i-butene            '&  
,'tr_2_butene         '&  
,'cis_2_butene        '&  
,'butadiene           '&  
,'n_butane            '&  
,'i-butane            '&  
,'1_pentene           '&  
,'2_pentene           '&  
,'n_pentane           '&  
,'2_Me_Butene         '&  
,'2_Me_butane         '&  
,'pentadienes         '&  
,'Isoprene            '&  
,'cyclopentene        '&   
,'cyclopentadiene     '&   
,'4_me_1_pentene      '&   
,'2_me_1_pentene      '&   
,'1_hexene            '&  
,'hexadienes          '&  
,'n_hexane            '&  
,'isohexanes          '&  
,'heptane             '&	 
,'octenes             '&	 
,'terpenes            '&	 
,'benzene             '&	 
,'toluene             '&	 
,'xylenes             '&	 
,'ethylbenzene        '&   
,'styrene             '&   
,'PAH                 '&   
,'Methanol            '&
,'Ethanol             '&
,'1_Propanol          '&
,'2_propanol          '&
,'Butanols            '&  
,'cyclopentanol       '&
,'phenol              '&
,'Formaldehyde        '& 
,'Acetald             '& 
,'Hydroxyacetaldehyde '& 
,'Acrolein            '& 
,'Propanal            '&    
,'Butanals            '&    
,'Hexanals            '&    
,'Heptanals           '&    
,'Acetone             '&
,'2_Butanone          '&
,'2_3_Butanedione     '&
,'Pentanones          '&
,'Hexanones           '&
,'Heptanones          '&
,'Octanones           '&
,'Benzaldehyde        '&
,'Furan               '&
,'2_Me_Furan          '&
,'3_Me_Furan          '&
,'2_ethylfuran        '&
,'2_4_dime_furan      '&
,'2_5_Dime_furan      '&
,'Tetrahydrofuran     '&
,'2_3_dihydrofuran    '&
,'benzofuran          '&
,'Furfural            '&
,'Me_format           '&
,'Me_Acetate          '&
,'Acetonitrile        '&
,'Acrylonitrile       '&
,'Propionitrile       '&
,'pyrrole             '&
,'trimethylpyrazole   '&
,'methylamine         '&
,'dimethylamine       '&
,'ethylamine          '&
,'trimethylamine      '&
,'n_pentylamine       '&
,'2_me_1_butylamine   '&
,'HFo                 '&
,'HAc                 '&
,'Propanoic           '&
,'H2                  '&
,'NOx                 '&
,'NOy                 '&
,'N2O                 '&
,'NH3                 '&
,'HCN                 '&
,'cyanogen            '&
,'N2                  '&
,'SO2                 '&
,'DMS                 '&
,'COS                 '&
,'CH3Cl               '&
,'CH3Br               '&
,'CH3I                '&
,'Hg                  '&
,'BBURN2              '&  !PM25
,'BBURN3              '&  !PM10
,'TC                  '&
,'OC                  '&
,'BC                  '&
! NUWRF EMK...Added dry mass burned (DM), which will be passed to WRF-Chem. 
,'DM                  '&
/) 	  	      

!---------------------------------------------------------------------------------------------------
real,    dimension(10,AeM_nspecies) :: emission_factor

data emission_factor/  &
!---------------------------------------------------------------------------------------------------
! Emission Factor
!---------------------------------------------------------------------------------------------------
! TropFor! Extra !Savanna! pasture | Biofuel! CharcMak ! CharcBurn ! AgResid ! Laboratory ! Molec. !
!	 !tropF  !	 ! cropland|	    !  	       !           !         !            ! Weight !
!---------------------------------------------------------------------------------------------------
  1580.,  1568., 1664.,    1664.,  1486.,     478.,       2543.,      1437.,   1476.,  44.01, & ! CO2
  103.,   106.,  63.,	   63.,    76.,       84.,        201.,       85.,     89.,    28.01, & ! CO
  6.8,    4.8,   2.2,	   2.2,    6.1,       12.2,       6.7,        4.6,     4.5,    16.04, & ! CH4
  8.1,    5.7,   3.4,	   3.4,    8.0,       17.5,       4.6,        11.2,    8.7,    15.00, & ! NHMC  									     
  0.402,  0.260, 0.269,    0.269,  0.967,     0.350,      0.245,      0.362,   0.308,  26.04, & ! C2H2
  1.938,  1.226, 0.845,    0.845,  1.604,     1.092,      0.584,      1.238,   1.095,  28.05, & ! C2H4
  1.202,  0.733, 0.325,    0.325,  0.937,     2.137,      0.656,      0.765,   1.065,  30.07, & ! C2H6
  0.013,  0.055, 0.026,    0.026,  0.027,     -9999,      0.071,      0.030,   -9999,  40.06, & ! C3H4
  1.359,  0.584, 0.342,    0.342,  1.080,     1.012,      0.399,      0.813,   0.671,  42.08, & ! C3H6
  1.042,  0.254, 0.087,    0.087,  0.389,     0.350,      0.207,      0.403,   0.151,  44.10, & ! C3H8
  0.125,  0.133, 0.081,    0.081,  0.252,     -9999,      0.113,      0.106,   -9999,  56.11, & ! 1-butene     
  0.109,  0.086, 0.038,    0.038,  0.197,     -9999,      0.087,      0.070,   0.410,  56.11, & ! i-butene
  0.050,  0.032, 0.024,    0.024,  0.155,     -9999,      0.037,      0.033,   0.040,  56.11, &! tr-2-butene
  0.042,  0.070, 0.019,    0.019,  0.093,     -9999,      0.023,      0.039,   0.026,  56.11, &! cis-2-butene
  0.094,  0.081, 0.066,    0.066,  0.197,     -9999,      0.062,      0.077,   0.136,  54.09, &! butadiene
  0.041,  0.073, 0.019,    0.019,  0.062,     -9999,      0.066,      0.047,   0.040,  58.12, &! n-butane
  0.015,  0.022, 0.006,    0.006,  0.027,     -9999,      0.010,      0.013,   -9999,  58.12, &! i-butane    
  0.056,  0.059, 0.022,    0.022,  0.316,     -9999,      0.028,      0.007,   0.011,  70.13, &! 1-pentene
  0.012,  0.012, 0.007,    0.007,  0.008,     -9999,      0.023,      0.009,   0.014,  70.13, &! 2-pentene (cis&trans)
  0.014,  0.059, 0.005,    0.005,  0.036,     -9999,      0.096,      0.022,   0.030,  72.15, &! n-pentane
  0.074,  0.033, 0.033,    0.033,  0.096,     -9999,      0.015,      0.006,   0.013,  70.13, &! 2-Me-Butene
  0.008,  0.030, 0.010,    0.010,  0.041,     -9999,      0.071,      0.015,   0.004,  72.15, &! 2-Me-butane
  -9999,  -9999, -9999,    -9999,  -9999,     -9999,      -9999,      -9999,   0.024,  -9999, &! pentadienes
  0.016,  0.103, 0.026,    0.026,  0.204,     -9999,      0.017,      0.043,   0.044,  68.12, &! Isoprene
  0.019,  0.019, 0.012,    0.012,  0.308,     -9999,      0.035,      0.016,   0.006,  68.12, &! cyclopentene
  0.027,  0.027, 0.016,    0.016,  0.020,     -9999,      0.052,      0.022,   0.023,  66.11, &! cyclopentadiene
  0.048,  0.052, 0.050,    0.050,  0.015,     -9999,      0.098,      0.041,   -9999,  84.16, &! 4-me-1-pentene
  -9999,  -9999, 0.004,    0.004,  -9999,     -9999,      -9999,      -9999,   -9999,  84.16, &! 2-me-1-pentene
  0.063,  0.094, 0.037,    0.037,  0.053,     -9999,      0.141,      0.059,   0.026,  84.16, &! 1-hexene
  -9999,  -9999, -9999,    -9999,  -9999,     -9999,      -9999,      -9999,   0.065,  -9999, &! hexadienes
  0.048,  0.052, 0.027,    0.027,  0.035,     -9999,      0.063,      0.039,   0.007,  86.17, &! n-hexane
  0.076,  0.079, 0.046,    0.046,  0.056,     -9999,      0.149,      0.063,   -9999,  86.17, &! isohexanes
  0.032,  0.033, 0.019,    0.019,  0.024,     -9999,      0.063,      0.026,   -9999,  98.17, &! heptane
  0.012,  0.005, 0.006,    0.006,  0.007,     -9999,      0.017,      0.007,   0.018,  112.22,&! octenes
  0.150,  0.223, 0.015,    0.015,  0.150,     -9999,      0.000,      0.015,   1.400,  136.24,&! terpenes
  0.402,  0.539, 0.229,    0.229,  1.331,     -9999,      1.124,      0.639,   0.355,  78.11, &! benzene
  0.249,  0.421, 0.129,    0.129,  0.713,     -9999,      0.374,      0.378,   0.319,  90.10, &! toluene
  0.059,  0.204, 0.043,    0.043,  0.370,     -9999,      0.148,      0.171,   0.087,  106.17,&! xylenes
  0.024,  0.048, 0.015,    0.015,  0.129,     -9999,      0.048,      0.056,   0.026,  106.17,&! ethylbenzene
  0.028,  0.128, 0.023,    0.023,  0.285,     -9999,      0.133,      0.028,   0.018,  104.15,&! styrene
  0.100,  0.100, 0.0024,   0.0024, 0.281,     -9999,      0.025,      0.100,   0.0640, -9999, &! PAH
  2.442,  1.866, 1.470,    1.470,  3.610,     12.300,     1.240,      2.007,   1.359,  32.04, &! Methanol
  0.017,  0.018, 0.010,    0.010,  0.013,     -9999,      0.033,      0.014,   0.097,  46.07, &! Ethanol
  0.041,  0.042, 0.025,    0.025,  0.030,     -9999,      0.080,      0.034,   0.248,  60.10, &! 1-Propanol
  -9999,  -9999, -9999,    -9999,  -9999,     -9999,      -9999,      -9999,   0.210,  -9999, &! 2-propanol
  0.009,  0.012, 0.008,    0.008,  0.008,     -9999,      0.022,      0.009,   0.025,  74.12, &! Butanols
  0.031,  0.044, 0.033,    0.033,  0.032,     -9999,      0.084,      0.035,   -9999,  86.13, &! cyclopentanol
  0.006,  0.876, 0.003,    0.003,  1.870,     2.750,      1.659,      0.699,   1.912,  94.11, &! phenol
  1.180,  1.876, 0.709,    0.709,  0.999,     1.060,      0.733,      0.970,   1.590,  30.03, &! Formaldehyde
  0.655,  0.501, 0.500,    0.500,  0.338,     -9999,      1.278,      0.538,   0.831,  45.06, &! Acetald
  0.961,  0.961, 0.961,    0.961,  0.660,     -9999,      0.961,      0.961,   0.961,  60.05, &! Hydroxyacetaldehyde
  0.178,  0.238, 0.075,    0.075,  0.047,     -9999,      0.347,      0.146,   0.372,  56.06, &! Acrolein
  0.074,  0.139, 0.008,    0.008,  0.061,     -9999,      0.145,      0.061,   0.545,  58.08, &! Propanal
  0.071,  0.206, 0.055,    0.055,  0.039,     -9999,      0.202,      0.085,   0.021,  72.11, &! Butanals
  0.031,  0.023, 0.014,    0.014,  0.082,     -9999,      0.091,      0.038,   0.003,  102.17,&! Hexanals
  0.003,  0.004, 0.003,    0.003,  0.003,     -9999,      0.008,      0.003,   0.004,  114.19,&! Heptanals
  0.694,  0.553, 0.517,    0.517,  0.210,     0.175,      1.355,      0.571,   0.726,  64.13, &! Acetone
  0.426,  0.455, 0.250,    0.250,  0.059,     -9999,      0.832,      0.350,   0.285,  72.11, &! 2-Butanone
  0.897,  0.925, 0.544,    0.544,  0.660,     -9999,      1.750,      0.737,   1.033,  86.08, &! 2,3-Butanedione
  0.028,  0.089, 0.016,    0.016,  0.035,     -9999,      0.092,      0.039,   0.038,  86.13, &! Pentanones
  -9999,  -9999, -9999,    -9999,  -9999,     -9999,      -9999,      -9999,   0.007,  -9999, &! Hexanones
  0.002,  0.006, 0.006,    0.006,  0.004,     -9999,      0.011,      0.005,   0.117,  114.19,&! Heptanones
  0.019,  0.023, 0.015,    0.015,  0.016,     -9999,      0.043,      0.018,   0.005,  128.22,&! Octanones
  0.027,  0.037, 0.030,    0.030,  0.023,     -9999,      0.070,      0.029,   0.058,  106.12,&! Benzaldehyde									      
  0.442,  0.427, 0.091,    0.091,  0.561,     1.150,      0.863,      0.363,   0.934,  68.08, &! Furan
  0.173,  0.472, 0.046,    0.046,  0.174,     -9999,      0.460,      0.194,   0.550,  82.10, &! 2-Me-Furan
  0.029,  0.052, 0.008,    0.008,  0.023,     -9999,      0.061,      0.026,   0.063,  82.10, &! 3-Me-Furan
  0.003,  0.006, 0.001,    0.001,  0.003,     -9999,      0.007,      0.003,   -9999,  94.10, &! 2-ethylfuran
  0.024,  0.019, 0.008,    0.008,  0.014,     -9999,      0.037,      0.015,   -9999,  96.13, &! 2,4-dime-furan      
  0.028,  0.053, 0.002,    0.002,  0.020,     -9999,      0.054,      0.023,   0.289,  96.13, &! 2,5-Dime-furan
  0.016,  0.022, 0.016,    0.016,  0.016,     -9999,      0.042,      0.018,   0.002,  72.12, &! Tetrahydrofuran
  0.013,  0.018, 0.013,    0.013,  0.013,     -9999,      0.033,      0.014,   -9999,  70.10, &! 2,3-dihydrofuran    
  0.015,  0.026, 0.014,    0.014,  0.016,     -9999,      0.042,      0.017,   0.014,  118.14,&! benzofuran
  0.384,  0.454, 0.233,    0.233,  0.240,     1.049,      0.749,      0.315,   2.054,  96.09, &! Furfural
  0.024,  0.025, 0.014,    0.014,  0.018,     -9999,      0.047,      0.020,   0.043,  60.05, &! Me-format
  0.095,  0.106, 0.053,    0.053,  0.070,     -9999,      0.186,      0.078,   0.112,  74.05, &! Me-Acetate
  0.180,  0.187, 0.198,    0.198,  0.180,     -9999,      0.180,      0.180,   0.173,  41.05, &! Acetonitrile
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.040,  53.10, &! Acrylonitrile
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.020,  55.07, &! Propionitrile
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.017,  67.09, &! pyrrole
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.124,  110.16,&! trimethylpyrazole
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.057,  31.06, &! methylamine  
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.062,  45.09, &! dimethylamine        
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.010,  45.09, &! ethylamine   
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.041,  59.11, &! trimethylamine       
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.439,  87.17, &! n-pentylamine        
  -9999., -9999.,-9999.,   -9999., -9999.,    -9999,      0.000,      -9999.,  0.140,  87.17, &! 2-me-1-butylamine    
  1.207,  2.433, 0.63,     0.63,   0.404,     0.450,      0.159,      0.220,   0.770,  47.02, &! HFo
  4.337,  3.500, 2.61,     2.61,   4.719,     8.920,      3.200,      0.800,   2.813,  60.05, &! HAc
  -9999,  -9999, -9999,    -9999,  -9999,     -9999,      -9999,      -9999,   0.560,  -9999, &! Propanoic acid
  3.80,   1.81,  0.99,     0.99,   1.762,     -9999,      4.670,      1.966,   -9999., 2.02,  &! H2
  1.85,   3.00,  2.35,     2.35,   1.19,      0.05,       2.68,       2.49,    1.87,   30.01, &! NOx
  -9999., 0.64,  2.20,     2.20,   -9999.,    0.00,       0.00,       -9999.,  0.00,   0.00,  &! NOy
  0.20,   0.26,  0.21,     0.21,   0.07,      0.02,       0.35,       0.10,    0.06,   44.01, &! N2O
  1.30,   1.71,  0.74,     0.74,   1.29,      0.23,       0.97,       1.30,    1.61,   18.02, &! NH3
  0.30,   0.81,  0.23,     0.23,   0.30,      0.30,       0.30,       0.30,    0.47,   27.02, &! HCN
  0.00,   -9999.,0.00,     0.00,   -9999.,    0.00,       0.00,       -9999.,  0.01,   52.04, &! cyanogen
  2.64,   2.64,  2.64,     2.64,   2.64,      -9999,      2.64,       2.64,    4.13,   28.01, &! N2
  0.57,   1.00,  0.37,     0.37,   0.36,      0.00,       0.40,       0.40,    0.87,   64.06, &! SO2
  -9999,  -9999, 0.0011,   0.0011, -9999,     -9999,      -9999,      -9999,   -9999,  0.00,  &! DMS
  0.040,  0.033, 0.015,    0.015,  0.040,     0.040,      0.040,      0.060,   0.032,  60.07, &! COS
  0.077,  0.049, 0.065,    0.065,  0.052,     0.010,      0.012,      0.228,   0.158,  50.49, &! CH3Cl
  0.0078, 0.0032,0.0017,   0.0017, 0.0030,    0.0030,     0.0030,     0.0030,  0.0008, 94.94, &! CH3Br
  0.0068, 0.0006,0.0006,   0.0006, 0.0010,    -9999,      0.0010,     0.0010,  -9999., 141.94,&! CH3I
  5.E-05, 5.E-05,0.E+0,    0.E+0,  5.E-05,    0.E-05,     5.E-05,     5.E-05,  3.7E-5, 200.59,&! Hg 
  9.1,    13.0,  4.9,	   4.9,    5.7,       2.1,        1.6,        3.9,     5.1,    -9999, &! PM2.5
  8.5,    17.6,  8.5,	   8.5,    7.1,       4.3,        2.4,        6.9,     7.2,    -9999, &! PM10 (old TPM)
  6.6,    8.3,   3.7,	   3.7,    5.1,       -9999,      6.3,        3.4,     3.6,    12.00, &! TC
  5.2,    9.1,   3.2,	   3.2,    3.6,       -9999,      4.8,        3.3,     2.5,    12.00, &! OC
! NUWRF EMK...Added dry mass burned (DM), which will be passed to WRF-Chem. 
!  0.66,   0.56,  0.46,     0.46,   0.41,      -9999,      1.50,       0.69,    1.08,   12.00 / ! BC
  0.66,   0.56,  0.46,     0.46,   0.41,      -9999,      1.50,       0.69,    1.08,   12.00, &! BC
 -9999,  -9999, -9999,    -9999,  -9999,      -9999,     -9999,      -9999,   -9999,   -9999 / ! DM


!index_bioma =2
!index_gas=2
!print*, emission_factor(index_bioma, index_gas)
!print*,spc_name(index_gas)

end module AeM_emission_factors
