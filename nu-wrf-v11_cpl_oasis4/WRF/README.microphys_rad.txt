
Replace all files in their corresponding directories (Registry, dyn_em, 
and phys) with files from our tar file. Then, you need to do "clean -a"
before you re-compile your WRF codes.

To use gsfcgce microphysics, set mp_physics = 7 in your namelist.input.
To use gsfc new lw and sw radiation , set 
          ra_lw_physics = 5 in your namelist.input.
          ra_sw_physics = 5 in your namelist.input.

By default, this program is pre-set for graupel and 3-ice scheme. If you wish 
to run for hail or 2-ice, please add "gsfcgce_hail" and/or "gsfcgce_2ice"
into physics section of namelist.input as the following

gsfcgce_hail = 0       for running gsfcgce microphysics with graupel
             = 1       for running gsfcgce microphysics with hail
        default value = 0
gsfcgce_2ice = 0       for running with snow, ice and graupel/hail
             = 1       for running with only ice and snow
             = 2       for running with only ice and graupel
                       (only used in very extreme situations)
             = 3       for running with warm rain only (no ice)
        default value = 0
gsfcgce_hail is always equal to 0 if gsfcgce_2ice is set to 1, 2 or 3.

If you choose to run gsfcgce microphysics with the default values, you 
don't need to add anything into namelist.input.

Adding gsfcgce microphysics and radiation into WRF doesn't affect the 
moving nest ability in WRF at all. If you wish to run WRF with moving 
nests, simply follow the instruction in configure.wrf to modfiy your 
configure.wrf.

This program will work only with WRF V3.2.1 released in 2010 but 
not any earlier versions. 

Good luck!
