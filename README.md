# beamingfacotr
A pre-computed beaming factor library is presented in BeamingFactor_Lib. (will soon uploaded in several days)

We present python based program BeamingFactor for users who want to compute factors themselves.
Installation(will be available soon in several days):
    pip install beamingfactor 
Usage:
    Users can specific their own flux calibrated spectra
    (in ascii format, first column being wavelength in Angstorm and second column being flux,seperated by comma.
    First line can be comments begin with '#'.)
    or synthetic spectra (direct download from PHOENIX, Atlas9, or TMAP).
    Supported spectype included: ATLAS9, PHOENIX, TMAP, and USER.
    
    Filter name should be given with band='Filter_Name'.
    Support Filter_Name: u,g,r,i,z,TESS,G,Gbp,Grp,U,B,V,R,I,kepler,NUV_CSST,u_c,g_c,r_c,i_c,z_c,y_c
    u,g,r,i,z are from SDSS;
    NUV_CSST,u_c,g_c,r_c,i_c,z_c,y_c are from CSST
    Example:
        factor(specfile='/path/to/file/file_name',spectype='ATLAS9',K=100) 
