# Beaming Facor Library
**A pre-computed beaming factor library is presented in directory *BeamingFactor_Lib***.  

# BeamingFactor  
## Installation:
&ensp;&ensp;We present python based program BeamingFactor for users who want to compute factors themselves.  

    pip install beamingfactor  
## Quick Start Usage:  
Users should specify at least three components for a run:  
- the spectrum file  
- the filter transmission curve  
- orbital parameter  

The program would return a list.  
If Constant_Factor is True:
- The first element is the beaming factor
- The second element is *D* index, default is *D*1, *D*2, and *D*5.
- The third and fourth elements are arrays about the radial velocity and corresponding beaming flux, respectively.  

If Constant_Factor is False:
- The first and second elements are arrays about the radial velocity and corresponding local beaming factor, respectively.
### Examples:

    from beamingfactor import BeamingFactor as bf
    bf.factor(specfile='/path/to/file/file_name',spectype='PHOENIX',band='V',K=100)  
### Input Spectra Files:  
Supported spectral type included: ATLAS9, PHOENIX, TMAP, and USER.  
Users can specify synthetic spectra (files direct download from [PHOENIX](https://phoenix.astro.physik.uni-goettingen.de/?page_id=15)[^1](Must use high resolution file), [ATLAS9](https://wwwuser.oats.inaf.it/fiorella.castelli/grids.html), or [TMAP](http://astro.uni-tuebingen.de/~rauch/TMAF/flux_H+He.html)), for example:  

    bf.factor(specfile='/path/to/file/file_name',spectype='ATLAS9',band='V',K=100)
Alternatively, users can specific their own flux calibrated spectra as follows:  

    bf.factor(specfile='/path/to/file/file_name',spectype='USER',band='V',K=100)  
- Files should be an **ascii** format **csv** file.  
- First column should be **wavelength** in Angstorm and second column being **flux**, seperated by comma.
- First line can be comments or colnames begin with **'#'**.
  
An example spectrum file:  

    #Wave,Flux  
    200,2.0  
    3000,2.1  
### Input Filter Passband:    
Filter name should be given with band='Filter_Name'.  
Support Filter_Name: u,g,r,i,z,TESS,G,Gbp,Grp,U,B,V,R,I,kepler,NUV_CSST,u_c,g_c,r_c,i_c,z_c,y_c.  
>    u,g,r,i,z are from SDSS;
> 
>    NUV_CSST,u_c,g_c,r_c,i_c,z_c,y_c are from CSST
>
> TESS for TESS
>
> kepler for *Kepler*
> 
For Example:  

    bf.factor(specfile='/path/to/file/file_name',spectype='ATLAS9',band='G',K=100)  
Alternatively, users can specify their own transmission curve as follows:

    bf.factor(specfile='/path/to/file/file_name',spectype='USER',band='USER',bandfile='/path/to/your/file',bandtype='E',K=100)  
- Files should be an **ascii** format **csv** file.  
- First column should be **wavelength** in Angstorm and second column being corresponding transmission efficiency, seperated by comma.
- First line can be comments or colnames begin with **'#'**.
  
An example filter transmission curve file:  

    #Wave,T  
    200,0.01  
    3000,0.5  
The default filter type is for energy counter detector (bandtype='E'). For photon counter detertor, please give bandtype='P'.
### Orbital Parameters:
Users can directly define the radial velocity speed range with RV_min (km/s) and RV_max (km/s).  
Alternatively, users can defined the orbaital parameters with K (km/s), V0 (km/s), e, omega (degree).  
For example:  

    bf.factor(specfile='/path/to/file/file_name',spectype='ATLAS9',band='G',K=100,V0=40,e=0.1,omega=20)
    bf.factor(specfile='/path/to/file/file_name',spectype='ATLAS9',band='G',RV_min=-100,RV_max=+120)
### Other Parameters:
#### Constant_Factor
By default, Constant_Factor is True. Program output would include an constant beaming factor.  
If Constant_Factor is set to be False, program output would include arrays of radial velocity and corresponding local beaming factor.  
#### Dindex
By default, Dindex is set to None and program output would include *D* indexes at 1, 2, and 5.  
Users can give a list to obtain different *D* indexes.  
For example: Dindex=[1.5,10] will give *D*1.5 and *D*10.  
[^1]: Must use high resolution (R=500,000) file.
