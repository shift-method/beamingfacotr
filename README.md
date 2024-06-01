# Beaming Facor Library
**A pre-computed beaming factor library is presented in directory BeamingFactor_Lib**.  
(Beaming factors fro SDSS(*ugriz*) systems under **PHOENIX** will soon uploaded in several days.)

# BeamingFactor  
## Installation:
&ensp;&ensp;We present python based program BeamingFactor for users who want to compute factors themselves.  
&ensp;&ensp;(will soon be available in several days)  

     pip install beamingfactor  
## Usage:  
### Examples:

    import beamingfactor as bf
    bf.factor(specfile='/path/to/file/file_name',spectype='PHOENIX',K=100)  
### Supported Parameters
### Supported Input Spectra Files:  
Users can specific synthetic spectra (files direct download from [PHOENIX](https://phoenix.astro.physik.uni-goettingen.de/?page_id=15)[^1](Must use high resolution file), [Atlas9](https://wwwuser.oats.inaf.it/fiorella.castelli/grids.html), or [TMAP](http://astro.uni-tuebingen.de/~rauch/TMAF/flux_H+He.html)), for example:  

    bf.factor(specfile='/path/to/file/file_name',spectype='ATLAS9')
Users can specific their own flux calibrated spectra as follows:  

    bf.factor(specfile='/path/to/file/file_name',spectype='USER')  
- Files should be an **ascii** format **csv** file.  
- First column should be **wavelength** in Angstorm and second column being **flux**, seperated by comma.
- First line can be comments or colnames begin with **'#'**.
  
For example:  

    #Wave,Flux  
    200,2.0  
    3000,2.1  
### Supported Input Filter Passband:  
Supported spectype included: ATLAS9, PHOENIX, TMAP, and USER.  
Filter name should be given with band='Filter_Name'.  
Support Filter_Name: u,g,r,i,z,TESS,G,Gbp,Grp,U,B,V,R,I,kepler,NUV_CSST,u_c,g_c,r_c,i_c,z_c,y_c  
>    u,g,r,i,z are from SDSS;
> 
>    NUV_CSST,u_c,g_c,r_c,i_c,z_c,y_c are from CSST
>
> TESS for TESS
>
> kepler for *Kepler*
> 
For Example:  

    bf.factor(specfile='/path/to/file/file_name',spectype='ATLAS9',K=100) 


[^1]: Must use high resolution (R=500,000) file.
