# This library contains beaming factors computed with four theory spectra library: [PHOENIX](https://phoenix.astro.physik.uni-goettingen.de/?page_id=15), [ATLAS9](https://wwwuser.oats.inaf.it/fiorella.castelli/grids.html), [TMAP](http://astro.uni-tuebingen.de/~rauch/TMAF/flux_H+He.html), and [TLUSTY](http://svo2.cab.inta-csic.es/theory/newov2/index.php?models=tlusty_mergedbin).  
# For more details on computation routine, please refer to: Zheng et al. submitted..
# File Description
Each file are archived in csv file (comma seperated ascii file).  
Files are named as *_{Filter_Name}.csv  
Filter_Name are name of filters as follows:  
>    u_SDSS,g_SDSS,r_SDSS,i_SDSS,z_SDSS are *u,g,r,i,z* from SDSS (Sloan Digital Sky Survey);
> 
>    NUVC,uC,gC,rC,iC,zC,yC are NUV,*u,g,r,i,z,y* from CSST
>
>    *U,B,V,R,I* are Johnson filter system
>
>    *G,Gbp,Grp* are *Gaia* filter system.
> 
> TESS for TESS satellite
>
> kepler for *Kepler* satellite
# PHOENIX Parameters Definition  
Factors for different alpha elements abundances are stored in eahc corresponding subdirectories.  
Each csv file contain 10 columns as follows:
- teff: Effective temperature.  
- logg: Surface gravity.  
- Z: Metallicity.  
- ebv: Stellar extinction E(B-V).  
- BF_FIT: Beaming factors calculated from linear fitting between beaming flux and radial velocities.  
- BF_ERR: Beaming factor uncertainties from linear fitting between beaming flux and radial velocities.
- BF_C: A reference index. Users can ignore this column.
- D1: D1 index. See definition in Zheng et al. submitted. 
- D2: D2 index. See definition in Zheng et al. submitted.
- D5: D5 index. See definition in Zheng et al. submitted.
# ATLAS9 Parameters Definition
Each csv file contain 13 columns as follows:
- teff: Effective temperature.  
- logg: Surface gravity. 
- [M/H]: Metallicity. **All values in csv files are multiplied by 10.**
- alpha: Alpha elements enhancement. False means solar abundance. True means alpha +0.4 enhanced.
- DY: Helium element enhancement. False for solar value. True for Helium enhanced by +0.1 dex. 
- vturb: Micro-turbulance velocity in km/s.
- ebv: Stellar extinction E(B-V).  
- BF_FIT: Beaming factors calculated from linear fitting between beaming flux and radial velocities.  
- BF_ERR: Beaming factor uncertainties from linear fitting between beaming flux and radial velocities.
- BF_C: A reference index. Users can ignore this column.
- D1: D1 index. See definition in Zheng et al. submitted. 
- D2: D2 index. See definition in Zheng et al. submitted.
- D5: D5 index. See definition in Zheng et al. submitted.
# TMAP Parameters Definition
Each csv file contain 11 columns as follows:
- teff: Effective temperature.  
- logg: Surface gravity.
- H: Hydrogen abundance in percentage by mass.
- He:	Helium abundance in percentage by mass.
- ebv: Stellar extinction E(B-V).  
- BF_FIT: Beaming factors calculated from linear fitting between beaming flux and radial velocities.  
- BF_ERR: Beaming factor uncertainties from linear fitting between beaming flux and radial velocities.
- BF_C: A reference index. Users can ignore this column.
- D1: D1 index. See definition in Zheng et al. submitted. 
- D2: D2 index. See definition in Zheng et al. submitted.
- D5: D5 index. See definition in Zheng et al. submitted.
# TLUSTY Parameters Definition
Each csv file contain 10 columns as follows:
- teff: Effective temperature.  
- logg: Surface gravity.
- Z: ratio of Z/Zsun. A value 0 refers to metal free.
- ebv: Stellar extinction E(B-V).  
- BF_FIT: Beaming factors calculated from linear fitting between beaming flux and radial velocities.  
- BF_ERR: Beaming factor uncertainties from linear fitting between beaming flux and radial velocities.
- BF_C: A reference index. Users can ignore this column.
- D1: D1 index. See definition in Zheng et al. submitted. 
- D2: D2 index. See definition in Zheng et al. submitted.
- D5: D5 index. See definition in Zheng et al. submitted.
