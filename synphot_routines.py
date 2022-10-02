import os
os.environ["PYSYN_CDBS"]
from os import listdir 
import numpy as np
import synphot as S
from synphot import  SourceSpectrum
from synphot.models import BlackBodyNorm1D
from synphot import SpectralElement, Observation, units
from synphot.models import Empirical1D
import astropy.units as u
from synphot import Observation, units
from matplotlib import pyplot as plt
import matplotlib 
import pandas as pd


from astropy.io import fits



plt.rcParams["figure.figsize"] = (8,5)

def normalize_model_spectrum(obs_spec, mod_spec, bin1, bin2):

	wl_range = range(bin1,bin2)
	obs_binflux = obs_spec.to_spectrum1d(wavelengths=wl_range)
	mod_binflux = mod_spec.to_spectrum1d(wavelengths=wl_range)

	r_flux = obs_binflux.flux/mod_binflux.flux

	med_r_flux = np.median(r_flux)

	norm_mod_spec = SourceSpectrum(Empirical1D,points=wl_range*u.angstrom,lookup_table=mod_binflux.flux*med_r_flux,keep_neg=True)
	return norm_mod_spec

def get_vega_synphot_table(filters,filter_file,filters_dir,vega,area,binrange1=1000,binrange2=50000,filter_plot=False):
	d =[]#used later for Pandas table
	vega_obsDict = {}
	for i in range(len(filters)):
	    bp=SpectralElement.from_file(filters_dir+'/'+filter_file[i])
	    
	    #convoluted way to change the wl units from nm to A
	    wl=bp.waveset*1E1
	    

	    th = bp.model.lookup_table
	    Filter=SpectralElement(Empirical1D,points=wl,lookup_table=th,keep_neg=True)
	                           
	    Filter.plot(left=1000, right=25000,title=filter_file[i])
	    binset=range(1000,50001)
	                                           
	    obs_vega_Filter=Observation(vega,Filter,binset=binset)
	    binflux=obs_vega_Filter.sample_binned(flux_unit='count',area=area)
	    #print('binflux = ',binflux)
	    
	    # Sample the "native" flux for comparison                                                                                     
	    flux_Vega=obs_vega_Filter(binset,flux_unit='count',area=area)
	    flux_Vega_sum=flux_Vega.sum()
	    #print('flux_Vega_sum = ',flux_Vega_sum)  
	    #print('flux_Vega count rate = ', obs_vega_Filter.countrate(area).value) 
	    #print('counts diffs = ', flux_Vega_sum.value/obs_vega_Filter.countrate(area).value)
	    #print(obs_vega_Filter.effstim(flux_unit=units.VEGAMAG,vegaspec=vega))
	    # Prepare Pandas table
	    d.append((filters[i],np.round(flux_Vega_sum,2),
	              np.round(obs_vega_Filter.effstim(u.ABmag),3),
	              np.round(obs_vega_Filter.pivot(None),2),
	              np.round(obs_vega_Filter.effective_wavelength(None),2)))

	    vega_obsDict[filter_file[i].strip(".txt")] = obs_vega_Filter                       
	#Create Pandas table, out of the loop
	obs_vega_Filter_table=pd.DataFrame(d,columns=('Filter','Counts_s-1','VegaMag','pivot_wl', 'effective_wl'))

	obs_vega_Filter_table = obs_vega_Filter_table.sort_values(by="effective_wl",ascending=True).reset_index(drop=True)           
	return obs_vega_Filter_table, vega_obsDict                           

def get_synphot_table(filters,filter_file,filters_dir,sp_STD,area,vega,binrange1=1000,binrange2=50000,filter_plot=False):

	obsDict = {}
	d =[]#used later for Pandas table

	obs_STD_Filters, obs_vega_Filters = [],[]
	for i in range(len(filters)):
		#print(filter_file[i])
		bp=SpectralElement.from_file(filters_dir+'/'+filter_file[i])
		#bp.waveset#.to(u.angstrom)
		#print(bp.waveset.shape)
		#convoluted way to change the wl units from micron to A. There must
		wl=bp.waveset*1E1
		print("{} wavelengths: {}".format(filter_file[i].strip(".txt"),wl))
		#print(wl.shape)
		th = bp.model.lookup_table
		#print(th.shape)
		Filter=SpectralElement(Empirical1D,points=wl,lookup_table=th,keep_neg=True)

		if filter_plot:
			Filter.plot(left=binrange1, right=binrange2,title=filter_file[i])
		
		binset=range(binrange1,binrange2+1)

		#print(binset,wl)
		try:
			obs_STD_Filter=Observation(sp_STD,Filter,binset=binset,force='taper')
		except S.exceptions.DisjointError:
			print("DisjointError for filter {}".format(filters[i]))
			continue
			
		binflux=obs_STD_Filter.sample_binned(flux_unit='count',area=area)
		# Sample the "native" flux for comparison                                           
		flux=obs_STD_Filter(binset,flux_unit='count',area=area)
		flux_sum=flux.sum()
		#print(flux_sum)

		obs_vega_Filter=Observation(vega,Filter,binset=binset)
		binflux=obs_vega_Filter.sample_binned(flux_unit='count',area=area)
		# Sample the "native" flux for comparison

		flux_Vega=obs_vega_Filter(binset,flux_unit='count',area=area)
		flux_Vega_sum=flux_Vega.sum()

		mag_STD_Filter=-2.5*np.log10(flux_sum/flux_Vega_sum)
		#print(obs_STD_Filter.effstim(flux_unit=units.VEGAMAG,vegaspec=vega))
		#print(obs_STD_Filter.effstim(u.ABmag))
		# print(filters[i],np.round(flux_sum,2),' VegaMag = ',np.round(mag_P

		#Other parameters
		# print("pivot wl:",np.round(obs_STD_Filter.pivot(),2)," effecti

		# Prepare Pandas table
		d.append((filters[i],np.round(flux_sum,2),
				  np.round(mag_STD_Filter,3),
				  np.round(obs_STD_Filter.pivot(),2),
				  np.round(obs_STD_Filter.effective_wavelength(),2)))
		

		obs_STD_Filters.append(obs_STD_Filter)
		obs_vega_Filters.append(obs_vega_Filter)

		obsDict[filter_file[i].strip(".txt")] = obs_STD_Filter
		#print(i)
	#Create Pandas table, out of the loop
	STD_table=pd.DataFrame(d,columns=('Filter','Counts_s-1','VegaMag','pivot_wl', 'effective_wl'))

	STD_table = STD_table.sort_values(by="effective_wl",ascending=True).reset_index(drop=True)

	return STD_table, obsDict #obs_STD_Filters, obs_vega_Filters

def get_synphot_spectra_with_diff_flux_units(inspec):



	print("from the original file")
	print(max(inspec(inspec.waveset)))

	#plot using the SourceSpectrum.plot method, specifying different flux u
	inspec.plot(left=1000,right=30000,flux_unit='fnu')
	inspec.plot(left=1000,right=30000,flux_unit='photlam')
	inspec.plot(left=1000,right=30000,flux_unit='flam')
	inspec.plot(left=1000,right=30000,flux_unit='Jy') #this is a Fnu-ty

	#Convert the flux *arrays* from photlam to fnu, flam, fJy
	spec_fnu=units.convert_flux(inspec.waveset,fluxes=inspec(inspec.waveset),out_flux_unit = units.FNU)
	spec_flam=units.convert_flux(inspec.waveset,fluxes=inspec(inspec.waveset),out_flux_unit = units.FLAM)
	spec_Jy=units.convert_flux(inspec.waveset,fluxes=inspec(inspec.waveset),out_flux_unit = u.Jy)
	#check:
	print("\n from converted fluxes")
	print(max(spec_fnu))
	print(max(spec_flam))
	print(max(spec_Jy))
	#ok
	         
	#then I rebuild the spectra with these these values
	spec_fnu=SourceSpectrum(Empirical1D,points=inspec.waveset,lookup_table=spec_fnu)
	spec_flam=SourceSpectrum(Empirical1D,points=inspec.waveset,lookup_table=spec_flam)
	spec_Jy=SourceSpectrum(Empirical1D,points=inspec.waveset,lookup_table=spec_Jy)

	# and create synphot spectra that have the right flux units.
	spec_fnu=spec_fnu(inspec.waveset,flux_unit=units.FNU)
	spec_flam=spec_flam(inspec.waveset,flux_unit=units.FLAM)
	spec_Jy=spec_Jy(inspec.waveset,flux_unit=u.Jy)
	#However, they are not "SourceSpectrum" so I think I have lost a good p
	#these objects, pretty bad.
	#check:
	print("\n from the rebuilt spectra")
	print(max(spec_fnu))
	print(max(spec_flam))
	print(max(spec_Jy))

	return spec_fnu, spec_flam, spec_Jy


def create_new_fits_bintab(lambda_col,flux_col,primary_hdu,err_col=None,new_fname=None):
	
	if err_col is not None:
		
		new_dhdu = fits.BinTableHDU.from_columns([lambda_col,flux_col,err_col])
		
	else:
		new_dhdu = fits.BinTableHDU.from_columns([lambda_col,flux_col])
	
	new_hdus = fits.HDUList([primary_hdu,new_dhdu])
	
	if new_fname is not None:
		
		new_hdus.writeto(new_fname,overwrite=True)
		
		
	return new_hdus

def create_new_fits_bintab2(lambda_angstrom,flx_data,primary_hdu,err_data=None,new_fname=None):

	lambda_col = fits.Column(name='WAVELENGTH',array=lambda_angstrom,format='D')
	flux_col = fits.Column(name='FLUX',array=flx_data,format='D')


	if err_data is not None:
		err_col = fits.Column(name='ERROR',array=err_data,format='D')
		new_dhdu = fits.BinTableHDU.from_columns([lambda_col,flux_col,err_col])
		
	else:
		new_dhdu = fits.BinTableHDU.from_columns([lambda_col,flux_col])
	
	new_hdus = fits.HDUList([primary_hdu,new_dhdu])
	
	if new_fname is not None:
		
		new_hdus.writeto(new_fname,overwrite=True)
		
		
	return new_hdus
