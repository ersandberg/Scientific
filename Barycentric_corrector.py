# Erik Sandberg, June 2016, Project MINERVA
# Goal: Calculate the barycentric correction for every spectra taken on a given night.


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import urllib2
from BeautifulSoup import BeautifulSoup
import bs4
from astropy.time import Time
import os

##### ---------------- Using Jason Eastman's barycentric velocity code ---------------- #####
# Query OSU page for barycentric correction
# jdutc can be a float of a list of floats. pmra, pmdec, rv, and parallax can be derived from the header. assume zmeas = 0, <-- WILL CAUSE LOW PRECISION

def barycorr(jdutc,ra,dec,pmra=0.0,pmdec=0.0,parallax=0.0,rv=0.0,zmeas=0.0,
             epoch=2451545.0,tbase=0.0,
             lon=-110.878977,lat=31.680407,elevation=2316.0):

    if type(jdutc) is list: jds = jdutc
    else: jds = [jdutc]

    url="http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.php?" +\
        "JDS=" + ','.join(map(str,jds)) + "&RA=" + str(ra) + "&DEC=" + str(dec) +\
        "&LAT=" + str(lat) + "&LON=" + str(lon) + "&ELEVATION=" + str(elevation) +\
        "&PMRA=" + str(pmra) + "&PMDEC=" + str(pmdec) + "&RV=" + str(rv) +\
        "&PARALLAX=" + str(parallax) + "&ZMEAS=" + str(zmeas) +\
        "&EPOCH=" + str(epoch) + "&TBASE=" + str(tbase)

    request = urllib2.Request(url)
    response = urllib2.urlopen(request)

    data = response.read()
    soup = BeautifulSoup(data) # Previously was (data,'html.parser')
    bc = soup.text.split() # Previously was soup.pre.text.split()

    if len(bc) == 1: return float(bc[0])
    return map(float,bc)


##### ---------------- End of imports/definitions, begin work ---------------- #####


print '-> Hello, User.'
print '-> Your expmeter.dat and spectra files should be in the night folder (i.e. foo/bar/n20160422).'
print '-> Your FAU files should be in folders corresponding to the telescope (i.e. foo/bar/n20160422/T1).'

# Example: directory = '/Users/Erik/Desktop/MINERVA/Barycentric_corrector/n20160523/'
directory = raw_input('-> Please input the directory to your night folder (i.e. foo/bar/n20160422: ')
if directory[-1] != '/' : directory = directory + '/' 


##### ---------------- Load in exposure meter data --------------- #####

text = np.loadtxt(directory + 'expmeter.dat', delimiter=',', dtype='str') # Exposure meter data
date = text[:,0]
fractional_hours = np.zeros(date.size) # Time data
expmeter_flux = text[:,1] # Flux data
for i in range(date.size):
    hours = float(date[i][11:13])
    minutes = float(date[i][14:16])
    seconds = float(date[i][17:])
    fractional_hours[i] = hours+minutes/60. + seconds/(60.**2.)


##### ---------------- Read in FAU files (different number of FAU files per telescope) --------------- #####
# Each telescope has different number of fits files, so create some dummy arrays that can be filled with 0s so all arrays match.
# Have the system perform (ls * > fau_names.txt) in each telescope directory
os.system('ls ' + directory + 'T1/ > ' + directory + 'T1/fau_names.txt') # File names of FAUs
os.system('ls ' + directory + 'T2/ > ' + directory + 'T2/fau_names.txt') # File names of FAUs
os.system('ls ' + directory + 'T3/ > ' + directory + 'T3/fau_names.txt') # File names of FAUs
os.system('ls ' + directory + 'T4/ > ' + directory + 'T4/fau_names.txt') # File names of FAUs
os.system('ls ' + directory + '*fits > ' + directory + 'data_names.txt') # File names of spectra

filesT1_badshape = np.loadtxt(directory + 'T1/' + 'fau_names.txt', dtype=str)
star_fau_filesT1 = np.zeros(2000).astype('S64')
star_fau_filesT1[:filesT1_badshape.size] = filesT1_badshape


filesT2_badshape = np.loadtxt(directory + 'T2/' + 'fau_names.txt', dtype=str)
star_fau_filesT2 = np.zeros(2000).astype('S64')
star_fau_filesT2[:filesT2_badshape.size] = filesT2_badshape


filesT3_badshape = np.loadtxt(directory + 'T3/' + 'fau_names.txt', dtype=str)
star_fau_filesT3 = np.zeros(2000).astype('S64')
star_fau_filesT3[:filesT3_badshape.size] = filesT3_badshape


filesT4_badshape = np.loadtxt(directory + 'T4/' + 'fau_names.txt', dtype=str)
star_fau_filesT4 = np.zeros(2000).astype('S64')
star_fau_filesT4[:filesT4_badshape.size] = filesT4_badshape

star_fau_files = np.vstack((star_fau_filesT1,star_fau_filesT2,star_fau_filesT3,star_fau_filesT4)) # Shape is 4 x 2000


#####----------------- PERFORM APERTURE PHOTOMETRY --------------- #####

xpix = np.arange(648) # np.arange(hd_FAU['NAXIS1'])
ypix = np.arange(486) # np.arange(hd_FAU['NAXIS2']) 
x,y = np.meshgrid(xpix,ypix)
aperture_radius = 10. # What is the fiber radius in pixels? That's what should go in here.
annulus_outer = 3.0*aperture_radius
annulus_inner = 2.0*aperture_radius
telescope_ID = np.array(['T1/','T2/','T3/','T4/'])
aperture_flux = np.zeros((star_fau_filesT1.size,telescope_ID.size)) # One per fau image, 4 sets, one for each telescope
fau_times = np.zeros((star_fau_filesT1.size,telescope_ID.size))
fau_jd = np.zeros((star_fau_filesT1.size,telescope_ID.size))
errors = []

print 'Currently performing aperture photometry on fibers, this may take several minutes'
for telescope in range(telescope_ID.size): # Run through every telescope
    if telescope==2: print '--> Halfway completed with aperture photometry.'
    for index,filename in enumerate(star_fau_files[telescope,:]): # Run through every spectrum
        try:
            star_FAU, hd_FAU = fits.getdata(directory + telescope_ID[telescope] + filename,header=True)
            fiber_x = hd_FAU['XFIBER'+ str(telescope+1)] # Fiber X position
            fiber_y = hd_FAU['YFIBER'+ str(telescope+1)] # Fiber Y position
            distance_from_fiber = ((x-fiber_x)**2+(y-fiber_y)**2)**0.5
            annulus_target = (distance_from_fiber<annulus_outer) & (distance_from_fiber>annulus_inner)
            pixels_in_aperture = (distance_from_fiber<(aperture_radius)).sum() # Total number of pixels in aperture
            sky_counts = np.median(star_FAU[annulus_target]) # Average flux from the background
            aperture_flux[index,telescope] = np.sum(star_FAU[(distance_from_fiber<aperture_radius)])-(sky_counts*pixels_in_aperture)
            hours = float(hd_FAU['DATE-OBS'][11:13]) # Hours
            minutes = float(hd_FAU['DATE-OBS'][14:16]) # Minutes
            seconds = float(hd_FAU['DATE-OBS'][17:]) # Seconds
            fau_times[index,telescope] = hours+minutes/60. + seconds/(60.**2) # Fractional hours, UTC
            fau_jd[index,telescope] = Time(hd_FAU['DATE-OBS'],format='isot',scale='utc').jd # Julian Date
        except Exception,e: 
            #print str(e) # There will always be errors associated with being unable to find 'directory/0.0' due to weird sizings. That error looks bad, but is NOT A PROBLEM.
            aperture_flux[index,telescope] = 0 # If there is an error, that flux is then = 0.
print 'Completed aperture photometry on fibers.'


##### ---------------- Get the starting and ending times of every spectra --------------- #####

data_names = np.loadtxt(directory + 'data_names.txt', delimiter=',', dtype='str') #Spectra data filenames
image_start = np.zeros(data_names.size) # UTC
image_end = np.zeros(data_names.size) # UTC
jd_start = np.zeros(data_names.size) # Julian date
jd_end = np.zeros(data_names.size) # Julian date
for index,spectrum_file in enumerate(data_names):
    try:
        image,hd = fits.getdata(spectrum_file, header=True) # Star spectra, just using to get the correct time of science images
        hours = float(hd['DATE-OBS'][11:13]) # Hours
        minutes = float(hd['DATE-OBS'][14:16]) # Minutes
        seconds = float(hd['DATE-OBS'][17:]) # Seconds
        exposure_time = float(hd['EXPTIME']) # Seconds
        jd_start[index] = hd['JD']
        jd_end[index] = hd['JD'] + exposure_time/86400. #86400 seconds per day
        image_start[index] = hours+minutes/60. + seconds/(60.**2.)
        image_end[index] = image_start[index]+(exposure_time)/(60.**2.)
    except Exception, e:
        print index, spectrum_file, str(e)




##### ---------------- Calculate BARYCENTRIC CORRECTION ---------------- #####
#times should be at every FAU exposure during one of the spectra

weighted_mean = np.zeros((data_names.size,telescope_ID.size)) # Number of spectra by number of telescopes
time_and_flux = np.zeros((data_names.size,250,2,telescope_ID.size)) # 250 is an arbitrary number. Each telescope has a different number of FAU exposures, so this array is a bit messy. This is the data that is attached as an extension to the spectra fits files.

for index, spectrum_file in enumerate(data_names):
    for telescope in range(telescope_ID.size):
        try:
            indices = np.logical_and(fau_jd[:,telescope] > jd_start[index] , fau_jd[:,telescope] < jd_end[index]) # Find the indices where the FAU times match up with when spectra is being taken
            bc_times = fau_jd[:,telescope][indices] # Get the FAU times at those indices
            bc_flux = aperture_flux[:,telescope][indices] # Get the FAU aperture flux at those indices
            time_and_flux[index,:sum(indices),0,telescope] = bc_times # Place the time stamps in the correct place (fits extension)
            time_and_flux[index,:sum(indices),1,telescope] = bc_flux  # Place the flux measurements in the correct place (fits extension)            
            barycentric_corr = barycorr(list(bc_times),hd['targra1'],hd['targdec1'],hd['pmra1'],hd['pmdec1'],hd['rv1']) # Calculate the barycentric correction at each of those times (i.e. indices)
            weighted_mean[index,telescope] = sum(barycentric_corr*bc_flux)/sum(bc_flux)
            print '--> Found barycentric correction for', spectrum_file,' T'+str(telescope+1)
            
        except Exception, e: # This occurs if there is NO FAU DATA during our spectrum
            print '--> Error: T' + str(telescope+1),' has failed at', spectrum_file, '. Most likely due to weather. Using the middle of the exposure to find a barycentric correction. '            
            time_midpoint = (jd_start[index]+jd_end[index])/2. # Time midpoint of spectra
            barycentric_corr = barycorr(time_midpoint,hd['targra1'],hd['targdec1'],hd['pmra1'],hd['pmdec1'],hd['rv1'])
            weighted_mean[index,telescope] = barycentric_corr # Not actually a "weighted" mean in this case. No FAU data means we can't get an actual weighted mean.
            


time_and_flux[np.where(time_and_flux==0)] = np.nan # Convert 0s to nans in 'time and flux' array to lessen confusion with extra numbers

##### ---------------- Append BARYCENTRIC CORRECTION to spectra ---------------- #####

# Attach barycentric correction to primary header (image[0].header), then attach Nx2x4 array (per spectrum) as a seconday header
for index,spectrum_file in enumerate(data_names):
   image = fits.open(spectrum_file)
   image[0].header['BARYCOR1'] = weighted_mean[index,0], 'T1 Flux weighted barycentric correction (m/s)'
   image[0].header['BARYCOR2'] = weighted_mean[index,1], 'T2 Flux weighted barycentric correction (m/s)'
   image[0].header['BARYCOR3'] = weighted_mean[index,2], 'T3 Flux weighted barycentric correction (m/s)'
   image[0].header['BARYCOR4'] = weighted_mean[index,3], 'T4 Flux weighted barycentric correction (m/s)'
   # Now for the Nx2x4 array of aperture photometries on the fibers.
   hdu = fits.PrimaryHDU(time_and_flux[index,:,:,:]) # Attach only ONE set of data per spectra
   image.append(hdu) # Append the extension ('hdu') to the spectra fits file ('image')
   image.writeto(spectrum_file, clobber=True) # Rewrite the fits file with new header
   print 'Appended aperture photometries and barycentric corrections to', spectrum_file, 'header'
print '--> Barycentric_corrector.py has been completed. '

