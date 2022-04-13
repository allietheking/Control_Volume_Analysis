
#################
### PACKAGES
#################

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import matplotlib.animation as manimation
import os,sys
import pandas as pd
import matplotlib as mpl
from scipy.interpolate import griddata
from shapely.geometry import Point
import matplotlib.pyplot as plt 

#################    
### USER INPUT
#################

# path to DWAQ run, hist file, and hist_bal file
path =  '../dwaq2019061300/'
histfn = os.path.join(path,'dwaq_hist.nc')
histbal_fn = os.path.join(path,'dwaq_hist_bal.nc') 

# path to shape functions defining control volumes and transects for DWAQ output
shpfn = '../shapefiles/Agg_exchange_lines.shp'
shpfn_poly = '../shapefiles/Agg_mod_contiguous.shp'

# names of variables to sum (in this case to get DIN) and stoichiometric multipliers
variables = ['no3', 'nh4']
multipliers = [1.00, 1.00]
#variables = ['no3', 'nh4', 'pon1', 'diat', 'zoopl_v', 'zoopl_e', 'zoopl_r']
#multipliers = [1.00, 1.00, 1.00, 0.15, 0.1818, 0.1818, 0.1818]

# start time string, for trimming data
start_time = np.datetime64('2012-10-01') 

# set filtering option: unfiltered (10-min), daily average, filtered over spring-neap cycle, or annual average
#filter_option = 'unfiltered' 
#filter_option = '3hr'
#filter_option = 'daily'
filter_option = 'spring-neap'
#filter_option = 'annual'

# set plotting increment in days (note for annual average, this won't be used)
dt_plot = 4. #3./24.

# set name of folder to store plots (must already exist) 
#plotfolder = 'DIN_Flux_3hr_Average_Streamlines'
plotfolder = 'DIN_Flux_Spring_Neap_Filtered_Streamlines'

# prefix for plots
plotname_prefix = 'din'

# plot title
#plot_title = 'DIN Flux, 3-Hour Average'
plot_title = 'DIN Flux, Filtered Over Spring-Neap Cycle'

# set a maximum radius defining the neighborhood for interpolating
# fluxes across transects to find flux vectors and also define power
# for inverse distance weithging
dmax = 6000.

# set max number of points to use in least squares approximation
Nmax = 4

# set maximum flux for coloring streamlines
qmax = 100  # good for din and tn annual average/spring-neap filtered
#qmax = 150 # use for daily averages
#qmax = 1000 # use for 3hr averages

# regular grid spanning South Bay
x = np.arange(553000.,593000.,250.)
y = np.arange(4138000.,4183000.,250.)
xg, yg = np.meshgrid(x,y)

# define extent of South Bay and Lower South Bay
x_SB = [551000.,595000.]
y_SB = [4138000.,4182000.]
x_LSB =[575800., 586300.]
y_LSB =[4144500., 4153400.]

#################
### FUNCTIONS
#################

def spring_neap_filter(y, dt_days = 1.0, low_bound=1/36.0, high_bound=1/48.0):
    low_bound = dt_days*low_bound
    high_bound = dt_days*high_bound
   
    if len(y) % 2:
        result = np.fft.rfft(y, len(y))
    else:
        result = np.fft.rfft(y)
    freq = np.fft.fftfreq(len(y))[:len(result)]
    factor = np.ones_like(result)
    factor[freq > low_bound] = 0.0

    sl = np.logical_and(high_bound < freq, freq < low_bound)

    a = factor[sl]
    # Create float array of required length and reverse
    a = np.arange(len(a) + 2).astype(float)[::-1]

    # Ramp from 1 to 0 exclusive
    a = (a/a[0])[1:-1]

    # Insert ramp into factor
    factor[sl] = a

    result = result * factor
    yf = np.fft.irfft(result, len(y))
   
    return yf
    
def Poly2Transect(left,right,pi='all'):
    # find the transects for each polygon and the signs based on the following
    # rule: transects with their 'from' segment left of the path are positive
    # otherwise negated; so left polygons are given a negative sign
    def p2t_i(i):
        indl = np.where(left==i)[0]
        indr = np.where(right==i)[0]
        signl = np.ones_like(indl)*-1
        signr = np.ones_like(indr)
        return {'transect':np.concatenate([indl,indr]),
                 'sign':np.concatenate([signl,signr])}
    if pi=='all':
        p2t = []
        for i in np.arange(pmax):
            p2t.append(p2t_i(i))
    elif isinstance(pi,int):
        p2t = p2t_i(pi)
    else:
        raise ValueError("The type of pi is not implemented")      
    return p2t

###################
### MAIN PROGRAM
###################

# load data from hist and hist_bal files
hdata = xr.open_dataset(histfn) # history file
hbdata = xr.open_dataset(histbal_fn) # history balance file

# find all transects from the output hist file
TransectBL = ['transect' in name for name in hdata.location_names.values[0]]
indT = np.where(TransectBL)[0]

# find all polygons from the output hist file 
PolygonBL = ['polygon' in name for name in hdata.location_names.values[0]]
indP = np.where(PolygonBL)[0] 

# find all polygons from the output hist blanace file
PolygonBL_bal = ['polygon' in name for name in hbdata.region.values]
indP_bal = np.where(PolygonBL_bal)[0] 

# sum variables (e.g. nh4 + no3), multiplying by their stoichiometric multipliers, for 
# both transects and polygons
varT = multipliers[0]*hdata.isel(nSegment=indT)[variables[0]].values
varP = multipliers[0]*hdata.isel(nSegment=indP)[variables[0]].values
if len(variables)>1:    
    for i in range(1,len(variables)):
        varT = varT + multipliers[i]*hdata.isel(nSegment=indT)[variables[i]].values
        varP = varP + multipliers[i]*hdata.isel(nSegment=indP)[variables[i]].values

# times
time = np.array([pd.to_datetime(t) for t in hdata.time.data])

# add averaging or filtering and trim dates before start time
ind = time>=start_time
time = time[ind]
if filter_option == 'unfiltered':
    varT = varT[ind,:]
    varP = varP[ind,:]
if filter_option == 'spring-neap':
    varTf = np.zeros(varT.shape)
    varPf = np.zeros(varP.shape)
    for i in range(len(varT[0,:])):
        varTf[:,i] = spring_neap_filter(varT[:,i], dt_days = 10/1440)
    for i in range(len(varP[0,:])):
        varPf[:,i] = spring_neap_filter(varP[:,i], dt_days = 10/1440)
    varT = varTf.copy()
    varP = varPf.copy()
    varT = varT[ind,:]
    varP = varP[ind,:] 
elif filter_option == '3hr':
    N = int(3*60/10)
    varTf = np.zeros(varT.shape)
    varPf = np.zeros(varP.shape)
    for i in range(len(varT[0,:])):
        varTf[:,i] = np.convolve(varT[:,i], np.ones((N,))/N, mode='same')
    for i in range(len(varP[0,:])):
        varPf[:,i] = np.convolve(varP[:,i], np.ones((N,))/N, mode='same')
    varT = varTf.copy() 
    varP = varPf.copy()
    varT = varT[ind,:]
    varP = varP[ind,:]
elif filter_option == 'daily':
    N = int(24*60/10)
    varTf = np.zeros(varT.shape)
    varPf = np.zeros(varP.shape)
    for i in range(len(varT[0,:])):
        varTf[:,i] = np.convolve(varT[:,i], np.ones((N,))/N, mode='same')
    for i in range(len(varP[0,:])):
        varPf[:,i] = np.convolve(varP[:,i], np.ones((N,))/N, mode='same')
    varT = varTf.copy() 
    varP = varPf.copy()
    varT = varT[ind,:]
    varP = varP[ind,:]
elif filter_option == 'annual':
    varT = varT[ind,:].mean(axis=0)
    varP = varP[ind,:].mean(axis=0)

# read control volume exchange lines and find control volume polygons to left and right
gdf = gpd.read_file(shpfn)
left = gdf.left
right = gdf.right
pmax = max(left.max(),right.max()) +1 # total number of polygons and they are zero-based

# map from polygons to transect -- not sure what this does           
p2t = Poly2Transect(left,right)    

# read control volume polygons from shapefile
poly_df = gpd.read_file(shpfn_poly)

# merge polygons into a single polygon representing the whole bay and
# extract the polygon from the geopandas dataframe as a shapely polygon
poly_allbay_gdf = poly_df.copy(deep=True)
poly_allbay_gdf['dummy']=1
poly_allbay_gdf = poly_allbay_gdf.dissolve(by='dummy')
poly_allbay = poly_allbay_gdf.iloc[0].geometry

# find centroids of the edges
xe = gdf.geometry.centroid.x.values
ye = gdf.geometry.centroid.y.values

# centroids of the polygons
poly_x = poly_df.geometry.centroid.x.values
poly_y = poly_df.geometry.centroid.y.values

# centroid of the polygon to the right of each edge
xright = poly_x[gdf.right.values]   
yright = poly_y[gdf.right.values]

# length of the line between the centroid of the edge 
# and the centroid of the polygon to the right
line_length = np.sqrt( (xe-xright)**2 + (ye-yright)**2 )

# this gives the lengths of the edges/transects/exchange lines
edge_length = gdf.geometry.length.values 

# loop through times, interpolate to centroids (using least squares approach), and plot
time0 = time[0]
while time0 < time[-1]:

    # find the time matching this time stamp
    if filter_option=='annual':
        varTmean = varT
    else:
        itime = np.argmax(time>=time0)
        varTmean = varT[itime,:]     # "mean" is a misnomer, carry over from earlier version of code
    
    # compute vector components of flux normal to each transect, 
    # normalizing the flux by edge length to get flux per unit length
    # (actually direction is not normal to the transect but points from one
    # polygon centroid to the other polygon centroid)
    qx =  (xright-xe)/line_length*varTmean/edge_length
    qy =  (yright-ye)/line_length*varTmean/edge_length
    
    # find the unit normal in the direction of the fluxes across the transects
    nx = qx/np.sqrt(qx**2 + qy**2)
    ny = qy/np.sqrt(qx**2 + qy**2)
    
    # filter out the fluxes across transects that are less than 800m long
    ind = edge_length > 800.
    nx1 = nx[ind]
    ny1 = ny[ind]
    xe1 = xe[ind]
    ye1 = ye[ind]
    qx1 = qx[ind]
    qy1 = qy[ind]
    
    # interpolate onto the regular grid using a least square inverse distance 
    # weighting method to find the flux vectors from the fluxes projected onto
    # the polygon transects
    
    # initialize flux vectors as nan at centroids
    # on regular grid
    qxg = np.nan*np.ones(np.shape(xg))
    qyg = np.nan*np.ones(np.shape(xg))
    
    # loop through the regular grid, evaluating flux vectors, using just the 5 closest 
    # points and the inverse distance cubed weighting matrix
    for i in range(len(x)):
        for j in range(len(y)):
            
            ## pick out coordinates on grid
            xc1 = x[i]
            yc1 = y[j]
            
            # check if point is inside the bay, and if it is, proceed with 
            # interpolation step
            if Point(xc1,yc1).within(poly_allbay):
            
                # find the distances between this grid point and points where fluxes
                # are defined
                d = np.sqrt((xe1-xc1)**2 + (ye1-yc1)**2)
                
                #############################################################
                # solve unweighted least square problem using only 4 closet
                # points in the neighborhood
                #############################################################
                
                # select points where d is less that 6000m 
                indn = d <= dmax
                d2 = d[indn]
                nx2 = nx1[indn]
                ny2 = ny1[indn]
                xe2 = xe1[indn]
                ye2 = ye1[indn]
                qx2 = qx1[indn]
                qy2 = qy1[indn]
                N = len(nx2)
                
                # use only the Nmax closest points
                if N > Nmax:
                    indn = np.argsort(d2)
                    indn = indn[0:Nmax]
                    N = Nmax
                    d2 = d2[indn]
                    nx2 = nx2[indn]
                    ny2 = ny2[indn]
                    xe2 = xe2[indn]
                    ye2 = ye2[indn]
                    qx2 = qx2[indn]
                    qy2 = qy2[indn]
                
                # only possible to solve matrix equation if have at least 3 points
                if N>=3:
                    A = np.zeros((N,2))
                    A[:,0] = nx2
                    A[:,1] = ny2
                    b = np.zeros(N)
                    b[:] = np.sqrt(qx2**2 + qy2**2)
                    W = np.zeros((N,N))
                    for n in range(N):
                        W[n,n] = 1/d2[n]**3
                    xhat = np.matmul(np.linalg.inv(np.matmul(np.transpose(A),np.matmul(W,A) )),
                                  np.matmul(np.transpose(A),np.matmul(W,b))) 
                    qxg[j,i] = xhat[0]
                    qyg[j,i] = xhat[1]
               
    # plot streamlines using the fine regular grid along with the fluxes 
    # across the polygon faces to make sure they look like they agree
    fig = plt.figure(figsize=(8,7))
    ax = fig.subplots()
    plt.axis((np.min(x_SB), np.max(x_SB), np.min(y_SB), np.max(y_SB)))
    norm = mpl.colors.Normalize(vmin = 0, vmax = qmax)
    plt.streamplot(xg,yg,qxg,qyg,color=np.sqrt(qxg**2+qyg**2),cmap='jet',norm=norm,density=4)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label('Flux Magnitude (g/m/d)')
    poly_allbay_gdf.plot(ax=ax,color='w',edgecolor='k')
    plt.xlabel('x (m, UTM Zone 10)')
    plt.ylabel('y (m)')
    plt.suptitle(plot_title)
    if filter_option=='annual':    # add time stamp unless this is an annual average
        fig.savefig(plotfolder + '/' + plotname_prefix + '_annual_avg.png')
    else:
        plt.title('%s' % time0)
        fig.savefig(plotfolder + '/' + plotname_prefix + '_%06d.png' % itime)
    plt.close()
    
    # advandce time counter, unless this in an annual average, in which case exit the loop
    if filter_option=='annual':
        break
    else:
        time0 = time0 + np.timedelta64(int(dt_plot*24),'h')
            
# to make into a movie, in terminal, navigate into folder with the plots and type
# ffmpeg -framerate 10 -start_number 61 -i din_%04d.png -f mp4 -vcodec h264 -pix_fmt yuv420p din_flux_streamlines_spring_neap_filtered.mp4
