import numpy as np
import os,sys,glob
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import rasterio,gdal
from pyproj import Proj, transform
import subprocess


def get_clim(ar):
    try:
        clim = np.percentile(ar.compressed(),(2,98))
    except:
        clim = np.percentile(ar,(2,98))
    return clim
def find_common_clim(im1,im2):
    perc1 = get_clim(im1)
    perc2 = get_clim(im2)
    perc = (np.min([perc1[0],perc2[0]]),np.max([perc1[1],perc2[1]]))
    abs_max = np.max(np.abs(perc))
    perc = (-abs_max,abs_max)
    return perc
def fn_2_ma(fn,b=1):
    ds = rasterio.open(fn)
    ar = ds.read(b)
    ndv = get_ndv(ds)
    ma_ar = np.ma.masked_less_equal(ar,ndv)
    return ma_ar
def get_ndv(ds):
    no_data = ds.nodatavals[0]
    if no_data == None:
        #this means no data is not set in tif tag, nead to cheat it from raster
        ndv = ds.read(1)[0,0]
    else:
        ndv = no_data
    return ndv

def plot_stereo_results(out_folder,ax):
    dem = glob.glob(out_folder+'*-DEM.tif')[0]
    l_img = glob.glob(out_folder+'*L.tif')[0]
    r_img = glob.glob(out_folder+'*R.tif')[0]
    disp = glob.glob(out_folder+'*F.tif')[0]
    error_fn = glob.glob(out_folder+'*In*.tif')[0]
    print("Found files {}\n {}\n {}\n {}\n {}\n".format(l_img,r_img,disp,error_fn,dem))
    dem_ma = fn_2_ma(dem,1)
    l_img_ma = fn_2_ma(l_img)
    r_img_ma = fn_2_ma(r_img)
    dx_ma = fn_2_ma(disp,1)
    dy_ma = fn_2_ma(disp,2)
    error_ma = fn_2_ma(error_fn,1)
    dem_ds = gdal.Open(dem)
    producttype = 'hillshade'
    hs_ds = gdal.DEMProcessing('',dem_ds,producttype,format='MEM')
    hs = hs_ds.ReadAsArray()
    #fig,ax = plt.subplots(3,2,figsize=(6,5))
    axa = ax.ravel()
    cmap_img = 'gray'
    cmap_disp = 'RdBu'
    cmap_error = 'inferno'
    plot_ar(l_img_ma,ax=axa[0],clim=get_clim(l_img_ma),cbar=False,cmap=cmap_img)
    print(get_clim(l_img_ma))
    plot_ar(l_img_ma,ax=axa[1],clim=get_clim(r_img_ma),cbar=False,cmap=cmap_img)
    disp_clim = find_common_clim(dx_ma,dy_ma)
    plot_ar(dx_ma,ax=axa[2],clim=disp_clim,cmap=cmap_disp,label='dx (px)')
    plot_ar(dy_ma,ax=axa[3],clim=disp_clim,cmap=cmap_disp,label='dy (px)')
    #plt.colorbar(dx_im,axa[2])
    #plt.colorbar(dy_im,axa[3])
    error_im = axa[4].imshow(error_ma,cmap=cmap_error,clim=get_clim(error_ma),interpolation='none')
    plot_ar(error_ma,ax=axa[4],clim=get_clim(error_ma),cmap=cmap_error,label='Intersection Error(m)')
    #plt.colorbar(error_im,ax=axa[4])
    axa[5].imshow(hs,cmap=cmap_img,clim=get_clim(hs),interpolation='none')
    plot_ar(dem_ma,ax=axa[5],clim=get_clim(dem_ma),label='HAE (m WGS84)',alpha=0.6)
    print(get_clim(dem_ma))
    #plt.colorbar(dem_im,axa[5])
    plt.tight_layout()
    print(type(dem_ma))
    
def plot_ar(im,ax,clim,cmap=None,label=None,cbar=True,alpha=1):
    if cmap:
        img = ax.imshow(im,cmap=cmap,clim=clim,alpha=alpha,interpolation='none')
    else:
        img = ax.imshow(im,clim=clim,alpha=alpha,interpolation='none')
    if cbar:
        divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #cax = divider.append_axes("right", size="5%", pad="2%")
        cax = divider.append_axes("right", size="2%", pad="1%")
        cb = plt.colorbar(img,cax=cax,ax=ax)
        cax.set_ylabel(label)
    ax.set_xticks([])
    ax.set_yticks([])
    
def subsetBBox(rast,proj_out):
    #originally written by Justing Pflug 
    # rasterio open data
    rB = rasterio.open(rast)
    proj_in = rB.crs
    # rasterio get bounding box
    [L,B,R,T] = rB.bounds

    if proj_in == proj_out:
        return L, R, T, B
    else:
        incord = Proj(init=proj_in)
        outcord = Proj(init=proj_out)

        [Left,Bottom] = transform(incord,outcord,L,B)
        [Right,Top] = transform(incord,outcord,R,T)
        return Left, Bottom, Right, Top
    
def download_srtm(raster,outfn):
    ### this is old
    import elevation
    minx,miny,maxx,maxy = rasterio.open(raster).bounds
    in_crs = rasterio.open(raster).crs
    minx,miny,maxx,maxy = subsetBBox(dem_file,'epsg:4326')
    elevation.clip(bounds=(minx, miny, maxx, maxy), output=outfn)
    #Clean cache
    elevation.clean()
    
    
def get_dem(demtype, bounds, api_key, out_fn=None, proj='EPSG:4326'):
    """
    download a DEM of choice from OpenTopography World DEM
    ## first written by David Shean
    Parameters
    ------------
    demtype: str
        type of DEM to fetch (e.g., SRTMGL1, SRTMGL1_E, SRTMGL3 etc)
    bounds: list
        geographic aoi extent in format (minx,miny,maxx,maxy)
    out_fn: str
        path to output filename
    t_srs: str
        output DEM projection
    Returns
    -----------
    out_DEM: str
        path to output DEM (useful if the downloaded DEM is reprojected to custom proj)
    """
    import requests
    from distutils.spawn import find_executable
    ### From David Shean
    base_url="https://portal.opentopography.org/API/globaldem?demtype={}&west={}&south={}&east={}&north={}&outputFormat=GTiff&API_Key={}"
    if out_fn is None:
        out_fn = '{}.tif'.format(demtype)
    if not os.path.exists(out_fn):
        #Prepare API request url
        #Bounds should be [minlon, minlat, maxlon, maxlat]
        url = base_url.format(demtype, *bounds, api_key)
        print(url)
        #Get
        response = requests.get(url)
        #Check for 200
        #Write to disk
        open(out_fn, 'wb').write(response.content)
    if proj != 'EPSG:4326':
        #Could avoid writing to disk and direclty reproject with rasterio, using gdalwarp for simplicity
        proj_fn = os.path.splitext(out_fn)[0]+'_proj.tif'
        if not os.path.exists(proj_fn):
            output_res = 30
            gdalwarp = find_executable('gdalwarp')
            gdalwarp_call = f"{gdalwarp} -r cubic -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -tr {output_res} {output_res} -t_srs '{proj}' {out_fn} {proj_fn}"
            print(gdalwarp_call)
            run_bash_command(gdalwarp_call)
        out_DEM = proj_fn
    else:
        out_DEM = out_fn
    return out_DEM

def run_bash_command(cmd):
    #written by Scott Henderson
    # move to asp_binder_utils
    """Call a system command through the subprocess python module."""
    print(cmd)
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal", -retcode, file=sys.stderr)
        else:
            print("Child returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

