from pathlib import Path
import argparse
import xarray as xr
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from metpy.interpolate import interpolate_to_grid
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Local functions:
#   - NS_format_func - format latitude axis labels
#   - EW_format_func - format longitude (E) axis labels
#

# 
# Latitude plotting fuction to get N/S labels right:
# see https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
#

def NS_format_func(value, tick_number):
    
    # add north/south latitude labels

    N = int(value)
    if N < 0:
        if N == -90:
            return ""
        else:
            return r"${0}S$".format(-N)
    elif N == 0:
        return "0"
    elif N > 0 :
        return r"${0}N$".format(N)
    else:
        return "oops!"
# 
# Latitude plotting fuction to get N/S labels right:
# see https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
#

def EW_format_func(value, tick_number):
    
    # add north/south latitude labels

    E = int(value)
    if E < 0:
        if E == -180:
            return "180"
        else:
            return r"${0}W$".format(-E)
    elif E == 0:
        return "0"
    elif E > 0 :
        return r"${0}E$".format(E)
    else:
        return "oops!"



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'fkessler_verify',
                    description = 'Produces simulated day 10 Kessler surface pressure and large scale precipitation plots for CAM-MPAS',
                    epilog = 'Comparison plots can be found in the  CESM documentation here: https://www2.cesm.ucar.edu/models/simpler-models/fkessler/index.html')
    parser.add_argument('filepath', help='path to FKESLER history file(s)')           # positional argument
    parser.add_argument('--saveplt',
                        action='store_true', help='save the plot to a file (default: display plot)')  # on/off flag
    parser.add_argument('--gblplt',
                        action='store_true', help='plot global field (default: plot prescribed region - (20:70N, 150W:30E)') # on/off flag
    parser.add_argument('--verbose',
                        action='store_true', help='enable verbose output') # on/off flag    
    args = parser.parse_args()
    print(args.filepath, args.saveplt, args.gblplt, args.verbose)
    
    print("Entering FKESSLER code.")
    # Revisit these inputs and control variables

#    gbl_plot = False  # Set the region to display
#    saveplt = False # Save plot (True) or display plot (False)?
    timname = "Time"
    samples_per_day = 1 # if time were a coordinate variable, this could be discerned directly
    base_date = "0001-01-01"

    print("Parameters set. Move to IO.")
    
    # Data loading and workflow (could make this CLI w/ argparse)
 #   hpath = Path("/Users/loft/Desktop/GLOBUS-MBP/FKESSLER/")
    hpath = Path(args.filepath)
    hfils = sorted(hpath.glob("*.h0.*"))

    # note: using combine/concat_dim wouldn't usually be necessary if the time coordinate were correct.
    if args.verbose:
        print(f"Found a total of {len(hfils)} files.")
    ds = xr.open_mfdataset(hfils, combine='nested', concat_dim=timname)
    # ds = xr.open_dataset(hfils[0]).load()
    if args.verbose:
        print("ds loaded (probably via dask)")
    PS = ds.PS
    PRECL=ds.PRECL

    Cell_y = ds.lat
    Cell_lat = Cell_y.squeeze()
    if args.verbose:
        print("\n\n======================Cell_lat====================\n")
        print(Cell_lat)

    Cell_x = ds.lon
#    Cell_lon = np.asarray(Cell_x).squeeze()
    Cell_lon = Cell_x.squeeze()
    if args.verbose:
        print("\n\n======================Cell_lon====================\n")
        print(Cell_lon)

    ps = PS.squeeze()
    precl=PRECL.squeeze()

    #####################
    # Get hPa on day 10. (look at timeslice 10*samples_per_day)
    # .01 here converts from Pa to hPa
    #####################

    ps10 = .01*ps[10*samples_per_day, :] # Convert from Pa to hPa
    precl10 = 8.64e7*precl[10*samples_per_day, :] #Convert from m/s to mm/day

    if args.verbose:
        print("\n\n======================ps====================\n")
        print(ps10)
    ncells=ps10.shape[0]

    if args.verbose:
        print("\n\n======================precl====================\n")
        print(precl10)
    ncells=precl10.shape[0]
    
    #####################
    # Set the resolution string according to what we find for the shape (number of cells) in the ps10 array.
    # Assumes quasi-uniform icosahedral mesh.
    #####################
    
    Gicos = int(math.log2((ncells-2)/10))
    if int(10.*pow(2,Gicos)+2) != ncells:
        print("ERROR: ncells ",ncells, "is an unrecognized mesh")
        exit()
        
    resolution = str(int(120/(pow(2,12-Gicos))))
    print("icosahedral resolution", resolution)
    print("\n ncells= ",ncells,"resolution=",resolution,"km\n")


    # define lat-lon mesh

    #
    # Define lat-lon mesh spacing
    # NB: hard wired to equally spaced mesh
    #
    
    dlon = 0.5
    dlat = dlon

    #
    # Define region to extract from image (img)
    #
    
    if args.gblplt:
        lon_min = -180.
        lon_max = 180.
        lat_min = -90.
        lat_max = 90.
    else:
        lon_min = -150.
        lon_max = 30.
        lat_min = 20.
        lat_max = 70.
    
    nlon = int((lon_max-lon_min)/dlon)
    nlat = int((lat_max-lat_min)/dlat)

    #
    # indices in the img where our region sits
    #
    
    lon_min_indx = -int((-180.-lon_min)/dlon)
    lon_max_indx = -int((-180.-lon_max)/dlon)
    lat_min_indx = -int((-90.-lat_min)/dlat)
    lat_max_indx = -int((-90.-lat_max)/dlat)
    
    lon = np.linspace(lon_min, lon_max, nlon)
    lat = np.linspace(lat_min, lat_max, nlat)
    
    gx, gy, ps10_img = interpolate_to_grid(np.asarray(Cell_lon), np.asarray(Cell_lat), np.asarray(ps10), interp_type='linear', hres=0.5, boundary_coords=None)

    gx, gy, precl10_img = interpolate_to_grid(np.asarray(Cell_lon), np.asarray(Cell_lat), np.asarray(precl10), interp_type='linear', hres=0.5, boundary_coords=None)
    
    
    ps10_img_reg = ps10_img[lat_min_indx:lat_max_indx, lon_min_indx:lon_max_indx]
    ps10_reg = xr.DataArray(
        data   = ps10_img_reg,
        dims   = ['latitude','longitude'],
        coords = {'latitude': lat, 'longitude': lon},
        attrs  = {'long_name': "Surface Pressure", 'units': "hPa"}
    )
    
    precl10_img_reg = precl10_img[lat_min_indx:lat_max_indx, lon_min_indx:lon_max_indx]
    precl10_reg = xr.DataArray(
        data   = precl10_img_reg,
        dims   = ['latitude','longitude'],
        coords = {'latitude': lat, 'longitude': lon},
        attrs  = {'long_name': "Large-scale precipitation rate", 'units': "mm/day"}
    )


    print("\n======================================= BEGIN PLOTTING  ==============================================\n")

    # Generate figure (set its size (width, height) in inches)

    fig = plt.figure(figsize=(8,6))

    gs = fig.add_gridspec(nrows=2, ncols = 1, hspace=0.5)

    axs = gs.subplots()

    #
    # Define font sizes for various parts of the plots:
    #

    fs_title=12
    fs_axlab=8
    fs_major=6
    fs_minor=4

    suptitle = "Day 10 (-compset FKESSLER -res " + resolution + "km)"
    fig.suptitle(suptitle, ha='center', fontsize=fs_title, fontweight="bold")
    
    axs[0].set_title("Surface pressure", loc='left', fontdict={'fontsize':fs_axlab})
    axs[0].set_title("hPa", loc='right', fontdict={'fontsize':fs_axlab})
    axs[1].set_title("Large-scale (stable) precipitation rate (liq + ice)", loc='left', fontdict={'fontsize':fs_axlab})
    axs[1].set_title("mm/day", loc='right', fontdict={'fontsize':fs_axlab})

    #
    # Generating surface pressure plot 
    #
    
    print("\nplotting PS\n")
    ax = axs[0]
    
    # x axis: set bottom major ticks with labels and minor ticks without labels
    ax.xaxis.set_major_locator(mticker.MultipleLocator(30))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(10))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(EW_format_func))
    ax.tick_params(which='both',right=True, top=True, labelright=False, labeltop=False, labelrotation=0)
    ax.tick_params(axis='both', which='major', labelsize=fs_major)
    
    # y axis
    ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(NS_format_func))
    # For the minor ticks, use no labels; default NullFormatter.
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(5))

    CS = ax.contourf(lon, lat, ps10_reg, 16, cmap="jet", alpha=0.8)
    
    axins1 = inset_axes(
        ax,
        width="70%",  # width: 70% of parent_bbox width works
        height="12%",  # height: 10% works
        loc="lower center",
        bbox_to_anchor=(.1, -0.12, .8, -0.80),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    #axins1.xaxis.set_ticks_position("bottom")
    cbar=fig.colorbar(CS, cax=axins1, orientation="horizontal")
    cbar.ax.tick_params(bottom = False, pad=16, labelsize=fs_major, grid_color='b', grid_alpha=1, grid_linewidth=2)
    cbar.ax.locator_params(nbins=16)
    
    print("\nplotting PRECL\n")

    #
    # Generating precip plot 
    #

    ax = axs[1]
    
    # x axis: set bottom major ticks with labels and minor ticks without labels
    ax.xaxis.set_major_locator(mticker.MultipleLocator(30))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(10))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(EW_format_func))
    ax.tick_params(which='both',right=True, top=True, labelright=False, labeltop=False, labelrotation=0)
    ax.tick_params(axis='both', which='major', labelsize=fs_major)
    
    # y axis
    ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(NS_format_func))
    # For the minor ticks, use no labels; default NullFormatter.
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(5))

    CS = ax.contourf(lon, lat, precl10_reg, levels=[0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64], cmap="jet", alpha=0.8)
    
    axins1 = inset_axes(
        ax,
        width="70%",  # width: 70% of parent_bbox width works
        height="12%",  # height: 10% works
        loc="lower center",
        bbox_to_anchor=(.1, -0.12, .8, -0.80),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    #axins1.xaxis.set_ticks_position("bottom")
    cbar=fig.colorbar(CS, cax=axins1, orientation="horizontal")
    cbar.ax.tick_params(bottom = False, pad=16, labelsize=fs_major, grid_color='b', grid_alpha=1, grid_linewidth=2)
    cbar.ax.locator_params(nbins=8)


    if args.saveplt:
        ofil = Path.home() / "fkessler_diag_plot.png"
        fig.savefig(ofil, bbox_inches='tight')
        print(f"Output saved to: {ofil}")
    else:
        plt.show()
    print('ALL DONE.')
