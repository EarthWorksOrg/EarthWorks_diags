from pathlib import Path

import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from metpy.interpolate import log_interpolate_1d

# Local functions:
#   - calc_xpyp : calculate time&zonal average eddy flux
#   - contour_zm : make the lat-height plot on an axis
#   - approx_pres : approximate pressure based on temperature
#   - construct_plev : makes the pressure grid to interpolate to
#   - convert_to_xr : utility to get back to plain Xarray object
#
def calc_xpyp(x, y=None, tname='time', lonname='lon'):
    """Calculate time mean, zonal mean covariance, <x'y'>.
    If y is not provided, calculate <x'x'>

    Average over time and/or lon if they are dimensions.
    * note: if neither time nor lon are dimensions, the
            result is xprime*yprime without averaging.
    * note 2: that shouldn't actually happen because the
              prime is a deviation from longitude;
              IF lon is not in x or y: ASSUME IT IS A DEVIATION ALREADY.
    Singleton dimensions are "squeezed".
    """
    xprime = x - x.mean(dim=lonname, keep_attrs=True).squeeze() if lonname in x.dims else x
    if y is None:
        y = x
        yprime = xprime
    else:
        yprime = y - y.mean(dim=lonname).squeeze() if lonname in y.dims else y

    avgdim = [d for d in x.dims if d in [tname, lonname]]
    return (xprime * yprime).mean(dim=avgdim, keep_attrs=True).squeeze()


def contour_zm(ax, data, ci, cmin, cmax, latname='lat', levname='lev'):
    """
    ax: an Axes object to use for the plot
    data: the zonal-averaged data, needs to be lev by lat
    ci: contour interval
    cmin: minimum contour
    cmax: maximum contour

    latname / lonname: strings to specify alterate dimension names
    """
    flevels = np.arange(cmin, cmax+ci, ci)
    colormap = 'PuOr_r'
    cnorm = mpl.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)
    if isinstance(data, xr.DataArray): 
        if (latname in data.coords) and (levname in data.coords):
            print("Identified lat and lev in data.")
            mlat, mlev = np.meshgrid(data[latname], data[levname])
            levaxlabel = getattr(data[levname], 'long_name', 'lev')
            lataxlabel = getattr(data[latname], 'long_name', 'lat')
            data = data.transpose(levname,latname)
        else:
            print("Did not know which coordinates to use, will try our best.")
            print(data.coords)
            dims = data.dims
            mlat, mlev = np.meshgrid(data[dims[1]], data[dims[0]])
            levaxlabel = data.dims[0]
            lataxlabel = data.dims[1]
            titleleft = ""
            titleright = ""
    else:
        raise ValueError("plotting function really needs xarray DataArray to work ok.")
    titleleft = getattr(data, 'long_name', '')
    titleright = getattr(data, 'units', '')
    IMG = ax.contourf(mlat, mlev, data, cmap=colormap, levels=flevels, norm=cnorm)
    CS = ax.contour(mlat, mlev, data, levels=flevels[::2], colors='k') #, linewidth=0.5)
    ax.clabel(CS, fontsize=9, inline=1)
    ax.set_ylabel(levaxlabel)
    ax.set_xlabel(lataxlabel)
    ax.set_title(titleleft, loc='left')
    ax.set_title(titleright, loc='right')
    ax.invert_yaxis()
    f = plt.gcf()
    f.colorbar(IMG, ax=ax, orientation='horizontal')


def approx_pres(density, temperature):
    Rd = 287.0             # Gas constant in (J kg^-1 K^-1))
    P0 = 100000.0          # Reference pressure (kg m^-1 s^-2)
    Cmb = 1000.0/P0        # Conversion from MKS to millibars
    Pt = Cmb*Rd*density*temperature  # pressure via equation of state
    Pt.name = "pressure"
    Pt.attrs['units'] = 'hPa'
    Pt.attrs['comment'] = 'pressure in hPa via equation of state'
    Pt.attrs["mdims"] = "1"
    # Pt.attrs["units"] = "J/m^3"
    Pt.attrs["long_name"] = "air pressure"
    return Pt


def construct_plev(nlev, verbose=False):
    """
    Set pressure levels: note these assume CAM ordering
    i.e. top of atmosphere (TOA) is level 1.
    """
    dsigma = 1.0/nlev
    lev = np.arange(0, nlev, 1) # CAM top of atmosphere lies on level 1.                                                                                                                 
    sigma_mid = dsigma*(lev+0.5)
    Psig = 1000.*sigma_mid
    Psig = xr.DataArray(Psig, dims="lev", name="lev", attrs={"units":"hPa"})
    if verbose:
        print("pressure on Sigma midpoints:")
        print(Psig)
    return Psig


def convert_to_xr(mpver, orig, intrpdim, intrpcoord):
    """
    Convert the MetPy object back to Xarray DataArray.
    Assumes we need to replace one coordinate variable because of interpolation.

    mpver : the MetPy output from interpolation
    orig : the original DataArray
    intrpdim : str, name of the dimension that was interpolated
    intrpcoord : the DataArray for interpolated coordinate
    """
    ocoords = orig.coords
    ocoords[intrpdim] = intrpcoord
    fshape = mpver.shape
    oshape = orig.shape
    assert len(fshape) == len(oshape), f'Ranks are different, compare: {fshape} to {oshape}'
    return xr.DataArray(mpver, dims=orig.dims, coords=ocoords, attrs=orig.attrs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'hs_diag',
        description = 'Produces zonally and temporally averaged plots of zonal wind, eddy temperature variance, Northward eddy momentum flux, and Northward eddy heat flux, as prescribed by Held & Suarez, 1994 tests for CAM-MPAS',
        epilog = 'Comparison plots can be found in the  CESM documentation here: https://www2.cesm.ucar.edu/models/simpler-models/held-suarez.html')
    parser.add_argument('filepath', help='path to FHS94 history file(s)')           # positional argument
    parser.add_argument('--saveplt',
                        action='store_true', help='save the plot to a file (default: display plot)')  # on/off flag
    parser.add_argument('--verbose',
                        action='store_true', help='enable verbose output') # on/off flag    
    args = parser.parse_args()
    if args.verbose:
        print(args.filepath, args.saveplt, args.verbose)

    print("Entering FHS94 analysis code.")
    interp_to_pressure = True  # whether to try to interpolate to pressure levels
    interp_early = False  # if F: try to interpolate only at the very end (smaller interpolation calc)
    time_sample = False  # TODO: allow specified temporal sampling
#    save_plot = True
    latname = "latitude"
    lonname = "longitude"
    levname = "lev"
    timname = "Time"
    samples_per_day = 4 # if time were a coordinate variable, this could be discerned directly
    base_date = "0001-01-01"

    print("Parameters set. Move to IO.")
    
    # Data loading and workflow (could make this CLI w/ argparse)
    # Some paths used for internal testing:
    #hpath = Path("/Users/loft/Desktop/Globus-MBP/")
    #hpath = Path("/glade/scratch/gdicker/val.FHS94.mpasa120.che.gnu/run/convertedOutputs_latlon")
    hpath = Path(args.filepath)
    hfils = sorted(hpath.glob("latlon_val.FHS94.mpasa120.che.gnu.cam.h1.0002-01*"))

    # hfils = [Path("/Users/brianpm/Documents/model_output/mpasa_hs94/latlon_val.FHS94.mpasa120.che.gnu.cam.h1.0001-01-01-00000.nc"),]

    # note: using combine/concat_dim wouldn't usually be necessary if the time coordinate were correct.
    print(f"Found a total of {len(hfils)} files.")
    ds = xr.open_mfdataset(hfils, combine='nested', concat_dim=timname)
    # ds = xr.open_dataset(hfils[0]).load()
    print("ds loaded (probably via dask)")

    # If there's not a proper time coordinate, make one:
    if timname not in ds.coords:
        print("Need to construct a time coordinate; will start at beginning of Year 0001")
        cft=xr.cftime_range(start="0001", periods=ds.dims[timname], freq=f"{samples_per_day}H", calendar="365_day")
        ds = ds.assign_coords({timname:cft})
    if timname != "time":
        ds = ds.rename({timname:"time"})
        print(f"renamed {timname} to time to avoid breaking xarray assumptions")
    
    if interp_to_pressure and interp_early:
        interp_axis = ds.T.dims.index(levname)
        plev = construct_plev(ds.dims[levname])
        pres = approx_pres(ds.rho, ds.T)
        # MetPy interpolation (slow)
        from dask.array.core import map_blocks
        print("interpolate T")
        # T = log_interpolate_1d(plev, pres, ds.T, axis=1)
        T = map_blocks(log_interpolate_1d, plev, pres, ds.T, dtype=ds.T.dtype, axis=interp_axis).compute()
        print("interpolate U")
        # U = log_interpolate_1d(plev, pres, ds.T, axis=1)
        U = map_blocks(log_interpolate_1d, plev, pres, ds.U, dtype=ds.U.dtype, axis=interp_axis).compute()
        print("interpolate V")
        # V = log_interpolate_1d(plev, pres, ds.V, axis=1)
        V = map_blocks(log_interpolate_1d, plev, pres, ds.V, dtype=ds.V.dtype, axis=interp_axis).compute()
        # I don't understand the MetPy object, revert to xarray:
        if (not isinstance(T, np.ndarray)) and (not isinstance(T, xr.DataArray)):
            T = convert_to_xr(T.m, ds.T, 'lev', plev)
            U = convert_to_xr(U.m, ds.U, 'lev', plev)
            V = convert_to_xr(V.m, ds.V, 'lev', plev)
    else:
        T = ds.T
        U = ds.U
        V = ds.V

    # Specify the time indices to use:
    #
    if time_sample:
        raise NotImplementedError("I have not gotten to the time sampling part.")
    else:
        print(f"The time sampling assumed to be 4/day, there are {len(T['time'])} time samples in the data.")

    # Calculations (assumes xarray used for I/O)
    #
    umean = U.mean(dim=("time",lonname), keep_attrs=True)
    vptpclim = calc_xpyp(V, T, lonname=lonname)
    upvpclim = calc_xpyp(U, V, lonname=lonname)
    tptpclim = calc_xpyp(T, lonname=lonname)

    if hasattr(umean,"compute"):
        umean = umean.compute()
    if hasattr(vptpclim, "compute"):
        vptpclim = vptpclim.compute()
    if hasattr(upvpclim, "compute"):
        upvpclim = upvpclim.compute()
    if hasattr(tptpclim, "compute"):
        tptpclim = tptpclim.compute()

    if args.verbose:
        print(f"{umean.shape = }")
        print(f"{vptpclim.shape = }")
        print(f"{upvpclim.shape = }")
        print(f"{tptpclim.shape = }")

    if interp_to_pressure and not interp_early:
        interp_axis = umean.dims.index(levname)
        print(f"{interp_axis = }, nlev = {ds.dims[levname]}")  # have to go back to dataset, not use umean
        plev = construct_plev(ds.dims[levname])
        pres = approx_pres(ds.rho.mean(dim=('time', lonname)).compute(), ds.T.mean(dim=('time',lonname)).compute())
        tmpu, tmpvt, tmpuv, tmptt = log_interpolate_1d(plev, pres, umean, vptpclim, upvpclim, tptpclim, axis=interp_axis)
        umean = tmpu
        # tmp = log_interpolate_1d(plev, pres, vptpclim, dtype=vptpclim.dtype, axis=interp_axis).compute()
        vptpclim = tmpvt
        # tmp = log_interpolate_1d(plev, pres, upvpclim, dtype=upvpclim.dtype, axis=interp_axis).compute()
        upvpclim = tmpuv
        # tmp = log_interpolate_1d(plev, pres, tptpclim, dtype=tptpclim.dtype, axis=interp_axis).compute()
        tptpclim = tmptt
        if (not isinstance(umean, np.ndarray)) and (not isinstance(umean, xr.DataArray)):
            umean = xr.DataArray(umean.m, dims=['lev',latname], coords={'lev':plev, latname:ds[latname]})
            vptpclim = xr.DataArray(vptpclim.m, dims=['lev',latname], coords={'lev':plev, latname:ds[latname]})
            upvpclim = xr.DataArray(upvpclim.m, dims=['lev',latname], coords={'lev':plev, latname:ds[latname]})
            tptpclim = xr.DataArray(tptpclim.m, dims=['lev',latname], coords={'lev':plev, latname:ds[latname]})

    umean = umean.assign_attrs(long_name="Mean zonal wind",units="m/s")
    vptpclim = vptpclim.assign_attrs(long_name="Northward eddy heat flux", units="K m s$^{-1}$")
    upvpclim = upvpclim.assign_attrs(long_name="Northward eddy momentum flux", units="m$^2$ s$^{-2}$")
    tptpclim = tptpclim.assign_attrs(long_name="Eddy temperature variance", units="K$^{2}$")

    if args.verbose:
        print(f"{umean.shape = }, {type(umean) =}")
        print(f"{vptpclim.shape = }, {type(vptpclim) =}")
        print(f"{upvpclim.shape = }, {type(upvpclim) =}")
        print(f"{tptpclim.shape = }, {type(tptpclim) =}")

    # at this point, resulting arrays must be 2D:
    for i, v in enumerate([umean, vptpclim, upvpclim, tptpclim]):
        if len(v.shape) != 2:
            raise ValueError(f"Something is wrong with array {i} shape: {v.shape}")
    
    # TODO: add option to dump these to a file.
    # make the multi-panel figure
    #
    fig, ax = plt.subplots(figsize=(9,9), ncols=2, nrows=2, constrained_layout=True)
    contour_zm(ax[0,0], umean, 2, -40, 40, latname=latname, levname=levname)
    contour_zm(ax[0,1], tptpclim, 5, -60, 60, latname=latname, levname=levname)
    contour_zm(ax[1,0], upvpclim, 4, -120,120, latname=latname, levname=levname)
    contour_zm(ax[1,1], vptpclim, 1, -40,40, latname=latname, levname=levname)

    if args.saveplt:
        ofil = Path.home() / "hs94_diag_plot.png"
        fig.savefig(ofil, bbox_inches='tight')
        print(f"Output saved to: {ofil}")
    else:
        plt.show()
    print('ALL DONE.')

    
