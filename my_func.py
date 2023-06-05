import numpy as np
from scipy.interpolate import interp1d
from math import sqrt, sin, cos, asin, pi, acos, log10, floor, ceil, exp, log
from math import gamma, tan
from scipy.special import gammainc
from mplchange import *
import cst
from datetime import datetime
from dateutil import tz


def Gamma(b):
    return (np.sqrt(1-b**2))**-1


def Beta(g):
    return np.sqrt(1-1/g**2)


def pri(x):
    print('%.8e' % x)
    return None


# luminosity distance (Gpc) at redshift z
def D_L(z):  # unit Gpc
    # Luminosity distance according to the concordance Lambda CDM model
    Omega_m = 0.308
    Omega_Lambda = 1. - Omega_m
    c = 3.e10
    H_0 = 6.78e6
    Nx = 100
    x = np.linspace(0, z, Nx)
    dx = x[2]-x[1]
    temp = 0
    for i in range(Nx):
        temp += dx/np.sqrt(Omega_m*(1+x[i])**3 + Omega_Lambda)
    return c*(1+z)/H_0 *temp/1e3


def D_A(z):   # angular diameter distance
    return D_L(z)/(1+z)**2


def dV_dz(z): # unit Gpc^3
    # comoving volume per redshift bin
    Omega_m = 0.308
    Omega_Lambda = 1. - Omega_m
    c = 3e10
    H_0 = 6.78e6
    Nx = 100
    x = np.linspace(0, z, Nx)
    dx = x[2]-x[1]
    temp = 0
    for i in range(Nx):
        temp += dx/sqrt(Omega_m*(1+x[i])**3 + Omega_Lambda)
    return 4*pi*c**3/H_0**3/1e9\
               *temp*temp/sqrt(Omega_m*(1+z)**3 + Omega_Lambda)


def z_DM(dm_igm):
    # redshift as a function of IGM DM
    zmin = 0.003
    zmax = 6.
    Nz = 300
    z_arr = np.logspace(log10(zmin), log10(zmax), Nz)
    z_ratio = z_arr[1]/z_arr[0]
    DM = np.zeros(Nz)  # dispersion measure, in pc cm^-3
    Omega_m = 0.308
    Omega_Lambda = 1. - Omega_m
    Omega_b = 0.0486
    c = 2.99792458e10
    H0 = 6.78e6/3.09e24
    figm = 0.83  # igm baryon fraction
    Ye = 0.875  # electron fraction
    zmin /= sqrt(z_ratio)   # left boundary
    temp = zmin + 0.5*(1-1.5*Omega_m)*zmin**2  # lowest-order integ
    for i in range(Nz):
        z = z_arr[i]
        dz = z*(sqrt(z_ratio) - 1./sqrt(z_ratio))
        temp += dz*(1 + z)/sqrt(Omega_m*(1+z)**3 + Omega_Lambda)
        DM[i] = temp
    z_arr *= sqrt(z_ratio)  # right boundary
    numeric_factor = 1./(8*pi*6.673e-8*1.675e-24*3.09e18)   # 1.152e11
    DM *= 3*c*H0*Omega_b*figm*Ye*numeric_factor
    y = interp1d(DM, z_arr)
    return y(dm_igm)


def DM_z(z):
    # redshift as a function of IGM DM
    Omega_m = 0.308
    Omega_Lambda = 1. - Omega_m
    Omega_b = 0.0486
    c = 2.99792458e10
    H0 = 6.78e6/3.09e24
    figm = 0.83  # igm baryon fraction
    Ye = 0.875  # electron fraction
    Nx = 300
    x_arr = np.logspace(log10(z/100.), log10(z), Nx, endpoint=True)
    dx_ratio = x_arr[2]/x_arr[1]
    temp = 0.
    for i in range(Nx):
        x = x_arr[i]
        dx = x*(sqrt(dx_ratio) - 1./sqrt(dx_ratio))
        if i == Nx-1:   # for the last step
            dx *= 0.5   # only integrate half a step
        temp += dx*(1 + x)/np.sqrt(Omega_m * (1+x)**3 + Omega_Lambda)
    numeric_factor = 1./(8*pi*6.673e-8*1.675e-24*3.09e18)   # 1.152e11
    DM = temp*3*c*H0*Omega_b*figm*Ye*numeric_factor
    return DM


def r_isco(a):    # ISCO radius in units of rg
    if a < -1 or a > 1:
        print('incorrect spin!')
        return 0
    x1 = (1+a)**(1./3)
    x2 = (1-a)**(1./3)
    z1 = 1 + x1*x2*(x1 + x2)
    z2 = sqrt(3*a*a + z1*z1)
    if a > 0:
        return 3 + z2 - sqrt((3-z1) * (3 + z1 + 2*z2))
    return 3 + z2 + sqrt((3-z1) * (3 + z1 + 2*z2))
    

def eta_disk(a):    # disk accretion efficiency at ISCO
    return 1 - sqrt(1 - 2./3/r_isco(a))

    
def Porb(m1, m2, a):   # orbital period of a binary system [in days]
    # a is the binary separation (or reduced semimajor axis) [in Rsun]
    # m1 and m2 are in Msun
    return 2*pi*sqrt((a*cst.Rsun)**3/(cst.G*(m1+m2)*cst.Msun))/86400


def tgw(m1, m2, a0, e0):
    # GW merger time [in yr] for given masses [in Msun]
    # e is eccentricity, a is the separation between the binary [in Rsun]
    print('tgw(m1 [Msun], m2 [Msun], a0 [Rsun], e0)= ### [yr]')
    print('initial orbital period [day]', Porb(m1, m2, a0))
    beta = (64./5)*cst.G**3 * m1 * m2 * (m1+m2) * cst.Msun**3 / cst.c**5
    Tc = a0**4*cst.Rsun**4/(4*beta)
    if e0 < 1e-5:   # circular orbit
        return Tc / cst.yr
    e = 1e-3*e0
    dlge = 1e-4
    integ = 0.
    while e < 1:
        de = min(e, 1-e)*dlge
        dinteg = de * e**(29/19) * (1 + (121./304)*e*e)**(1181./2299) * (1 - e*e)**(-1.5)
        integ += dinteg
        e += de
        if e > e0:
            # linear interpolation
            integ += dinteg/de * (e0-e)
            break
    T_over_Tc = integ * 48./19 * (1-e0*e0)**4/e0**(48./19) \
                * (1 + 121./304*e0*e0)**(-4*870./2299)
    return T_over_Tc * Tc/cst.yr


def tgw_fit(m1, m2, a0, e0):
    # GW merger time [in yr] for given masses [in Msun]
    # e is eccentricity, a is the separation between the binary [in Rsun]
    # print('tgw(m1 [Msun], m2 [Msun], a0 [Rsun], e0)= ### [Gyr]')
    # print('initial orbital period [day]', Porb(m1, m2, a0))
    beta = (64./5)*cst.G**3 * m1 * m2 * (m1+m2) * cst.Msun**3 / cst.c**5
    Tc = a0**4*cst.Rsun**4/(4*beta)
    T0 = Tc * (1-e0**2)**3.5
    p = 2.66   # a fitting parameter
    return T0 * (768./425 - p*(1-e0*e0)**0.5 + (p-343./425)*(1-e0*e0)**0.8) / cst.yr


def absM_Lnu(Lnu):  # absolute magnitude for specific luminosity [erg/s/Hz]
    fnu = Lnu/(4*pi*(10*cst.pc)**2)*1e23  # flux density in Jy
    return -2.5*log10(fnu/3630.780547701)


def Lnu_absM(absM):
    fnu = 3630.780547701 * 1e-23 * 10**(-absM/2.5)
    return fnu * (4*pi*(10*cst.pc)**2)


def pltimg(ax, xarr, yarr, zarr, xlabel, ylabel, zlabel, cmap,
           CB_levels, CB_ticklabels, flag_contour):
    min_val, max_val = np.amin(zarr), np.amax(zarr)
    ax.set_xlabel(xlabel, labelpad=-2)
    ax.set_ylabel(ylabel)
    im = ax.imshow(zarr.transpose(),
                   interpolation='bicubic', origin='lower',
                   cmap=cmap, aspect='auto', alpha=0.7,
                   extent=(min(xarr), max(xarr),
                           min(yarr), max(yarr)))
    im.set_clim(vmin=min_val, vmax=max_val)

    CB = pl.colorbar(im, ax=ax, ticks=CB_levels)
    CB.ax.set_yticklabels(CB_ticklabels)
    CB.ax.set_ylabel(zlabel, labelpad=3)
    CB.ax.minorticks_off()

    if flag_contour:
        X, Y = np.meshgrid(xarr, yarr)
        CS = ax.contour(X, Y, zarr.transpose(),
                        CB_levels, linestyles='solid',
                        colors='k', linewidths=2, alpha=0.5)
        fmt = {}
        for l, s in zip(CS.levels, CB_ticklabels):
            fmt[l] = s
        pl.clabel(CS, CS.levels, inline=True, fmt=fmt,
                  fontsize=30, colors=None)


def time(from_zone):
    zones = [from_zone]
    names = []
    for zone in zones:
        if zone in ['ET', 'Easter Time', 'Easter_Time',
                    'easter_time', 'Easter', 'easter',
                    'edt', 'est', 'et', 'US/Eastern',
                    'EDT', 'EST', 'America/New_York']:
            name = 'America/New_York'
        elif zone in ['US/Central', 'CT', 'CDT', 'CST',
                      'ct', 'cdt', 'cst', 'America/Chicago']:
            name = 'US/Central'
        elif zone in ['US/Mountain', 'America/Denver',
                       'MST', 'MDT', 'MT', 'mst', 'mdt',
                       'mt', 'America/Phoenix', 'US/Arizona']:
            name = 'US/Mountain'
        elif zone in ['UC/Pacific', 'PST', 'PDT', 'PT',
                       'pst', 'pdt', 'pt', 'America/Los_Angeles']:
            name = 'US/Pacific'
        elif zone in ['Asia/Shanghai', 'Beijing', 'China',
                      'Beijing_Time', 'beijing', 'beijing_time',
                      'Beijing Time', 'beijing time', 'BJT', 'bjt']:
            name = 'Asia/Shanghai'
        elif zone in ['Israel', 'Asia/Jerusalem', 'IST', 'IDT',
                      'ist', 'idt']:
            name = 'Israel'
        elif zone in ['UTC', 'utc', 'ut']:
            name = 'UTC'
        else:
            name = 'UTC'
        names += [name]
    print(names[0])
    return datetime.now(tz.gettz(names[0]))



def round_sig(x, sig=2):
    # round a number to a given significant digit
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def ximax_polytrope(n):
    return 0.05583*n**4 - 0.2165*n**3 + 0.5656*n**2 + 0.1914*n + 2.544


def phimax_polytrope(n):
    return 0.02393*n**4 - 0.2363*n**3 + 0.989*n**2 - 2.402*n + 4.767
