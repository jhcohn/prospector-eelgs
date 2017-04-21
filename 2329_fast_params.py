import numpy as np
import os
from prospect.models import priors, sedmodel
from prospect.sources import CSPBasis
from sedpy import observate
from astropy.cosmology import WMAP9
# from td_io import load_zp_offsets

tophat = priors.tophat
logarithmic = priors.logarithmic
APPS = os.getenv('APPS')

#############
# RUN_PARAMS
#############

'''
    # do Photmask if -99: False
    # email Themiya
    tau = 10**fast['ltau']/1e9
    tage = 10**fast['lage']/1e9
    logmass = fast['lmass']
    dust2 = fast['Av']/ 1.086

    # obj_maggies.append((obj[step] / (10 ** ((23.9 - 25) / 2.5))))  # according to catalog
    # errs_maggies.append((obj[step + 1] / (10 ** ((23.9 - 25) / 2.5))))
    # obj_maggies.append((obj[step] / (10 ** ((23.9 - 25) / 2.5))) * 10 ** - (25 / 2.5))  # old, wrong
    # errs_maggies.append((obj[step + 1] / (10 ** ((23.9 - 25) / 2.5))) * 10 ** - (25 / 2.5))  # old, wrong
'''

run_params = {'verbose': True,
              'debug': False,
              'outfile': 'output/2329',
              # 'outfile': os.getenv('APPS') + '/threedhst_bsfh/results/td_massive/td_massive',
              'nofork': True,
              # Optimizer params
              'ftol': 0.5e-5,
              'maxfev': 5000,
              # MCMC params
              'nwalkers': 128,  # 620,
              'nburn': [50, 100],  # [150, 200, 400, 600],
              'niter': 500,  # 2500,
              'interval': 0.2,
              # Model info
              'zcontinuous': 2,
              'compute_vega_mags': False,
              'initial_disp': 0.1,
              'interp_type': 'logarithmic',
              'agelims': [0.0, 8.0, 8.5, 9.0, 9.5, 9.8, 10.0],
              # Data info
              'photname': 'cosmos_2329_maggies_phot_direct.dat',  # 'cosmos_2329_maggies_phot_direct.dat',
              'errname': 'cosmos_2329_maggies_err_direct.dat',  # 'cosmos_2329_maggies_err_direct.dat',
              'fastname': '/home/jonathan/cosmos/cosmos-2329-fout',  # '/home/jonathan/cosmos/cosmos.v1.3.6.fout',
              'objname': '2329',
#              'photname': APPS + '/threedhst_bsfh/data/COSMOS_td_massive.cat',
#              'datname': APPS + '/threedhst_bsfh/data/COSMOS_td_massive.dat',
#              'fastname': APPS + '/threedhst_bsfh/data/COSMOS_td_massive.fout',
              # 'objname': '579',
              'objid': 0
              }
run_params['outfile'] = run_params['outfile'] + '_' + run_params['objname']


#############
# OBS
#############
'''
filts = [n+'_filter' for n in ['B', 'G', 'I', 'IA427', 'IA484', 'IA505', 'IA527', 'IA624', 'IA709', 'IA738', 'R', 'U',
                               'V', 'Rp', 'Z', 'Zp', 'Hl', 'Hs', 'J1', 'J2', 'J3', 'Ks', 'NB118', 'NB209', 'F125W',
                               'F140W', 'F160W', 'F606W', 'F814W', 'uvista_J', 'uvista_H', 'uvista_Ks', 'uvista_Y',
                               'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']]
'''

folder = 'cat_filters/'
mega_prefix = 'megaprime-cfht_mega_'
mega = [folder+mega_prefix+n for n in ['u_cfh9301', 'g_cfh9401', 'r_cfh9601', 'i_cfh9701', 'z_cfh9801']]  # U G R I Z
mega_1 = [mega[1], mega[3]]  # G, I
mega_2 = [mega[2], mega[0]]  # R, U
mega_3 = [mega[4]]  # Z
capak_prefix = 'CAPAK_v2-'
capak_suffix = '_subaru'
capak = [folder+capak_prefix+n+capak_suffix for n in ['V', 'B', 'r', 'z']]  # V, B, Rp, Zp
capak_0 = [capak[0], capak[2]]  # V, Rp
capak_1 = [capak[1]]  # B
capak_2 = [capak[3]]  # Zp
subaru_prefix = 'Subaru_MB-IA'
subaru = [folder+subaru_prefix+n for n in ['427', '484', '505', '527', '624', '709', '738']]
fourstar_pre = 'FOURSTAR-'
fourstar_suffix = '_cam_optics_sky'
fstar = [folder+fourstar_pre+n+fourstar_suffix for n in ['Hlong', 'Hshort', 'J1', 'J2', 'J3', 'Ks']]
nb_pre = 'LCO_FourStar_NB'
# nb_suff = '_filter'
nb = [folder+nb_pre+n for n in ['118', '209']]
hst_prefix = 'hst-wfc3-'
hst = [folder+hst_prefix+n for n in ['IR-f125w', 'IR-f140w', 'IR-f160w', 'UVIS-f606w', 'UVIS-f814w']]
# hst = [folder+hst_prefix+n for n in ['IR-f125w', 'IR-f160w', 'UVIS-f606w', 'UVIS-f814w']]  # remove f140W for 2329
vista_prefix = 'VISTA-'
vista_suff = '_system+atmos'
vista = [folder+vista_prefix+n+vista_suff for n in ['J', 'H', 'Ks', 'Y']]
irac_pre = 'IRAC-irac_tr'
irac_suffix = '_2004-08-09'
irac = [folder+irac_pre+n+irac_suffix for n in ['1', '2', '3', '4']]

# filts = [folder+n for n in [capak_0, mega_1, subaru, mega_2, capak_1, fstar, nb, hst, vista, irac]
# filts = [capak[0], mega[1], mega[3], subaru, mega[2], mega[0], mega[4], capak[1], capak[2], fstar, nb, hst, vista, irac]
filts = capak_1 + mega_1 + subaru + mega_2 + capak_0 + mega_3 + capak_2 + fstar + nb + hst + vista + irac
filtersets = (filts, filts)

def load_obs(objid, phottable='cosmos_2329_maggies_phot_direct.dat',  # 'cosmos_2329_maggies_phot_direct.dat',  # cosmos_1824_maggies_phot.dat
             errtable='cosmos_2329_maggies_err_direct.dat', **kwargs):  # cosmos_1824_maggies_err.dat  'cosmos_2329_maggies_err_direct.dat',

    """Load photometry from an ascii file.  Assumes the following columns:
    `objid`, `filterset`, [`mag0`,....,`magN`] where N >= 11.  The User should
    modify this function (including adding keyword arguments) to read in their
    particular data format and put it in the required dictionary.

    :param objid:
        The object id for the row of the photomotery file to use.  Integer.
        Requires that there be an `objid` column in the ascii file.

    :param phottable:
        Name (and path) of the ascii file containing the photometry.

    :returns obs:
        Dictionary of observational data.
    """
    # Writes your code here to read data.  Can use FITS, h5py, astropy.table,
    # sqlite, whatever.
    # e.g.:
    # import astropy.io.fits as pyfits
    # catalog = pyfits.getdata(phottable)

    # Here we will read in an ascii catalog of magnitudes as a numpy structured
    # array
    with open(phottable, 'r') as f:
        # drop the comment hash
        header = f.readline().split()[1:]
    catalog = np.genfromtxt(phottable, comments='#',
                            dtype=np.dtype([(n, np.float) for n in header]))

    # Find the right row
    ind = catalog['objid'] == float(objid)  # = True
    # Here we are dynamically choosing which filters to use based on the object
    # and a flag in the catalog.  Feel free to make this logic more (or less)
    # complicated.
    filternames = filtersets[int(catalog[ind]['filterset']) ]  # catalog[True]['filterset'] = 0.0

    # GET INFO FROM PHOTTABLE
    maggies = []
    i = 2  # skipping 'objid' and 'filterset' indices (i = 0, 1)
    while i < len(catalog[ind][0]):
        maggies.append(catalog[ind][0][i])  # append photometry in units of maggies
        i += 1
    # mags = np.array(mags)
    maggies = np.array(maggies)
    print(maggies)

    # Build output dictionary.
    obs = {}
    # This is a list of sedpy filter objects.    See the
    # sedpy.observate.load_filters command for more details on its syntax.
    obs['filters'] = observate.load_filters(filternames)
    obs['wave_effective'] = np.array([filt.wave_effective for filt in obs['filters']])

    # This is a list of maggies, converted from mags.  It should have the same
    # order as `filters` above.
    # obs['maggies'] = np.squeeze(10**(-mags/2.5))
    # INPUT IS NOW IN MAGGIES YAY
    obs['maggies'] = maggies

    # FLUX UNCERTAINTIES
    with open(errtable, 'r') as f_err:
        err_header = f_err.readline().split()[1:]
    err_catalog = np.genfromtxt(errtable, comments='#',
                            dtype=np.dtype([(n, np.float) for n in err_header]))
    errs = []
    i = 2  # skipping objid and filterset indices (i = 0, 1)
    '''
    while i < len(err_catalog[ind][0]):
        if err_catalog[ind][0][i] / catalog[ind][0][i] < 0.05:  # BUCKET ERROR FLOORS
            errs.append(catalog[ind][0][i] * 0.05)  # append errors in units of 0.05 * flux (abs mags) if errs too small
        else:  # BUCKET ERROR FLOORS
            errs.append(err_catalog[ind][0][i])  # append errors in units of abs mags (NOW IN MAGGIES)
        i += 1
    '''
    while i < len(err_catalog[ind][0]):  # error floor implemented in maggies_phot.py that writes phot & err .dat files
        errs.append(err_catalog[ind][0][i])  # append errors in units of abs mags (NOW IN MAGGIES)
        i += 1

    errs = np.array(errs)
    print('err', errs)
    # obs['maggies_unc'] = np.squeeze(10**(-errs/2.5)) ALREADY ENTERED AS MAGGIES NOW!
    phot_mask = (maggies != errs) & (maggies != -99.0)
    obs['phot_mask'] = phot_mask
    '''
    obs['phot_mask'] = [True] * len(obs['maggies'])
    for i in range(len(obs['maggies'])):
        if obs['maggies'][i] == -2.72652161939e-08:  # -99 * 10 ** -6 / 3631:
            print('-99')
            obs['phot_mask'][i] = False
    '''

    obs['maggies_unc'] = errs

    # Here we mask out any NaNs or infs
    # obs['phot_mask'] = np.isfinite(np.squeeze(mags))  # BUCKET WHAT'S THIS?
    obs['phot_mask'] = np.isfinite(np.squeeze(maggies))  # BUCKET WHAT'S THIS?
    # We have no spectrum.
    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['logify_spectrum'] = False

    # Add unessential bonus info.  This will be stored in output
    # obs['dmod'] = catalog[ind]['dmod']
    obs['objid'] = objid
    f.close()  # closing file

    return obs


##########################
# TRANSFORMATION FUNCTIONS
##########################
def transform_logmass_to_mass(mass=None, logmass=None, **extras):
    return 10 ** logmass


def load_gp(**extras):
    return None, None


def add_dust1(dust2=None, **extras):
    return 0.86 * dust2


def tie_gas_logz(logzsol=None, **extras):
    return logzsol


#############
# MODEL_PARAMS
#############

model_params = []

###### BASIC PARAMETERS ##########
model_params.append({'name': 'zred', 'N': 1,
                     'isfree': False,
                     'init': 3.0815,  # 3.077,  # 0.0,
                     'units': '',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0, 'maxi': 4.0}})

model_params.append({'name': 'add_igm_absorption', 'N': 1,
                     'isfree': False,
                     'init': 1,
                     'units': None,
                     'prior_function': None,
                     'prior_args': None})

model_params.append({'name': 'add_agb_dust_model', 'N': 1,
                     'isfree': False,
                     'init': False,
                     'units': None,
                     'prior_function': None,
                     'prior_args': None})

model_params.append({'name': 'logmass', 'N': 1,
                     'isfree': True,
                     'init': 11.23,
                     'units': 'Msun',
                     'prior_function': priors.tophat,
                     # 'prior_args': {'mini': 9.0, 'maxi': 10.5}})
                     'prior_args': {'mini': 1.0, 'maxi': 14.0}})

model_params.append({'name': 'mass', 'N': 1,
                     'isfree': False,
                     'init': 1e11,
                     'depends_on': transform_logmass_to_mass,
                     'units': 'Msun',
                     'prior_function': priors.tophat,
                     # 'prior_args': {'mini': 1e9, 'maxi': 10 ** 10.5}})
                     'prior_args': {'mini': 1e1, 'maxi': 1e14}})

model_params.append({'name': 'logzsol', 'N': 1,
                     'isfree': False,
                     'init': 0.0,  # 0.02
                     'init_disp': 0.4,
                     'log_param': True,
                     'units': r'$\log (Z/Z_\odot)$',
                     'prior_function': tophat,
                     'prior_args': {'mini': -2.0, 'maxi': 0.19}})

###### SFH   ########
model_params.append({'name': 'sfh', 'N': 1,
                     'isfree': False,
                     'init': 1,
                     'units': 'type',
                     'prior_function_name': None,
                     'prior_args': None})

model_params.append({'name': 'tau', 'N': 1,
                     'isfree': True,
                     'init': 10,
                     'init_disp': 0.5,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     # 'prior_args': {'mini': 0.1, 'maxi': 1.0}})
                     'prior_args': {'mini': 0.1, 'maxi': 100.0}})

model_params.append({'name': 'tage', 'N': 1,
                     'isfree': True,
                     'init': 1.0,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     # 'prior_args': {'mini': 0.10001, 'maxi': 1.0}})
                     'prior_args': {'mini': 0.01, 'maxi': 14.0}})

model_params.append({'name': 'tburst', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
                     'init_disp': 1.0,
                     'units': '',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0, 'maxi': 10.0}})

model_params.append({'name': 'fburst', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
                     'init_disp': 0.5,
                     'units': '',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0, 'maxi': 0.2}})

model_params.append({'name': 'fconst', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
                     'units': '',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0, 'maxi': 1.0}})

model_params.append({'name': 'sf_start', 'N': 1,
                     'isfree': False,
                     'reinit': True,
                     'init': 0.0,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0,
                                    'maxi': 14.0}})

########    IMF  ##############
model_params.append({'name': 'imf_type', 'N': 1,
                     'isfree': False,
                     'init': 1,  # 1 = chabrier
                     'units': None,
                     'prior_function_name': None,
                     'prior_args': None})

######## Dust Absorption ##############
model_params.append({'name': 'dust_type', 'N': 1,
                     'isfree': False,
                     'init': 2,
                     'units': 'index',
                     'prior_function_name': None,
                     'prior_args': None})

model_params.append({'name': 'dust2', 'N': 1,
                     'isfree': True,
                     'init': 0.0,
                     'init_disp': 0.2,
                     'units': '',
                     'prior_function': tophat,
                     # 'prior_args': {'mini': 0.0, 'maxi': 0.5}})
                     'prior_args': {'mini': 0.0, 'maxi': 4.0}})

###### Dust Emission ##############
model_params.append({'name': 'add_dust_emission', 'N': 1,
                     'isfree': False,
                     'init': 0,
                     'units': None,
                     'prior_function': None,
                     'prior_args': None})

###### Nebular Emission ###########
model_params.append({'name': 'add_neb_emission', 'N': 1,
                     'isfree': False,
                     'init': False,
                     'units': r'log Z/Z_\odot',
                     'prior_function_name': None,
                     'prior_args': None})

model_params.append({'name': 'gas_logz', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
                     'units': r'log Z/Z_\odot',
                     'prior_function': tophat,
                     'prior_args': {'mini': -2.0, 'maxi': 0.5}})

model_params.append({'name': 'gas_logu', 'N': 1,
                     'isfree': False,
                     'init': -2.0,
                     'units': '',
                     'prior_function': tophat,
                     'prior_args': {'mini': -4, 'maxi': -1}})

####### Calibration ##########
model_params.append({'name': 'phot_jitter', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
                     'units': 'mags',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0, 'maxi': 0.2}})
model_params.append({'name': 'peraa', 'N': 1,
                     'isfree': False,
                     'init': False})

model_params.append({'name': 'mass_units', 'N': 1,
                     'isfree': False,
                     'init': 'mstar'})

#### resort list of parameters

#### resort list of parameters
#### so that major ones are fit first
parnames = [m['name'] for m in model_params]
fit_order = ['logmass', 'tage', 'tau', 'dust2']
tparams = [model_params[parnames.index(i)] for i in fit_order]
for param in model_params:
    if param['name'] not in fit_order:
        tparams.append(param)
model_params = tparams


###### REDEFINE MODEL FOR MY OWN NEFARIOUS PURPOSES ######
class BurstyModel(sedmodel.SedModel):
    def prior_product(self, theta):
        """
        Return a scalar which is the ln of the product of the prior
        probabilities for each element of theta.  Requires that the
        prior functions are defined in the theta descriptor.

        :param theta:
            Iterable containing the free model parameter values.

        :returns lnp_prior:
            The log of the product of the prior probabilities for
            these parameter values.
        """
        lnp_prior = 0

        for k, v in self.theta_index.iteritems():
            start, end = v
            this_prior = np.sum(self._config_dict[k]['prior_function']
                                (theta[start:end], **self._config_dict[k]['prior_args']))

            if (not np.isfinite(this_prior)):
                print('WARNING: ' + k + ' is out of bounds')
            lnp_prior += this_prior
        return lnp_prior


def load_sps(**extras):
    sps = CSPBasis(**extras)
    return sps


def load_model(objname='', datname='', fastname='', agelims=[], **extras):
    ###### REDSHIFT ######
    ### open file, load data
    '''
    with open(datname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1],'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(datname, comments = '#', delimiter=' ',
                     dtype = dtype)
    '''

    with open(fastname, 'r') as f:
        hdr = f.readline().split()
    # dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    fast = np.loadtxt(fastname, comments='#', delimiter=' ')  # , dtype=dtype)

    '''
    with open(fastname, 'r') as f:
        # drop the comment hash
        header = f.readline().split()
    fast = np.genfromtxt(fastname, comments='#', dtype=np.dtype([(n, np.float) for n in header]))
    f.close()
    '''
    idx = 0  # fast['id'] == objname
    print(fast)

    zred = fast[1]  # 3.0815
    print(zred)

    #### INITIAL VALUES
    tau = 10 ** fast[2] / 1e9  # 10 ** 7.80 / 1e9  #
    tage = 10 ** fast[4] / 1e9  # 10 ** 8.40 / 1e9  #
    logmass = fast[6]  # 11.23  #
    dust2 = fast[5] / 1.086  # 1.70 / 1.086  #
    # logzsol = 0  # np.log10(fast['metal'])

    n = [p['name'] for p in model_params]
    model_params[n.index('tau')]['init'] = tau
    model_params[n.index('tage')]['init'] = tage
    model_params[n.index('logmass')]['init'] = logmass
    model_params[n.index('dust2')]['init'] = dust2
    # model_params[n.index('logzsol')]['init'] = logzsol


    #### CALCULATE TUNIV #####
    tuniv = WMAP9.age(zred).value
    model_params[n.index('tage')]['prior_args']['maxi'] = tuniv

    #### INSERT REDSHIFT INTO MODEL PARAMETER DICTIONARY ####
    zind = n.index('zred')
    model_params[zind]['init'] = zred

    #### CREATE MODEL
    model = BurstyModel(model_params)

    return model


model_type = BurstyModel