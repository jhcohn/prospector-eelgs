import numpy as np
import os
from prospect.models import priors, sedmodel
from prospect.sources import FastStepBasis  # NEW, replaced CSPBasis
from sedpy import observate
from astropy.cosmology import WMAP9

tophat = priors.tophat
logarithmic = priors.logarithmic

#############
# RUN_PARAMS
#############

id = str(1824)  # 1824 (eelg), 2329 (normal), 3921 (normal-er)

run_params = {'verbose': True,
              'debug': False,
              'outfile': 'output/' + id,
              'nofork': True,
              # Optimizer params
              'ftol': 0.5e-5,
              'maxfev': 5000,
              # MCMC params
              'nwalkers': 140,
              'nburn': [50, 100],
              'niter': 500,  # 900
              'interval': 0.2,
              # Model info
              'zcontinuous': 2,
              'compute_vega_mags': False,
              'initial_disp': 0.1,
              'interp_type': 'logarithmic',
              'agelims': [0.0, 8.0, 8.5, 9.0, 9.5, 9.8, 10.0],  # NEW reset these? (see load_model)
              # Data info
              'photname': '/home/jonathan/cosmos/cosmos.v1.3.8.cat',
              'datname': '/home/jonathan/cosmos/cosmos.v1.3.8.cat',
              # 'fastname': '/home/jonathan/cosmos/cosmos.v1.3.6.awk.fout',  # .fout edited to correct format using awk
              'zname': '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout',  # main cat has z_spec, but not z_phot
              'objname': id,
              'convergence_check_interval': 100,  # Fix convergence test problem
              'convergence_kl_threshold': 0.0  # Fix convergenve test problem
              }
run_params['outfile'] = run_params['outfile'] + '_' + run_params['objname']

############
# FILTERS
#############
# Names of filter columns in the datafile (datname)
filternames = ['B', 'G', 'I', 'IA427', 'IA484', 'IA505', 'IA527', 'IA624', 'IA709', 'IA738', 'R', 'U', 'V', 'Rp', 'Z',
               'Zp', 'Hl', 'Hs', 'J1', 'J2', 'J3', 'Ks', 'NB118', 'NB209', 'F125W', 'F140W', 'F160W', 'F606W', 'F814W',
               'UVISTA_J', 'UVISTA_H', 'UVISTA_Ks', 'UVISTA_Y', 'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']

# Names of filter response curve files (in the Filters folder)
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
nb = [folder+nb_pre+n for n in ['118', '209']]
hst_prefix = 'hst-wfc3-'
hst = [folder+hst_prefix+n for n in ['IR-f125w', 'IR-f140w', 'IR-f160w', 'UVIS-f606w', 'UVIS-f814w']]
vista_prefix = 'VISTA-'
vista_suff = '_system+atmos'
vista = [folder+vista_prefix+n+vista_suff for n in ['J', 'H', 'Ks', 'Y']]
irac_pre = 'IRAC-irac_tr'
irac_suffix = '_2004-08-09'
irac = [folder+irac_pre+n+irac_suffix for n in ['1', '2', '3', '4']]

filts = capak_1 + mega_1 + subaru + mega_2 + capak_0 + mega_3 + capak_2 + fstar + nb + hst + vista + irac
filtersets = (filts, filts)

############
# OBS
#############
def load_obs(photname, objname, err_floor=0.05, zperr=True, **extras):
    '''
    photname: photometric file location
    objname: number of object in the 3D-HST COSMOS photometric catalog
    err_floor: the fractional error floor (0.05 = 5% floor)
    zp_err: inflate the errors by the zeropoint offsets from Skelton+14
    '''

    ### open file, load data
    with open(photname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(photname, comments='#', delimiter=' ',
                     dtype=dtype)

    ### extract filters, fluxes, errors for object
    # from ReadMe: "All fluxes are normalized to an AB zeropoint of 25, such that: magAB = 25.0-2.5*log10(flux)
    obj_idx = (dat['id'] == objname)
    # print(dat[obj_idx]['id'], 'idx')
    filters = np.array(filts)  # [f[2:] for f in dat.dtype.names if f[0:2] == 'f_'])
    flux = np.squeeze([dat[obj_idx]['f_' + f] for f in filternames])
    unc = np.squeeze([dat[obj_idx]['e_' + f] for f in filternames])

    ### define photometric mask, convert to maggies
    phot_mask = (flux != unc) & (flux != -99.0)
    maggies = flux * 10**-6 / 3631  # flux [uJy] * 1e-6 [Jy / uJy] * 1 [maggy] / 3631 [Jy]
    maggies_unc = unc * 10**-6 / 3631
    # print(maggies, 'maggies')
    # print(flux, 'flux')
    # print(maggies_unc, 'maggies_unc')
    # print(unc, 'unc')

    ### implement error floor
    maggies_unc = np.clip(maggies_unc, maggies * err_floor, np.inf)

    ### build output dictionary
    obs = {}
    obs['filters'] = observate.load_filters(filters)
    obs['wave_effective'] = np.array([filt.wave_effective for filt in obs['filters']])
    obs['phot_mask'] = phot_mask
    obs['maggies'] = maggies
    obs['maggies_unc'] = maggies_unc
    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['logify_spectrum'] = False

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


def transform_logtau_to_tau(tau=None, logtau=None, **extras):
    return 10**logtau

#############
# MODEL_PARAMS
#############

model_params = []

###### BASIC PARAMETERS ##########
model_params.append({'name': 'zred', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
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
                     'init': 10.0,
                     'units': 'Msun',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini': 1.0, 'maxi': 14.0}})

model_params.append({'name': 'mass', 'N': 1,
                     'isfree': False,
                     'init': 1e10,
                     'depends_on': transform_logmass_to_mass,
                     'units': 'Msun',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini': 1e1, 'maxi': 1e14}})

model_params.append({'name': 'logzsol', 'N': 1,
                     'isfree': False,
                     'init': 0.0,
                     'init_disp': 0.4,
                     'log_param': True,
                     'units': r'$\log (Z/Z_\odot)$',
                     'prior_function': tophat,
                     'prior_args': {'mini': -2.0, 'maxi': 0.19}})

###### SFH   ########
model_params.append({'name': 'sfh', 'N': 1,
                     'isfree': False,
                     'init': 0,  # NEW: 0 = non-parametric SFH; 1 = parametric SFH
                     'units': 'type',
                     'prior_function_name': None,
                     'prior_args': None})

model_params.append({'name': 'tau', 'N': 1,
                     'isfree': False,
                     'init': 10,
                     'depends on': transform_logtau_to_tau,
                     'init_disp': 0.5,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.1,
                                    'maxi': 100.0}})

model_params.append({'name': 'logtau', 'N': 1,
                        'isfree': False,  # NEW turn off
                        'init': 1,
                        'init_disp': 0.5,
                        'units': 'Gyr',
                        'prior_function':tophat,
                        'prior_args': {'mini':-1, 'maxi':2}})

model_params.append({'name': 'tage', 'N': 1,
                     'isfree': False,  # NEW turn off
                     'init': 1.0,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.1, 'maxi': 14.0}})  # 0.01, 'maxi': 14.0}})

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
                     'init': 0.0,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     'prior_args': {'mini': 0.0,
                                    'maxi': 14.0}})

model_params.append({'name': 'agebins', 'N': 1,  # NEW
                        'isfree': False,
                        'init': [],
                        'units': 'log(yr)',
                        'prior_function': priors.tophat,
                        'prior_args':{'mini':0.1, 'maxi':15.0}})

model_params.append({'name': 'sfr_fraction', 'N': 1,  # NEW
                        'isfree': True,
                        'init': [],
                        'units': 'Msun',
                        'prior_function': priors.tophat,
                        'prior_args':{'mini':0.0, 'maxi':1.0}})

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
                     'prior_args': {'mini': 0.0,
                                    'maxi': 4.0}})

###### Dust Emission ##############
model_params.append({'name': 'add_dust_emission', 'N': 1,
                     'isfree': False,
                     'init': 0,
                     'units': None,
                     'prior_function': None,
                     'prior_args': None})

###### Nebular Emission ###########
model_params.append({'name': 'add_neb_emission', 'N': 1,  # BUCKET NEED TO TURN ANYTHING ELSE ON TOO OR NO?
                     'isfree': False,
                     'init': True,  # BUCKET1 emission lines make init: True instead of isfree?
                     'units': r'log Z/Z_\odot',
                     'prior_function_name': None,  # BUCKET1 prior_function_name: None throws error when isfree: True
                     'prior_args': None})

model_params.append({'name': 'gas_logz', 'N': 1,
                     'isfree': False,  # BUCKET need to turn this on too or no?
                     'init': 0.0,
                     'depends_on': tie_gas_logz,  # BUCKET1 tie_gas_logz added
                     'units': r'log Z/Z_\odot',
                     'prior_function': tophat,
                     'prior_args': {'mini': -2.0, 'maxi': 0.5}})

model_params.append({'name': 'gas_logu', 'N': 1,
                     'isfree': True,  # BUCKET1 emission lines
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


#### resort list of parameters ####
#### so that major ones are fit first ####
parnames = [m['name'] for m in model_params]
# fit_order = ['logmass', 'tage', 'logtau', 'dust2']  # for FAST mimic
fit_order = ['logmass', 'sfr_fraction', 'dust2']
tparams = [model_params[parnames.index(i)] for i in fit_order]
for param in model_params:
    if param['name'] not in fit_order:
        tparams.append(param)
model_params = tparams


###### REDEFINE MODEL FOR MY OWN NEFARIOUS PURPOSES ######
class BurstyModel(sedmodel.SedModel):  # NEW, replacing below

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

        # sum of SFH fractional bins <= 1.0
        print(self.theta_index['sfr_fraction'], 'theta_index')
        print(theta, 'theta')
        # print(theta[0], 'sfr_theta_index')
        if 'sfr_fraction' in self.theta_index:
            sfr_fraction = theta[self.theta_index['sfr_fraction']]
            if np.sum(sfr_fraction) > 1.0:
                return -np.inf

        for k, v in self.theta_index.iteritems():
            this_prior = np.sum(self._config_dict[k]['prior_function']
                                (theta[v], **self._config_dict[k]['prior_args']))

            if (not np.isfinite(this_prior)):
                print('WARNING: ' + k + ' is out of bounds')
            lnp_prior += this_prior
        return lnp_prior

'''
###### REDEFINE MODEL FOR MY OWN NEFARIOUS PURPOSES ######
class BurstyModel(sedmodel.SedModel):  #old, replaced by above
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
'''


class FracSFH(FastStepBasis):  # NEW CLASS
    def get_galaxy_spectrum(self, **params):
        self.update(**params)

        #### here's the custom fractional stuff
        fractions = np.array(self.params['sfr_fraction'])
        bin_fractions = np.append(fractions, (1 - np.sum(fractions)))
        time_per_bin = []
        for (t1, t2) in self.params['agebins']: time_per_bin.append(10 ** t2 - 10 ** t1)
        bin_fractions *= np.array(time_per_bin)
        bin_fractions /= bin_fractions.sum()

        mass = bin_fractions * self.params['mass']
        mtot = self.params['mass'].sum()

        time, sfr, tmax = self.convert_sfh(self.params['agebins'], mass)
        self.ssp.params["sfh"] = 3  # Hack to avoid rewriting the superclass
        self.ssp.set_tabular_sfh(time, sfr)
        wave, spec = self.ssp.get_spectrum(tage=tmax, peraa=False)

        return wave, spec / mtot, self.ssp.stellar_mass / mtot

def load_sps(**extras):
    sps = FracSFH(**extras)  # NEW, replaced CSPBasis
    return sps


def load_model(objname='', datname='', zname='', agelims=[], **extras):
    ###### REDSHIFT ######
    ### open file, load data

    with open(datname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1],'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(datname, comments = '#', delimiter=' ',
                     dtype = dtype)

    with open(zname, 'r') as fz:
        hdr_z = fz.readline().split()
    dtype_z = np.dtype([(hdr_z[1],'S20')] + [(n, np.float) for n in hdr_z[2:]])
    zout = np.loadtxt(zname, comments = '#', delimiter=' ', dtype = dtype_z)

    idx = dat['id'] == objname
    zred = zout['z_spec'][idx][0]  # BUCKET FIX THE INDEX HERE
    if zred == -99:
        zred = zout['z_peak'][idx][0]

    print(zred, 'zred')

    #### CALCULATE TUNIV #####
    tuniv = WMAP9.age(zred).value
    print(tuniv, 'tuniv')

    n = [p['name'] for p in model_params]
    model_params[n.index('tage')]['prior_args']['maxi'] = tuniv

    #### NONPARAMETRIC SFH ######  # NEW
#     agelims[-1] = np.log10(tuniv*1e9)
    ncomp = len(agelims) - 1
    agelims = [0.0, 7.0, 8.0, (8.0 + (np.log10(tuniv*1e9)-8.0)/4), (8.0 + 2*(np.log10(tuniv*1e9)-8.0)/4),
               (8.0 + 3*(np.log10(tuniv*1e9)-8.0)/4), np.log10(tuniv*1e9)]
    agebins = np.array([agelims[:-1], agelims[1:]])
    # calculate the somethings: [0, a, b, b + (f-b)/4, b + 2*(f-b)/4, b + 3*(f-b)/4, b + 4*(f-b)/4 = f]

    #### INSERT REDSHIFT INTO MODEL PARAMETER DICTIONARY ####
    zind = n.index('zred')
    model_params[zind]['init'] = zred

    #### SET UP AGEBINS
    model_params[n.index('agebins')]['N'] = ncomp
    model_params[n.index('agebins')]['init'] = agebins.T

    #### FRACTIONAL MASS INITIALIZATION  # NEW
    # N-1 bins, last is set by x = 1 - np.sum(sfr_fraction)
    model_params[n.index('sfr_fraction')]['N'] = ncomp-1
    model_params[n.index('sfr_fraction')]['prior_args'] = {
                                                           'maxi':np.full(ncomp-1,1.0),
                                                           'mini':np.full(ncomp-1,0.0),
                                                           # NOTE: ncomp instead of ncomp-1 makes the prior take into
                                                           # account the implicit Nth variable too
                                                          }
    model_params[n.index('sfr_fraction')]['init'] =  np.zeros(ncomp-1)+1./ncomp
    model_params[n.index('sfr_fraction')]['init_disp'] = 0.02

    #### CREATE MODEL
    model = BurstyModel(model_params)

    return model


model_type = BurstyModel
