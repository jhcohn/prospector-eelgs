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

run_params = {'verbose': True,
              'debug': False,
              'outfile': 'output/' + '',
              'nofork': True,
              # Optimizer params
              'ftol': 0.5e-5,
              'maxfev': 5000,
              # MCMC params
              'nwalkers': 140,
              'nburn': [50, 100],
              'niter': 1200,
              'interval': 0.2,
              # Model info
              'zcontinuous': 2,
              'compute_vega_mags': False,
              'initial_disp': 0.1,
              'interp_type': 'logarithmic',
              'agelims': [0.0, 8.0, 8.5, 9.0, 9.5, 9.8, 10.0],  # NEW (see load_model)
              # Data info
              'field': 'cosmos',
              'objname': '1824',
              'photname': '',
              'zname': '',
              'convergence_check_interval': 100,  # Fix convergence test problem
              'convergence_kl_threshold': 0.04  # Fix convergence test problem
              }
run_params['outfile'] = run_params['outfile'] + '_' + run_params['objname']  # + '_' +run_params['field']

############
# FILTERS
#############
folder = 'cat_filters/'  # folder in /sedpy-0.1.0-py2.7.egg/sedpy/data/filters containing all filters for all fields

irac_pre = 'IRAC-irac_tr'  # SAME for all fields
irac_suf = '_2004-08-09'  # SAME for all fields
irac = [folder + irac_pre + n + irac_suf for n in ['1', '2', '3', '4']]  # SAME for all fields
KsHI = ['VLT-hawki_k_ETC']  # SAME for CDFS, UDS (see: cdfs/eazy/cdfs.v1.6.9.translate, uds/eazy/uds.v1.5.8.translate)
fourstar_pre = 'FOURSTAR-'  # SAME for all fields
fourstar_suf = '_cam_optics_sky'  # SAME for all fields
subaru_pre = 'Subaru_MB-IA'  # SAME for all fields
nb_pre = 'LCO_FourStar_NB'  # SAME for COSMOS, CDFS
nb = [folder + nb_pre + n for n in ['118', '209']]  # SAME for COSMOS, CDFS


# CDFS
# Names of filter columns in the datafile (datname)
cdfs_filternames = ['B', 'I', 'R', 'U', 'V', 'Z', 'Hs', 'Hl', 'J1', 'J2', 'J3', 'Ks', 'KsHI', 'NB118', 'NB209', 'F098M',
                    'F105W', 'F125W', 'F140W', 'F160W', 'F814W', 'IA484', 'IA527', 'IA574', 'IA598', 'IA624', 'IA651',
                    'IA679', 'IA738', 'IA767', 'IA797', 'IA856', 'WFI_V', 'WFI_Rc', 'WFI_U38', 'tenisK', 'IRAC_36',
                    'IRAC_45', 'IRAC_58', 'IRAC_80']

# Names of filter response curve files (in the Filters folder)
acs_pre = 'hst-ACS_update_sep07-wfc_'
acs_suf = '_t77'
acs = [folder+acs_pre+n+acs_suf for n in ['f435w', 'f775w', 'f606w', 'f850lp', 'f814w']]
acs_1 = [acs[0], acs[1]]  # f435w, f775w = B, I
acs_2 = [acs[2], acs[3]]  # f606w, f850lp = V, Z
acs_3 = [acs[4]]  # F814W
vimos_pre = 'ESO-'
vimos = [vimos_pre+n for n in ['VIMOS-R', 'vimos_u']]  # R, U
fstar = [folder+fourstar_pre+n+fourstar_suf for n in ['Hshort', 'Hlong', 'J1', 'J2', 'J3', 'Ks']]  # DIFF ORDER
hstIR_pre = 'hst-wfc3-IR-f'
hstIR = [folder+hstIR_pre+n for n in ['098m', '105w', '125w', '140w', '160w']]
sub = [folder+subaru_pre+n for n in ['484', '527', '574', '598', '624', '651', '679', '738', '768', '797', '856']]
wfi_pre = 'ESO-'
wfi = [folder+wfi_pre+n for n in ['WFI-V89_843', 'WFI-Rc162_844', 'wfi_BB_U38_ESO841']]  # WFI_V, WFI_Rc, WFI_U38
tenisK = ['WIRCam-cfh8302_Ks']  # (according to: cdfs/eazy/cdfs.v1.6.9.translate)

cdfs_filts = acs_1 + vimos + acs_2 + fstar + KsHI + nb + hstIR + acs_3 + sub + wfi + tenisK + irac


# COSMOS
# Names of filter columns in the datafile (datname)
cos_filternames = ['B', 'G', 'I', 'IA427', 'IA484', 'IA505', 'IA527', 'IA624', 'IA709', 'IA738', 'R', 'U', 'V', 'Rp',
                   'Z', 'Zp', 'Hl', 'Hs', 'J1', 'J2', 'J3', 'Ks', 'NB118', 'NB209', 'F125W', 'F140W', 'F160W', 'F606W',
                   'F814W', 'UVISTA_J', 'UVISTA_H', 'UVISTA_Ks', 'UVISTA_Y', 'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']

# Names of filter response curve files (in the Filters folder)
mega_pre = 'megaprime-cfht_mega_'
mega = [folder+mega_pre+n for n in ['u_cfh9301', 'g_cfh9401', 'r_cfh9601', 'i_cfh9701', 'z_cfh9801']]  # U G R I Z
mega_1 = [mega[1], mega[3]]  # G, I
mega_2 = [mega[2], mega[0]]  # R, U
mega_3 = [mega[4]]  # Z
capak_pre = 'CAPAK_v2-'
capak_suf = '_subaru'
capak = [folder+capak_pre+n+capak_suf for n in ['V', 'B', 'r', 'z']]  # V, B, Rp, Zp
capak_0 = [capak[0], capak[2]]  # V, Rp
capak_1 = [capak[1]]  # B
capak_2 = [capak[3]]  # Zp
sub = [folder+subaru_pre+n for n in ['427', '484', '505', '527', '624', '709', '738']]
fstar = [folder+fourstar_pre+n+fourstar_suf for n in ['Hlong', 'Hshort', 'J1', 'J2', 'J3', 'Ks']]
hst_pre = 'hst-wfc3-'
hst = [folder+hst_pre+n for n in ['IR-f125w', 'IR-f140w', 'IR-f160w', 'UVIS-f606w', 'UVIS-f814w']]
vista_pre = 'VISTA-'
vista_suf = '_system+atmos'
vista = [folder+vista_pre+n+vista_suf for n in ['J', 'H', 'Ks', 'Y']]

cos_filts = capak_1 + mega_1 + sub + mega_2 + capak_0 + mega_3 + capak_2 + fstar + nb + hst + vista + irac


# UDS
# Names of filter columns in the datafile (datname)
uds_filternames = ['u', 'B', 'V', 'R', 'i', 'z', 'J1', 'J2', 'J3', 'Hs', 'Hl', 'Ks', 'J', 'H', 'K', 'KsHI', 'F125W',
                   'F140W', 'F160W', 'F606W', 'F814W', 'Y', 'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']

# Names of filter response curve files (in the Filters folder)
mega_u = ['megaprime-cfht_mega_u_cfh9301']
ukidss_pre = 'UKIDSS-'
ukidss_suf = '_qe'
ukidss = [ukidss_pre+n+ukidss_suf for n in ['B', 'R', 'i', 'z']]
ukidss_1 = [ukidss[0]]  # B
ukidss_2 = [ukidss[1], ukidss[2], ukidss[3]]  # R, i, z
subaru_V = ['COSMOS-SUBARU_filter_V']
fstar = [folder+fourstar_pre+n+fourstar_suf for n in ['J1', 'J2', 'J3', 'Hshort', 'Hlong', 'Ks']]  # DIFF ORDER
jhk_pre = 'UKIDSS-Table0'
jhk_suf = '_online'
jhk = [jhk_pre+n+jhk_suf for n in ['4', '5', '6']]  # J, H, K
hst_pre = 'hst-wfc3-'
hst = [folder+hst_pre+n for n in ['IR-f125w', 'IR-f140w', 'IR-f160w', 'UVIS-f606w']]
f814 = ['hst-ACS_update_sep07-wfc_f814w_t81']
vlt_y = ['VLT-hawki_y_ETC']

uds_filts = mega_u + ukidss_1 + subaru_V + ukidss_2 + fstar + jhk + KsHI + hst + f814 + vlt_y + irac


############
# OBS
#############


def get_names(field):
    if field == 'cosmos':
        photname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'
        zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'
        fname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.fout'
        # cosmos.v1.3.6_EL_Z004.awk.fout, cosmos.v1.3.6_EL_Z02.awk.fout
        filternames = cos_filternames
        filts = cos_filts

    elif field == 'cdfs':
        photname = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
        zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
        fname = '/home/jonathan/cdfs/cdfs.v1.6.9_EL_Z004.awk.fout'  # cdfs.v1.6.9.awk.fout'
        # cdfs.v1.6.9_EL_Z004.awk.fout, cosmos.v1.6.9_EL_Z02.awk.fout
        filternames = cdfs_filternames
        filts = cdfs_filts

    elif field == 'uds':
        photname = '/home/jonathan/uds/uds.v1.5.10.cat'
        zname = '/home/jonathan/uds/uds.v1.5.8.awk.zout'
        fname = '/home/jonathan/uds/uds.v1.5.8.awk.fout'
        # uds.v1.5.8_EL_Z004.awk.fout, uds.v1.5.8_EL_Z02.awk.fout
        filternames = uds_filternames
        filts = uds_filts

    return photname, zname, fname, filternames, filts


def load_obs(field, objname, err_floor=0.05, zperr=True, **extras):
    """
    photname: photometric file location
    objname: number of object in the 3D-HST COSMOS photometric catalog
    err_floor: the fractional error floor (0.05 = 5% floor)
    zp_err: inflate the errors by the zeropoint offsets from Skelton+14
    """

    photname, zname, fname, filternames, filts = get_names(field)

    # OPEN FILE, LOAD DATA
    with open(photname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(photname, comments='#', delimiter=' ', dtype=dtype)

    # EXTRACT FILTERS, FLUXES, ERRORS FOR OBJECT
    obj_idx = (dat['id'] == objname)
    # print(dat[obj_idx]['id'], 'idx')

    filters = np.array(filts)  # [f[2:] for f in dat.dtype.names if f[0:2] == 'f_'])
    # print(filters)
    flux = np.squeeze([dat[obj_idx]['f_' + f] for f in filternames])
    unc = np.squeeze([dat[obj_idx]['e_' + f] for f in filternames])

    # DEFINE PHOTOMETRIC MASK< CONVERT TO MAGGIES
    phot_mask = (flux != -99.0)
    # print('mask', phot_mask)
    maggies = flux * 10**-6 / 3631  # flux [uJy] * 1e-6 [Jy / uJy] * 1 [maggy] / 3631 [Jy]
    maggies_unc = unc * 10**-6 / 3631
    # print(maggies, 'maggies')
    # print(flux, 'flux')
    # print(maggies_unc, 'maggies_unc')
    # print(unc, 'unc')

    # ERROR FLOOR
    maggies_unc = np.clip(maggies_unc, maggies * err_floor, np.inf)  # for any unc < err_floor, replace with err_floor

    # BUILD OUTPUT DICTIONARY
    obs = {}
    obs['filters'] = observate.load_filters(filters)

    # NEW, MASK POINTS AROUND LY-ALPHA
    with open(zname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype)
    zred = zout['z_spec'][obj_idx][0]  # use z_spec
    if zred == -99:  # if z_spec doesn't exist
        zred = zout['z_peak'][obj_idx][0]  # use z_phot
    wave_eff = np.array([filt.wave_effective for filt in obs['filters']])

    width = np.array([filt.effective_width for filt in obs['filters']])  # get effective width of each filter
    wave_range = [None] * len(wave_eff)
    for i in range(len(wave_eff)):
        wave_range[i] = wave_eff[i] - width[i]/2, wave_eff[i] + width[i]/2  # store filter ranges

    range_rest = [None] * len(wave_eff)  # convert filter ranges to rest frame
    wave_rest = []  # convert normal filter effective center points to rest frame (this is not new)
    for i in range(len(wave_eff)):
        range_rest[i] = (wave_range[i][0] / (1 + zred)), (wave_range[i][1] / (1+zred))
        wave_rest.append(wave_eff[i] / (1 + zred))

    ly_mask = []
    for i in range(len(range_rest)):
        ly_mask.append((range_rest[i][0] < 1216) or ((range_rest[i][0] < 1180) & (range_rest[i][1] > 1180))
                       or ((range_rest[i][0] < 1260) & (range_rest[i][1] > 1260)))
    # print('mask', ly_mask)  # mask out ly_alpha and all points blueward of 1216
    for i in range(len(phot_mask)):
        if ly_mask[i]:
            phot_mask[i] = False

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
                     'isfree': False,  # fixedmet
                     'init': 0.0,  # noelg, trying -0.7, was 0.0. Maybe try -0.4? -0.5? -0.3?
                     'init_disp': 0.4,
                     'log_param': True,
                     'units': r'$\log (Z/Z_\odot)$',
                     'prior_function': tophat,
                     'prior_args': {'mini': -2.0, 'maxi': 0.19}})  # 'mini': -1.5

###### SFH ########
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
                     'prior_args': {'mini': 0.1, 'maxi': 100.0}})

model_params.append({'name': 'logtau', 'N': 1,
                     'isfree': False,  # NEW turn off
                     'init': 1,
                     'init_disp': 0.5,
                     'units': 'Gyr',
                     'prior_function': tophat,
                     'prior_args': {'mini': -1, 'maxi': 2}})

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
                     'prior_args': {'mini': 0.0, 'maxi': 14.0}})

model_params.append({'name': 'agebins', 'N': 1,  # NEW
                     'isfree': False,
                     'init': [],
                     'units': 'log(yr)',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini': 0.1, 'maxi': 15.0}})

model_params.append({'name': 'sfr_fraction', 'N': 1,  # NEW
                     'isfree': True,
                     'init': [],
                     'units': 'Msun',
                     'prior_function': priors.tophat,
                     'prior_args': {'mini': 0.0, 'maxi': 1.0}})

######## IMF ##############
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
                     'prior_args': {'mini': 0.0, 'maxi': 4.0}})

###### Dust Emission ##############
model_params.append({'name': 'add_dust_emission', 'N': 1,
                     'isfree': False,
                     'init': 1,  # init: 1 = add dust; init: 0 = no dust
                     'units': None,
                     'prior_function': None,
                     'prior_args': None})

###### Nebular Emission ###########
model_params.append({'name': 'add_neb_emission', 'N': 1,
                     'isfree': False,
                     'init': True,  # BUCKET1 emission lines --> init: True
                     'units': r'log Z/Z_\odot',
                     'prior_function_name': None,
                     'prior_args': None})

model_params.append({'name': 'gas_logz', 'N': 1,
                     'isfree': False,  # DECOUPLE: True (False when coupled) (False for non-EELGs)
                     'init': 0.0,
                     # 'depends_on': tie_gas_logz,  # BUCKET1 emission lines --> tie_gas_logz  # DECOUPLE --> remove
                     'units': r'log Z/Z_\odot',
                     'prior_function': tophat,
                     'prior_args': {'mini': -2.0, 'maxi': 0.5}})

model_params.append({'name': 'gas_logu', 'N': 1,
                     'isfree': False,  # BUCKET1 emission lines --> isfree: True (OR False because not important)
                     'init': -2.0,
                     'units': '',
                     'prior_function': tophat,
                     'prior_args': {'mini': -4, 'maxi': -1}})

model_params.append({'name': 'add_neb_continuum', 'N': 1,  # BUCKET1
                     'isfree': False,
                     'init': True,
                     'units': '',
                     'prior_function_name': None,
                     'prior_args': None})

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
# TURN ON DUST EMISSION FOR LOW-z
# mstar = stellar mass; to convert to SFH, need mass formed (requires calling fsps to convert)


# RE-SORT LIST OF PARAMETERS SO THAT MAJOR ONES ARE FIT FIRST
parnames = [m['name'] for m in model_params]
# fit_order = ['logmass', 'tage', 'logtau', 'dust2']  # for FAST mimic
fit_order = ['logmass', 'sfr_fraction', 'dust2']
tparams = [model_params[parnames.index(i)] for i in fit_order]
for param in model_params:
    if param['name'] not in fit_order:
        tparams.append(param)
model_params = tparams


# REDEFINE MODEL FOR MY OWN NEFARIOUS PURPOSES
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
        # print(self.theta_index['sfr_fraction'], 'theta_index')
        # print(theta, 'theta')
        # print(theta[0], 'sfr_theta_index')
        if 'sfr_fraction' in self.theta_index:
            sfr_fraction = theta[self.theta_index['sfr_fraction']]
            if np.sum(sfr_fraction) > 1.0:
                return -np.inf

        for k, v in self.theta_index.iteritems():
            this_prior = np.sum(self._config_dict[k]['prior_function']
                                (theta[v], **self._config_dict[k]['prior_args']))

            if not np.isfinite(this_prior):
                print('WARNING: ' + k + ' is out of bounds')
            lnp_prior += this_prior
        return lnp_prior


class FracSFH(FastStepBasis):  # NEW CLASS
    def get_galaxy_spectrum(self, **params):
        self.update(**params)

        # here's the custom fractional stuff
        fractions = np.array(self.params['sfr_fraction'])
        bin_fractions = np.append(fractions, (1 - np.sum(fractions)))
        time_per_bin = []
        for (t1, t2) in self.params['agebins']:
            time_per_bin.append(10 ** t2 - 10 ** t1)
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


# def load_model(objname='', datname='', zname='', agelims=[], **extras):
def load_model(objname, field, agelims=[], **extras):
    # REDSHIFT
    # open file, load data

    photname, zname, fname, filtername, filts = get_names(field)

    with open(photname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(photname, comments='#', delimiter=' ', dtype=dtype)

    with open(zname, 'r') as fz:
        hdr_z = fz.readline().split()
    dtype_z = np.dtype([(hdr_z[1], 'S20')] + [(n, np.float) for n in hdr_z[2:]])
    zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype_z)

    idx = dat['id'] == objname  # creates T/F array: True when dat[id] = objname
    zred = zout['z_spec'][idx][0]  # use z_spec
    if zred == -99:  # if z_spec doesn't exist
        zred = zout['z_peak'][idx][0]  # use z_phot
    print(zred, 'zred')

    # NEW FOUT (MASS-METALLICITY)
    with open(fname, 'r') as ff:
        hdr_f = ff.readline().split()
    dtype_f = np.dtype([(hdr_f[1], 'S20')] + [(n, np.float) for n in hdr_f[2:]])
    fout = np.loadtxt(fname, comments='#', delimiter=' ', dtype=dtype_f)
    idx_f = fout['id'] == objname  # creates T/F array: True when fout[id] = objname
    lmass = fout['lmass'][idx_f][0]
    print(lmass, 'lmass')
    # met = [-1.5, 0.5, -0.3]
    if lmass <= 9.7:
        met = [-1.0, 0.0, -0.6]
    elif 9.7 < lmass <= 10.0:
        met = [-0.9, 0.1, -0.4]
    elif 10.0 < lmass <= 10.3:
        met = [-0.8, 0.15, -0.2]
    elif lmass > 10.3:
        met = [-0.2, 0.3, 0.0]

    # CALCULATE AGE OF THE UNIVERSE (TUNIV) AT REDSHIFT ZRED
    tuniv = WMAP9.age(zred).value
    print(tuniv, 'tuniv')

    n = [p['name'] for p in model_params]
    # model_params[n.index('tage')]['prior_args']['maxi'] = tuniv

    # NEW FOUT MASS-METALLICITY
    model_params[n.index('logzsol')]['prior_args'] = {'mini': met[0], 'maxi': met[1]}  # {'mini': -0.8, 'maxi': 0.15}
    model_params[n.index('logzsol')]['init'] = met[2]  # -0.25
    # 4942_cosmos_noelgmet was run with met[2] = -0.4
    # 4942_cosmos_noelg25 was run with met[2] = -0.25

    # NONPARAMETRIC SFH
    agelims = [0.0, 7.48, 8.0, 8.48, 9.0, (9.0 + (np.log10(tuniv*1e9) - 9.0)/2), np.log10(tuniv*1e9)]
    # 0, 30 Myr, 100 Myr, 300 Myr, 1 Gyr, two evenly spaced to tuniv
    ncomp = len(agelims) - 1
    agebins = np.array([agelims[:-1], agelims[1:]])  # why agelims[1:] instead of agelims[0:]?

    # INSERT REDSHIFT INTO MODEL PARAMETER DICTIONARY
    zind = n.index('zred')
    model_params[zind]['init'] = zred

    # SET UP AGEBINS
    model_params[n.index('agebins')]['N'] = ncomp
    model_params[n.index('agebins')]['init'] = agebins.T

    # FRACTIONAL MASS INITIALIZATION  # NEW
    # N-1 bins, last is set by x = 1 - np.sum(sfr_fraction)
    model_params[n.index('sfr_fraction')]['N'] = ncomp-1
    model_params[n.index('sfr_fraction')]['prior_args'] = {
                                                           'maxi': np.full(ncomp-1, 1.0),
                                                           'mini': np.full(ncomp-1, 0.0),
                                                           # NOTE: ncomp instead of ncomp-1 makes the prior take into
                                                           # account the implicit Nth variable too
                                                          }
    model_params[n.index('sfr_fraction')]['init'] = np.zeros(ncomp-1)+1./ncomp
    model_params[n.index('sfr_fraction')]['init_disp'] = 0.02

    # CREATE MODEL
    model = BurstyModel(model_params)

    return model


model_type = BurstyModel

'''
RUNNING WITH
mpirun -n 4 python prospector.py --param_file=noelg_otherbins_params.py --outfile=1969_uds_nother --niter=1200
--field=uds --objname=1969
'''