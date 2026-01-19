import numpy as np

def spherical_center(ra_deg, dec_deg):
    """
    Compute a centre on the celestial sphere by averaging unit vectors.

    Parameters
    ----------
    ra_deg, dec_deg : float, float
        The RA and Dec in degrees

    Returns
    -------
    ra0, dec0 : float, float
        The centre RA and Dec in degrees
    """
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    xm, ym, zm = np.mean(x), np.mean(y), np.mean(z)
    r = np.sqrt(xm**2 + ym**2 + zm**2)
    xm, ym, zm = xm/r, ym/r, zm/r

    dec0 = np.arcsin(zm)
    ra0 = np.arctan2(ym, xm) % (2*np.pi)

    return np.rad2deg(ra0), np.rad2deg(dec0)

def angular_separation(ra_deg, dec_deg, ra0_deg, dec0_deg):
    """
    Compute the Great-circle angular separation in radians.

    Parameters
    ----------
    ra_deg, dec_deg : iterable, iterable
        RA and Dec values in degrees of an iterable of stars

    ra0_deg, dec0_deg : float, float
        RA and Dec values of the cluster centre

    Returns
    -------
    ang_sep : float
        The angular separation in radians
    """
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    ra0 = np.deg2rad(ra0_deg)
    dec0 = np.deg2rad(dec0_deg)

    cosang = np.sin(dec)*np.sin(dec0) + np.cos(dec)*np.cos(dec0)*np.cos(ra - ra0)
    cosang = np.clip(cosang, -1.0, 1.0)
    ang_sep = np.arccos(cosang)
    return ang_sep


### CLUSTER MODELS #####################################################


def sigma_king_cluster(R, Sigma0, rc, rt):
    """
    A King profile for the surface density of the cluster.

    Parameters
    ----------
    R : numpy.ndarray
        A array of radius values for which to compute the surface profile

    Sigma0 : float
        The central surface density

    r_c, r_t, float, float
        The core radius and tidal radius

    Returns
    -------
    Sigma : numpy.ndarray
        An array with the surface density at the radii R
    """
    term_t = 1.0 / np.sqrt(1.0 + (rt/rc)**2)
    term_R = 1.0 / np.sqrt(1.0 + (R/rc)**2)
    Sigma = Sigma0 * (term_R - term_t)**2
    Sigma = np.where(R <= rt, Sigma, 0.0)
    return Sigma

def sigma_exp_cluster(R, Sigma0, Rs):
    """
    Expomential model for the surface density 
    of the cluster

    Parameters
    -----------
    R : numpy.ndarray
        A array of radius values for which to compute the surface profile

    Sigma0 : float
        The central surface density

    r_s : float
        The exponential scale of the model

     Returns
    -------
    Sigma : numpy.ndarray
        An array with the surface density at the radii R   
    """
    return Sigma0 * np.exp(-R / Rs)

### PRIORS ####################################################################

def log_halfnormal(x, s):
    """
    Computes the log-pdf of a half-normal 
    probability distribution, i.e. a distribution 
    that's 0 below 0 and Gaussian above.

    Parameters
    ----------
    x : np.ndarray or float
        An array or scalar with the value(s) for which to 
        compute the prior pdf

    s : float
        The scale parameter of the distribution

    Returns
    -------
    pdf : float or np.ndarray
        the log-pdf for the values in x
    """
    if x < 0:
        return -np.inf
    # HalfNormal density: sqrt(2/pi)/s * exp(-x^2/(2s^2))
    return float(np.log(np.sqrt(2/np.pi)) - np.log(s) - 0.5*(x/s)**2)

def log_uniform(x, a, b):
    """
    Computes the log-pdf for a uniform distribution
    between bounds a and b, outside of which the 
    pdf will be zero (or the log-pdf -inf)

    Parameters
    ----------
    x : np.ndarray or float
        An array or scalar with the value(s) for which to 
        compute the prior pdf

    a, b : float, float
        The lower and upper boundary for the distribution

    Returns
    -------
    pdf : float or np.ndarray
        the log-pdf for the values in x    
    """
    if (x < a) or (x > b):
        return -np.inf
    return float(np.log(b - a))

def log_gaussian(x, mu, sigma):
    """
    Log-pdf of a Normal (Gaussian) distribution.

    Parameters
    ----------
    x : float or np.ndarray
        Point(s) at which to evaluate the log-density.
    mu : float
        Mean of the Normal distribution.
    sigma : float
        Standard deviation of the Normal distribution. Must be > 0.

    Returns
    -------
    logp : float or np.ndarray
        Log-density evaluated at `x`. Returns `-np.inf` where parameters
        are invalid (e.g., sigma <= 0) or where evaluation is undefined.
    """
    x = np.asarray(x)
    if not np.isfinite(mu) or not np.isfinite(sigma) or sigma <= 0:
        return np.full_like(x, -np.inf, dtype=float) if x.ndim else -np.inf

    z = (x - mu) / sigma
    logp = -0.5 * (np.log(2.0 * np.pi) + 2.0 * np.log(sigma) + z**2)
    return logp


def log_student_t(x, nu, mu=0.0, sigma=1.0):
    """
    Log-pdf of a Student's t distribution.

    This is the distribution of:
        X = mu + sigma * T,
    where T ~ StudentT(nu) (standard t with df=nu).

    Parameters
    ----------
    x : float or np.ndarray
        Point(s) at which to evaluate the log-density.
    nu : float
        Degrees of freedom. Must be > 0.
    mu : float, optional
        Location parameter (default is 0.0).
    sigma : float, optional
        Scale parameter. Must be > 0 (default is 1.0).

    Returns
    -------
    logp : float or np.ndarray
        Log-density evaluated at `x`. Returns `-np.inf` where parameters
        are invalid (e.g., nu <= 0, sigma <= 0).
    """
    x = np.asarray(x)
    if (not np.isfinite(nu)) or (not np.isfinite(mu)) or (not np.isfinite(sigma)) or nu <= 0 or sigma <= 0:
        return np.full_like(x, -np.inf, dtype=float) if x.ndim else -np.inf

    z = (x - mu) / sigma
    # log C where C = Gamma((nu+1)/2) / (Gamma(nu/2) * sqrt(nu*pi) * sigma)
    logC = (
        gammaln((nu + 1.0) / 2.0)
        - gammaln(nu / 2.0)
        - 0.5 * (np.log(nu) + np.log(np.pi))
        - np.log(sigma)
    )
    logp = logC - 0.5 * (nu + 1.0) * np.log1p((z**2) / nu)
    return logp

def log_half_student_t(x, nu, sigma=1.0):
    """
    Log-pdf of a half-Student's t distribution (support x >= 0).

    Definition
    ----------
    If T ~ StudentT(nu, loc=0, scale=sigma), then X = |T| has the
    half-Student's t distribution. Its PDF is:
        f_X(x) = 2 * f_T(x) for x >= 0,
        f_X(x) = 0           for x < 0.

    Parameters
    ----------
    x : float or np.ndarray
        Point(s) at which to evaluate the log-density.
    nu : float
        Degrees of freedom. Must be > 0.
    sigma : float, optional
        Scale parameter. Must be > 0 (default is 1.0).

    Returns
    -------
    logp : float or np.ndarray
        Log-density evaluated at `x`. Returns `-np.inf` for x < 0 or
        invalid parameters.
    """
    x = np.asarray(x)

    if (not np.isfinite(nu)) or (not np.isfinite(sigma)) or nu <= 0 or sigma <= 0:
        return np.full_like(x, -np.inf, dtype=float) if x.ndim else -np.inf

    logp = np.full_like(x, -np.inf, dtype=float)
    mask = x >= 0
    if np.any(mask):
        logp_t = logpdf_student_t(x[mask], nu=nu, mu=0.0, sigma=sigma)
        logp[mask] = np.log(2.0) + logp_t

    return logp
