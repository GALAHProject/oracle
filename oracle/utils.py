# coding: utf-8

""" General utilities """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["atomic_number", "element", "reflect_about", "latexify", "readable_dict"]

import collections
import logging
from functools import update_wrapper
from threading import RLock

import numpy as np
from scipy.special import wofz

logger = logging.getLogger("oracle")


def unpack_atomic_transition(transition, **defaults):

    _defaults = {
        "synthesise_surrounding": 1.0,
        "opacity_contribution": 1.0,
        "damp1": 0,
        "damp2": 0
    }
    _defaults.update(defaults)

    if isinstance(transition, dict):
        wavelength, species, excitation_potential, loggf = [transition[k] \
            for k in ("wavelength", "species", "excitation_potential", "loggf")]

        damp1 = transition.get("damp1", _defaults["damp1"])
        damp2 = transition.get("damp2", _defaults["damp2"])
        opacity_contribution = transition.get("opacity_contribution",
            _defaults["opacity_contribution"])
        synthesise_surrounding = transition.get("synthesise_surrounding",
            _defaults["synthesise_surrounding"])

    else:

        damp1, damp2 = _defaults["damp1"], _defaults["damp2"]
        opacity_contribution = _defaults["opacity_contribution"]
        synthesise_surrounding = _defaults["synthesise_surrounding"]

        wavelength, species, excitation_potential, loggf = transition[:4]
        if len(transition) > 4:
            damp1 = transition[4]

            if len(transition) > 5:
                damp2 = transition[5]

                if len(transition) > 6:
                    synthesise_surrounding = transition[6]

                    if len(transition) > 7:
                        synthesise_surrounding = transition[7]

    return (wavelength, species, excitation_potential, loggf, damp1, damp2,
        synthesise_surrounding, opacity_contribution)



# This shouldn't be here, but I don't know where to put it yet.
def estimate_microturbulence(effective_temperature, surface_gravity):
    """
    Estimate microtubulence from relations between effective temperature and
    surface gravity. For giants (logg < 3.5) the relationship employed is from
    Kirby et al. (2008, ) and for dwarfs (logg >= 3.5) the Reddy et al. (2003, 
        ) relation is used.

    :param effective_temperature:
        The effective temperature of the star in Kelvin.

    :type effective_temperature:
        float

    :param surface_gravity:
        The surface gravity of the star.

    :type surface_gravity:
        float

    :returns:
        The estimated microturbulence (in km/s) from the given stellar parameters.

    :rtype:
        float
    """

    if surface_gravity >= 3.5:
        return 1.28 + 3.3e-4 * (effective_temperature - 6000) \
            - 0.64 * (surface_gravity - 4.5)
    else:
        return 2.70 - 0.509 * surface_gravity


def readable_dict(keys, values):
    return ", ".join(["{0} = {1:.2f}".format(k, v) for k, v in zip(keys, values)])


def update_recursively(original, new):
    """
    Recursively update a nested dictionary.
    
    :param original:
        The original nested dictionary to update.

    :type original:
        dict

    :param new:
        The nested dictionary to use to update the original.

    :type new:
        dict

    :returns:
        The updated original dictionary.

    :rtype:
        dict
    """

    for k, v in new.iteritems():
        if isinstance(v, collections.Mapping) \
        and isinstance(original.get(k, None), collections.Mapping):
            r = update_recursively(original.get(k, {}), v)
            original[k] = r
        else:
            original[k] = new[k]
    return original


def invert_mask(mask, data=None, limits=(0, 10e5), padding=0):
    """
    Invert a mask to return regions that should be considered. If ``data`` are
    provided then they will instruct on the wavelenth limits to be considered.
    Otherwise, ``limits`` will be used.

    :param mask:
        The regions to mask.

    :type mask:
        :class:`numpy.ndarray`

    :param data: [optional]
        The observed spectra.

    :type data:
        list of :class:`oracle.specutils.Spectrum1D` objects

    :param limits: [optional]
        The wavelength limits to consider.

    :type limits:
        tuple
    """

    if data is not None:
        limits = np.sort(np.array([[s.disp[0], s.disp[-1]] for s in data]).flatten())
        limits = [np.min(limits), np.max(limits)]

    points = []
    points.extend(limits)
    points.extend(mask.flatten())
    points = np.sort(points)

    # Find the start of the limts
    i = points.searchsorted(limits)
    inverted_mask = points[i[0] + (i[0] % 2):i[1] + (i[1] % 2)].reshape(-1, 2)
    inverted_mask[:, 0] -= padding
    inverted_mask[:, 1] += padding
    return inverted_mask



_CacheInfo = collections.namedtuple("CacheInfo",
    ["hits", "misses", "maxsize", "currsize"])

class _HashedSeq(list):
    __slots__ = 'hashvalue'

    def __init__(self, tup, hash=hash):
        self[:] = tup
        self.hashvalue = hash(tup)

    def __hash__(self):
        return self.hashvalue


def _make_key(args, kwds, typed,
             kwd_mark = (object(),),
             fasttypes = {int, str, frozenset, type(None)},
             sorted=sorted, tuple=tuple, type=type, len=len):
    'Make a cache key from optionally typed positional and keyword arguments'
    key = args
    if kwds:
        sorted_items = sorted(kwds.items())
        key += kwd_mark
        for item in sorted_items:
            key += item
    if typed:
        key += tuple(type(v) for v in args)
        if kwds:
            key += tuple(type(v) for k, v in sorted_items)
    elif len(key) == 1 and type(key[0]) in fasttypes:
        return key[0]
    return _HashedSeq(key)


def lru_cache(maxsize=128, typed=False):
    """Least-recently-used cache decorator.

    If *maxsize* is set to None, the LRU features are disabled and the cache
    can grow without bound.

    If *typed* is True, arguments of different types will be cached separately.
    For example, f(3.0) and f(3) will be treated as distinct calls with
    distinct results.

    Arguments to the cached function must be hashable.

    View the cache statistics named tuple (hits, misses, maxsize, currsize) with
    f.cache_info().  Clear the cache and statistics with f.cache_clear().
    Access the underlying function with f.__wrapped__.

    See:  http://en.wikipedia.org/wiki/Cache_algorithms#Least_Recently_Used

    """

    # Users should only access the lru_cache through its public API:
    #       cache_info, cache_clear, and f.__wrapped__
    # The internals of the lru_cache are encapsulated for thread safety and
    # to allow the implementation to change (including a possible C version).

    def decorating_function(user_function):

        cache = dict()
        stats = [0, 0]                  # make statistics updateable non-locally
        HITS, MISSES = 0, 1             # names for the stats fields
        make_key = _make_key
        cache_get = cache.get           # bound method to lookup key or return None
        _len = len                      # localize the global len() function
        lock = RLock()                  # because linkedlist updates aren't threadsafe
        root = []                       # root of the circular doubly linked list
        root[:] = [root, root, None, None]      # initialize by pointing to self
        nonlocal_root = [root]                  # make updateable non-locally
        PREV, NEXT, KEY, RESULT = 0, 1, 2, 3    # names for the link fields

        if maxsize == 0:

            def wrapper(*args, **kwds):
                # no caching, just do a statistics update after a successful call
                result = user_function(*args, **kwds)
                stats[MISSES] += 1
                return result

        elif maxsize is None:

            def wrapper(*args, **kwds):
                # simple caching without ordering or size limit
                key = make_key(args, kwds, typed)
                result = cache_get(key, root)   # root used here as a unique not-found sentinel
                if result is not root:
                    stats[HITS] += 1
                    return result
                result = user_function(*args, **kwds)
                cache[key] = result
                stats[MISSES] += 1
                return result

        else:

            def wrapper(*args, **kwds):
                # size limited caching that tracks accesses by recency
                key = make_key(args, kwds, typed) if kwds or typed else args
                with lock:
                    link = cache_get(key)
                    if link is not None:
                        # record recent use of the key by moving it to the front of the list
                        root, = nonlocal_root
                        link_prev, link_next, key, result = link
                        link_prev[NEXT] = link_next
                        link_next[PREV] = link_prev
                        last = root[PREV]
                        last[NEXT] = root[PREV] = link
                        link[PREV] = last
                        link[NEXT] = root
                        stats[HITS] += 1
                        logger.debug("Using cache")
                        return result
                result = user_function(*args, **kwds)
                with lock:
                    root, = nonlocal_root
                    if key in cache:
                        # getting here means that this same key was added to the
                        # cache while the lock was released.  since the link
                        # update is already done, we need only return the
                        # computed result and update the count of misses.
                        pass
                    elif _len(cache) >= maxsize:
                        # use the old root to store the new key and result
                        oldroot = root
                        oldroot[KEY] = key
                        oldroot[RESULT] = result
                        # empty the oldest link and make it the new root
                        root = nonlocal_root[0] = oldroot[NEXT]
                        oldkey = root[KEY]
                        oldvalue = root[RESULT]
                        root[KEY] = root[RESULT] = None
                        # now update the cache dictionary for the new links
                        del cache[oldkey]
                        cache[key] = oldroot
                    else:
                        # put result in a new link at the front of the list
                        last = root[PREV]
                        link = [last, root, key, result]
                        last[NEXT] = root[PREV] = cache[key] = link
                    stats[MISSES] += 1
                return result

        def cache_info():
            """Report cache statistics"""
            with lock:
                return _CacheInfo(stats[HITS], stats[MISSES], maxsize, len(cache))

        def cache_clear():
            """Clear the cache and cache statistics"""
            with lock:
                cache.clear()
                root = nonlocal_root[0]
                root[:] = [root, root, None, None]
                stats[:] = [0, 0]

        wrapper.__wrapped__ = user_function
        wrapper.cache_info = cache_info
        wrapper.cache_clear = cache_clear
        return update_wrapper(wrapper, user_function)

    return decorating_function


def rounder(*decimals, **decimal_kwargs):
    def decorator(function):
        def wrapper(*args, **kwargs):
            
            rounded_args = []
            for arg, decimal in zip(args, decimals):
                if decimal is None:
                    if isinstance(arg, (np.core.records.record, )):
                        rounded_args.append(tuple(arg))
                    else:
                        rounded_args.append(arg)
                else:
                    if isinstance(arg, (float, int)):
                        rounded_args.append(np.round(arg, decimal))
                    else:
                        rounded_args.append(tuple(np.round(arg, decimal)))


            # Extend with arguments that don't have precision requirements
            missing_args = len(args) - len(rounded_args)
            if missing_args > 0:
                rounded_args.extend(args[-missing_args:])
            return function(*rounded_args, **kwargs)
        return wrapper
    return decorator



#        return ("teff", "xi", "logg", "[M/H]")
# teff, logg, xi, [M/H]

def stellar_jacobian(stellar_parameters, *args):
    """ Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations from the Sun """

    logger.info("Updating approximation of Jacobian matrix for stellar parameter determination")

    """"
    teff, vt, logg, feh = stellar_parameters[:4]
    original_jacobian = np.array([
        [ 5.4393e-08*teff - 4.8623e-04, -7.2560e-02*vt + 1.2853e-01,  1.6258e-02*logg - 8.2654e-02,  1.0897e-02*feh - 2.3837e-02],
        [ 4.2613e-08*teff - 4.2039e-04, -4.3985e-01*vt + 8.0592e-02, -5.7948e-02*logg - 1.2402e-01, -1.1533e-01*feh - 9.2341e-02],
        [-3.2710e-07*teff + 2.8178e-03,  3.8185e-02*vt - 1.6601e-01, -1.2006e-01*logg - 3.5816e-02, -2.8592e-04*feh + 1.4257e-02],
        [-1.7822e-07*teff + 1.8250e-03,  3.5564e-01*vt - 1.1024e-00, -1.2114e-01*logg + 4.1779e-01, -1.8847e-01*feh - 1.0949e-00]
    ])
    """

    teff, logg, feh, vt = stellar_parameters[:4]
    # This is total black magic. Like, wizard style.
    return np.array([
        [+5.4393e-08*teff - 4.8623e-04, +1.6258e-02*logg - 8.2654e-02, +1.0897e-02*feh - 2.3837e-02, -7.2560e-02*vt + 1.2853e-01],
        [+4.2613e-08*teff - 4.2039e-04, -5.7948e-02*logg - 1.2402e-01, -1.1533e-01*feh - 9.2341e-02, -4.3985e-01*vt + 8.0592e-02],
        [-3.2710e-07*teff + 2.8178e-03, -1.2006e-01*logg - 3.5816e-02, -2.8592e-04*feh + 1.4257e-02, +3.8185e-02*vt - 1.6601e-01],
        [-1.7822e-07*teff + 1.8250e-03, -1.2114e-01*logg + 4.1779e-01, -1.8847e-01*feh - 1.0949e-00, +3.5564e-01*vt - 1.1024e-00]
    ]).T


def reflect_about(a, limits):
    """
    Similar to :func:`numpy.clip`, except it just reflects about some limiting axes.

    :param a:
        The array of values to reflect.

    :type a:
        :class:`numpy.array`

    :param limits:
        The upper and lower limits to reflect about. Use ``None`` for no limit.

    :type limits:
        A two length tuple or list-type.

    :returns:
        The reflected array.

    :rtype:
        :class:`numpy.array`
    """

    lower, upper = limits
    if lower is not None:
        a[a < lower] = lower + (lower - a[a < lower])
    if upper is not None:
        a[a > upper] = upper - (a[a > upper] - upper)
    return a


def latexify(labels, default_latex_labels=None):
    """
    Return a LaTeX-ified label.

    Args:
        labels (str or list-type of str objects): The label(s) to latexify. 
        default_latex_labels (dict): Dictionary of common labels to use.

    Returns:
        LaTeX-ified label.
    """

    common_labels = {
        "teff": "$T_{\\rm eff}$ (K)",
        "feh": "[Fe/H]",
        "logg": "$\log{g}$",
        "alpha": "[$\\alpha$/Fe]",
        "xi": "$\\xi$ (km s$^{-1}$)"
    }

    if default_latex_labels is not None:
        common_labels.update(default_latex_labels)
    
    listify = True
    if isinstance(labels, str):
        listify = False
        labels = [labels]

    latex_labels = []
    for label in labels:

        if label.startswith("doppler_sigma_"):
            color = ["blue", "green", "red", "ir"][int(label.split("_")[-1])]
            latex_labels.append("$\sigma_{\\rm doppler," + color + "}$ ($\\AA{}$)")

        else:
            latex_labels.append(common_labels.get(label, label))

    if not listify:
        return latex_labels[0]

    return latex_labels


def atomic_number(element):
    """
    Return the atomic number of a given element.

    :param element:
        The short-hand notation for the element (e.g., Fe).

    :type element:
        str

    :returns:
        The atomic number for a given element.

    :rtype:
        int
    """
    
    if not isinstance(element, (unicode, str)):
        raise TypeError("element must be represented by a string-type")

    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    
    if element not in periodic_table:
        return ValueError("element '{0}' is not known".format(element))

    return periodic_table.index(element) + 1


def element(atomic_number):
    """
    Return the element of a given atomic number.

    :param atomic_number:
        The atomic number for the element in question (e.g., 26).

    :type atomic_number:
        int-like

    :returns:
        The short-hand element for a given atomic number.

    :rtype:
        str
    """

    atomic_number = int(atomic_number)
    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    return periodic_table[atomic_number - 1]