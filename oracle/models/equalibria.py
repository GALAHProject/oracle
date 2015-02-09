# coding: utf-8

from __future__ import absolute_import, print_function

""" Equialibria (excitation & ionisation balance) model for stellar spectra. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from scipy import stats, optimize as op

from oracle import atmospheres, specutils, synthesis, utils
from oracle.models.model import Model
from oracle.models import transitions
from oracle.models.jacobian \
    import approximate as stellar_parameter_jacobian_approximation

logger = logging.getLogger("oracle")

class EqualibriaModel(Model):

    """
    A class to forward model stellar spectra. This class performs excitation and
    ionization balances in order to determine stellar parameters.
    """

    # Default configuration for a EqualibriaModel class
    config = {
        "model": {
            "redshift": True,
            "instrumental_resolution": True,
            "continuum": False,
            "profile_function": "gaussian"
        },
        "settings": {
            "threads": 1
        }
    }
    _line_parameter_format = "A({species:.1f}@{wavelength:.2f})"


    def __init__(self, configuration):
        """
        Initialise an Equialibria Model class.
        
        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string 
            containing the configuration in a YAML-dump format. Only YAML-style
            formats are accepted.

        :type configuration:
            str or dict
        """

        super(EqualibriaModel, self).__init__(configuration)
        return None


    def parameters(self, data):
        """
        Return the model parameters for some data. The model configuration is
        applicable for a large number of observed channels, so the number of 
        observed channels is required to determine the exact number of model
        parameters.

        :param data:
            The observed data.

        :type num_data_channels:
            list of :class:`oracle.specutils.Spectrum1D` objects

        :returns:
            A list of model parameters.

        :rtype:
            tuple
        """

        parameters = ["effective_temperature", "surface_gravity", "[M/H]",
            "microturbulence"]

        # Single radial velocity for all channels
        if self.config["model"]["redshift"] == True:
            parameters.append("v_rad")

        # Different radial velocity for each channel?
        elif isinstance(self.config["model"]["redshift"], (tuple, list, )):
            parameters.extend(["v_rad.{}".format(i) for i, channel_v_rad in \
                zip(range(len(data)), self.config["model"]["redshift"])\
                if channel_v_rad])

        # Instrumental broadening
        if self.config["model"]["instrumental_resolution"]:
            parameters.extend(["instrumental_resolution.{}".format(i) \
                for i in range(len(data))])

        # Continuum treatment
        if isinstance(self.config["model"]["continuum"], (tuple, list)):
            # List contains order for each channel
            for i, order in \
                zip(range(len(data)), self.config["model"]["continuum"]):
                parameters.extend(["continuum.{0}.{1}".format(i, j) \
                    for j in range(order + 1)])

        # Atomic line abundances
        num_atomic_transitions = 0
        wavelength_ranges = [(each.disp[0], each.disp[-1]) for each in data]

        atomic_transitions = self.config["model"].get("atomic_transitions", None)
        if atomic_transitions is None:
            # If not specified, try and load simplistic information from a file
            atomic_transitions = np.loadtxt(self.config["model"]["atomic_transitions_filename"])

        for atomic_transition in atomic_transitions:

            if isinstance(atomic_transition, dict):
                wavelength = atomic_transition["wavelength"]
                species = atomic_transition["species"]

            else:
                # Assume list
                wavelength, species = atomic_transition[:2]

            # Ensure the wavelength is within the observed region
            for wavelength_start, wavelength_end in wavelength_ranges:
                if wavelength_end >= wavelength and wavelength >= wavelength_start:
                    parameters.append(self._line_parameter_format.format(
                        wavelength=wavelength, species=species))
                    num_atomic_transitions += 1

        if 2 > num_atomic_transitions:
            raise ValueError("less than two atomic lines provided for equalibria")

        return tuple(parameters)


    def initial_theta(self, data, full_output=False, **kwargs):
        """
        Return an initial guess of the model parameters theta using no prior
        information.

        :param data:
            The observed data.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects
        """

        # Get initial theta from the base model class, which will estimate the
        # radial velocities, continuum coefficients, and stellar parameters
        theta, r_chi_sq, model_dispersion, model_fluxes = super(EqualibriaModel,
            self).initial_theta(data, full_output=True, **kwargs)

        if full_output:
            return (theta, r_chi_sq, model_dispersion, model_fluxes)
        return theta


    def estimate_stellar_parameters(self, data, initial_theta=None,
        fitting_frequency=0, **kwargs):
        """
        Return point estimates for the stellar parameters (effective temperature,
        surface gravity, metallicity, microturbulence) using an excitation and
        ionisation balance approach.

        :param data:
            The observed data.

        :type data:
            list of `oracle.specutils.Spectrum1D` objects

        :param initial_theta: [optional]
            Initial estimates of the stellar parameters, radial velocities, and
            continuum parameters (where appropriate).

        :type initial_theta:
            dict

        :param fitting_frequency: [optional]
            This parameter specifies how frequently the profiles be re-fit
            and/or synthesised. If set to zero (default), then fitting is only
            performed once with the initial theta parameters.

        :type fitting_frequency:
            int

        :returns:
            Point estimates of the effective temperature, surface gravity,
            metallicity, and microturbulence.
        """

        if initial_theta is None:
            initial_theta = self.initial_theta(data)

        else:
            logger.warn("assert that all parameters are provided by initial theta")

        if fitting_frequency != 0:
            raise NotImplementedError("sorry")
            # [TODO]

        # Get the initial stellar parameters
        initial_stellar_parameters = [initial_theta[label] for label in \
            ("effective_temperature", "surface_gravity", "[M/H]", 
                "microturbulence")]

        kwds = kwargs.copy()
        kwds.update({
            "full_output": False,
            "initial_theta": initial_theta
        })

        # Get the initial state
        initial_state, transitions_table, atomic_transitions = \
            self.equalibrium_state(data, *initial_stellar_parameters,
                initial_theta=initial_theta, full_output=True)

        # Get usable_transitions
        use = transitions_table["is_filtered"].astype(bool) \
            * ~transitions_table["is_outlier"].astype(bool)
        neutral = (transitions_table["species"][use] % 1) == 0
        ionised = (transitions_table["species"][use] % 1) > 0
        reduced_equivalent_width = np.log(transitions_table["equivalent_width"]\
            /transitions_table["wavelength"])
        #transitions = transitions_table[use].view(float).reshape(use.sum(), -1)[:, :7]

        def state_function(x):
            print("STATE IN", x)

            # Calculate abundances given the stellar parameters
            try:
                abundances = synthesis.atomic_abundances(transitions_table[use],
                    x[:3], microturbulence=x[3])

            except ValueError:
                return np.array([np.nan, np.nan, np.nan, np.nan])

            # Calculate the state
            excitation_slope = stats.linregress(
                x=transitions_table["excitation_potential"][use][neutral],
                y=abundances[neutral])

            # Calculate the ionisation state
            ionisation_state = abundances[neutral].mean() \
                - abundances[ionised].mean()

            # Calculate the abundance state
            metallicity = x[2]
            #metallicity = mh
            abundance_state = (abundances \
                - (atmospheres.solar_abundance(transitions_table["species"][use])
                    + metallicity)).mean()

            # Slope with reduced equivalent width and line abundance
            line_strength_slope = stats.linregress(
                x=reduced_equivalent_width[use][neutral], y=abundances[neutral])

            # Calculate the state
            #state = np.array([
            #    excitation_slope[0],
            #    ionisation_state,
            #    abundance_state,
            #    line_strength_slope[0]
            #])


            # Re-arrange the state
            #teff, vt, logg, feh = stellar_parameters[:4]
            state = np.array([
                excitation_slope[0],
                ionisation_state,
                0.1 * abundance_state,
                line_strength_slope[0],
            ])
            
            if np.any(~np.isfinite(state)):
                raise WTFError()
            print("STATE OUT", x, state, (state**2).sum())
            return state

        # Re-arrange the initial stellar parameters
        initial_stellar_parameters = [initial_theta[label] for label in \
            ("effective_temperature",
                 "surface_gravity", "[M/H]","microturbulence", )]

        """
        # This is for a fitting_frequency == 1
        state_function = lambda x: self.equalibrium_state(data, *x, **kwds)

        # Optimise the state function
        result = op.fsolve(state_function, initial_stellar_parameters,
            fprime=stellar_parameter_jacobian_approximation, col_deriv=True,
            epsfcn=0, xtol=xtol, maxfev=maxfev, full_output=True)
        """

        # Optimisation
        xtol = kwargs.pop("xtol", 1e-10)
        maxfev = kwargs.pop("maxfev", 100)

        result = op.fsolve(state_function, initial_stellar_parameters,
            fprime=stellar_parameter_jacobian_approximation, col_deriv=1,
            epsfcn=0, xtol=xtol, maxfev=maxfev, full_output=True)

        raise a


        raise NotImplementedError

    def equalibrium_state(self, data, effective_temperature, surface_gravity,
        metallicity, microturbulence, initial_theta=None, full_output=False,
        **kwargs):
        """
        Return the equilibrium state information for a given set of stellar
        parameters. The equilibrium state information includes the:

        (1) slope of line abundance with the excitation potential
        (2) the mean difference between neutral and ionised lines
        (4) the difference of the input metallicty and mean output metallicity
        (3) slope of line abundance with reduced equivalent width
        
        :param data:
            The observed spectra.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects

        :param effective_temperature:
            The effective temperature to calculate the equilibrium state.

        :type effective_temperature:
            float

        :param surface_gravity:
            The surface gravity to calculate the equilibrium state.

        :type surface_gravity:
            float

        :param metallicity:
            The overall scaled-solar metallicity to calculate the equilibrium
            state.

        :type metallicity:
            float

        :param microturbulence:
            The photospheric microturbulence to calculate the equilibrium state.

        :type microturbulence:
            float

        :param full_output: [optional]
            Return all optional parameters.

        :type full_output:
            bool
        """

        logger.debug("--------------------------------------------------")
        logger.debug("effective temperature = {0:.0f}".format(effective_temperature))
        logger.debug("surface_gravity = {0:.2f}".format(surface_gravity))
        logger.debug("metallicity = {0:.2f}".format(metallicity))
        logger.debug("microturbulence = {0:.2f}".format(microturbulence))
        logger.debug("--------------------------------------------------")

        if initial_theta is None:
            initial_theta = {}

        interpolator = kwargs.pop("_interpolator", None)
        if interpolator is None:
            # Interpolate a model atmospheres
            atmosphere_kwargs = kwargs.pop("atmosphere_kwargs", {})
            interpolator = atmospheres.Interpolator(**atmosphere_kwargs)

        # Interpolate the photospheric quantities
        interpolator_kwargs = kwargs.pop("atmosphere_interpolator_kwargs", {})
        photospheric_structure = interpolator(effective_temperature,
            surface_gravity, metallicity, **interpolator_kwargs)        

        # We want to save the profile information and the details of the
        # measured atomic transitions
        profile_initial_theta = kwargs.pop("profile_initial_theta", {})
        profile_initial_theta.update({
            "effective_temperature": effective_temperature,
            "surface_gravity": surface_gravity,
            "metallicity": metallicity,
            "microturbulence": microturbulence,
        })

        # Initialise all the atomic transitions
        atomic_transitions = [transitions.AtomicTransition(**each) \
            for each in self.config["model"]["atomic_transitions"]]
        #logger.warn("ONLY USING FE FOR THE MOMENT #TODO #YOLO")

        for atomic_transition in atomic_transitions:

            # Find which channel the wavelength is in
            ci = minimum_pixel_sampling(data, atomic_transition.wavelength)[1][0]

            # Get the continuum order and limit the continuum order in this
            # region to 2, otherwise we will be overfitting
            co = self._continuum_order(ci)
            co = co if co >= 0 else None
            _profile_initial_theta = profile_initial_theta.copy()
            if co > -1:

                for i in range(co + 1):
                    k = "continuum.{0}.{1}".format(ci, i)
                    if k in initial_theta:
                        _profile_initial_theta["continuum.{}".format(i)] = initial_theta[k]

            fit_profile_kwargs = atomic_transition._init_kwargs.get(
                "_fit_profile_kwargs", {})
            fit_profile_kwargs.update({"full_output": True})
            fit_profile_kwargs.setdefault("continuum_order", co)

            try:
                result = atomic_transition.fit_profile(data[ci],
                    initial_theta=_profile_initial_theta, **fit_profile_kwargs)

            except ValueError:
                continue

            # Save the information
            #(optimal_theta, oc, orc, omf, op_info)

        # Filter out unacceptable measurements
        filtered = kwargs.pop("filter_transition",
            lambda x: (200 > x.equivalent_width > 5))

        is_filtered = np.zeros(len(atomic_transitions), dtype=bool)
        atomic_transitions_arr = np.zeros((len(atomic_transitions), 7))
        for i, atomic_transition in enumerate(atomic_transitions):
            is_filtered[i] = filtered(atomic_transition)
            atomic_transitions_arr[i, :] = [
                atomic_transition.wavelength,
                atomic_transition.species,
                atomic_transition.excitation_potential,
                atomic_transition.loggf,
                atomic_transition.van_der_waals_broadening,
                0,
                atomic_transition.equivalent_width
            ]

        fuck = np.core.records.fromarrays(atomic_transitions_arr.T,
            names=("wavelength", "species", "excitation_potential", "loggf",
                "van_der_waals_broadening", "damp2", "equivalent_width"))


        returned_abundances = synthesis.atomic_abundances(fuck[is_filtered],
            [effective_temperature, surface_gravity, metallicity],
            microturbulence=microturbulence)
        

        # Create a filler array so that filtered transitions will have an
        # abundance of 'nan'
        all_abundances = np.array([np.nan] * len(atomic_transitions_arr))
        all_abundances[is_filtered] = returned_abundances

        is_outlier = np.zeros(len(all_abundances), dtype=bool) # Assume none
        atomic_transitions_rec = np.core.records.fromarrays(
            np.vstack([np.array(atomic_transitions_arr).T,
                all_abundances, is_filtered, is_outlier]),
            names=(
                "wavelength", "species", "excitation_potential", "loggf",
                "van_der_waals_broadening", "damp2", "equivalent_width",
                "abundance", "is_filtered", "is_outlier"))

        # Calculate the excitation and ionisation state
        logger.warn("Assuming only one kind of atomic species")
        neutral_lines = (atomic_transitions_rec["species"] % 1) == 0
        ionised_lines = ~neutral_lines

        outliers_removed = False
        while True:

            # Apply filters
            use_lines = atomic_transitions_rec["is_filtered"].astype(bool) * \
                ~atomic_transitions_rec["is_outlier"].astype(bool)
            neutral_lines *= use_lines
            ionised_lines *= use_lines

            # Excitation slope 
            excitation_slope = stats.linregress(
                x=atomic_transitions_rec["excitation_potential"][neutral_lines],
                y=atomic_transitions_rec["abundance"][neutral_lines])

            # Calculate the ionisation state
            ionisation_state = \
                atomic_transitions_rec["abundance"][neutral_lines].mean() \
                    - atomic_transitions_rec["abundance"][ionised_lines].mean()

            # Calculate the abundance state
            abundance_state = (atomic_transitions_rec["abundance"][use_lines] \
                - (atmospheres.solar_abundance(atomic_transitions_rec["species"][use_lines])
                    + metallicity)).mean()

            # Slope with reduced equivalent width and line abundance
            reduced_equivalent_width = np.log(
                atomic_transitions_rec["equivalent_width"]\
                    /atomic_transitions_rec["wavelength"])
            line_strength_slope = stats.linregress(
                x=reduced_equivalent_width[neutral_lines],
                y=atomic_transitions_rec["abundance"][neutral_lines])

            # Calculate the state
            state = np.array([
                excitation_slope[0],
                ionisation_state,
                abundance_state,
                line_strength_slope[0]
            ])

            if outliers_removed:
                final_state = (excitation_slope, ionisation_state,
                    abundance_state, line_strength_slope)
                break

            else:
                # Remove the outliers?
                initial_state = (excitation_slope, ionisation_state,
                    abundance_state, line_strength_slope)

                outlier_limit = kwargs.pop("outlier_sigma_clip", 3)
                if outlier_limit is None or not np.isfinite(outlier_limit) \
                or 0 >= outlier_limit:
                    # Don't remove any outliers
                    outliers_removed, final_state = False, initial_state
                    break

                # We don't want to remove stars that are simply on the edge of
                # distributions because if our initial guess was very wrong, we
                # could just be removing the weakest or strongest lines.

                # Instead we will remove lines that are discrepant from the line
                # fits.

                # Outliers in excitation potential vs abundance:
                x = atomic_transitions_rec["excitation_potential"][use_lines]
                y = atomic_transitions_rec["abundance"][use_lines]
                line = excitation_slope[0] * x + excitation_slope[1]

                differences = np.abs(line - y)
                excitation_sigma = differences/np.std(differences)

                # Outliers in reduced equivalent width vs abundance:
                x = reduced_equivalent_width[use_lines]
                y = atomic_transitions_rec["abundance"][use_lines]
                line = line_strength_slope[0] * x + line_strength_slope[1]

                differences = np.abs(line - y)
                line_strength_sigma = differences/np.std(differences)

                # Update the is_outliers mask
                atomic_transitions_rec["is_outlier"][use_lines] = \
                    (excitation_sigma > outlier_limit) | \
                    (line_strength_sigma > outlier_limit)

                outliers_removed = True
                continue # to re-fit the lines

        """
        import matplotlib.pyplot as plt

        # Show which was used and which wasn't
        for profile, atomic_transition in zip(profiles, atomic_transitions_rec):
            if profile[3] is None: continue

            # Find this in the atomic_transitions

            # Show  whether is an outlier/filtered.

            fig, ax = plt.subplots()
            used = atomic_transition["is_filtered"].astype(bool) \
                * ~atomic_transition["is_outlier"].astype(bool)
        
            profile_spectra = profile[3]
            ax.plot(profile_spectra["data_spectrum"].disp, profile_spectra["data_spectrum"].flux, "k")

            # fitted_spectrum, data_spectrum, continuum_spectrum, blending_spectrum_unsmoothed, blending_spectrum_smoothed
            ax.plot(profile_spectra["continuum_spectrum"].disp,
                profile_spectra["continuum_spectrum"].flux, c="#666666")
            if profile_spectra["blending_spectrum_unsmoothed"] is not None:
                ax.plot(profile_spectra["blending_spectrum_unsmoothed"].disp,
                    profile_spectra["blending_spectrum_unsmoothed"].flux, 'b:')

            if profile_spectra["blending_spectrum_smoothed"] is not None:
                ax.plot(profile_spectra["blending_spectrum_smoothed"].disp,
                    profile_spectra["blending_spectrum_smoothed"].flux, 'b')

            ax.plot(profile_spectra["fitted_spectrum"].disp, profile_spectra["fitted_spectrum"].flux, "rg"[used])

            ax.axvline(profile[0].get("wavelength", atomic_transition.wavelength), linestyle="-", c="#666666")
            ax.set_xlim(profile_spectra["data_spectrum"].disp[0], profile_spectra["data_spectrum"].disp[-1])
            ax.set_xticklabels(["{0:.2f}".format(e) for e in ax.get_xticks()])

            ax.set_title("EW {0:.1f} mA, filtered: {1}, outlier: {2}, used: {3}"
                .format(atomic_transition["equivalent_width"], atomic_transition["is_filtered"],
                    atomic_transition["is_outlier"], ["No", "Yes"][used]))

            filename = "fig-{0:.1f}.png".format(atomic_transition["wavelength"])
            fig.savefig(filename)

            plt.close("all")
            print("created {0}".format(filename))


        fig, ax = plt.subplots(2)

        # Show the data
        outlier = atomic_transitions_rec["is_outlier"].astype(bool)
        ax[0].scatter(atomic_transitions_rec["excitation_potential"][~outlier],
            atomic_transitions_rec["abundance"][~outlier], facecolor="k")
        ax[1].scatter(reduced_equivalent_width[~outlier],
            atomic_transitions_rec["abundance"][~outlier], facecolor="k")

        ax[0].scatter(atomic_transitions_rec["excitation_potential"][outlier],
            atomic_transitions_rec["abundance"][outlier], facecolor="r")
        ax[1].scatter(reduced_equivalent_width[outlier],
            atomic_transitions_rec["abundance"][outlier], facecolor="r")

        # Show the initial state before outliers were removed
        if outliers_removed:
            x = np.array(ax[0].get_xlim())
            m, b = initial_state[0][:2]
            ax[0].plot(x, x * m + b, 'r')

            x = np.array(ax[1].get_xlim())
            m, b = initial_state[3][:2]
            ax[1].plot(x, x * m + b, 'r')

        # Show the final state
        x = np.array(ax[0].get_xlim())
        m, b = final_state[0][:2]
        ax[0].plot(x, x * m + b, 'k')

        x = np.array(ax[1].get_xlim())
        m, b = final_state[3][:2]
        ax[1].plot(x, x * m + b, 'k')

        print("state is ", state)
        fig.savefig("state.png")

        print("Outliers are")
        for each in atomic_transitions_rec[outlier]:
            print(each["wavelength"], each["species"])

        plt.close("all")
        """

        if full_output:
            return (state, atomic_transitions_rec, atomic_transitions)
        return state




def minimum_pixel_sampling(data, wavelength):
    """
    Find the minimum pixel size in the data at the given wavelength.
    """

    pixel_size = np.nan
    channel_indices = []
    for i, spectrum in enumerate(data):
        if spectrum.disp[-1] >= wavelength >= spectrum.disp[0]:
            channel_indices.append(i)

            index = spectrum.disp.searchsorted(wavelength)

            diff = np.diff(spectrum.disp)
            if pixel_size > diff[index] or not np.isfinite(pixel_size):
                pixel_size = diff[index]

    if len(channel_indices) == 0 and not np.isfinite(pixel_size):
        raise ValueError("cannot find wavelength {0:.2f} in any data channel"\
            .format(wavelength))

    return (pixel_size, channel_indices)

