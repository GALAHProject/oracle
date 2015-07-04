 def fit_all_atomic_transitions(self, data, effective_temperature=None,
        surface_gravity=None, metallicity=None, microturbulence=None,
        wavelength_region=2.01, oversampling_rate=4, **kwargs):
        """
        Fit all clean transitions within a single data channel.
        """

        photosphere = None
        # Interpolate a photosphere if we have the information to do so.
        if None not in (effective_temperature, surface_gravity, metallicity):
            photosphere = self._interpolator(
                effective_temperature, surface_gravity, metallicity)

        # Identify clean transitions that are actually in the data channels.
        # If data_indices is -1 it means the line was not found in the data

        clean = self.atomic_transitions["clean"]
        data_indices = \
            wavelengths_in_data(self.atomic_transitions["wavelength"], data)
        measurable = clean * (data_indices > -1)

        # Create the subset containing the measurable lines.
        transition_indices = np.where(measurable)[0]
        transitions = self.atomic_transitions[measurable]
        data_indices = data_indices[measurable]


        # Join by data indices because we will be doing one full channel at a
        # time.

        if "profile_stddev" not in self.atomic_transitions.dtype.names:
            self.atomic_transitions.add_column(table.Column(
                name="profile_stddev",
                data=[np.nan] * len(self.atomic_transitions)))

        import matplotlib.pyplot as plt
        fitted_profiles = []
        for data_index in np.unique(data_indices):

            data_channel = data[data_index]
            pixels = np.zeros(data_channel.disp.size, dtype=bool)
            is_blended = np.zeros(data_channel.disp.size, dtype=bool)

            comparison_regions = []
            _ = np.where(data_indices == data_index)[0]
            wavelengths = np.sort(transitions["wavelength"][_])
            for wavelength in wavelengths:
                start, end = (wavelength - wavelength_region, \
                    wavelength + wavelength_region)

                # Does this overlap with the last region?
                if len(comparison_regions) == 0 \
                or start > comparison_regions[-1][1]:
                    # If not add this region.
                    comparison_regions.append([start, end])

                else:
                    # If so just extend the last region to cover this one.
                    comparison_regions[-1][1] = end

            # is_blended will be a masked array indicating whether a pixel needs
            # to be interpolated from the synthetic spectra.

            synthetic_disp, synthetic_flux = [], []
            for start, end in comparison_regions:

                # Save the data pixels for this region.
                indices = data_channel.disp.searchsorted([start, end]) + [0, 1]
                pixels.__setslice__(indices[0], indices[1], True)

                relevant_lines = \
                    (end >= self.atomic_transitions["wavelength"]) * \
                    (self.atomic_transitions["wavelength"] >= start)

                # Apply custom masks.
                for relevant_line in self.atomic_transitions[relevant_lines]:
                    if not isinstance(relevant_line["custom_mask"],
                        (str, unicode)): continue

                    wavelength = relevant_line["wavelength"]
                    mask_name = relevant_line["custom_mask"]
                    mask = np.clip(self.config["custom_mask"][mask_name.strip()],
                        max([start, wavelength - wavelength_region]),
                        min([end, wavelength + wavelength_region]))
                    pixels[indices[0]:indices[1]] *= \
                        self._mask(data_channel.disp[indices[0]:indices[1]], mask)

                # Synthesise spectra as required, and save the data pixels.
                blending_lines = ~clean * relevant_lines

                # If there are no blending lines in this region, we should go on
                if not np.any(blending_lines):
                    continue

                # Check that we can actually do the synthesis
                if photosphere is None:
                    raise ValueError(
                        "there are blending lines in the region between {0:.0f}"
                        "to {1:.0f} (so synthesis will be needed) but not all "
                        "the stellar parameters were given".format(start, end))

                # Mark these x-pixels as being blended with a background spec.
                is_blended.__setslice__(indices[0], indices[1], True)

                # Synthesise an over-sampled spectrum.
                x = data_channel.disp[indices[0]:indices[1]]
                synth_pixel_size = np.diff(x).mean()/oversampling_rate
                disp, flux = synthesis.moog.synthesise(
                    self.atomic_transitions[blending_lines], photosphere,
                    microturbulence=microturbulence,
                    wavelength_region=[x.min(), x.max()],
                    wavelength_step=synth_pixel_size)

                # Store the synthetic spectrum
                synthetic_disp.extend(disp)
                synthetic_flux.extend(flux)


            # Build a model.
            # Needs to have the disp/flux data stored internally.
            # Needs to add in all the absorption profiles at each point.

            # Local continuum at each region?


            # Prepare the data
            x = data_channel.disp[pixels]
            y = data_channel.flux[pixels]
            is_blended = is_blended[pixels]

            # Prepare the synthetic data.
            synthetic_disp = np.array(synthetic_disp)
            synthetic_flux = np.array(synthetic_flux)

            # Define the model.
            class StellarSpectrum(modeling.Fittable1DModel):

                resolution = modeling.Parameter(default=20000)

                stddev_m = modeling.Parameter(default=1e-5)
                stddev_c = modeling.Parameter(default=1e-2)

                @staticmethod
                def evaluate(x, resolution, stddev_m, stddev_c):

                    # Non-blended regions will be 1.
                    model = np.ones(x.size, dtype=float)

                    # Convolve the synthetic spectrum
                    if synthetic_disp.size > 0:
                        binner = _rebinner(synthetic_disp, x, resolution)
                        model[is_blended] = \
                            np.array(synthetic_flux * binner)[is_blended]
                    return model

            # Create the model, then compound it with absorption profiles at
            # each of the wavelengths.
            model = StellarSpectrum(resolution=10000)
            for wavelength in wavelengths:

                idx = x.searchsorted(wavelength)
                if is_blended[idx]:
                    depth = \
                        synthetic_flux[synthetic_disp.searchsorted(wavelength)]
                else:
                    depth = 1

                # TODO this is a bad approximation
                amplitude = np.clip(depth - y[idx], 0, 1)

                model *= modeling.models.GaussianAbsorption1D(
                    mean=wavelength, amplitude=amplitude,
                    stddev=0.10)

            # Tie the std. devs. to the spectral resolution, and set boundaries
            # for the amplitude.
            for n in range(1, len(wavelengths) + 1):
                #model.tied["stddev_{}".format(n)] = lambda _: _.stddev_1 + \
                #    _.resolution_0 * getattr(_, "mean_{}".format(n))
                """

                if np.any(0.01 >= np.abs(np.array([4890.759, 4891.494]) - wavelengths[n-1])):
                    print("skipping on {}".format(n))

                else:
                    model.tied["stddev_{}".format(n)] = lambda _: \
                        _.stddev_c_0 + _.stddev_m_0 * (getattr(_, "mean_{}".format(n)) \
                            - getattr(_, "mean_1"))

                    init = model.stddev_m_0 * wavelengths[n - 1] + model.stddev_c_0
                    setattr(model, "stddev_{}".format(n), init)


                """
                model.bounds["stddev_{}".format(n)] = [0, 1]
                model.bounds["amplitude_{}".format(n)] = [0, 1]
                model.fixed["mean_{}".format(n)] = True



            fit = modeling.fitting.LevMarLSQFitter()
            model.bounds["stddev_m_0"] = [0, None]
            model.bounds["resolution_0"] = [7500, 40000]

            fitted = fit(model, x, y)

            np.sort(transitions["wavelength"][_])

            fig, ax = plt.subplots()
            ax.plot(x, y, c='k')
            ax.scatter(x, y, facecolor='k')
            ax.plot(synthetic_disp, synthetic_flux, "r:")
            ax.plot(x, model(x), 'r-.')

            """
            model2 = np.ones(x.size, dtype=float)

            # Convolve the synthetic spectrum
            binner = _rebinner(synthetic_disp, x, 20000)
            model2[is_blended] = \
                np.array(synthetic_flux * binner)[is_blended]

            ax.plot(x, model2, c='g')
            """


            for wavelength in wavelengths:
                ax.axvline(wavelength, c="#666666", zorder=-1)


            ax.plot(x, fitted(x), c='r')


            for j in range(1, len(wavelengths) + 1):

                clean = self.atomic_transitions["clean"]
                index = \
                    np.argmin(np.abs(self.atomic_transitions["wavelength"][clean] \
                        - getattr(fitted, "mean_{}".format(j))))
                index = np.where(clean)[0][index]

                stddev = getattr(fitted, "stddev_{}".format(j))
                amplitude = getattr(fitted, "amplitude_{}".format(j))

                self.atomic_transitions["equivalent_width"][index] = \
                    1000 * np.sqrt(2 * np.pi) * amplitude * stddev
                self.atomic_transitions["profile_stddev"][index] = stddev.value


        sensible = np.isfinite(self.atomic_transitions["equivalent_width"]) \
            * (self.atomic_transitions["equivalent_width"] > 0)
        # Eliminate all other abundances, then re-calculate them.
        self.atomic_transitions["abundance"] = np.nan
        self.atomic_transitions["abundance"][sensible] = \
            synthesis.moog.atomic_abundances(self.atomic_transitions[sensible],
                photosphere, microturbulence=microturbulence)


        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)

        # Add profile_* columns to the atomic_transitions table:
        # profile_stddev, profile_amplitude
        ok = np.isfinite(self.atomic_transitions["abundance"])

        ax[0].scatter(self.atomic_transitions["excitation_potential"][ok],
            self.atomic_transitions["abundance"][ok], facecolor='k')
        ax[1].scatter(
            np.log(self.atomic_transitions["equivalent_width"]/self.atomic_transitions["wavelength"])[ok],
            self.atomic_transitions["abundance"][ok], facecolor="k")

        ax[2].scatter(self.atomic_transitions["wavelength"][ok],
            self.atomic_transitions["profile_stddev"][ok], facecolor="k")

        state, info = utils.equalibrium_state(
            self.atomic_transitions, metallicity=metallicity, full_output=True)
        coefficients = info["coefficients"]

        x = self.atomic_transitions["excitation_potential"][ok]
        ax[0].plot(x, np.polyval(coefficients[0], x), c='b')
        x = np.log(self.atomic_transitions["equivalent_width"] \
            / self.atomic_transitions["wavelength"])[ok]
        ax[1].plot(x, np.polyval(coefficients[1], x), c='b')

        plt.show()

        raise a





    def fit_atomic_transitions2(self, data, effective_temperature=None,
        surface_gravity=None, metallicity=None, microturbulence=None,
        wavelength_region=2.0, opacity_contributes=1.0, **kwargs):
        """
        Fit absorption profiles to the atomic transitions in this model and
        account for the nearby blends.

        :param data:
            The observed spectra.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D`

        :param effective_temperature: [sometimes optional]
            The effective temperature for the photosphere.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type effective_temperature:
            float

        :param surface_gravity: [sometimes optional]
            The surface gravity for the photosphere.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type surface_gravity:
            float

        :param metallicity: [sometimes optional]
            The scaled-solar metallicity for the photosphere.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type metallicity:
            float

        :param microturbulence: [sometimes optional]
            The microturbulence for the photosphere. This is not required for
            <3D> models.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type microturbulence:
            float

        :param wavelength_region: [optional]
            The +/- region (in Angstroms) around each atomic transition to
            consider.

        :type wavelength_region:
            float

        :param outlier_modeling: [optional]
            Detect inaccurate fits due to blending lines and account for them.
            If enabled, then a poor fit is detected when the centroid of the
            atomic transition is not within 5 per cent of the data, or if the
            profile FWHM has hit an upper boundary. When this happens the code
            will detect regions that are most discrepant, and attempt to fit
            them with additional profiles (using the same profile FWHM).

        :type outlier_modeling:
            bool

        :param max_outlier_profiles: [optional]
            The maximum number of outlier profiles to add for each atomic
            transition. This keyword argument is ignored if `outlier_modeling`
            is set to `False`.
        """

        oversampling_rate = int(kwargs.pop("oversampling_rate", 4))
        if 1 > oversampling_rate:
            raise ValueError("oversampling rate must be a positive integer")


        photosphere = None
        # Interpolate a photosphere if we have the information to do so.
        if None not in (effective_temperature, surface_gravity, metallicity):
            photosphere = self._interpolator(
                effective_temperature, surface_gravity, metallicity)

        # Identify clean transitions that are actually in the data.
        # If data_indices is -1 it means the line was not found in the data
        data_indices = \
            wavelengths_in_data(self.atomic_transitions["wavelength"], data)
        measurable = self.atomic_transitions["clean"] * (data_indices > -1) \
            * ((self.atomic_transitions["species"] == 26.0) | \
                (self.atomic_transitions["species"] == 26.1))

        # Create the subset containing the measurable lines.
        transition_indices = np.where(measurable)[0]
        transitions = self.atomic_transitions[measurable]
        data_indices = data_indices[measurable]


        plotting = True
        if plotting:
            import matplotlib.pyplot as plt

        stopper = -1

        if "profile_stddev" not in self.atomic_transitions.dtype.names:
            self.atomic_transitions.add_column(table.Column(
                name="profile_stddev",
                data=[np.nan] * len(self.atomic_transitions)))

        fitted_profiles = []
        fit = modeling.fitting.LevMarLSQFitter()
        for i, (transition_index, transition, data_index) \
        in enumerate(zip(transition_indices, transitions, data_indices)):

            if i < stopper:
                continue

            wavelength = transition["wavelength"]

            # Look for nearby transitions within the wavelength region and
            # ignore this line.
            blending = wavelength_region >= \
                np.abs(self.atomic_transitions["wavelength"] - wavelength)
            blending[transition_index] = False

            # Slice the data +/- some region.
            spectrum = data[data_index]
            disp_indices = spectrum.disp.searchsorted([
                wavelength - wavelength_region,
                wavelength + wavelength_region
                ]) + [0, 1]
            x = spectrum.disp.__getslice__(*disp_indices)
            y = spectrum.flux.__getslice__(*disp_indices)


            # Continuum!
            model = modeling.models.Polynomial1D(degree=1, c0=1, c1=0)
            #model = modeling.models.Const1D(amplitude=1)

            bounds = {}
            #fixed = ["c0_0", "c1_0"]
            fixed = []

            # Anything to synthesise?
            if np.any(blending):
                if photosphere is None:
                    raise ValueError("transition at {0:.3f} has blending lines "
                        "within {1:.0f} (so a synthesis approach is needed) but"
                        " not all stellar parameters were given".format(
                            wavelength, wavelength_region))

                # Synthesise a spectrum (with oversampling)
                synth_pixel_size = np.diff(x).mean()/oversampling_rate
                synth_disp, synth_flux = synthesis.moog.synthesise(
                    self.atomic_transitions[blending], photosphere,
                    microturbulence=microturbulence,
                    wavelength_region=[
                        wavelength - (wavelength_region + opacity_contributes),
                        wavelength + (wavelength_region + opacity_contributes)
                    ],
                    wavelength_step=synth_pixel_size)


                class StellarSpectrum(modeling.Fittable1DModel):

                    resolution = modeling.Parameter(default=20000)

                    @staticmethod
                    def evaluate(x, resolution):
                        binner = _rebinner(synth_disp, x, resolution)
                        return np.array(synth_flux * binner)

                model *= StellarSpectrum()
                bounds["resolution_1"] = [10000, 40000]

                bounds["amplitude_2"] = [0, 1]
                bounds["stddev_2"] = [0, 0.3]
                fixed.append("mean_2")

            else:

                bounds["amplitude_1"] = [0, 1]
                bounds["stddev_1"] = [0, 0.3]
                fixed.append("mean_1")

            fixed.append("amplitude_0")

            # Add Gaussian profile.
            idx = x.searchsorted(wavelength)
            model *= modeling.models.GaussianAbsorption1D(
                mean=wavelength, amplitude=1.0 - y[idx],
                stddev=0.10)

            # Tie and bounds.
            for k in fixed:
                model.fixed[k] = True

            for k, v in bounds.items():
                model.bounds[k] = v

            if np.any(blending):
                model.tied["resolution_2"] = lambda _: transition/(_.stddev_2 * 2.65)
            else:
                model.tied["resolution_1"] = lambda _: transition/(_.stddev_1 * 2.65)


            fitted = fit(model, x, y)


            # Any masks to apply?
            if plotting:

                fig, ax = plt.subplots()
                ax.plot(x, y, c='k')

                if isinstance(transition["custom_mask"], (str, unicode)) \
                and len(transition["custom_mask"].strip()) > 0:
                    mask_name = transition["custom_mask"].strip()
                    mask = self._mask(x, self.config["custom_mask"][mask_name])
                    print(len(x), len(y))
                    x = x[mask]
                    y = y[mask]
                    print(len(x), len(y))
                    print(mask)

                    ax.scatter(x, y, c='k', lw=2)

                ax.plot(x, model(x), 'r:')
                ax.plot(x, fitted(x), c='r')
                ax.axvline(wavelength, c="#666666", zorder=-1)

            # Save the equivalent widths.

            if i == stopper:
                if plotting:
                    plt.show()
                raise a

            amplitude = getattr(fitted, "amplitude_1",
                getattr(fitted, "amplitude_2", None)).value
            stddev = getattr(fitted, "stddev_1",
                getattr(fitted, "stddev_2", None)).value
            self.atomic_transitions["equivalent_width"][transition_index] = \
                1000 * np.sqrt(2 * np.pi) * amplitude * stddev
            self.atomic_transitions["profile_stddev"][transition_index] = stddev


        sensible = np.isfinite(self.atomic_transitions["equivalent_width"]) \
            * (self.atomic_transitions["equivalent_width"] > 0)
        # Eliminate all other abundances, then re-calculate them.
        self.atomic_transitions["abundance"] = np.nan
        self.atomic_transitions["abundance"][sensible] = \
            synthesis.moog.atomic_abundances(self.atomic_transitions[sensible],
                photosphere, microturbulence=microturbulence)


        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)

        # Add profile_* columns to the atomic_transitions table:
        # profile_stddev, profile_amplitude
        ok = np.isfinite(self.atomic_transitions["abundance"])

        ax[0].scatter(self.atomic_transitions["excitation_potential"][ok],
            self.atomic_transitions["abundance"][ok], facecolor='k')
        ax[1].scatter(
            np.log(self.atomic_transitions["equivalent_width"]/self.atomic_transitions["wavelength"])[ok],
            self.atomic_transitions["abundance"][ok], facecolor="k")

        ax[2].scatter(self.atomic_transitions["wavelength"][ok],
            self.atomic_transitions["profile_stddev"][ok], facecolor="k")

        state, info = utils.equalibrium_state(
            self.atomic_transitions, metallicity=metallicity, full_output=True)
        coefficients = info["coefficients"]

        x = self.atomic_transitions["excitation_potential"][ok]
        ax[0].plot(x, np.polyval(coefficients[0], x), c='b')
        x = np.log(self.atomic_transitions["equivalent_width"] \
            / self.atomic_transitions["wavelength"])[ok]
        ax[1].plot(x, np.polyval(coefficients[1], x), c='b')

        plt.show()

        return None





    def fit_atomic_transitions(self, data, effective_temperature=None,
        surface_gravity=None, metallicity=None, microturbulence=None,
        wavelength_region=2.5, outlier_modeling=False, max_outlier_profiles=5,
        **kwargs):
        """
        Fit absorption profiles to the atomic transitions in this model and
        account for the nearby blends.

        :param data:
            The observed spectra.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D`

        :param effective_temperature: [sometimes optional]
            The effective temperature for the photosphere.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type effective_temperature:
            float

        :param surface_gravity: [sometimes optional]
            The surface gravity for the photosphere.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type surface_gravity:
            float

        :param metallicity: [sometimes optional]
            The scaled-solar metallicity for the photosphere.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type metallicity:
            float

        :param microturbulence: [sometimes optional]
            The microturbulence for the photosphere. This is not required for
            <3D> models.

            When there are nearby blending transitions, these lines need to be
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type microturbulence:
            float

        :param wavelength_region: [optional]
            The +/- region (in Angstroms) around each atomic transition to
            consider.

        :type wavelength_region:
            float

        :param outlier_modeling: [optional]
            Detect inaccurate fits due to blending lines and account for them.
            If enabled, then a poor fit is detected when the centroid of the
            atomic transition is not within 5 per cent of the data, or if the
            profile FWHM has hit an upper boundary. When this happens the code
            will detect regions that are most discrepant, and attempt to fit
            them with additional profiles (using the same profile FWHM).

        :type outlier_modeling:
            bool

        :param max_outlier_profiles: [optional]
            The maximum number of outlier profiles to add for each atomic
            transition. This keyword argument is ignored if `outlier_modeling`
            is set to `False`.
        """

        oversampling_rate = int(kwargs.pop("oversampling_rate", 4))
        if 1 > oversampling_rate:
            raise ValueError("oversampling rate must be a positive integer")

        # TODO this should really be changed to a limit based on radial velocity.
        wavelength_tolerance = abs(kwargs.pop("wavelength_tolerance", 0))

        # Allow the user to specify bounds on the data.
        # (And here we will set some sensible ones.)
        common_bounds = kwargs.pop("bounds", {
            "stddev": (0, 0.3),
            "amplitude": (0, 1)
        })
        if "wavelength" in common_bounds or "mean" in common_bounds:
            raise ValueError("apply bounds on the profile location through the "
                "wavelength_tolerance keyword argument")

        # Create a handy function to update the compound model properties.
        def _update_compound_model(compound_model, wavelength):

            # Deal with the line we care about first.
            single_model = hasattr(compound_model, "mean")
            transition_mean_key = ["mean_0", "mean"][single_model]
            if wavelength_tolerance > 0:
                compound_model.bounds[transition_mean_key] = (
                    wavelength - wavelength_tolerance,
                    wavelength + wavelength_tolerance
                )
            else:
                compound_model.fixed[transition_mean_key] = True

            for param_name in compound_model.param_names:
                if param_name.startswith("mean_") \
                and param_name != transition_mean_key:
                    # It's an outlier line. Fix the wavelength.
                    compound_model.fixed[param_name] = True

                # Apply common bounds.
                try:
                    prefix, num = param_name.split("_")

                except ValueError:
                    prefix, num = param_name, "0"

                if prefix in common_bounds.keys():
                    compound_model.bounds[param_name] = common_bounds[prefix]

                # Tie the stddevs to the original absorption profile.
                # The stddevs between the absorption transition we care
                # about and the outlier transition are related by:
                # R = lambda_1/delta_lambda_1 = lambda_2/delta_labmda_2
                if prefix == "stddev" and num != "0":
                    compound_model.tied[param_name] = lambda _: _.stddev_0 * \
                        getattr(compound_model, "mean_{}".format(num))/_.mean_0
            return True


        photosphere = None
        # Interpolate a photosphere if we have the information to do so.
        if None not in (effective_temperature, surface_gravity, metallicity):
            photosphere = self._interpolator(
                effective_temperature, surface_gravity, metallicity)

        # Identify clean transitions that are actually in the data.
        # If data_indices is -1 it means the line was not found in the data
        data_indices = \
            wavelengths_in_data(self.atomic_transitions["wavelength"], data)
        measurable = self.atomic_transitions["clean"] * (data_indices > -1)

        # Create additional columns in the atomic_transitions table if needed.
        if "profile_amplitude" not in self.atomic_transitions.dtype.names:
            self.atomic_transitions.add_column(table.Column(
                name="profile_amplitude",
                data=[np.nan] * len(self.atomic_transitions)))
        if "profile_stddev" not in self.atomic_transitions.dtype.names:
            self.atomic_transitions.add_column(table.Column(
                name="profile_stddev",
                data=[np.nan] * len(self.atomic_transitions)))

        # Create the subset containing the measurable lines.
        transition_indices = np.where(measurable)[0]
        transitions = self.atomic_transitions[measurable]
        data_indices = data_indices[measurable]

        fitted_profiles = []
        for transition_index, transition, data_index \
        in zip(transition_indices, transitions, data_indices):

            wavelength = transition["wavelength"]

            # Look for nearby transitions within the wavelength region and
            # ignore this line.
            blending = wavelength_region >= \
                np.abs(self.atomic_transitions["wavelength"] - wavelength)
            blending[transition_index] = False

            # Slice the data +/- some region.
            spectrum = data[data_index]
            disp_indices = spectrum.disp.searchsorted([
                wavelength - wavelength_region,
                wavelength + wavelength_region
                ]) + [0, 1]
            x = spectrum.disp.__getslice__(*disp_indices)
            y = spectrum.flux.__getslice__(*disp_indices)

            # Calculate an initial stddev value based on the x spacing.
            initial_stddev = 2 * 5 * np.diff(x).mean()

            # Apply any custom mask.
            # TODO This should just remove the masked pixels from x and y,
            #      because setting them to NaN will break the fitter.


            # TODO the amplitude initial guess will have to be udpated in the
            #      presence of continuum.

            synthesised_spectra = {}
            _ = x.searchsorted(wavelength)
            initial_amplitude = 1.0 - y[_]

            # Any continuum?
            continuum_degree = self._continuum_degree(data_index)
            if continuum_degree > 0:
                # TODO I specify order and astropy uses degree. Switch to degree!
                #profile_init *= modeling.models.Polynomial1D(continuum_degree + 1)

                # Set initial estimates of continuum.
                # TODO

                # This will probably fuck up the parameter names.
                raise NotImplementedError

            # Any synthesis?
            if np.any(blending):
                if photosphere is None:
                    raise ValueError("transition at {0:.3f} has blending lines "
                        "within {1:.0f} (so a synthesis approach is needed) but"
                        " not all stellar parameters were given".format(
                            wavelength, wavelength_region))

                # Synthesise a spectrum (with oversampling)
                synth_pixel_size = np.diff(x).mean()/oversampling_rate
                synth_disp, synth_flux = synthesis.moog.synthesise(
                    self.atomic_transitions[blending], photosphere,
                    microturbulence=microturbulence,
                    wavelength_region=[x.min(), x.max()],
                    wavelength_step=synth_pixel_size)

                # Create a custom class that uses the synthesised spectrum.
                class GaussianAbsorption1D(modeling.Fittable1DModel):

                    amplitude = modeling.Parameter(default=1)
                    mean = modeling.Parameter(default=0)
                    stddev = modeling.Parameter(default=1)
                    # TODO see issue 18 on GitHub

                    @staticmethod
                    def evaluate(x, amplitude, mean, stddev):
                        convolved = ndimage.gaussian_filter1d(
                            synth_flux, stddev/synth_pixel_size)
                        sampled = np.interp(x, synth_disp, convolved, 1, 1)
                        return sampled * (1.0 - \
                            modeling.models.Gaussian1D.evaluate(x, amplitude,
                                mean, stddev))

                profile = GaussianAbsorption1D(
                    mean=wavelength, amplitude=initial_amplitude,
                    stddev=initial_stddev)

                # Save the initial profile + synthesised spectra
                # TODO multiply by some continuum
                synth_only_model = GaussianAbsorption1D(
                    mean=wavelength, amplitude=0, stddev=initial_stddev)
                synthesised_spectra["initial_synthesis"] = synth_only_model(x)
                profile_only_model = modeling.models.GaussianAbsorption1D(
                    mean=wavelength, amplitude=initial_amplitude,
                    stddev=initial_stddev)
                synthesised_spectra["initial_profile"] = profile_only_model(x)

            else:
                profile = modeling.models.GaussianAbsorption1D(
                    mean=wavelength, amplitude=initial_amplitude,
                    stddev=initial_stddev)

                # TODO multiply by some continuum
                synthesised_spectra["initial_synthesis"] = np.ones(len(x))
                synthesised_spectra["initial_profile"] = profile(x)

            # Update the bounds, fixed, and tied properties.
            _update_compound_model(profile, wavelength)

            j, fitter = 1, modeling.fitting.LevMarLSQFitter()
            synthesised_spectra["initial_composite"] = profile(x)

            while True:

                # Fit the profile.
                fitted = fitter(profile, x, y)

                # Break here if we have no more outlier modeling to do.
                if not outlier_modeling or j > max_outlier_profiles: break

                # Limitingly-high stddev values are good indicators of nearby
                # lines that have not been accounted for.

                # Having a large % difference between the profile and the data
                # at the transition point is another good indicator of nearby
                # lines that have not been accounted for.
                if (j == 1 and fitted.stddev == profile.bounds["stddev"][1])   \
                or (j > 1 and fitted.stddev_0 == profile.bounds["stddev_0"][1])\
                or not (1.05 > y[_]/fitted(x[_]) > 0.95): #absorption is 5% off

                    # Add a(nother) outlier absorption profile to this model at
                    # the location where a blending line is most likely to be.

                    # But ignore locations near existing lines so we don't just
                    # pile up absorption profiles in the same place.
                    existing_means = [getattr(fitted, key) \
                        for key in fitted.param_names if key[:4] == "mean"]

                    difference = abs(fitted(x) - y)
                    for mean in existing_means:
                        __ = np.clip(x.searchsorted([
                            mean - 3 * initial_stddev,
                            mean + 3 * initial_stddev
                        ]) + [0, 1], 0, len(difference))
                        difference.__setslice__(__[0], __[1], 0)

                    most_discrepant = difference.argmax()

                    # TODO continuum will fuck this up too.
                    profile *= modeling.models.GaussianAbsorption1D(
                        mean=x[most_discrepant], stddev=initial_stddev,
                        amplitude=1.0 - y[most_discrepant])

                    # Update the bounds, fixed, and tied properties.
                    _update_compound_model(profile, wavelength)
                    j += 1

                else:
                    # No outlier treatment required, apparently.
                    break

            # Create final copies of things
            synthesised_spectra["fitted_composite"] = fitted(x)

            # Save the final model fit, and the model parameters.
            parameters = dict(zip(fitted.param_names, fitted.parameters))
            fitted_profiles.append((x, y, synthesised_spectra, parameters))

            # Update the equivalent width
            amplitude = \
                parameters.get("amplitude", parameters.get("amplitude_0", None))
            stddev = parameters.get("stddev", parameters.get("stddev_0", None))

            # Integral of Gaussian = amplitude * sigma * sqrt(2 * pi)
            #              (in mA) *= 1000
            self.atomic_transitions["profile_stddev"][transition_index] = stddev
            self.atomic_transitions["profile_amplitude"][transition_index] = \
                amplitude
            self.atomic_transitions["equivalent_width"][transition_index] = \
                1000 * np.sqrt(2 * np.pi) * amplitude * stddev

        sensible = np.isfinite(self.atomic_transitions["equivalent_width"]) \
            * (self.atomic_transitions["equivalent_width"] > 0)
        # Eliminate all other abundances, then re-calculate them.
        self.atomic_transitions["abundance"] = np.nan
        self.atomic_transitions["abundance"][sensible] = \
            synthesis.moog.atomic_abundances(self.atomic_transitions[sensible],
                photosphere, microturbulence=microturbulence)

        # Supply some metadata to the atomic_transitions table
        self.atomic_transitions.meta["profiles_given_stellar_parameters"] = \
            [effective_temperature, surface_gravity, metallicity, microturbulence]
        self.atomic_transitions.meta["abundances_given_stellar_parameters"] = \
            [effective_temperature, surface_gravity, metallicity, microturbulence]

        #x return fitted_profiles

        # At this point should we consider re-fitting lines that are deviant
        # from the wavelength vs stddev plot
        for i, (x, y, synthesised_spectra, parameters) in enumerate(fitted_profiles):

            fig, ax = plt.subplots()
            ax.plot(x,y,c='k')
            ax.plot(x, synthesised_spectra["initial_profile"], "r", lw=1.5, label="Initial profile")
            ax.plot(x, synthesised_spectra["initial_synthesis"], "b", lw=1.5, label="Initial synthesis")
            ax.plot(x, synthesised_spectra["initial_composite"], "m", label="Initial composite")
            ax.plot(x, synthesised_spectra["fitted_composite"], "g", label="Fitted composite")

            ax.legend()
            for p, v in parameters.items():
                if p[:4] == "mean" and p not in ("mean", "mean_0"):
                    ax.axvline(v, c="r")

        # Only update those with good quality constraints.
        fig, ax = plt.subplots(3)

        # Add profile_* columns to the atomic_transitions table:
        # profile_stddev, profile_amplitude

        ax[0].scatter(self.atomic_transitions["excitation_potential"],
            self.atomic_transitions["abundance"], facecolor='k')
        ax[1].scatter(
            np.log(self.atomic_transitions["equivalent_width"]/self.atomic_transitions["wavelength"]),
            self.atomic_transitions["abundance"], facecolor="k")

        ax[2].scatter(self.atomic_transitions["wavelength"],
            self.atomic_transitions["profile_stddev"], facecolor="k")

        state, info = utils.equalibrium_state(
            self.atomic_transitions, full_output=True)
        coefficients = info["coefficients"]
        ok = np.isfinite(self.atomic_transitions["abundance"])
        x = self.atomic_transitions["excitation_potential"][ok]
        ax[0].plot(x, np.polyval(coefficients[0], x), c='b')
        x = np.log(self.atomic_transitions["equivalent_width"] \
            / self.atomic_transitions["wavelength"])[ok]
        ax[1].plot(x, np.polyval(coefficients[1], x), c='b')
        raise a
