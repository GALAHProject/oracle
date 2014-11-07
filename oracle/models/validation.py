

def validate_configuration(configuration):
    return True

def validate(model):

    return True
"""
        # Check the configuration is valid.
        self._validate()
        return None


    def _validate(self):

        # Default things that we should have.
        self.config.setdefault("mask", [])
        self.config["mask"] = np.array(self.config["mask"])

        if not self.config.has_key("model"):
            raise KeyError("no model information specified")

        validation_functions = {
            "continuum": self._validate_continuum,
            "elements": self._validate_elements
        }
        for item, state in self.config["model"].iteritems():
            if not state or item not in validation_functions: continue
            validation_functions[item]()

        return True
"""