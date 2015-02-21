#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" oracle, the suppository of all wisdom """ 

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import argparse
import cPickle as pickle
import logging
import os
import sys
from time import time

# Third-party.
import matplotlib.pyplot as plt

# Module-specific.
import oracle

logger = logging.getLogger("oracle")

# Usage: oracle estimate model.yaml <filenames>
#        oracle estimate model.yaml -r <read_from_filename>


def common_basename(filenames):
    if isinstance(filenames, (str, )):
        filenames = [filenames]
    common_prefix, ext = os.path.splitext(os.path.commonprefix(
        map(os.path.basename, filenames)))
    common_prefix = common_prefix.rstrip("_-")
    return common_prefix if len(common_prefix) > 0 else "source"

def estimate(args):
    """ Estimate model parameters by cross-correlation against a grid. """

    # Do we have to read from file?
    if args.read_from_filename:
        with open(args.spectrum_filenames[0], "r") as fp:
            all_sources = [r.split() for r in map(str.strip, fp.readlines())]
    else:
        all_sources = [args.spectrum_filenames]

    # Create the model.
    # TODO this is just a basic GALAH model. We should really read from a
    # model filename instead.
    model = oracle.models.Model({ "model": {
            "redshift": True,
            "continuum": 3,
            "continuum_mask": [
    #            [4899, 4905],
                [7592, 7730]
            ],
            "cross_correlation_mask": [
                [7500, 7730]
            ]
        },
        "settings": {
            "threads": 4
        }
        })

    successful, exceptions = 0, 0
    for i, filenames in enumerate(all_sources):
        
        basename = common_basename(filenames)
        logger.info("Sources for #{0} (basename {1}): {2}".format(i + 1,
            basename, ", ".join(filenames)))

        t_init = time()
        try:
            data = map(oracle.specutils.Spectrum1D.load, filenames)
            initial_theta, expected_dispersion, expected_flux = \
                model.initial_theta(data, full_output=True)

        except:
            exceptions += 1
            logger.exception("Exception raised when trying to analyse source #"\
                "{0}".format(i + 1))
            if args.debug: raise

        else:
            successful += 1
            logger.info("Completed successfully in {0:.1f} seconds".format(
                time() - t_init))

            # Try and add v_helio?
            try:
                initial_theta["v_helio"] = data[0].v_helio
            except KeyError:
                logger.exception("Could not calculate heliocentric velocity "
                    "correction for source #{}".format(i + 1))
                initial_theta["v_helio"] = np.nan

            logger.info("Initial model parameters for source #{0} is {1}".format(
                i + 1, initial_theta))

            # Save the initial theta information to somewhere.
            output_filename = "initial-{}.pkl".format(basename)
            with open(output_filename, "wb") as fp:
                pickle.dump((initial_theta, expected_dispersion, expected_flux),
                    fp, -1)
            logger.info("Saved output to {0}".format(output_filename))

            if args.plotting:
                plot_filename = "source-{}-initial.png".format(basename)
                fig, axes = plt.subplots(len(data))
                for j, (ax, channel) in enumerate(zip(axes, data)):
                    ax.plot(channel.disp, channel.flux, c="k")
                    ylim = ax.get_ylim()
                    ax.plot(expected_dispersion, expected_flux, c="r", zorder=-1)
                    ax.set_xlim(channel.disp.min(), channel.disp.max())
                    ax.set_ylim(ylim)
                    ax.set_ylabel("Counts")

                ax.set_xlabel("Wavelength")
                fig.tight_layout()
                fig.savefig(plot_filename)
                logger.info("Saved figure to {0}".format(plot_filename))
                plt.close("all")

    logger.info("{0} successful, {1} exceptions".format(successful, exceptions))
    


def parser(input_args=None):

    parser = argparse.ArgumentParser(
        description="oracle, the suppository of all wisdom",
        epilog="See 'oracle COMMAND -h' for help on a specific command.")

    # Create subparsers
    subparsers = parser.add_subparsers(title="command", dest="command",
        description="Specify the command to perform.")

    # Create a parent subparser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true", default=False, 
        help="Vebose mode")
    parent_parser.add_argument(
        "--overwrite", dest="overwrite",action="store_true", default=False,
        help="Overwrite existing files if they already exist")
    parent_parser.add_argument(
        "--debug", dest="debug", action="store_true", default=False,
        help="Enable debug mode. Any suppressed exception during run-time will "
            "be re-raised")

    # Create parser for the estimate command
    estimate_parser = subparsers.add_parser(
        "estimate", parents=[parent_parser],
        help="Estimate model parameters quickly.")
    #estimate_parser.add_argument(
    #    "model", type=str,
    #    help="The path for the YAML-formatted model filename")
    estimate_parser.add_argument(
        "-r", action="store_true", dest="read_from_filename", default=False,
        help="Read input spectra from a single filename")
    estimate_parser.add_argument(
        "--no-plots", dest="plotting", action="store_false", default=True,
        help="Disable plotting")
    estimate_parser.add_argument(
        "spectrum_filenames", nargs="+",
        help="Filenames of (observed) spectroscopic data")
    estimate_parser.set_defaults(func=estimate)

    args = parser.parse_args(input_args)
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    return args

    
def main():
    args = parser()
    return args.func(args)

if __name__ == "__main__":
    args = parser(sys.argv[1:])
    args.func(args)