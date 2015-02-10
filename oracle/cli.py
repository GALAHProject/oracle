#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" oracle, the suppository of all wisdom """ 

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import argparse
import logging
import sys

import oracle

logger = logging.getLogger("oracle")

# Usage: oracle estimate model.yaml <filenames>
#        oracle estimate model.yaml -r <read_from_filename>


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
                [4899, 4905],
                [7592, 7730]
            ],
            "cross_correlation_mask": [
                [7500, 7730]
            ]
        }})

    successful, exceptions = 0, 0
    for i, filenames in enumerate(all_sources):
        try:
            data = map(oracle.specutils.Spectrum1D.load, filenames)
            initial_theta, r_chi_sq, expected_dispersion, expected_flux = \
                model.initial_theta(data, full_output=True)

        except:
            exceptions += 1
            logger.exception("Exception raised when trying to analyse source #"\
                "{0}".format(i + 1))
            if args.debug: raise

        else:
            successful += 1
            print("Initial chi-sq is {0:.2f} for theta: {1}".format(r_chi_sq,
                initial_theta))

    print("{0} successful, {1} exceptions".format(successful, exceptions))


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