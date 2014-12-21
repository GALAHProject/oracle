# coding: utf-8

from __future__ import absolute_import

""" oracle, the suppository of all wisdom """ 

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"
__version__ = "0.01"

import logging
from . import (models, synthesis, specutils, utils)

logger = logging.getLogger("oracle")
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s")
