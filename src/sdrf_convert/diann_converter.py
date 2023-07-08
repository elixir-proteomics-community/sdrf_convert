"""
Module for converting SDRF filres to DIA-NN cfg file and CLI calls
"""

from copy import deepcopy
from io import IOBase
from pathlib import Path
import re
from typing import ClassVar, Dict, Iterator, List, Set, Tuple, Type, Union

# 3rd party imports
import pandas as pd

from sdrf_convert.abstract_converter import AbstractConverter

class DiannConverter(AbstractConverter):
    """
    SDRF converter for DIA-NN search engine
    """
    def __init__(self, diann_params: str):
            super().__init__()
            self.diann_params = diann_params


    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[modification parameters]": [pd.StringDtype()],
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
        "comment[data file]": [pd.StringDtype()],
        "comment[cleavage agent details]": [pd.StringDtype()],
        "comment[number of missed cleavages]": [pd.Int64Dtype()]
    }

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {

    }