"""
Module for converting SDRF filres to DIA-NN CLI call
"""
import os
from typing import ClassVar, Dict, List, Type, Union, Any
from io import IOBase
from pathlib import Path

# 3rd party imports
import pandas as pd
from abstract_converter import AbstractConverter

class DiannConverter(AbstractConverter):
    """
    SDRF converter for DIA-NN search engine
    """
    def __init__(self, diann_path: Path, diann_params: str, raw_files_path: Path) -> None:
        super().__init__()
        self.diann_path = diann_path
        self.diann_params = diann_params
        self.raw_files_path = raw_files_path
        

    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[modification parameters]": [pd.StringDtype()],
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
        "comment[data file]": [pd.StringDtype()],
        "comment[cleavage agent details]": [pd.StringDtype()]
    }

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[number of missed cleavages]": [pd.Int64Dtype()]
    }

    SDRF_COL_NAME_DIANN_PARAM_MAP: ClassVar[Dict[str, List[str]]] = {
        "comment[number of missed cleavages]": ["comment[modification parameters]*"],
        "comment[modification parameters]": ["--fixed-mod", "--var-mod"],
        "comment[precursor mass tolerance]": ["--mass-acc-ms1"],
        "comment[fragment mass tolerance]": ["--mass-acc"],
        "comment[data file]": ["--f"],
        "comment[cleavage agent details]": ["--cut"]

    }

    PROTESES_MAP: ClassVar[Dict[str, str]] = {
        "Trypsin/P": "K*,R*",
        "Trypsin": "K*,R*,!*P",
        "Lys-C": "(?<=K)(?!P)",
        "Chymotrypsin": "",
        "Asp-N": "",
        "Glu-C": ""

    }

    def parse_file_names(self, sdrf_df: pd.DataFrame) -> str:
        """
        Takes a path to the MS file(s) and their names from SDRF file
        and parse name(s) to the DIA-NN's parameter format.

        Returns
        -------
        str
            String with parsed names of MS file(s)
            e.g. --f folder/subfolder/file_1.raw
        """
        strinf = ''

        file_names = sdrf_df['comment[data file]'].values
        #apply(lambda x: strinf + '--f ' + path_to_files + '/' + x).to_string()

        strinf = ''
        for name in file_names:
            strinf = strinf + ' --f ' + self.raw_files_path + '/' + name

        return strinf

    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Any:
        """
        Convert the SDRF file to the desired format.
        The output depends on the targeted software, it may be a config file only, 
        a CLI string or a Tuple containing both.

        Returns
        -------
        Any
            The converted SDRF file
        """
        self.init_converter(sdrf)
        return str(self.diann_path) + ' ' + self.parse_file_names(self.sdrf_df) + ' ' + self.diann_params
