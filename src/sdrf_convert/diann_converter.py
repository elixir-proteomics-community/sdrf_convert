"""
Module for converting SDRF filres to DIA-NN cfg file and CLI calls
"""

from copy import deepcopy
from io import IOBase, StringIO
from pathlib import Path
import re
from typing import ClassVar, Dict, Iterator, List, Set, Tuple, Type, Union

# 3rd party imports
import pandas as pd

from abstract_converter import AbstractConverter

SDRF_CELL_SEPARATOR = '\t'

class DiannConverter(AbstractConverter):
    """
    SDRF converter for DIA-NN search engine
    """
    def __init__(self, sdrf: pd.DataFrame, diann_params: StringIO):
            super().__init__()
            self.sdrf = sdrf
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

    @classmethod
    def check_sdrf_path(self, sdrf_path: Union[str, StringIO, Path]) -> StringIO:
        """
        Checks the type of sdrf path and convert i to the StringIO if needed.

        Returns
        -------
        StringIO
            Sdrf path string in StringIO type
        """
        if not isinstance(sdrf_path, StringIO):
            sdrf_path = StringIO(sdrf_path)
        else:
            return sdrf_path
        return sdrf_path

    @classmethod
    def read_sdrf(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> pd.DataFrame:
        """
        Read SDRF file into a pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            The SDRF file as pandas DataFrame

        Raises
        ------
        TypeError
            If sdrf is not a pandas DataFrame or a StringIO object
        """
        # Check input type
        if not isinstance(sdrf, (pd.DataFrame, IOBase, Path)):
            raise TypeError((
                "Expected sdrf to be either a pandas DataFrame "
                f"or a StringIO object but got {type(sdrf)}"
            ))
        if isinstance(sdrf, pd.DataFrame):
            # Return SDRF DataFrame
            return sdrf
        # Convert to DataFrame
        else:
            return pd.read_csv(sdrf, sep=SDRF_CELL_SEPARATOR, dtype={
                col: dtypes[0]
                for col, dtypes in self.COLUMN_PROPERTIES.items()
            } # Use the first accepted dtype for each column
        )

    @classmethod
    def check_necessary_columns_and_types(self, sdrf_df: pd.DataFrame) -> None:
        """
        Check if all necessary columns are present in the SDRF 
        file and if they have the correct type.

        Raises
        ------
        ValueError
            If column is not present in SDRF file
        TypeError
            If column dtype is not correct
        """
        for col_name, col_types in self.COLUMN_PROPERTIES.items():
            col_names: List[str] = self.find_columns(sdrf_df, col_name)
            for col_name in col_names:
                if not col_name in sdrf_df.columns:
                    raise ValueError(f"Column `{col_name}` not found in SDRF file")
                if not sdrf_df[col_name].dtype in col_types:
                    raise TypeError((
                        f"Column `{col_name}` has wrong type. "
                        f"Expected one of {col_types} but got `{sdrf_df[col_name].dtype}`"
                    ))
                
    @classmethod
    def ontology_str_to_dict(cls, ontology_str: str) -> Dict[str, str]:
        """
        Convert ontology string to a dict.
        (Maybe ontology string is not the correct term.)

        Parameters
        ----------
        ontology_str : str
            NT=Oxidation;TA=M;AC=UNIMOD:35;MT=variable

        Returns
        -------
        Dict[str, str]
            {
                "NT": "Oxidation",
                "TA": "M",
                "AC": "UNIMOD:35",
                "MT": "variable"
            }
        """
        ontology_dict: Dict[str, str] = {}
        ontology_split: str = ontology_str.split(";")
        for elem in ontology_split:
            elem = elem.strip()
            if elem == "":
                continue
            elem_split = elem.split("=")
            if len(elem_split) != 2:
                raise ValueError(f"Invalid ontology string: {ontology_str}")
            ontology_dict[elem_split[0].strip()] = elem_split[1].strip()

        return ontology_dict
    
    

PATH_TO_SDRF = '/home/uladzislau/Documents/Projects/O2TD_2022_project/SDRF_annotation/sdrf_final_test.sdrf.tsv'

test_df = pd.read_csv(PATH_TO_SDRF)

diann_params = '--out --fixed-mod'

diann_convert = DiannConverter(test_df, diann_params)

PATH_TO_SDRF_TYPE = diann_convert.check_sdrf_path(PATH_TO_SDRF)
sdrf = diann_convert.read_sdrf(PATH_TO_SDRF_TYPE)