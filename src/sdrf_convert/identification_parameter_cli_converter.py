"""
Module for converting SDRF files to Comet params files and CLI calls
"""

# std imports
import argparse
from io import IOBase
from pathlib import Path
import re
from typing import Any, ClassVar, Dict, Iterator, List, Type, Union, Tuple
import warnings

# 3rd party imports
import pandas as pd
from pyteomics.mass.unimod import Unimod

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter

class IdentificationParameterCliConverter(AbstractConverter):
    """
    SDRF convert for [IdentificationParametersCLI](http://compomics.github.io/projects/denovogui/wiki/IdentificationParametersCLI)
    """
    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
    }

    AMINO_ACID_PATTERN: ClassVar[re.Pattern] = re.compile(r"[A-Z]")
    """Pattern to find amino acids in TA attribute of modification
    """

    def __init__(self):
        """
        Creates a new instance of the FlashLFQConverter class
        """
        # pylint: disable=useless-parent-delegation
        # Subclassing is needed for automatic CLI generation
        super().__init__()

    @classmethod
    def parse_tolerance(cls, tolerance: pd.Series) -> List[Tuple[float, str]]:
        """
        Returns list of tuples containing the tolerance and the unit

        Arguments
        ---------
        tolerance : pd.Series
            Series containing the tolerance and the unit

        Returns
        -------
        List[Tuple(float, str)]
            List of tuples containing the tolerance and the unit
        """
        return list(map(
            lambda x: (float(x[0]), x[1]),
            map(
                lambda x: x.strip().split(" "),
                tolerance
            )
        ))
    
    def convert_modifications(self) -> Iterator[Tuple[bool, str]]:
        """
        Converts modifications from SDRF to IdentificationParametersCLI.


        Yields
        ------
        Iterator[Tuple[bool, str]]
            Iterator over tuples containing a boolean (true if fixed) and the modification string for IdentificationParametersCLI
        """
        mod_columns: List[str] = self.find_columns(self.sdrf_df, 'comment[modification parameters]*')
        if len(mod_columns) > 0 and self.unimod_db is None:
            # Ignore warnings from pyteomics when initializing the Unimod database
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.unimod_db = Unimod()
        for col_name in mod_columns:
            for mod_ontology in self.sdrf_df[col_name]:
                modification_dict = self.ontology_str_to_dict(mod_ontology)
                # TODO: Accession should include the onthology, e.g, UNIMOD:1
                mod_accession = modification_dict['AC']
                mod: Any = None
                try:
                    mod = self.unimod_db.get(int(mod_accession))
                except KeyError as err:
                    raise KeyError(f"Unknown Unimod accession: {mod_accession}") from err

                is_fixed: bool = modification_dict['MT'].lower() == "fixed"

                for amino_acid_match in self.AMINO_ACID_PATTERN.finditer(modification_dict['TA']):
                    yield (
                        is_fixed,
                        f"{mod.code_name} of {amino_acid_match.group().upper()}"
                    )

    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[str]:
        """
        Creates a IdentificationParametersCLI call from the SDRF file with is used with the PeptideShakerCLI.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path
        
        Returns
        -------
        Iterator[float]
            Iterator over the largest precursor mass tolerance in the SDRF file
        """
        # init
        tool_params: str = ""
        self.init_converter(sdrf)
        prec_tolerance = max(
            self.parse_tolerance(self.sdrf_df["comment[precursor mass tolerance]"]),
            key=lambda x: x[0]
        )
        prec_tol_unit: int  = -1
        match prec_tolerance[1]:
            case "ppm":
                prec_tol_unit = 0
            case "Da":
                prec_tol_unit = 1
            case _:
                raise ValueError(f"Unknown precursor mass tolerance unit {prec_tolerance[1]}")
        tool_params += f"-prec_tol {prec_tolerance[0]} -prec_ppm {prec_tol_unit} "
        frag_tolerance = max(
            self.parse_tolerance(self.sdrf_df["comment[fragment mass tolerance]"]),
            key=lambda x: x[0]
        )
        frag_tol_unit: int  = -1
        match frag_tolerance[1]:
            case "ppm":
                frag_tol_unit = 0
            case "Da":
                frag_tol_unit = 1
            case _:
                raise ValueError(f"Unknown fragment mass tolerance unit {frag_tolerance[1]}")
        tool_params += f"-frag_tol {frag_tolerance[0]} -frag_ppm {frag_tol_unit} "

        fixed_mods = set()
        variabel_mods = set()

        for is_fixed, mod in self.convert_modifications():
            if is_fixed:
                fixed_mods.add(mod)
            else:
                variabel_mods.add(mod)

        tool_params += f"-fixed_mods \"{','.join(fixed_mods)}\" -variable_mods \"{','.join(variabel_mods)}\" "

        yield tool_params

    @classmethod
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        # Read the comet params
        converter = cls()
        param_iter = converter.convert(Path(cli_args.sdrf_file))
        print(next(param_iter))

    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("identificationparameterscli", help="SDRF to config converter for FlashLFQ")
        tool_parser.set_defaults(func=cls.convert_via_cli)
