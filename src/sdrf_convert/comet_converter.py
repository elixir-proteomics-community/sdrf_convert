"""
Module for converting SDRF files to Comet params files and CLI calls
"""

# std imports
import argparse
from copy import deepcopy
from io import IOBase
from pathlib import Path
import re
from typing import ClassVar, Dict, Iterator, List, Set, Tuple, Type, Union

# 3rd party imports
import pandas as pd

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter

class CometConverter(AbstractConverter):
    """
    SDRF convert for Comet search engine
    """
    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "source name": [pd.StringDtype()],
        "assay name": [pd.StringDtype()],
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
        "comment[data file]": [pd.StringDtype()],
        "comment[cleavage agent details]": [pd.StringDtype()]
    }

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[modification parameters]": [pd.StringDtype()],
    }

    COMET_UNITS: ClassVar[Dict[str, int]] = {
        "amu": 0,
        "mmu": 1,
        "ppm": 2,
        "da": 0
    }
    """Comet representation of mass tolerance units
    """

    COMET_ENZYM_INFO_MAP: ClassVar[Dict[str, int]] = {
        "MS:1001045": 0,    # unspecified cleavage / Cut_everywhere TODO: check if this is correct
        "MS:1001251": 1,    # Trypsin
        "MS:1001313": 2,    # Trypsin/P
        "MS:1001309": 3,    # LysC
        "MS:1003093": 4,    # LysN
        "MS:1001303": 5,    # ArgC
        "MS:1001304": 6,    # AspN
        "MS:1001307": 7,    # CNBr
        "MS:1001917": 8,    # GluC
        "MS:1001311": 9,    # PepsinA
        "MS:1001306": 10,   # Chymotrypsin
        "MS:1001955": 11    # No cleavage / No_cut
    }
    """Cleavage enzymes mapped to COMET_ENZYM_INFO
    """


    def __init__(
        self,
        comet_params: str,
        group_similar_searches: bool = False,

    ):
        """
        Creates a new instance of the CometConverter class

        Parameters
        ----------
        sdrf_file : Union[pd.DataFrame, IOBase, Path]
            _description_
        comet_params : str
            _description_
        group_similar_searches : bool, optional
            If True samples with equal search parameters will be grouped 
            and only one Comet params file and CLI call per group is created.
            If False each it will be grouped by assay name. By default False
        """
        super().__init__()
        self.comet_params = comet_params
        # self.comet_params = self.cleanup_params(self.comet_params) # TODO implement cleanup
        self.group_similar_searches = group_similar_searches



    def get_escaped_basename_of_data_file(self, data_file: str) -> str:
        """
        Returns the name of the data file with escaped spaces.

        Parameters
        ----------
        data_file : str
            Data file as stated in SDRF

        Returns
        -------
        str
            Whitespace escaped data file name
        """
        sample_file_path = Path(data_file)
        sample_file_name = sample_file_path.name
        return sample_file_name.replace(" ", "\\ ")

    def cleanup_params(self, params: str) -> str:
        """
        Some parameters are preconfigured and need should be cleaned up
        before values from SDRF are inserted.

        Parameters
        ----------
        params : str
            Comet parameter file content

        Returns
        -------
        str
            Cleaned up Comet parameter file content
        """
        raise NotImplementedError("Not implemented yet")

    def convert_modifications(self, sample: pd.DataFrame) -> Iterator[Tuple[re.Pattern, str]]:
        """
        Converts modifications from SDRF to Comet format.

        Parameters
        ----------
        sample : pd.DataFrame
            Row from SDRF

        Yields
        ------
        Iterator[Tuple[re.Pattern, str]]
            Pattern for replacing the modification in the params file and th
        """
        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):
            modification_dict = self.ontology_str_to_dict(sample[col_name])
            _mod_accession = modification_dict['AC']
            # TODO: Include pyteomics for unimod lookup
            # TODO: Add converter parameter for custom modifications
            raise NotImplementedError("Not implemented yet")

    def convert_sample(self, row_idx: int) -> Tuple[str, str]:
        """
        Creates a Comet CLI call and params file for the given sample.

        Parameters
        ----------
        row_idx : int

        Returns
        -------
        Tuple[str, str]
            _description_
        """
        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]
        sample_config: str = deepcopy(self.comet_params)

        # set cleavage enzyme
        cleavage_agent_details_dict = self.ontology_str_to_dict(
            sample['comment[cleavage agent details]']
        )
        cleavage_agent_num: int = self.COMET_ENZYM_INFO_MAP[
            cleavage_agent_details_dict['AC']
        ]
        sample_config = re.sub(
            r"search_enzyme_number = \d",
            f"search_enzyme_number = {cleavage_agent_num}", sample_config
        )

        # set peptide mass tolerance and unit
        peptide_mass_tolerance_split: str = sample['comment[precursor mass tolerance]'].split()
        peptide_mass_tolerance = float(peptide_mass_tolerance_split[0])
        peptide_mass_tolerance_unit = self.COMET_UNITS[peptide_mass_tolerance_split[1]]
        sample_config = re.sub(
            r"peptide_mass_tolerance = \d+(.\d+){0,1}",
            f"peptide_mass_tolerance = {peptide_mass_tolerance}", sample_config
        )
        sample_config = re.sub(
            r"peptide_mass_units = \d",
            f"peptide_mass_units = {peptide_mass_tolerance_unit}", sample_config
        )

        # set fragment mass tolerance and unit
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        fragment_mass_tolerance = float(fragment_mass_tolerance_split[0]) / 20
        fragment_mass_tolerance_unit = fragment_mass_tolerance_split[1]
        if fragment_mass_tolerance_unit.lower()== "amu":
            fragment_mass_tolerance = fragment_mass_tolerance / 1000
        sample_config = re.sub(
            r"peptide_mass_tolerance = \d.\d",
            f"peptide_mass_tolerance = {fragment_mass_tolerance}", sample_config
        )

        sample_file_name = self.get_escaped_basename_of_data_file(sample['comment[data file]'])
        cli_string: str = f"-P<PARAMS> -D<FASTA> {sample_file_name}"

        return cli_string, sample_config


    def get_grouping(self) -> pd.core.groupby.DataFrameGroupBy:
        """
        Groups samples by assay nae or by similar search parameters.

        Returns
        -------
        pd.core.groupby.DataFrameGroupBy
            Grouped samples
        """
        if not self.group_similar_searches:
            return self.sdrf_df.groupby(["assay name"])
        # Collect columns with search params
        search_param_cols: Set[str] = set(self.COLUMN_PROPERTIES.keys())
        search_param_cols.update(self.present_optional_columns)
        search_param_cols.remove('assay name')
        search_param_cols.remove('source name')
        search_param_cols.remove('comment[data file]')
        return self.sdrf_df.groupby(list(search_param_cols))


    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[Tuple[str, str]]:
        """
        Convert SDRF file to Comet CLI string and corresponding input file.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path
        
        Returns
        -------
        (str, str)
            Tuple containing the CLI call and content of the Comet input files
            (cli call, comet.params content)
        """
        # init
        self.init_converter(sdrf)

        # Group rows by search params
        grouping: pd.core.groupby.DataFrameGroupBy = self.get_grouping()
        for _, rows in grouping.groups.items():
            grouped_df = self.sdrf_df.iloc[rows]
            # create CLI call and params based on first row of group...
            cli_str, params = self.convert_sample(grouped_df.index[0])
            # ... and append the data files of the other rows
            for row_idx in grouped_df.index[1:]:
                cli_str += " " + self.get_escaped_basename_of_data_file(
                    grouped_df.loc[row_idx, "comment[data file]"]
                )
            yield cli_str, params

    @classmethod
    def convert_via_cli(cli_args: argparse.Namespace):
        # Read the comet params
        comet_params: str = Path(cli_args.comet_params).read_text()
        converter = CometConverter(
            comet_params,
            cli_args.group_similar_searches
        )
        converter.convert(cli_args.sdrf_file)
    
    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("comet", help="SDRF to config converter for Comet")
        tool_parser.add_argument("comet_params", help="Comet params file")
        tool_parser.add_argument(
            "-g", "--group-similar-searches", default=False, action="store_true", help="Group samples with equal search parameters in one config file and CLI command"
        )
        tool_parser.set_defaults(func=cls.convert_via_cli)