"""
Module for converting SDRF files to Comet params files and CLI calls
"""

# std imports
import argparse
from copy import deepcopy
from enum import IntEnum, unique
from io import IOBase
from pathlib import Path
import re
from typing import Any, ClassVar, Dict, Iterator, List, Set, Tuple, Type, Union

# 3rd party imports
import pandas as pd
from pyteomics.mass.unimod import Unimod

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter

@unique
class CometVarModTerminal(IntEnum):
    """
    Enum for Comet variable modification terminal constraints
    """
    PROTEIN_N = 0
    PROTEIN_C = 1
    PEPTIDE_N = 2
    PEPTIDE_C = 3

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

    ONTOLOGY_ID_COMET_ENZYM_MAP: ClassVar[Dict[str, int]] = {
        "MS:1001956": 0,    # unspecified cleavage / Cut_everywhere
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
        "MS:1001955": 11,    # No cleavage / No_cut
    }
    """Cleavage enzymes ontology mapped to COMET_ENZYM_INFO
    """

    ONTOLOGY_TERM_COMET_ENZYM_MAP: ClassVar[Dict[str, int]] = {
        "unspecific cleavage": 0,       # unspecified cleavage / Cut_everywhere
        "trypsin": 1,                   # Trypsin
        "trypsin/p": 2,                 # Trypsin/P
        "lys-c": 3,                     # LysC
        "trypsin/K": 3,                 # LysC
        "lys-n": 4,                     # LysN
        "arg-c": 5,                     # ArgC
        "clostripain": 5,               # ArgC
        "trypsin/r": 5,                 # ArgC
        "asp-n": 6,                     # AspN
        "cnbr": 7,                      # CNBr
        "glutamyl endopeptidase": 8,    # GluC
        "gluc": 8,                      # GluC
        "staphylococcal protease": 8,   # GluC
        "pepsina": 9,                   # PepsinA
        "chymotrypsin": 10,             # Chymotrypsin
        "no cleavage": 11,              # No cleavage / No_cut
    }
    """Cleavage enzymes names and synonymes mapped to COMET_ENZYM_INFO
    """

    MAX_VAR_MODS: ClassVar[int] = 9
    """Maximum number of variable modifications per sample
    """

    def __init__(
        self,
        comet_params: str,
        max_variable_modification: int,
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
        self.comet_params = self.cleanup_params(comet_params)
        self.max_variable_modification = max_variable_modification
        self.group_similar_searches = group_similar_searches
        self.unimod = Unimod()



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
    
    @classmethod
    def set_param(cls, config: str, param_name: str, value: str) -> str:
        """
        Set a parameter in a Comet parameter file.
        Keeps the comment if present.

        Parameters
        ----------
        config : str
            Comet parameter file content
        param_name : str
            Name of the parameter to set
        value : str
            Value to set
        
        Returns
        -------
        str
        """

        pattern: re.Pattern = re.compile(fr"^{param_name} = .+?$", re.MULTILINE)
        return re.sub(pattern, fr"{param_name} = {value}", config)

    def cleanup_params(self, config: str) -> str:
        """
        Some parameters are preconfigured and should be cleaned up
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

        # Unset default modifications
        cleaned_config: str = self.set_param(config, "variable_mod01", "0.0 X 0 3 -1 0 0 0.0")
        return self.set_param(cleaned_config, "add_C_cysteine", "0.0000")


    def convert_modifications(self, sample: pd.DataFrame) -> Iterator[Tuple[re.Pattern, str]]:
        """
        Converts modifications from SDRF to Comet format.
        It is necessary to to now the monoisotopic mass (MM) of the modification and if it is fixed or variable (MT).
        If MT is missing an error is raised. If MM is missing, the Unimod database is queried. If the modification is not found an error is raised.

        Parameters
        ----------
        sample : pd.DataFrame
            Row from SDRF

        Yields
        ------
        Iterator[Tuple[re.Pattern, str]]
            Pattern for finding the param and replacement including the value
        """
        # Convert modifications
        fix_mod_table: List[Dict[str, Any]] = []
        var_mod_table: List[Dict[str, Any]] = []
        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):
            mod = self.ontology_str_to_dict(sample[col_name])
            # if MT is missing, raise error as the information can not be found in Unimod
            if 'MT' not in mod:
                raise ValueError((
                    f"Invalid modification parameters in sample {sample['source name']}: {sample[col_name]}. "
                    "MT attribute is not mandatory for SDRF but for Comet we need to now."
                )) 
            # if MM is missing, query Unimod and check if the modification is found
            if 'MM' not in mod:
                umod = self.unimod.get(mod['NT'])
                if umod is not None:
                    mod['MM'] = umod.monoisotopic_mass
                else:
                    raise ValueError((
                        f"Invalid modification parameters in sample {sample['source name']}: {sample[col_name]}. "
                        "MM attribute is not mandatory for SDRF but for Comet we need to now. "
                        f"Could not find Unimod entry for {mod['NT']}."
                    ))
            # Sort into fixed and variable modifications
            if mod['MT'].lower() == 'fixed':
                fix_mod_table.append(mod)
            else:
                var_mod_table.append(mod)

        # Convert modifications to Comet format
        var_mod_counter: int = 0
        for mod in var_mod_table:
            # Comma separated list of target amino acids
            targets = ",".join(self.get_plain_modification_targets(mod["TA"]))
            # Distance constraint to terminal, default -1 for anywhere
            dist_constraint = -1
            # List if terminal constraints. For any N/C-term, two separate entries are needed in the params file
            terminals = [CometVarModTerminal.PROTEIN_N] # only protein N-term, ignored when dist_constraint is -1
            if 'PP' in mod:
                match mod['PP'].lower():
                    case "protein n-term":
                        dist_constraint = 0
                        terminals = [CometVarModTerminal.PROTEIN_N]
                    case "protein c-term":
                        dist_constraint = 0
                        terminals = [CometVarModTerminal.PROTEIN_C]
                    case "any n-term":
                        dist_constraint = 0
                        terminals = [CometVarModTerminal.PROTEIN_N, CometVarModTerminal.PEPTIDE_N]
                    case "any c-term":
                        dist_constraint = 0
                        terminals = [CometVarModTerminal.PROTEIN_C, CometVarModTerminal.PEPTIDE_C]

            for terminal in terminals:
                var_mod_counter += 1
                # Check if we exceed the maximum number of variable modifications
                if var_mod_counter > self.MAX_VAR_MODS:
                    raise ValueError(f"Too many variable modifications in sample {sample['source name']}")
                # Yield match pattern and replacement for variable modifications
                yield (
                    re.compile(fr"^(variable_mod0{var_mod_counter} =) .+?$", re.MULTILINE),
                    fr"\1 {mod['MM']} {targets} 0 {self.max_variable_modification} {dist_constraint} {terminal.value} 0 0.0"
                )

        # Yield match patterns and replacements for fixed modifications 
        for mod in fix_mod_table:
            for target in self.get_plain_modification_targets(mod["TA"]):
                yield (
                    re.compile(fr"^(add_{target}_[a-z]+ =) .+?$", re.MULTILINE),
                    fr"\1 {mod['MM']}"
                )


    def convert_sample(self, row_idx: int) -> Tuple[str, str, str]:
        """
        Creates a Comet CLI call and params file for the given sample.

        Parameters
        ----------
        row_idx : int

        Returns
        -------
        Tuple[str, str, str]
            Tuple containing the CLI call, content of the Comet param file and the base name of the data file
        """
        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]
        sample_config: str = deepcopy(self.comet_params)

        # set cleavage enzyme
        cleavage_agent_details_dict = self.ontology_str_to_dict(
            sample['comment[cleavage agent details]']
        )
        cleavage_agent_lookup_errors: List[str] = []
        cleavage_agent_num: int = -1
        # Check AC first
        if 'AC' in cleavage_agent_details_dict:
            try:
                cleavage_agent_num = self.ONTOLOGY_ID_COMET_ENZYM_MAP[
                    cleavage_agent_details_dict['AC']
                ]
            except KeyError:
                cleavage_agent_lookup_errors.append(f"AC '{cleavage_agent_details_dict['AC']}' is unknown for Comet")
        # if cleavage agent is not found by AC, try NT
        if cleavage_agent_num == -1 and 'NT' in cleavage_agent_details_dict:
            try:
                cleavage_agent_num = self.ONTOLOGY_TERM_COMET_ENZYM_MAP[
                    cleavage_agent_details_dict['NT'].lower()
                ]
            except KeyError:
                cleavage_agent_lookup_errors.append(f"NT '{cleavage_agent_details_dict['NT']}' is unknown for Comet")
        # if cleavage agent is not found by AC or NT, raise error
        if cleavage_agent_num == -1:
            cleavage_agent_lookup_errors_message = "\n\t- ".join(cleavage_agent_lookup_errors)
            raise ValueError(f"Invalid cleavage agent: {sample['comment[cleavage agent details]']}.\n\t- {cleavage_agent_lookup_errors_message}")
        sample_config = self.set_param(
            sample_config,
            "search_enzyme_number",
            str(cleavage_agent_num)
        )

        # set peptide mass tolerance and unit
        peptide_mass_tolerance_split: str = sample['comment[precursor mass tolerance]'].split()
        peptide_mass_tolerance = float(peptide_mass_tolerance_split[0])
        peptide_mass_tolerance_unit = self.COMET_UNITS[peptide_mass_tolerance_split[1]]
        sample_config = self.set_param(
            sample_config,
            "peptide_mass_tolerance",
            str(peptide_mass_tolerance)
        )
        sample_config = self.set_param(
            sample_config,
            "peptide_mass_units",
            str(peptide_mass_tolerance_unit)
        )

        # set fragment mass tolerance and unit
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        fragment_mass_tolerance = float(fragment_mass_tolerance_split[0]) / 20
        fragment_mass_tolerance_unit = fragment_mass_tolerance_split[1]
        if fragment_mass_tolerance_unit.lower()== "amu":
            fragment_mass_tolerance = fragment_mass_tolerance / 1000
        sample_config = self.set_param(
            sample_config,
            "peptide_mass_tolerance",
            str(fragment_mass_tolerance)
        )

        # set modifications
        sample_config = self.set_param(
            sample_config,
            "max_variable_mods_in_peptide",
            str(self.max_variable_modification)
        )
        for (pattern, replacement) in self.convert_modifications(sample):
            print(pattern, replacement)
            sample_config = re.sub(pattern, replacement, sample_config)

        sample_file_name = self.get_escaped_basename_of_data_file(sample['comment[data file]'])
        cli_string: str = f"-P<PARAMS> -D<FASTA> {sample_file_name}"

        return cli_string, sample_config, sample_file_name


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
    
    @classmethod
    def get_clean_group_label(cls, group_label: str) -> str:
        """
        Cleans group label from special characters

        Parameters
        ----------
        group_label : str
            Group label to clean

        Returns
        -------
        str
            Cleaned group label
        """
        return re.sub(r"[^a-zA-Z0-9]", "_", group_label)


    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[Tuple[str, str, str]]:
        """
        Convert SDRF file to Comet CLI string and corresponding input file.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path
        
        Returns
        -------
        (str, str, str)
            Tuple containing the CLI call, content of the Comet input files and the and ID for the params file
            if group_similar_searches is True, the ID is the group index, otherwise it is the sample file name.
            (cli call, comet.params content, params file ID)
        """
        # init
        self.init_converter(sdrf)

        # Group rows by search params
        grouping: pd.core.groupby.DataFrameGroupBy = self.get_grouping()
        for group_idx, (_, rows) in enumerate(grouping.groups.items()):
            grouped_df = self.sdrf_df.iloc[rows]
            # # Create an identifier. If not grouping by similar searches, use assay name
            # identifier = grouped_df.iloc[0]["assay name"] if not self.group_similar_searches else group_counter
            # create CLI call and params based on first row of group...
            cli_str, params, sample_file_name = self.convert_sample(grouped_df.index[0])
            # ... and append the data files of the other rows
            for row_idx in grouped_df.index[1:]:
                cli_str += " " + self.get_escaped_basename_of_data_file(
                    grouped_df.loc[row_idx, "comment[data file]"]
                )
            params_file_id: str = self.get_clean_group_label(f"{group_idx}") if self.group_similar_searches else f"{sample_file_name}"
            yield cli_str, params, params_file_id

    @classmethod
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        # Read the comet params
        comet_params: str = Path(cli_args.comet_params).read_text()
        converter = CometConverter(
            comet_params,
            cli_args.max_variable_modification,
            cli_args.group_similar_searches,
        )
        config_folder = Path(cli_args.config_folder)
        for cli_str, params, params_file_id in converter.convert(Path(cli_args.sdrf_file)):
            print(cli_str)
            config_folder.joinpath(f"{params_file_id}.comet.params").write_text(params, encoding="utf-8")
    
    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("comet", help="SDRF to config converter for Comet. Monoisotopic mass and information about fixed or variable is necessary for PTMs.")
        tool_parser.add_argument("max_variable_modification", type=int, help="Max. number of variable modifications per sample (will be applied to all var. modification)")
        tool_parser.add_argument("comet_params", help="Comet params file")
        tool_parser.add_argument("config_folder", help="Folder to store the generated config files")
        tool_parser.add_argument(
            "-g", "--group-similar-searches", default=False, action="store_true", help="Group samples with equal search parameters in one config file and CLI command"
        )
        tool_parser.set_defaults(func=cls.convert_via_cli)