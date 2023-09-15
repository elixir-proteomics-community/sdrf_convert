"""
Abstract class for SDRF file converters.
"""

# std imports
import argparse
from io import IOBase
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Set, Type, Union

# 3rd party imports
import pandas as pd
from pyteomics import mass

SDRF_CELL_SEPARATOR: str = "\t"
"""Separator used in SDRF file
"""

class AbstractConverter:
    """
    Abstract class for SDRF file converters.
    """

    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {}
    """Column names (key) and dtypes (value).
    If a column may occurs multiple times, add a star at the end of the column name
    (e.g. "comment[modification parameters]*")
    """

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {}
    """Columns which are optional.
    If a column may occurs multiple times, add a star at the end of the column name
    (e.g. "comment[modification parameters]*")
    """

    def __init__(self):
        """
        Initialize class.
        """

        self.sdrf_df: pd.DataFrame = pd.DataFrame()
        self.present_optional_columns: Set[str] = set()
        self.unimod_db: mass.Unimod = None

    @classmethod
    def find_columns(cls, sdrf_df: pd.DataFrame, col_name: str) -> List[str]:
        """
        Find duplicate columns in SDRF file.

        Returns
        -------
        List[str]
            List of duplicate columns pandas will mark them with a suffix 
            (e.g. "comment[modification parameters].1")
        """
        if not col_name.endswith("*"):
            return [col_name]
        plain_col_name: str = col_name[:-1]
        return list(filter(
            lambda col: col.startswith(plain_col_name), sdrf_df.columns
        ))

    @classmethod
    def check_necessary_columns_and_types(cls, sdrf_df: pd.DataFrame) -> None:
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
        for col_name, col_types in cls.COLUMN_PROPERTIES.items():
            col_names: List[str] = cls.find_columns(sdrf_df, col_name)
            for col_name in col_names:
                if not col_name in sdrf_df.columns:
                    raise ValueError(f"Column `{col_name}` not found in SDRF file")
                if not sdrf_df[col_name].dtype in col_types:
                    raise TypeError((
                        f"Column `{col_name}` has wrong type. "
                        f"Expected on of {col_types} but got `{sdrf_df[col_name].dtype}`"
                    ))

    @classmethod
    def check_optional_columns_and_types(cls, sdrf_df: pd.DataFrame) -> Set[str]:
        """
        Checks if optional columns are present in the SDRF file and if they have the correct type.

        Returns
        -------
        Set[str]
            Set of optional columns which are present in the SDRF file

        Raises
        ------
        TypeError
            If column dtype is not correct
        """
        present_optional_columns: Set[str] = set()
        for org_col_name, col_types in cls.COLUMN_PROPERTIES.items():
            col_names: List[str] = cls.find_columns(sdrf_df, org_col_name)
            for col_name in col_names:
                if col_name in sdrf_df.columns:
                    present_optional_columns.add(org_col_name)
                    # If column is present but with incorrect type
                    if not sdrf_df[col_name].dtype in col_types:
                        raise TypeError((
                            f"Column {col_name} has wrong type. "
                            f"Expected {col_types} but got {sdrf_df[col_name].dtype}"
                        ))
        return present_optional_columns

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
            elem_split = elem.split("=", maxsplit=1)
            if len(elem_split) != 2:
                raise ValueError(f"Invalid ontology string: {ontology_str}")
            ontology_dict[elem_split[0].strip()] = elem_split[1].strip()

        return ontology_dict

    @classmethod
    def read_sdrf(cls, sdrf: Union[pd.DataFrame, IOBase, Path]) -> pd.DataFrame:
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
        return pd.read_csv(
            sdrf,
            sep=SDRF_CELL_SEPARATOR,
            dtype={
                col: dtypes[0]
                for col, dtypes in cls.COLUMN_PROPERTIES.items()
            } # Use the first accepted dtype for each column
        )
    
    def get_unimod_from_NT(self, nt_entry: str) -> Dict[str, Any]:
        """
        Uses the given nt_entry for a query to Unimod (by pyteomics) and
        returns the resulting entry, or an empty record, if no corresponding
        entry was found in Unimod

        Parameters
        ----------
        nt_entry : str
            The title of the modification, like Oxidation (given by NT in the SDRF columns)

        Returns
        -------
        Dict[str, Any]:
            an Unimod entry
        """
        if self.unimod_db == None:
            self.unimod_db = mass.Unimod()
        
        mod = self.unimod_db.by_title(nt_entry)
        if (len(mod) == 0) :
            mod = self.unimod_db.by_name(nt_entry)
        
        return mod

    def init_converter(self, sdrf: Union[pd.DataFrame, IOBase, Path]):
        """
        1. Reads the SDRF and sets it to self.sdrf_df
        2. Validates the SDRF
        3. Stores the present optional columns in self.present_optional_columns
        (It is intentionally not the constructor so it is easier
        to call the convert methods multiple time
        for different SDRF with the same converter settings. This could also be a set-property)
        
        Raises
        ------
        ValueError
            If column is not present in SDRF file
        TypeError
            If column dtype is not correct or if sdrf is not a pandas DataFrame,
            file handler or path
        """
        self.present_optional_columns = set()
        self.sdrf_df = self.read_sdrf(sdrf)
        self.check_necessary_columns_and_types(self.sdrf_df)
        self.present_optional_columns = self.check_optional_columns_and_types(self.sdrf_df)


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
        raise NotImplementedError(f"Method convert() not implemented for {self.__class__.__name__}")
    
    @classmethod
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        """Uses the CLI arguments convert a SDRF file to the tool config.
        Example for an implemented subclass:
        ```python
        
        def convert_via_cli(cli_args: argparse.Namespace):
            sdrf_path = Path(cli_args.sdrf)
            positional_arg1 = cli_args.positional_arg1
            converter = ToolConverter(positional_arg1)
            converter.convert(sdrf_path)
        ```

        Parameters
        ----------
        cli_args : argparse.Namespace
            Collected CLI arguments
        """
        raise NotImplementedError("Method convert_via_cli() not implemented for AbstractConverter")
    
    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        """
        Add additional CLI arguments for the subclass. You do not need to add a SDRF-file parameter, as this is already done in the main CLI.
        Example for an implemented subclass:
        ```python
        def add_cli_args(subparsers: argparse._SubParsersAction):
            # Create a parser for the subclass
            tool_parser = subparsers.add_parser("<TOOLNAME>", help="SDRF to config converter for <TOOLNAME>")
            # Adding parameters for the specific tool
            tool_parser.add_argument("positional-arg1", type=str, help="Positional arument 1")
            tool_parser.add_argument("--additional-arg1", "-a", required=True, type=int, default=1, help="An additional argument")
            # Add the convert_via_cli method as callback
            tool_parser.set_defaults(func=cls.convert_via_cli)
        ```
        `convert_via_cli` will automatically be called with the parsed arguments.
        Example for calling the CLI: `python -m sdrf_convert <TOOLNAME> <positional-arg1> --additional-arg1 2`

        Parameters
        ----------
        subparsers : argparse._SubParsersAction
            Subparsers for the subclass
        """
        raise NotImplementedError(f"Method add_cli_args() not implemented for {cls.__name__}")
