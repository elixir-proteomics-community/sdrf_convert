"""
Module for converting SDRF files to Comet params files and CLI calls
"""

# std imports
import argparse
from io import IOBase
from pathlib import Path
from typing import ClassVar, Dict, Iterator, List, Type, Union

# 3rd party imports
import pandas as pd

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter

class FlashLFQConverter(AbstractConverter):
    """
    SDRF convert for FlashLFQ
    """
    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[precursor mass tolerance]": [pd.StringDtype()],
    }

    def __init__(self):
        """
        Creates a new instance of the FlashLFQConverter class
        """
        super().__init__()

    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[float]:
        """
        Returns the largest precursor mass tolerance in the SDRF file as FlashLFQ get all PSMs in one files.

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
        self.init_converter(sdrf)
        tolerances = list(map(
            lambda x: x.strip().split(" ")[0],
                filter(
                    lambda x: x.strip().endswith("ppm"),
                    self.sdrf_df["comment[precursor mass tolerance]"]
                )
        ))
        if len(tolerances) == 0:
            raise ValueError("No precursor mass tolerance with unit PPM found in SDRF file")

        yield max(tolerances)

    @classmethod
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        # Read the comet params
        converter = cls()
        param_iter = converter.convert(Path(cli_args.sdrf_file))
        print(f"--ppm {next(param_iter)}")

    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("flashlfq", help="SDRF to config converter for FlashLFQ")
        tool_parser.set_defaults(func=cls.convert_via_cli)
