# std imports
import argparse
import importlib
from pathlib import Path

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter

# Import all subclassed converters automatically so they are available in AbstractConverter.__subclasses__()
# and can be automatically added to the CLI
for converter_path in Path(__file__).parent.glob('*_converter.py'):
    if converter_path.stem == 'abstract_converter':
        continue
    importlib.import_module(f"sdrf_convert.{converter_path.stem}")

class CommandLineInterface:
    """
    Command line interface for the sdrf_convert package
    """

    def __init__(self):
        self.__parser = argparse.ArgumentParser(description='Convert an SDRF to your favorite tool\'s config')
        self.__parser.add_argument('sdrf_file', help='Path to the SDRF file')
        #self.__parser.set_defaults(func=lambda: self.__parser.print_help())
        subparsers = self.__parser.add_subparsers()

        # Iterate over all subclasses of AbstractConverter and add their CLI arguments
        for converter in AbstractConverter.__subclasses__():
            converter.add_cli_args(subparsers)

    def start(self):
        """
        Converts command line arguments and calls callback function of the selected converter
        """
        args = self.__parser.parse_args()
        args.func(args)
