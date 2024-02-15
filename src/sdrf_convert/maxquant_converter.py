"""
Module for converting SDRF files to MaxQuant xml param files
"""

# std imports
import argparse
from io import IOBase
from pathlib import Path
from typing import ClassVar, Dict, Iterator, List, Tuple, Type, Union
import xml.etree.ElementTree as ET

# 3rd party imports
import pandas as pd
from pyteomics import mass

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter



class maxquantConverter(AbstractConverter):
    """
    SDRF convert for MaxQuant search engine
    """

    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "source name": [pd.StringDtype()],
        "assay name": [pd.StringDtype()],
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
        "comment[data file]": [pd.StringDtype()],
        "comment[cleavage agent details]": [pd.StringDtype()],
        "characteristics[organism]": [pd.StringDtype()]
    }

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[modification parameters]": [pd.StringDtype()],
    }

    


    def __init__(self,
                 raw_files_path: str = "./",
                 fasta_file: str = "/path/to/FASTA",
                 maxquant_files_path: str = "./",
                 output_path: str = "./maxquant_identification_",
                 max_missed_cleavages: int = 2,
                 max_threads: int = 4,
                 group_similar_searches: bool = True):
        """
        Creates a new instance of the CometConverter class

        Parameters
        ----------
        group_similar_searches : bool, optional
            If True samples with equal search parameters will be grouped 
            and only one Comet params file and CLI call per group is created.
            If False each it will be grouped by assay name. By default False
        """
        super().__init__()

        self.raw_files_path : str = raw_files_path
        self.fasta_file : str = fasta_file
        self.maxquant_files_path : str = maxquant_files_path
        self.output_path : str = output_path
        self.max_missed_cleavages: int = max_missed_cleavages
        self.max_threads: int = max_threads
        self.group_similar_searches: bool = group_similar_searches




    def convert_row(self, row_idx: int, raw_file_rows: int) -> str:
        """
        Creates an maxquant params file the sample at the given row of the SDRF file.
        

        Parameters
        ----------
        row_idx : int

        Returns
        -------
        str
            representation of the maxquant params file
        """

        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]

        # initiate the XML tree for the params file
        xml_read = ET.parse(Path(self.maxquant_files_path))
        xml_params = xml_read.getroot()
        
        # the spectrum file
        try:
            raw_files_path:Path = Path(self.raw_files_path)
            raw_files_path_str = str(raw_files_path.absolute())
            xml_filePaths = xml_params.find('filePaths')
        except:
            raise TypeError("path to MzMLs doesn't exist")
        for file in raw_file_rows:
            related_row: pd.DataFrame = self.sdrf_df.iloc[file]
            print(related_row['comment[data file]'])
            ET.SubElement(xml_filePaths, "string").text = str(raw_files_path_str + "/" + related_row['comment[data file]'])


        #fasta file details
        fasta_path : Path = Path(self.fasta_file)
        xml_params.find('fastaFiles').find('FastaFileInfo').find('fastaFilePath').text = str(fasta_path)
    
        
        # set peptide mass tolerance and check unit
        parameterGroup = xml_params.find('parameterGroups').find('parameterGroup')
        """
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        try:
            fragment_mass_tolerance = float(fragment_mass_tolerance_split[0])
        except:
            raise TypeError("Fragment mass tolerance must be numeric, given was '{tol}".format(tol=fragment_mass_tolerance_split[0]))
        
        fragment_mass_tolerance_unit = fragment_mass_tolerance_split[1].strip()
        if fragment_mass_tolerance_unit.casefold() == "ppm":
            parameterGroup.find('firstSearchTol').text = fragment_mass_tolerance
        else:
            raise SyntaxError("Fragment mass tolerance must be ppm or Daltons, but is {unit}".format(unit=fragment_mass_tolerance_unit))
        """
        

        # set number of threads
        xml_params.find('numThreads').text = str(self.max_threads)

        # modifications
        self.__parse_modifications(sample, xml_params)

        # missed cleavages
        parameterGroup.find('maxMissedCleavages').text = str(self.max_missed_cleavages)

        # cleavage agent details
        cleavage_agent_details_dict: dict = self.ontology_str_to_dict(
            sample['comment[cleavage agent details]']
        )
        if cleavage_agent_details_dict['NT']:
            ET.SubElement(parameterGroup.find('enzymes'), "string").text = cleavage_agent_details_dict['NT']
    


        # get the XML tree and format with indentations
        tree = ET.ElementTree(xml_params)
        ET.indent(tree)
        
        return """<?xml version="1.0" encoding="utf-8"?>\n""" + ET.tostring(xml_params, encoding='unicode')

    
    def __parse_modifications(self, sample: pd.DataFrame, xml_params: ET) -> None:
        """
        Parses the modifications of the given sample (row) and adds the
        information
        """

        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):            
            modification_dict = self.ontology_str_to_dict(sample[col_name])
            # Mods falls die Masse direkt in der SDRF gegeben ist

            parameterGroup = xml_params.find('parameterGroups').find('parameterGroup')
            type = "variable"
            if modification_dict['MT'] and modification_dict['MT'].casefold() == "fixed":
                type = "fixed"
            targets_in_row = modification_dict['TA'].split(",")
            for target in targets_in_row:
                if len(target) != 1:
                    target = target.strip("[']")
                if type == "variable":
                    ET.SubElement(parameterGroup.find('variableModifications'), "string").text = str(modification_dict['NT']+" ("+target+")")
                if type == "fixed":
                    ET.SubElement(parameterGroup.find('fixedModifications'), "string").text  = str(modification_dict['NT']+" ("+target+")")
            

        
            


    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[Tuple[str, str]]:
        """
        Convert SDRF file to MaxQuant input file.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path object
        
        Returns
        -------
        (str)
            The Maxquant parameter files for each row/ each groupting
        """

        # init
        self.init_converter(sdrf)

        #group the sdrf entrys by their search Parameters
        grouping: pd.core.groupby.DataFrameGroupBy = self.get_grouping()
        
        #create and return config.json for each grouping
        for _, rows in grouping.groups.items():
            grouped_df = self.sdrf_df.iloc[rows]
            yield grouped_df['comment[data file]'] + ".maxquant_input.xml", self.convert_row(grouped_df.index[0], grouped_df.index)


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
        search_param_cols.remove("assay name")
        search_param_cols.remove('source name')
        search_param_cols.remove('comment[data file]')
        return self.sdrf_df.groupby(list(search_param_cols))

    @classmethod
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        # Read the MaxQuant params
        converter = maxquantConverter(
            raw_files_path=cli_args.raw_files_path,
            fasta_file=cli_args.fasta,
            maxquant_files_path=cli_args.maxquant_files_path,
            output_path=cli_args.output_path,
            max_missed_cleavages=cli_args.max_missed_cleavages,
            max_threads=cli_args.max_threads,
            group_similar_searches=cli_args.group_similar_searches
        )


        results = converter.convert(Path(cli_args.sdrf_file))
        i = 0
        
        for params in results:
            file_path : Path = Path(converter.output_path + "/" + params[0][0])
            with open(file_path, mode='wt') as file:
                file.write(params[1])
                if (i > 0):
                    print(str(file_path.absolute()))
            i += 1
        


    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("maxquant", help="SDRF parser to config parameter files for Maxquant")

        tool_parser.add_argument(
            "-r", "--raw-files-path", default="./", type=str, action="store", metavar="/path/to/RAWs", help="Path to the place, where the RAW (or mzML) files are"
        )
        tool_parser.add_argument(
            "-f", "--fasta", default="/path/to/FASTA", type=str, action="store", metavar="/path/to/FASTA", help="Path to the used FASTA file for identification"
        )
        tool_parser.add_argument(
            "-x", "--maxquant-files-path", default="./", type=str, action="store", metavar="/path/to/maxquant_files", help="Path to the position to the original Maxquant XML file"
        )
        tool_parser.add_argument(
            "-o", "--output-path", default="./maxquant_identification_", type=str, action="store", metavar="/path/to/output/prefix", help="Output path for the identification results, including prefix for the files"
        )
        tool_parser.add_argument(
            "-m", "--max-missed-cleavages", default=2, type=int, action="store", metavar="X", help="Maximum missed cleavages for the search (not provided by SDRF)"
        )
        tool_parser.add_argument(
            "-t", "--max-threads", default=4, type=int, action="store", metavar="X", help="Maximum number of used threads for teh identification"
        )
        tool_parser.add_argument(
            "-g", "--group-similar-searches", default=True, type=bool, action="store", metavar="X", help="Group equal search settings in one Config.json"
        )
        tool_parser.set_defaults(func=cls.convert_via_cli)
