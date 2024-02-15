"""
Module for converting SDRF files to Fragger param files and CLI calls
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



class msfraggerConverter(AbstractConverter):
    """
    SDRF convert for msfragger search engine
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
                 msfragger_files_path: str = "./",
                 output_path: str = "./msfragger_identification_",
                 max_missed_cleavages: int = 2,
                 max_threads: int = 11,
                 group_similar_searches: bool = True,
                 decoy_tag: str = "_rev"):
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
        self.msfragger_files_path : str = msfragger_files_path
        self.output_path : str = output_path
        self.max_missed_cleavages: int = max_missed_cleavages
        self.max_threads: int = max_threads
        self.group_similar_searches: bool = group_similar_searches
        self.decoy_tag: str = decoy_tag




    def convert_row(self, row_idx: int, raw_file_rows: int) -> str:
        """
        Creates an msfragger params file the sample at the given row of the SDRF file.
        

        Parameters
        ----------
        row_idx : int

        Returns
        -------
        str
            representation of the msfragger params file
        """

        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]

        # initiate Array for the params file
        params = self.get_parameter(str(self.msfragger_files_path))

        #fasta file details
        fasta_path : Path = Path(self.fasta_file)
        params['database_name'] = str(fasta_path)
        
        
        # set peptide mass tolerance and check unit
        peptide_mass_tolerance_split: str = sample['comment[precursor mass tolerance]'].split()
        try:
            peptide_mass_tolerance = float(peptide_mass_tolerance_split[0])
        except:
            raise TypeError("Peptide mass tolerance must be numeric, given was '{tol}".format(tol=peptide_mass_tolerance_split[0]))

        params['precursor_mass_lower'] = "-"+str(peptide_mass_tolerance)			
        params['precursor_mass_upper'] = str(peptide_mass_tolerance)
        params['precursor_true_tolerance'] = str(peptide_mass_tolerance)
        peptide_mass_tolerance_unit = peptide_mass_tolerance_split[1].strip()
        if peptide_mass_tolerance_unit.casefold() == "ppm":
            params['precursor_mass_units'] = "1"
            params['precursor_true_units'] = "1"	
        elif peptide_mass_tolerance_unit.casefold() == "da":
            params['precursor_mass_units'] = "0"
            params['precursor_true_units'] = "0"	
        else:
            raise SyntaxError("Peptide mass tolerance must be ppm or Daltons, but is {unit}".format(unit=fragment_mass_tolerance_unit))
        
        # set fragment mass tolarance
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        try:
            fragment_mass_tolerance = float(fragment_mass_tolerance_split[0])
        except:
            raise TypeError("Fragment mass tolerance must be numeric, given was '{tol}".format(tol=fragment_mass_tolerance_split[0]))
        params["fragment_mass_tolerance"] = str(fragment_mass_tolerance)
        fragment_mass_tolerance_unit = fragment_mass_tolerance_split[1].strip()
        if fragment_mass_tolerance_unit.casefold() == "ppm":
            params['fragment_mass_units'] = "1"	
        elif fragment_mass_tolerance_unit.casefold() == "da":
            params['fragment_mass_units'] = "0"
        else:
            raise SyntaxError("Fragment mass tolerance must be ppm or Daltons, but is {unit}".format(unit=fragment_mass_tolerance_unit))
        			
        # set decoy prefix
        params['decoy_prefix'] = self.decoy_tag

        # set number of threads
        params['num_threads'] = str(self.max_threads)

        # modifications
        self.__parse_modifications(sample, params)

        # missed cleavages
        params['allowed_missed_cleavage_1'] = str(self.max_missed_cleavages)
        # cleavage agent details
        cleavage_agent_details_dict: dict = self.ontology_str_to_dict(
            sample['comment[cleavage agent details]']
        )
        if (cleavage_agent_details_dict['NT'] == 'Trypsin'):
            params['search_enzyme_name_1'] = 'Trypsin'
            params['search_enzyme_cut_1'] = 'KR'
            params['search_enzyme_nocut_1'] = 'P'
            params['search_enzyme_sense_1'] = ''
        else:
            raise SyntaxError("Could not parse cleavage agent from {cad}".format(cad=str(cleavage_agent_details_dict)))

        
        return params

    
    def __parse_modifications(self, sample: pd.DataFrame, params: str) -> None:
        """
        Parses the modifications of the given sample (row) and adds the
        information
        """
        i = 0
        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):            
            modification_dict = self.ontology_str_to_dict(sample[col_name])
            i = i + 1
            # Mods falls die Masse direkt in der SDRF gegeben ist
            type = "variable"
            if modification_dict['MT'] and modification_dict['MT'].casefold() == "fixed":
                type = "fixed"
            mod = self.get_unimod_from_NT(modification_dict['NT'])
            if type == "variable":
                mod_space = str("variable_mod_0"+str(i))
                if (modification_dict['NT'] == "Oxidation"):
                    params[mod_space] = str(mod['mono_mass']) + " M 2"
            

        
            


    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[Tuple[str, str]]:
        """
        Convert SDRF file to msfragger input file.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path object
        
        Returns
        -------
        (str)
            The msfragger parameter files for each row/ each groupting
        """

        # init
        self.init_converter(sdrf)

        #group the sdrf entrys by their search Parameters
        grouping: pd.core.groupby.DataFrameGroupBy = self.get_grouping()
        
        #create and return config.json for each grouping
        for _, rows in grouping.groups.items():
            grouped_df = self.sdrf_df.iloc[rows]
            yield grouped_df['comment[data file]'] + ".msfragger.params", self.convert_row(grouped_df.index[0], grouped_df.index)

    def get_parameter(self, file):
        result = {}
        with open(file) as fd:
            for line in fd:
                i=0
                key = ''
                value = ''
                try:
                    if line != "\n":
                        if "#" in line:
                            array = line.split('#')
                        else:
                            array = [line,""]
                        try:
                            if "=" in array[0]:
                                key, value = array[0].split('=')
                            key, value = key.strip(), value.strip()
                            result[key] = value
                        except Exception as m:
                            print(m)
                except Exception as e:
                    print(e)
        return result


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
        # Read the MSFragger params
        converter = msfraggerConverter(
            raw_files_path=cli_args.raw_files_path,
            fasta_file=cli_args.fasta,
            msfragger_files_path=cli_args.msfragger_files_path,
            output_path=cli_args.output_path,
            max_missed_cleavages=cli_args.max_missed_cleavages,
            max_threads=cli_args.max_threads,
            group_similar_searches=cli_args.group_similar_searches,
            decoy_tag = cli_args.decoy_tag
        )


        results = converter.convert(Path(cli_args.sdrf_file))
        
        i=0
        for params in results:
            for data in params[1]:
                file_path : Path = Path(converter.output_path + "/" + params[0][0])
                with open(file_path, mode='a') as file:  
                    file.write(data)
                    file.write(" = ")
                    file.write(params[1][data])
                    file.write("\n")
        


    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("msfragger", help="SDRF parser to config parameter files for msfragger")

        tool_parser.add_argument(
            "-r", "--raw-files-path", default="./", type=str, action="store", metavar="/path/to/RAWs", help="Path to the place, where the RAW (or mzML) files are"
        )
        tool_parser.add_argument(
            "-f", "--fasta", default="/path/to/FASTA", type=str, action="store", metavar="/path/to/FASTA", help="Path to the used FASTA file for identification"
        )
        tool_parser.add_argument(
            "-x", "--msfragger-files-path", default="./", type=str, action="store", metavar="/path/to/msfragger_files", help="Path to the position to the original msfragger XML file"
        )
        tool_parser.add_argument(
            "-o", "--output-path", default="./msfragger_identification_", type=str, action="store", metavar="/path/to/output/prefix", help="Output path for the identification results, including prefix for the files"
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
        tool_parser.add_argument(
            "-d", "--decoy-tag", default="rev_", type=str, action="store", metavar="X", help="decoy tag"
        )
        tool_parser.set_defaults(func=cls.convert_via_cli)
