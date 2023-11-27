"""
Module for converting SDRF files to Sage config file
"""

# std imports
import argparse
from io import IOBase
from pathlib import Path
from typing import ClassVar, Dict, List, Type, Union
import copy

# 3rd party imports
import json as json
import pandas as pd

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter



class SageConverter(AbstractConverter):
    """
    SDRF convert for Sage search engine
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


    '''vvv  Config Datei wie Sie als beispiel Datei auf in der SAGE Dokumentation angegeben ist.  vvv'''
    SAGE_ORIGINAL_CONFIG = {
        "database": {
            "bucket_size": 8192,
            "enzyme": {
                "missed_cleavages": 2,
                "min_len": 7,
                "max_len": 50,
                "cleave_at": "KR",
                "restrict": "P"
            },
            "fragment_min_mz": 150.0,
            "fragment_max_mz": 1500.0,
            "peptide_min_mass": 500.0,
            "peptide_max_mass": 5000.0,
            "ion_kinds": [
                "b",
                "y"
            ],
            "min_ion_index": 2,
            "max_variable_mods": 2,
            "static_mods": {
            },
            "variable_mods": {
            },
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": "human_contam.fasta"
        },
        "precursor_tol": {
        },
        "fragment_tol": {
        },
        "deisotope": True,
        "chimera": True,
        "min_peaks": 15,
        "max_peaks": 250,
        "min_matched_peaks": 4,
        "max_fragment_charge": 1,
        "report_psms": 2,
        "predict_rt": True,
        "output_directory": "s3://bucket/prefix",
        "mzml_paths": []
    }

    def __init__(self,
                 raw_files_path: str = "./",
                 fasta_file: str = "/path/to/FASTA",
                 output_path: str = "./Sage_identification",
                 allowed_miscleavages: int = 2, 
                 max_parent_charge: int = 4,
                 group_similar_searches: bool = True):
        """Creates a new instance of the Sage_Converter Class
        
        Parameters
        ------------------
        raw_files_path: Path to mzML Data
        fasta_file: Path to Fasta File
        output_path: Path to the Folder in wich the config.json and the Sage results will be saved
        allowed_miscleavages: Number of maximum missed cleaveges
        max_parent_charge: Maximum Parent charge
        group_similar_searches: If true mzMLs with equal Search Parameters will only produce one config.json
        
        """

        super().__init__()
        self.raw_files_path : str = raw_files_path
        self.fasta_file : str = fasta_file
        self.output_path : str = output_path
        self.allowed_miscleavages : int = allowed_miscleavages
        self.max_parent_charge :int = max_parent_charge
        self.group_similar_searches: bool = group_similar_searches

    def convert_row(self, row_idx, mzML_rows):
        """
        Creates the config.json for the given row index and in cludes Paths to similat mzMLs
        
        Parameters:
        
        row_idx: index of the row, that is used as reference for the Search Parameters given in mzML

        mzML_rows: index of all Rows containing mzML with equal Search Parameters
        
        """
        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]

        """copying the dictonary to preserve the original"""         
        SAGE_CONFIG = copy.deepcopy(self.SAGE_ORIGINAL_CONFIG)
        
       
        # Set several search Parameters not given by the sdrf
        SAGE_CONFIG["output_directory"] = self.output_path
        SAGE_CONFIG["database"]["fasta"] = self.fasta_file
        SAGE_CONFIG["database"]["enzyme"]["missed_cleavages"] = self.allowed_miscleavages
        SAGE_CONFIG["max_fragment_charge"]= self.max_parent_charge

        # add all mzML Paths with the Search Parametes of row_idx
        for mzML in mzML_rows:
            related_mzML_row: pd.DataFrame = self.sdrf_df.iloc[mzML]
            SAGE_CONFIG["mzml_paths"].append(self.raw_files_path + "/" + related_mzML_row['comment[data file]'])

        # set Peptide Mass tolarance and unit
        precursor_mass_tolerance_split: str = sample['comment[precursor mass tolerance]'].split()
        try:
            precursor_mass_tolerance = float(precursor_mass_tolerance_split[0])
        except:
            raise TypeError("Peptide mass tolerance must be numeric, given was '{tol}".format(tol=precursor_mass_tolerance_split[0]))
        
        if(precursor_mass_tolerance_split[1]=='ppm'):
            SAGE_CONFIG["precursor_tol"]["ppm"] = [-precursor_mass_tolerance, precursor_mass_tolerance]
        elif(precursor_mass_tolerance_split[1]=='Da'):
            SAGE_CONFIG["precursor_tol"]["da"] = [-precursor_mass_tolerance, precursor_mass_tolerance]
        else:
            raise TypeError("Peptide mass tolarance must be given in da or ppm")

        # set fragment mass tolerance and unit
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        try:
            fragment_mass_tolerance = float(fragment_mass_tolerance_split[0])
        except:
            raise TypeError("Fragment mass tolerance must be numeric, given was '{tol}".format(tol=fragment_mass_tolerance_split[0]))
        
        if(fragment_mass_tolerance_split[1]=='ppm'):
            SAGE_CONFIG["fragment_tol"]["ppm"] = [-fragment_mass_tolerance, fragment_mass_tolerance]
        elif(fragment_mass_tolerance_split[1]=='Da'):
            SAGE_CONFIG["fragment_tol"]["da"] = [-fragment_mass_tolerance, fragment_mass_tolerance]
        else:
            raise TypeError("Fragment mass tolarance must be given in da or ppm")
        
        # set cleavage agent details
        cleavage_agent_details_dict: dict = self.ontology_str_to_dict(
            sample['comment[cleavage agent details]']
        )
        if (cleavage_agent_details_dict['NT'] == 'Trypsin'):
            SAGE_CONFIG['database']['enzyme']['cleave_at'] = "KR"
            SAGE_CONFIG['database']['enzyme']['restrict'] = "P"
        else:
            raise SyntaxError("Could not parse cleavage agent from {cad}".format(cad=str(cleavage_agent_details_dict)))


        # Set Modification Parameters if they are listet in mzML
        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):
            modification_dict = self.ontology_str_to_dict(sample[col_name])
            
            try:
                modification_value = float(modification_dict['MM'])
            except:
                raise TypeError("MM of modification must be numeric '{tol}".format(tol=fragment_mass_tolerance_split[0]))

            type = "variable"
            if modification_dict['MT'] and modification_dict['MT'].casefold() == "fixed":
                type = "fixed"
            if type == "variable":
                SAGE_CONFIG['database']['variable_mods'][modification_dict['TA'][2]] = [modification_value]
            if type == "fixed":
                SAGE_CONFIG['database']['static_mods'][modification_dict['TA'][2]] = modification_value

        # write and return the json
        jsonpath = Path(self.output_path + "/Config" + str(row_idx) + ".json")
        jsonpath.write_text(json.dumps(SAGE_CONFIG))
        return jsonpath
   


    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> json:
        """
        Convert the SDRF file to the a json per Search Parameter Set.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path object

        Returns
        -------
        Paths to config.jsons
            The Path to the converted SDRF file
        """

        # init
        self.init_converter(sdrf)
        
        #group the sdrf entrys by their search Parameters
        grouping: pd.core.groupby.DataFrameGroupBy = self.get_grouping()
        
        #create and return config.json for each grouping
        for _, rows in grouping.groups.items():
            grouped_df = self.sdrf_df.iloc[rows]
            self.convert_row(grouped_df.index[0], grouped_df.index)
            print ("To run Sage in Commandline use sage " + self.output_path + "/Config" + str(grouped_df.index[0]) + ".json")
            return Path(self.output_path+"/Config"+str(grouped_df.index[0])+".json")
            




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
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        # Read the X!Tandem params
        converter = SageConverter(
            raw_files_path=cli_args.raw_files_path,
            fasta_file=cli_args.fasta,
            output_path=cli_args.output_path,
            allowed_miscleavages=cli_args.allowed_miscleavages,
            max_parent_charge=cli_args.max_parent_charge, 
            group_similar_searches=cli_args.group_similar_searches
        )
        results = converter.convert(Path(cli_args.sdrf_file))
        
        
   
    @classmethod
    def add_cli_args(cls, subparsers: argparse._SubParsersAction):
        tool_parser = subparsers.add_parser("sage", help="SDRF parser to config jason file for sage")

        tool_parser.add_argument(
            "-r", "--raw-files-path", default="./", type=str, action="store", metavar="/path/to/RAWs", help="Path to the place, where the RAW mzML files are"
        )
        tool_parser.add_argument(
            "-f", "--fasta", default="/path/to/FASTA", type=str, action="store", metavar="/path/to/FASTA", help="Path to the used FASTA file for identification"
        )
        tool_parser.add_argument(
            "-o", "--output-path", default="./SAGE_identification", type=str, action="store", metavar="/path/to/output/prefix", help="Output path for the identification results, including prefix for the files"
        )
        tool_parser.add_argument(
            "-m", "--allowed-miscleavages", default=2, type=int, action="store", metavar="X", help="Maximum missed cleavages for the search (not provided by SDRF)"
        )
        tool_parser.add_argument(
            "-c", "--max-parent-charge", default=4, type=int, action="store", metavar="X", help="Maximum charge state of the parent mass"
        )
        tool_parser.add_argument(
            "-g", "--group-similar_searches", default=True, type=bool, action="store", metavar="X", help="Group equal search settings in one Config.json"
        )
        tool_parser.set_defaults(func=cls.convert_via_cli)
