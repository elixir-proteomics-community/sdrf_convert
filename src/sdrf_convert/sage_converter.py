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


    

    def __init__(self,
                 raw_files_path: str = "./",
                 fasta_file: str = "/path/to/FASTA",
                 output_path: str = "./Sage_identification",
                 allowed_miscleavages: int = 2, 
                 max_parent_charge: int = 4,
                 group_similar_searches: bool = True,
                 sage_params: str = "/path/to/config",
                 predict_rt: bool = False,
                 generate_decoy: bool = False,
                 decoy_tag: str = "_rev"
                 ):
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
        self.sage_params: str = sage_params
        self.predict_rt: bool = predict_rt
        self.generate_decoy: bool = generate_decoy
        self.decoy_tag: str = decoy_tag

    def convert_row(self, row_idx, mzML_rows):
        """
        Creates the config.json for the given row index and in cludes Paths to similat mzMLs
        
        Parameters:
        
        row_idx: index of the row, that is used as reference for the Search Parameters given in mzML

        mzML_rows: index of all Rows containing mzML with equal Search Parameters
        
        """
        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]

        """copying the dictonary to preserve the original"""
        SAGE_ORIGINAL_CONFIG = json.load(open(Path(self.sage_params)))         
        SAGE_CONFIG = copy.deepcopy(SAGE_ORIGINAL_CONFIG)
        
       
        # Set several search Parameters not given by the sdrf
        try:
            output_path:Path = Path(self.output_path)
            output_path_str:str = str(output_path.absolute())
            SAGE_CONFIG["output_directory"] = output_path_str
        except:
            raise TypeError("Output Path doesnt exist")

        try:
            fasta_path:Path = Path(self.fasta_file)
            fasta_path_str = str(fasta_path.absolute())
            SAGE_CONFIG["database"]["fasta"] = fasta_path_str
        except:
            raise TypeError("path to fasta doesnt exist")

        SAGE_CONFIG["database"]["enzyme"]["missed_cleavages"] = self.allowed_miscleavages
        SAGE_CONFIG["max_fragment_charge"]= self.max_parent_charge
        print(self.predict_rt, self.generate_decoy, self.decoy_tag)
        SAGE_CONFIG["predict_rt"] = self.predict_rt
        SAGE_CONFIG["database"]["generate_decoys"] = self.generate_decoy
        SAGE_CONFIG["database"]["decoy_tag"] = self.decoy_tag

        # add all mzML Paths with the Search Parametes of row_idx
        try:
            raw_files_path:Path = Path(self.raw_files_path)
            raw_files_path_str = str(raw_files_path.absolute())
        except:
            raise TypeError("path to MzMLs doesn't exist")
        for mzML in mzML_rows:
            related_mzML_row: pd.DataFrame = self.sdrf_df.iloc[mzML]
            SAGE_CONFIG["mzml_paths"].append(raw_files_path_str + "/" + related_mzML_row['comment[data file]'].split('.')[0] + ".mzML")

        # set Peptide Mass tolarance and unit
        precursor_mass_tolerance_split: str = sample['comment[precursor mass tolerance]'].split()
        try:
            precursor_mass_tolerance = float(precursor_mass_tolerance_split[0])
        except:
            raise TypeError("Peptide mass tolerance must be numeric, given was '{tol}".format(tol=precursor_mass_tolerance_split[0]))
        
        if(precursor_mass_tolerance_split[1]=='ppm'):
            SAGE_CONFIG["precursor_tol"].clear()
            SAGE_CONFIG["precursor_tol"]["ppm"] = [-precursor_mass_tolerance, precursor_mass_tolerance]
        elif(precursor_mass_tolerance_split[1]=='Da'):
            SAGE_CONFIG["precursor_tol"].clear()
            SAGE_CONFIG["precursor_tol"]["da"] = [-precursor_mass_tolerance, precursor_mass_tolerance]
        else:
            print("used pre cursor tolarance of the original config")

        # set fragment mass tolerance and unit
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        try:
            fragment_mass_tolerance = float(fragment_mass_tolerance_split[0])
        except:
            raise TypeError("Fragment mass tolerance must be numeric, given was '{tol}".format(tol=fragment_mass_tolerance_split[0]))
        
        if(fragment_mass_tolerance_split[1]=='ppm'):
            SAGE_CONFIG["fragment_tol"].clear()
            SAGE_CONFIG["fragment_tol"]["ppm"] = [-fragment_mass_tolerance, fragment_mass_tolerance]
        elif(fragment_mass_tolerance_split[1]=='Da'):
            SAGE_CONFIG["fragment_tol"].clear()
            SAGE_CONFIG["fragment_tol"]["da"] = [-fragment_mass_tolerance, fragment_mass_tolerance]
        else:
            print("Used the Fragment tolarance of Sage Config file")
        
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
        self.__parse_modifications(sample, SAGE_CONFIG)

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


    

    def __parse_modifications(self, sample: pd.DataFrame, SAGE_CONFIG: dict) -> None:
        """
        Parses the modifications of the given sample (row) and adds the
        information
        """

        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):            
            modification_dict = self.ontology_str_to_dict(sample[col_name])
            # Mods falls die Masse direkt in der SDRF gegeben ist
            if 'MM' in modification_dict:
                try:
                    modification_value:float = float(modification_dict['MM'])
                except:
                    raise TypeError("MM of modification must be numeric '{tol}".format(tol=fragment_mass_tolerance_split[0]))
            elif 'NT' in modification_dict:
                mod = self.get_unimod_from_NT(modification_dict['NT'])
                modification_value = mod['mono_mass']
            else:
                modification_value = 0
                print ("modification mass of ", sample[col_name], "could not be found")


            type = "variable"
            if modification_dict['MT'] and modification_dict['MT'].casefold() == "fixed":
                type = "fixed"
            targets_in_row = modification_dict['TA'].split(",")
            for target in targets_in_row:
                if len(target) != 1:
                    target = target.strip("[']")
                if type == "variable":
                    SAGE_CONFIG['database']['variable_mods'][target] = [modification_value]
                if type == "fixed":
                    SAGE_CONFIG['database']['static_mods'][target] = modification_value

    @classmethod
    def convert_via_cli(cls, cli_args: argparse.Namespace):
        # Read the X!Tandem params
        converter = SageConverter(
            raw_files_path=cli_args.raw_files_path,
            fasta_file=cli_args.fasta,
            output_path=cli_args.output_path,
            allowed_miscleavages=cli_args.allowed_miscleavages,
            max_parent_charge=cli_args.max_parent_charge, 
            group_similar_searches=cli_args.group_similar_searches,
            sage_params = cli_args.sage_params,
            predict_rt = cli_args.predict_rt,
            generate_decoy = cli_args.generate_decoy,
            decoy_tag = cli_args.decoy_tag,
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
            "-g", "--group-similar-searches", default=True, type=bool, action="store", metavar="X", help="Group equal search settings in one Config.json"
        )
        tool_parser.add_argument(
            "-p", "--sage-params", default="Sage/config/original", type=str, action="store", metavar="X", help="Base/ Example Config.json that serves als a reference for default settings"
        )
        tool_parser.add_argument(
            "--predict-rt", default=False, type=bool, action="store", metavar="X", help="Bool value determaning if predict rt should be true, default false"
        )
        tool_parser.add_argument(
            "--generate-decoy", default=False, type=bool, action="store", metavar="X", help="Bool value determining if decoys should be generated by sage"
        )
        tool_parser.add_argument(
            "-t", "--decoy-tag", default="rev_", type=str, action="store", metavar="X", help="decoy tag"
        )
        tool_parser.set_defaults(func=cls.convert_via_cli)
