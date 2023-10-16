"""
Module for converting SDRF files to DIA-NN CLI call and config file (for GUI version DIA-NN)
"""
from typing import ClassVar, Dict, List, Type, Union, Any
from io import IOBase
from pathlib import Path

# 3rd party imports
import pandas as pd
from abstract_converter import AbstractConverter


class DiannConverter(AbstractConverter):
    """
    SDRF converter for DIA-NN search engine
    """

    def __init__(self, diann_path: Path, diann_params: str, raw_files_path: Path) -> None:
        super().__init__()
        self.diann_path = diann_path
        self.diann_params = diann_params
        self.raw_files_path = raw_files_path

    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[modification parameters]": [pd.StringDtype()],
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
        "comment[data file]": [pd.StringDtype()],
        "comment[cleavage agent details]": [pd.StringDtype()]
    }

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[number of missed cleavages]": [pd.Int64Dtype()]
    }

    SDRF_COL_NAME_DIANN_PARAM_MAP: ClassVar[Dict[str, Any]] = {
        "comment[number of missed cleavages]": "--missed-cleavages",
        "comment[precursor mass tolerance]": "--mass-acc-ms1",
        "comment[fragment mass tolerance]": "--mass-acc",
        "comment[data file]": "--f",
        "comment[cleavage agent details]": "--cut"

    }

    CLEAVAGE_SITES_MAP: ClassVar[Dict[str, str]] = {
        "Trypsin/P": "K*,R*",
        "Trypsin": "K*,R*,!*P",
        "Lys-C": "K*",
        "Chymotrypsin": "F*,L*,M*,W*,Y*,!*P",
        "Asp-N": "*D",
        "Glu-C": "E*"

    }

    def parse_file_names(self, sdrf_df: pd.DataFrame) -> str:
        """
        Takes a path to the MS file(s) and their names from SDRF file
        and parse name(s) to the DIA-NN's parameter format.

        Returns
        -------
        str
            String with parsed names of MS file(s)
            e.g. "--f folder/subfolder/file_1.raw"
        """
        strinf = ''

        file_names = sdrf_df['comment[data file]'].values
        #apply(lambda x: strinf + '--f ' + path_to_files + '/' + x).to_string()

        for name in file_names:
            strinf = strinf + ' --f ' + self.raw_files_path + '/' + name

        return strinf

    def parse_diann_params(self, sdrf: pd.DataFrame, opt_columns: set) -> str:
        """
        Takes values from OPTIONAL and REQUIRED columns
        in SDRF (except data file) and parse them into
        DIA-NN format according to SDRF_COL_NAME_DIANN_PARAM_MAP.

        Returns
        -------
        str
            String with parsed parameters for DIA-NN
            taken from SDRF REQUIRED and OPTIONAL columns
            e.g. "--var-mod UniMod:1,42.010565,*n --cut K*,R*,!*P --mass-acc-ms1 0"
        """
        # parse precursor and fragment mass tolerance (ppm)
        prec_mass_tol = sdrf['comment[precursor mass tolerance]'].values[0]
        frag_mass_tol = sdrf['comment[fragment mass tolerance]'].values[0]
        mass_tol_str = self.SDRF_COL_NAME_DIANN_PARAM_MAP['comment[precursor mass tolerance]'] + \
            ' ' + str(prec_mass_tol) + ' ' + self.SDRF_COL_NAME_DIANN_PARAM_MAP[
                'comment[fragment mass tolerance]'] + ' ' + str(frag_mass_tol) + ' '

        # parse modifications
        # parser will check hom many columns with modification details are in SDRF and if it is > 1
        # pars ethem separately
        mod_str = ''
        for col_name in sdrf.columns:
            if 'comment[modification parameters]' in col_name:
                mod_details_dict = self.ontology_str_to_dict(sdrf[col_name].values[0])
                n_t = mod_details_dict['NT']
                t_a = mod_details_dict['TA']
                m_t = mod_details_dict['MT']
                unimod_id = self.get_unimod_from_NT(n_t)
                if m_t == 'Fixed':
                    cur_mod_str = '--fixed-mod' + ' ' + 'UniMod:' + str(unimod_id['record_id']) + \
                        ',' + str(unimod_id['mono_mass']) + ',' + t_a + ' '
                elif m_t == 'Variable':
                    cur_mod_str = '--var-mod' + ' ' + 'UniMod:' + str(unimod_id['record_id']) + \
                        ',' + str(unimod_id['mono_mass']) + ',' + t_a + ' '
                mod_str = mod_str + cur_mod_str

        # parse proteases
        protease_details_dict = self.ontology_str_to_dict(
            sdrf['comment[cleavage agent details]'].values[0])
        protease_name = protease_details_dict['NT']
        cleavage_site_regex = self.CLEAVAGE_SITES_MAP[protease_name]
        clvg_ag_str = self.SDRF_COL_NAME_DIANN_PARAM_MAP[
            'comment[cleavage agent details]'] + ' ' + cleavage_site_regex + ' '

        # parse missed cleavage
        if len(opt_columns) > 0:
            missed_cleavages = sdrf['comment[number of missed cleavages]'].values[0]
            missed_cleavages_str = self.SDRF_COL_NAME_DIANN_PARAM_MAP[
                'comment[number of missed cleavages]'] + ' ' + str(missed_cleavages)

        else:
            missed_cleavages_str = ''

        return mass_tol_str + mod_str + clvg_ag_str + missed_cleavages_str

    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> tuple:
        """
        Convert the SDRF file to the desired format.
        The output depends on the targeted software, it may be a config file only, 
        a CLI string or a Tuple containing both.

        Returns
        -------
        tuple
            CLI command (for terminal execution of DIA-NN) and config file (any extension, e.g .cfg)
            containing command line paramenters that can be passed in GUI version
            of DIA-NN via --cfg parameter in Additional options window.
        """
        self.init_converter(sdrf)
        cli = str(self.diann_path) + ' ' \
            + self.parse_file_names(self.sdrf_df) + ' '\
            + self.diann_params + ' ' +\
            self.parse_diann_params(self.sdrf_df, self.present_optional_columns)

        with open(self.raw_files_path + '/config_diann.cfg', 'wt', encoding='UTF-8') as cfg:
            cfg.write(cli)

        return cli
    