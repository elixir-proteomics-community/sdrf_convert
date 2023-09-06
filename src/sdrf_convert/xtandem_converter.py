
"""
Module for converting SDRF files to XTandem param files and CLI calls
"""

# std imports
from copy import deepcopy
from io import IOBase
from pathlib import Path
import re
from typing import ClassVar, Dict, Iterator, List, Set, Tuple, Type, Union

import xml.etree.ElementTree as ET
from pyteomics import mass

# 3rd party imports
import pandas as pd

# internal imports
from sdrf_convert.abstract_converter import AbstractConverter



class XTandemConverter(AbstractConverter):
    """
    SDRF convert for XTandem search engine
    """

    COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[precursor mass tolerance]": [pd.StringDtype()],
        "comment[fragment mass tolerance]": [pd.StringDtype()],
        "comment[data file]": [pd.StringDtype()],
        "comment[cleavage agent details]": [pd.StringDtype()],
        "characteristics[organism]": [pd.StringDtype()]
    }

    OPTIONAL_COLUMN_PROPERTIES: ClassVar[Dict[str, List[Type]]] = {
        "comment[modification parameters]*": [pd.StringDtype()],
    }

    CLEAVAGE_REGEX_MAP: ClassVar[Dict[str, str]] = {
        "MS:1001956"          : "[X]|[X]",
        "unspecific cleavage" : "[X]|[X]",
        "MS:1001251"          : "[RK]|{P}",
        "Trypsin"             : "[KR]|{P}",
        "MS:1001313"          : "[KR]",
        "Trypsin/P"           : "[KR]",
        "MS:1001309"          : "[K]|{P}",
        "Lys-C"               : "[K]|{P}",
        "MS:1003093"          : "[K]",
        "Lys-N"               : "[K]",
        "MS:1001303"          : "[R]|{P}",
        "Arg-C"               : "[R]|{P}",
        "MS:1001304"          : "[BD]",
        "Asp-N"               : "[BD]",
        "MS:1001307"          : "[M]",
        "CNBr"                : "[M]",
    }
    """Cleavage enzymes mapped to X!Tandem regex
    """

    def __init__(self, group_similar_searches: bool = False, ):
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
        self.group_similar_searches = group_similar_searches
        self.unimod_db = mass.Unimod()



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

    def __convert_modifications(self, sample: pd.DataFrame) -> Iterator[Tuple[re.Pattern, str]]:
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

    def __convert_row(self, row_idx: int) -> str:
        """
        Creates an XTandem params file the sample at the given row of the SDRF file.
        
        As X!Tandem cannot process RAW files, the actual file path in each XML file 
        will be given as "/path/to/converted/DATA_FILE_NAME" for each data file in the SDRF.

        Parameters
        ----------
        row_idx : int

        Returns
        -------
        str
            representation of the XTandem params file
        """

        sample: pd.DataFrame = self.sdrf_df.iloc[row_idx]
        sample_config: str = ""

        # initiate the XML tree for the params file
        xml_params = ET.Element("bioml")
        
        file_path = "/path/to/converted/" + sample['comment[data file]']
        ET.SubElement(xml_params, "note", type="input", label="spectrum, path").text = file_path

        # TODO: set the output path via parameter
        ET.SubElement(xml_params, "note", type="input", label="output, path").text = "/tmp/xtandem_identification_" + sample['comment[data file]'] + ".xml"
        ET.SubElement(xml_params, "note", type="input", label="output, path hashing").text = "no"

        # TODO: add taxonomy information!!!
        ET.SubElement(xml_params, "note", type="input", label="list path, taxonomy information").text = "sdrf_convert_taxonomy.xml"
        ET.SubElement(xml_params, "note", type="input", label="protein, taxon").text = sample['characteristics[organism]']

        ###
        # set fragment mass tolerance and unit
        fragment_mass_tolerance_split: str = sample['comment[fragment mass tolerance]'].split()
        try:
            fragment_mass_tolerance = float(fragment_mass_tolerance_split[0])
        except:
            raise TypeError("Fragment mass tolerance must be numeric, given was '{tol}".format(tol=fragment_mass_tolerance_split[0]))
        ET.SubElement(xml_params, "note", type="input", label="spectrum, fragment monoisotopic mass error").text = str(fragment_mass_tolerance)

        fragment_mass_tolerance_unit = fragment_mass_tolerance_split[1].strip()
        if fragment_mass_tolerance_unit.casefold() == "ppm":
            fragment_mass_tolerance_unit = "ppm"
        elif fragment_mass_tolerance_unit.casefold() == "da" or fragment_mass_tolerance_unit.casefold() == "daltons":
            fragment_mass_tolerance_unit = "Daltons"
        else:
            raise SyntaxError("Fragment mass tolerance must be ppm or Daltons, but is {unit}".format(unit=fragment_mass_tolerance_unit))
        ET.SubElement(xml_params, "note", type="input", label="spectrum, fragment monoisotopic mass error units").text = fragment_mass_tolerance_unit

        # TODO: check or parameter for monoisotopic or average errors (see also above)
        ET.SubElement(xml_params, "note", type="input", label="spectrum, fragment mass type").text = "monoisotopic"

        ###
        # set peptide mass tolerance and unit
        peptide_mass_tolerance_split: str = sample['comment[precursor mass tolerance]'].split()
        try:
            peptide_mass_tolerance = float(peptide_mass_tolerance_split[0])
        except:
            raise TypeError("Peptide mass tolerance must be numeric, given was '{tol}".format(tol=peptide_mass_tolerance_split[0]))
        
        # TODO: should the parent mass error be plus/minus the value or half the value?
        ET.SubElement(xml_params, "note", type="input", label="spectrum, parent monoisotopic mass error plus").text = str(peptide_mass_tolerance)
        ET.SubElement(xml_params, "note", type="input", label="spectrum, parent monoisotopic mass error plus").text = str(peptide_mass_tolerance)
        
        # only ppm and Daltons (Da) is allowed as unit
        peptide_mass_tolerance_unit = peptide_mass_tolerance_split[1].strip()
        if peptide_mass_tolerance_unit.casefold() == "ppm":
            peptide_mass_tolerance_unit = "ppm"
        elif peptide_mass_tolerance_unit.casefold() == "da" or peptide_mass_tolerance_unit.casefold() == "daltons":
            peptide_mass_tolerance_unit = "Daltons"
        else:
            raise SyntaxError("Peptide mass tolerance must be ppm or Daltons, but is {unit}".format(unit=peptide_mass_tolerance_unit))
        ET.SubElement(xml_params, "note", type="input", label="spectrum, parent monoisotopic mass error units").text = peptide_mass_tolerance_unit
        
        # TODO: check or parameter for monoisotopic or average errors
        ET.SubElement(xml_params, "note", type="input", label="spectrum, parent monoisotopic mass isotope error").text = "yes"
        
        # TODO: set max parent charge with param
        ET.SubElement(xml_params, "note", type="input", label="spectrum, maximum parent charge").text = "4"

        # modifications
        self.__parse_modifications(sample, xml_params)

        # missed cleavages
        # TODO: set missed cleavages by param
        max_missed_cleavages: int = 2
        ET.SubElement(xml_params, "note", type="input", label="scoring, maximum missed cleavage sites").text = str(max_missed_cleavages)

        # cleavage agent details
        cleavage_agent_details_dict: dict = self.ontology_str_to_dict(
            sample['comment[cleavage agent details]']
        )
        cleavage_site: str = None
        # CS param does not work, as X!Tandem uses non-standard reular expressions
        if ('AC' in self.CLEAVAGE_REGEX_MAP) and (cleavage_agent_details_dict['AC'] in self.CLEAVAGE_REGEX_MAP):
            cleavage_site = self.CLEAVAGE_REGEX_MAP[cleavage_agent_details_dict['AC']]
        else:
            if cleavage_agent_details_dict['NT'] in self.CLEAVAGE_REGEX_MAP:
                cleavage_site = self.CLEAVAGE_REGEX_MAP[cleavage_agent_details_dict['NT']]
            else:
                raise SyntaxError("Could not parse cleavage agent from {cad}".format(cad=str(cleavage_agent_details_dict)))
        
        ET.SubElement(xml_params, "note", type="input", label="protein, cleavage site").text = cleavage_site
        ET.SubElement(xml_params, "note", type="input", label="protein, cleavage semi").text = "no"

        # TODO: set threads with param
        ET.SubElement(xml_params, "note", type="input", label="spectrum, threads").text = "4"

        ###
        # for now, we assume we want all results (e.g. for post-processing by Percolator, PIA, OpenMS)
        ET.SubElement(xml_params, "note", type="input", label="output, maximum valid expectation value").text = "0.1"
        ET.SubElement(xml_params, "note", type="input", label="output, parameters").text = "yes"
        ET.SubElement(xml_params, "note", type="input", label="output, proteins").text = "yes"
        ET.SubElement(xml_params, "note", type="input", label="output, results").text = "all"
        ET.SubElement(xml_params, "note", type="input", label="output, sort results by").text = "spectrum"
        ET.SubElement(xml_params, "note", type="input", label="output, spectra").text = "yes"
        ET.SubElement(xml_params, "note", type="input", label="output, xsl path").text = "tandem-style.xsl"

        ET.SubElement(xml_params, "note", type="input", label="protein, use minimal annotations").text = "yes"
        ET.SubElement(xml_params, "note", type="input", label="scoring, pluggable scoring").text = "no"

        # get the XML tree and format with indentations
        tree = ET.ElementTree(xml_params)
        ET.indent(tree)

        return """<?xml version="1.0"?>
<?xml-stylesheet type="text/xsl" href="tandem-input-style.xsl"?>\n""" + ET.tostring(xml_params, encoding='unicode')
    

    def __parse_modifications(self, sample: pd.DataFrame, xml_params: ET.Element) -> None:
        """
        Parses the modifications of the given sample (row) and adds the
        information to teh XML document.
        """

        mods = dict()
        mods["fixed_res"] = list()
        mods["variable_res"] = list()
        mods["quick_acetyl"] = False
        mods["quick_pyrolidone"] = False

        for col_name in self.find_columns(self.sdrf_df, 'comment[modification parameters]*'):
            modification_dict = self.ontology_str_to_dict(sample[col_name])
            mod = self.get_unimod_from_NT(modification_dict['NT'])
            type = "variable"
            if modification_dict['MT'] and modification_dict['MT'].casefold() == "fixed":
                type = "fixed"
            mod_str = self.__create_modification_string(modification_dict['TA'], mod['mono_mass'])

            if type == "variable":
                if (mod['record_id'] == 27) or (mod['record_id'] == 28):
                    mods["quick_pyrolidone"] = True
                elif (mod['record_id'] == 1):
                    mods["quick_acetyl"] = True
                else:
                    mods["variable_res"].append(mod_str)
            else:
                mods["fixed_res"].append(mod_str)
        
        # add fixed modifications
        fixed_mods_str = ""
        for mod in mods["fixed_res"]:
            if len(fixed_mods_str) > 0:
                fixed_mods_str += ","
            fixed_mods_str += mod
        ET.SubElement(xml_params, "note", type="input", label="residue, modification mass").text = fixed_mods_str
        
        # add variable modifications
        variable_mods_str = ""
        for mod in mods["variable_res"]:
            if len(variable_mods_str) > 0:
                variable_mods_str += ","
            variable_mods_str += mod
        ET.SubElement(xml_params, "note", type="input", label="potential modification mass").text = variable_mods_str

        # add quick modifications
        ET.SubElement(xml_params, "note", type="input", label="protein, quick acetyl").text = "yes" if mods["quick_acetyl"] else "no"
        ET.SubElement(xml_params, "note", type="input", label="protein, quick pyrolidone").text = "yes" if mods["quick_pyrolidone"] else "no"
        return mods
    

    def __create_modification_string(self, residues_str: str, mass_shift: float) -> str:
        """
        Create the modifications in a string for the X!Tandem XML
        file, given the residues like in the SDRF file and the
        specific mass shift.

        Parameters
        ----------
        residues_str : str
            comma separated list of residues (amino acids) for teh modification
        
        mass_shift: float
            the mass shift of the modification

        Returns
        -------
        str
            representation of the modification in XTandem style
        """
        residues: List[str] = residues_str.split(",")

        mod_str = ""
        for res in residues:
            if len(mod_str) > 0:
                mod_str += ","
            mod_str += str(mass_shift) + "@" + res

        return mod_str
    

    def __create_taxonomy_stump(self) -> str:
        """
        Creates a stum XML for the taxonomy definitions, which also needs to
        be created for the X!Tandem search (one for all rows together)

        Returns
        -------
        str
            taxonomy XML for all taxonomies in the SDRF file
        """
        # initiate the XML tree for the taxonomy file
        xml_params = ET.Element("bioml", label="x! taxon-to-file matching list")


        for taxonomy in set(self.sdrf_df['characteristics[organism]'].tolist()):
            taxon = ET.SubElement(xml_params, "taxon", label=taxonomy)
            ET.SubElement(taxon, "file", format="peptide", URL="/path/to/FASTA/for/"+taxonomy)

         # get the XML tree and format with indentations
        tree = ET.ElementTree(xml_params)
        ET.indent(tree)

        return """<?xml version="1.0"?>\n""" + ET.tostring(xml_params, encoding='unicode') 


    def convert(self, sdrf: Union[pd.DataFrame, IOBase, Path]) -> Iterator[str]:
        """
        Convert SDRF file to XTandem CLI string and corresponding input file.

        Parameters
        ----------
        sdrf : Union[pd.DataFrame, IOBase, Path]
            SDRF file as pandas DataFrame, file handle or path object
        
        Returns
        -------
        (str)
            The X!Tandem parameter files for each row plus one additional for
            the taxonomy.
        """
        # init
        self.init_converter(sdrf)

        # get the first row for now
        yield self.__convert_row(0)
        # finally, return a taxonomy XML file‚
        yield self.__create_taxonomy_stump()


