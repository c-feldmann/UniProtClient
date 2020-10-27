import pandas as pd
from tqdm import tqdm
import urllib.parse
import urllib.request
import numpy as np
from typing import *
from time import sleep


class _UniProtClient:
    def __init__(self, base_url):
        self._base_url = base_url

    @staticmethod
    def _query(query_string) -> str:
        with urllib.request.urlopen(query_string) as f:
            response = f.read()
        return response.decode('utf-8')

    @staticmethod
    def _response2dictlist(response_string) -> List[dict]:
        header_row = response_string.split("\n")[0]
        header_items = header_row.split("\t")
        r_dict_list = []
        for line in response_string.split("\n")[1:]:
            if not line:
                continue
            line_items = line.split("\t")
            assert len(header_items) == len(line_items), (header_items, line_items)
            r_dict_list.append(dict(zip(header_items, line_items)))
        return r_dict_list

    @staticmethod
    def _chunkwise(iterables, chunk_size):
        for i in range(0, len(iterables), chunk_size):
            chunk = iterables[i:i + chunk_size]
            yield chunk


class UniProtMapper(_UniProtClient):
    def __init__(self, from_id: str, to_id: str):
        """For mapping of protein IDs to another ID. Uses UniProt API.
        for valid parameters see: https://www.uniprot.org/help/api_idmapping

        :param to_id:

        Parameters
        ----------
        from_id: origin ID string
        to_id: target ID string

        Examples
        ________
        gi2uniprotmapping =  UniProtMapper("P_GI", "ACC") # This class mapps form GI-number to Uniprot IDs

        """
        super().__init__("https://www.uniprot.org/uploadlists/")
        self._from_id = from_id
        self._to_id = to_id
        self._data_format = "tab"

    def map_protein_ids(self, protein_list: List[str], chunk_size: int = 500) -> pd.DataFrame:
        final_dict_list = []

        pbar = tqdm(total=len(protein_list))
        try:
            for chunk in self._chunkwise(protein_list, chunk_size):
                chunklist = "+".join(chunk)
                server_query = f"?from={self._from_id}&to={self._to_id}&format={self._data_format}&query={chunklist}"
                req = "".join([self._base_url, server_query])
                server_response = self._query(req)
                server_response_formatted = self._response2dictlist(server_response)
                final_dict_list.extend(server_response_formatted)
                pbar.update(len(chunk))
        finally:
            pbar.close()
        valid_mappings = pd.DataFrame(final_dict_list)
        invalid_ids = set(protein_list) - set(valid_mappings["From"].unique())
        invalid_mapping = pd.DataFrame()
        invalid_mapping["From"] = sorted(invalid_ids)
        invalid_mapping["To"] = np.nan
        return pd.concat([valid_mappings, invalid_mapping])


def simle_name_from(long_name):
    """ Extracts primary name from uniprot string containing all names.

    Additional names are given in brackets and parantheses
    """
    out = []
    buffer = []
    in_bracket = 0
    in_square_bracket = 0
    for letter in long_name:
        if letter == "(":
            in_bracket += 1
            buffer.append(letter)
        elif letter == ")":
            in_bracket -= 1
            buffer.append(letter)
        elif letter == "[":
            in_square_bracket += 1
            buffer.append(letter)
        elif letter == "]":
            in_square_bracket -= 1
            buffer.append(letter)
        else:
            # If not in bracket
            if in_bracket == 0 and in_square_bracket == 0:
                if letter == " ":
                    buffer.append(letter)
                elif buffer:
                    out.extend(buffer)
                    buffer = []
                    out.append(letter)
                else:
                    out.append(letter)
            else:
                buffer.append(letter)
    assert in_bracket == 0
    assert in_square_bracket == 0
    return "".join(out)


class UniProtProteinInfo(_UniProtClient):
    """
    TODO: Find out why UniProt returns duplicate rows.
    """
    def __init__(self, column_list: Optional[List[str]] = None):
        """

        Parameters
        ----------
        column_list: strings of column identifiers. [1]


        ..References
        ..__________
        [1] https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
        """
        super().__init__("https://www.uniprot.org/uniprot/")
        if column_list is None:
            column_list = ["id", "entry_name", "protein_names", "families", "organism", "ec", "genes(PREFERRED)",
                           "go(molecular_function)"]
        self.columns = ",".join(column_list)

    def load_protein_info(self, protein_list: List[str], chunk_size: int = 200, sleeptime=0):
        final_dict_list = []
        with tqdm(total=len(protein_list)) as p_bar:
            for protein_chunk in self._chunkwise(protein_list, chunk_size):
                joined_proteins = "+OR+".join(protein_chunk)
                server_query = f"?query=accession:{joined_proteins}&format=tab&columns={self.columns}"
                req = "".join([self._base_url, server_query])
                server_response = self._query(req)
                server_response_formatted = self._response2dictlist(server_response)
                final_dict_list.extend(server_response_formatted)
                p_bar.update(len(protein_chunk))
                sleep(sleeptime)

        valid_mappings = pd.DataFrame(final_dict_list)
        if "protein_names" in self.columns:
            valid_mappings["prefered_name"] = valid_mappings["Protein names"].apply(simle_name_from)
        invalid_ids = set(protein_list) - set(valid_mappings["Entry"].unique())
        invalid_mapping = pd.DataFrame()
        invalid_mapping["Entry"] = sorted(invalid_ids)

        # For some reason the data returned from UniProt contains duplicates. Potential Error in query?
        return pd.concat([valid_mappings, invalid_mapping]).drop_duplicates()
