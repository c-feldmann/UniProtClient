import pandas as pd
from tqdm import tqdm
import urllib.parse
import urllib.request
import numpy as np
from typing import *


class UniProtMapper:
    def __init__(self, from_id: str, to_id: str):
        """For mapping of protein IDs to another ID. Uses UniProt API.
        for valid parameters see: https://www.uniprot.org/help/api_idmapping

        :param from_id:
        :param to_id:
        """
        self._from_id = from_id
        self._to_id = to_id
        self._data_format = "tab"

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

    def map_protein_ids(self, protein_list: List[str], chunk_size: int = 500):
        final_dict_list = []

        pbar = tqdm(total=len(protein_list))
        for i in range(0, len(protein_list), chunk_size):
            chunk = protein_list[i:i + chunk_size]
            chunklist = "+".join(chunk)
            base_url = "https://www.uniprot.org/uploadlists/"
            server_query = f"?from={self._from_id}&to={self._to_id}&format={self._data_format}&query={chunklist}"
            req = "".join([base_url, server_query])
            server_response = self._query(req)
            server_response_formatted = self._response2dictlist(server_response)
            final_dict_list.extend(server_response_formatted)
            pbar.update(len(chunk))
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
        elif letter == ")":
            in_bracket -= 1
        elif letter == "[":
            in_square_bracket += 1
        elif letter == "]":
            in_square_bracket -= 1
        else:
            if in_bracket == 0 and in_square_bracket == 0:
                if buffer and letter != " ":
                    out.extend(buffer)
                    buffer = []
                out.append(letter)
    assert in_bracket == 0
    assert in_square_bracket == 0
    return "".join(out).rstrip(" ")


class UniProtProteinInfo:
    def __init__(self, column_list: Optional[List[str]] = None):
        if column_list is None:
            column_list = ["id", "entry_name", "protein_names", "families", "organism", "ec", "genes(PREFERRED)",
                           "go(molecular_function)"]
        self.columns = ",".join(column_list)

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
            # Replacing "" with np.nan
            line_items = [i if i != "" else np.nan for i in line_items]
            assert len(header_items) == len(line_items), (header_items, line_items)
            r_dict_list.append(dict(zip(header_items, line_items)))
        return r_dict_list

    def load_protein_info(self, protein_list: List[str]):
        final_dict_list = []
        pbar = tqdm(total=len(protein_list))
        try:
            for protein in protein_list:
                base_url = "https://www.uniprot.org/uniprot/"
                server_query = f"?query=accession:{protein}&format=tab&columns={self.columns}"
                req = "".join([base_url, server_query])
                server_response = self._query(req)
                server_response_formatted = self._response2dictlist(server_response)
                final_dict_list.extend(server_response_formatted)
                pbar.update(1)
        finally:
            pbar.close()
        valid_mappings = pd.DataFrame(final_dict_list)
        if "protein_names" in self.columns:
            valid_mappings["prefered_name"] = valid_mappings["Protein names"].apply(simle_name_from)
        invalid_ids = set(protein_list) - set(valid_mappings["Entry"].unique())
        invalid_mapping = pd.DataFrame()
        invalid_mapping["Entry"] = sorted(invalid_ids)

        return pd.concat([valid_mappings, invalid_mapping])

