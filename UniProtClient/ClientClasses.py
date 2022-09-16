import pandas as pd
from tqdm.auto import tqdm
import requests
import numpy as np
import sys
from typing import *
from time import sleep
from UniProtClient.parsing import simple_name_from
from UniProtClient.parsing import extract_families


class _UniProtClient:
    def __init__(self, base_url):
        self._base_url = base_url

    @staticmethod
    def _query(query_string) -> str:
        for i in range(10):
            try:
                response = requests.get(query_string)
                return response.text
            except ConnectionResetError:
                sleep(i*10)

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
        --------
        gi2uniprotmapping =  UniProtMapper("P_GI", "ACC") # This class mapps form GI-number to Uniprot IDs

        """
        super().__init__("https://rest.uniprot.org/id-mapping")
        self._from_id = from_id
        self._to_id = to_id
        self._data_format = "tsv"

    def map_protein_ids(self, protein_list: List[str], chunk_size: int = 500) -> pd.DataFrame:
        final_dict_list = []

        pbar = tqdm(total=len(protein_list))
        try:
            for chunk in self._chunkwise(protein_list, chunk_size):
                chunklist = "+".join(chunk)
                server_query = f"?from={self._from_id}&to={self._to_id}&format={self._data_format}&query={chunklist}"
                req = "".join([self._base_url, server_query])
                print(req)
                server_response = self._query(req)
                print(server_response)
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


class UniProtProteinInfo(_UniProtClient):
    """ A class for information retrieval about proteins form UniProt.
    """
    def __init__(self, column_list: Optional[List[str]] = None, merge_multi_fam_associations: Optional[str] = "string",
                 tqdm: bool = True):
        """

        Parameters
        ----------
        column_list: List[str]
            list of strings with column identifiers. [1]
        merge_multi_fam_strings: str, default = "string"
            'sting': when a protein belongs to multiple families: merge names.
            'list': families will be stored as list. Including proteins from only one family.
            None: Protein row is cloned and each row contains one family association.
        tqdm: bool
            show tqdm progress bar

        References
        ----------
        [1] https://www.uniprot.org/help/return_fields
        """
        super().__init__("https://rest.uniprot.org/uniprotkb/")
        if column_list is None:
            column_list = ["accession", "id", "protein_name", "protein_families", "organism_name", "organism_id",
                           "ec", "gene_primary", "go_f"]
        column_list = [self._reformat_column_string(col_id, lower=False) for col_id in column_list]
        if "id" not in column_list:
            column_list.append("id")
        self.columns = ",".join(column_list)
        if merge_multi_fam_associations in ["string", "list", None, False]:
            self.merge_multi_fam_associations = merge_multi_fam_associations
        else:
            raise NotImplementedError(f"Invalid choice for 'merge_multi_fam_strings': "
                                      f"{self.merge_multi_fam_associations}")
        self.tqdm = tqdm

    @staticmethod
    def _reformat_column_string(column_name: str, lower=True) -> str:
        """ A cheap and hacky string formatting procedure replacing spaces with underscores"""
        column_name_reformat = column_name
        while "  " in column_name_reformat:
            column_name_reformat = column_name_reformat.replace("  ", " ")
        column_name_reformat = column_name_reformat.replace(" (", "(")
        column_name_reformat = column_name_reformat.replace(" )", ")")
        column_name_reformat = column_name_reformat.replace(" ", "_")
        if lower:
            column_name_reformat = column_name_reformat.lower()
        return column_name_reformat

    def load_protein_info(self, protein_list: List[str], chunk_size: int = 200, sleeptime=0) -> pd.DataFrame:
        final_dict_list = []
        with tqdm(total=len(protein_list), disable=not self.tqdm) as p_bar:
            for protein_chunk in self._chunkwise(protein_list, chunk_size):
                joined_proteins = "+OR+accession:".join(protein_chunk)
                server_query = f"search?query=accession:{joined_proteins}&format=tsv&fields={self.columns}"
                req = "".join([self._base_url, server_query])
                server_response = self._query(req)
                server_response_formatted = self._response2dictlist(server_response)
                final_dict_list.extend(server_response_formatted)
                p_bar.update(len(protein_chunk))
                sleep(sleeptime)

        valid_entry_df = pd.DataFrame(final_dict_list)
        if "protein_name" in self.columns:
            name_list = []
            for idx, row in valid_entry_df.iterrows():
                try:
                    primary_name = simple_name_from(row["Protein names"])
                except AssertionError as err:
                    primary_name = None
                    print("{}: {}".format(row["Entry"], err), file=sys.stderr)
                name_list.append(primary_name)
            valid_entry_df["primary_name"] = name_list
        if valid_entry_df.shape[0]:
            invalid_ids = set(protein_list) - set(valid_entry_df["Entry"].unique())
        else:
            invalid_ids = set(protein_list)
        invalid_entry_df = pd.DataFrame()
        invalid_entry_df["Entry"] = sorted(invalid_ids)
        if valid_entry_df.shape[0]:
            combined_df = pd.concat([valid_entry_df, invalid_entry_df])
        else:
            combined_df = invalid_entry_df
        column_name_mapping = {old_name: self._reformat_column_string(old_name) for old_name in combined_df.columns}
        combined_df.rename(columns=column_name_mapping, inplace=True)
        combined_df.set_index("entry", inplace=True)

        # Handling protein-families:
        if "protein_families" in combined_df.columns:
            combined_df.loc[combined_df.protein_families.isna(), "protein_families"] = ""
            extracted_family_list = []
            for entry, row in combined_df.iterrows():
                fam_dict_list = extract_families(row["protein_families"])

                if self.merge_multi_fam_associations:
                    fam_dict = {"entry": str(entry)}
                    for category in ["subfamily", "family", "superfamily"]:
                        if self.merge_multi_fam_associations == "string":
                            merged_name = [x[category] if x[category] else "-" for x in fam_dict_list]
                            merged_name = "; ".join(merged_name)
                            if merged_name == "-":
                                merged_name = None
                        elif self.merge_multi_fam_associations == "list":
                            merged_name = [x[category] if x[category] else None for x in fam_dict_list]
                        else:
                            raise NotImplementedError(f"Invalid choice for 'merge_multi_fam_strings':"
                                                      f" {self.merge_multi_fam_associations}")

                        fam_dict[category] = merged_name
                    extracted_family_list.append(fam_dict)
                else:
                    for fam_dict in fam_dict_list:
                        fam_dict["entry"] = str(entry)
                        extracted_family_list.append(fam_dict)
            if extracted_family_list:
                extracted_family_df = pd.DataFrame(extracted_family_list)
                extracted_family_df.set_index("entry", inplace=True)
                combined_df = combined_df.merge(extracted_family_df, left_index=True, right_index=True, how="left")
        return combined_df
