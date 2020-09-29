import unittest
import pandas as pd
from UniProtClient import UniProtMapper


class TestMapping(unittest.TestCase):
    def test_gi2uniprot_mapping(self):
        """Testing correct mapping from gi number to uniprot ids on 3 exemplary proteins"""
        gi_numbers = ["224586929", "224586929", "4758208"]
        uniprot_ids = ["B4DZW8", "Q9Y2R2", "P51452"]
        mapper = UniProtMapper("P_GI", "ACC")
        mapped_ids = mapper.map_protein_ids(gi_numbers).To.tolist()
        self.assertEqual(uniprot_ids, mapped_ids)


if __name__ == '__main__':
    unittest.main()
