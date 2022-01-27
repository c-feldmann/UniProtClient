from typing import *
from UniProtClient.errors import ParsingError
import warnings as warn


def simple_name_from(long_name: str) -> str:
    """ Extracts primary name from uniprot string containing all names.

        Additional names are given in brackets or parentheses. These names are removed.

        Parameters
        __________
        long_name: str
            name given by UniProt

        Returns
        _______
        str
            shortened name
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
    if in_bracket != 0 or in_square_bracket != 0:
        warn.warn(f"Error processing: {long_name}\n Returning input name!")
        return long_name
    return "".join(out)


def extract_families(fam_string: str) -> List[Dict[str, Optional[str]]]:
    """Function to parse a string with family-associations of a protein from UniProt.


    Parameters
    ----------
    fam_string: str
        contains the family-associations of a protein given by UniProt. Multiple associations are seperated by ';'
        Different levels (e.g sub- or super-families are seperated by ',')

    Returns
    -------
        List[Dict[str, Optional[str]]]
            Each item in the list corresponds to a family association. In most cases each protein has only one family.
            The items (dictionaries) contain information about subfamily, family, and superfamily.
    """
    if not fam_string:
        return [{"subfamily": None, "family": None, "superfamily": None}]

    # Splitting by "family;" instead of ";" avoids splitting families with ";" in the name. Not perfect but it works.
    fam_group_list = fam_string.split("family;")
    # Reattaching 'family' without ';'. Last element does not end with "family;" therefore no reattachment.
    fam_group_list = [f"{x}family" for x in fam_group_list[:-1]] + [fam_group_list[-1]]
    fam_group_list = [x.strip() for x in fam_group_list] # Removing excess spaces at start and end.

    out_dict_list = []
    for fam_group in fam_group_list:
        family_list = fam_group.split("family,")
        family_list = [f"{x}family" for x in family_list[:-1]] + [family_list[-1]]

        fam_dict = {}
        for category in ["subfamily", "family", "superfamily"]:
            entries = [f for f in family_list if f" {category}" in f]
            if len(entries) == 0:
                fam_dict[category] = None
            elif len(entries) == 1:
                fam_dict[category] = entries[0]
            else:
                raise ParsingError(f"sting contains multiple {category}s: {fam_group}")
        if fam_dict:
            out_dict_list.append(fam_dict)
    if out_dict_list:
        return out_dict_list
    else:
        return [{"subfamily": None, "family": None, "superfamily": None}]
