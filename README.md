# UniProtClient
Python classes in this package allow convenient access to [UniProt](https://www.uniprot.org/) for protein ID mapping and information retrieval.

## Installation in Conda
If not already installed, install **pip** and **git**:  
```
conda install git
conda install pip
```
Then install via pip:
```
pip install git+git://github.com/c-feldmann/UniProtClient
```

## Usage
### Mapping
Protein IDs differ from database to database. The class *UniProtMapper* can be utilized for mapping of protein IDs from one database to corresponding IDs of another database, specified by [letter codes](https://www.uniprot.org/help/api_idmapping).  


```python
from UniProtClient import UniProtMapper
origin_database = 'P_GI'  # PubChem Gene ID
target_database = 'ACC'  # UniProt Accession
gi_2_acc_mappig = UniProtMapper(origin_database, target_database)
```

The obtained object has a function called `map_protein_ids`, which takes a list of strings with protein IDs as input, returning a pandas DataFrame. The DataFrame has two columns: "From" and "To" referring to the origin and target ID, respectively.


```python
gi_numbers = ['224586929', '224586929', '4758208'] # IDs should be represented as a list of strings
# a pandas DataFrame is returned containing the columns "From" and "To"
mapping_df = gi_2_acc_mappig.map_protein_ids(gi_numbers)
uniprot_accessions = mapping_df['To'].tolist()
```

    100%|██████████| 3/3 [00:10<00:00,  3.45s/it]



```python
mapping_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>From</th>
      <th>To</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>224586929</td>
      <td>B4DZW8</td>
    </tr>
    <tr>
      <th>1</th>
      <td>224586929</td>
      <td>Q9Y2R2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4758208</td>
      <td>P51452</td>
    </tr>
  </tbody>
</table>
</div>



### Protein information
UniProt provides a varity of protein specific information, such as protein family, organism, function, EC-number, and many more.
The class *UniProtProteinInfo* is initialized with [column identifier](https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames) specifing the requested information. Spaces in column names should be substituted by underscores.  
If no columns are specified the default is used:

| Column-ID |
|:------:|
| id |
| entry_name |
| protein_names |
| families |
| organism |
| ec |
| genes(PREFERRED) |
| go(molecular_function) |

The column "protein_names" contains all protein names, where secondary names are given in brackets or parenthesis. If this column is requested, the primary name is extracted and added as a new column, called "primary_name".


```python
from UniProtClient import UniProtProteinInfo
info = UniProtProteinInfo()
```


```python
info.load_protein_info(["B4DZW8", "Q9Y2R2", "P51452"])
```

    100%|██████████| 3/3 [00:01<00:00,  1.51it/s]





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>entry_name</th>
      <th>protein_names</th>
      <th>protein_families</th>
      <th>organism</th>
      <th>ec_number</th>
      <th>gene_names(primary)</th>
      <th>gene_ontology(molecular_function)</th>
      <th>primary_name</th>
    </tr>
    <tr>
      <th>entry</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Q9Y2R2</th>
      <td>PTN22_HUMAN</td>
      <td>Tyrosine-protein phosphatase non-receptor type...</td>
      <td>Protein-tyrosine phosphatase family, Non-recep...</td>
      <td>Homo sapiens (Human)</td>
      <td>3.1.3.48</td>
      <td>PTPN22</td>
      <td>kinase binding [GO:0019900]; non-membrane span...</td>
      <td>Tyrosine-protein phosphatase non-receptor type 22</td>
    </tr>
    <tr>
      <th>P51452</th>
      <td>DUS3_HUMAN</td>
      <td>Dual specificity protein phosphatase 3 (EC 3.1...</td>
      <td>Protein-tyrosine phosphatase family, Non-recep...</td>
      <td>Homo sapiens (Human)</td>
      <td>3.1.3.16; 3.1.3.48</td>
      <td>DUSP3</td>
      <td>cytoskeletal protein binding [GO:0008092]; MAP...</td>
      <td>Dual specificity protein phosphatase 3</td>
    </tr>
    <tr>
      <th>B4DZW8</th>
      <td>B4DZW8_HUMAN</td>
      <td>cDNA FLJ55436, highly similar to Tyrosine-prot...</td>
      <td></td>
      <td>Homo sapiens (Human)</td>
      <td></td>
      <td></td>
      <td>protein tyrosine phosphatase activity [GO:0004...</td>
      <td>cDNA FLJ55436, highly similar to Tyrosine-prot...</td>
    </tr>
  </tbody>
</table>
</div>




```python

```
