# AF2 structures

AlphaFold2-predicted structures for ProteinGym proteins.

## Source

Files originate from the **ProteinGym v1.3** substitution-benchmark structure
bundle (`ProteinGym_AF2_structures.zip`), downloaded from the Marks Lab:

```bash
VERSION="v1.3"
FILENAME="ProteinGym_AF2_structures.zip"
curl -o ${FILENAME} https://marks.hms.harvard.edu/proteingym/ProteinGym_${VERSION}/${FILENAME}
unzip ${FILENAME} && rm ${FILENAME}
```

The unzipped archive is kept in `ProteinGym_AF2_structures/`. Each file is named
`<UniProt_ID>.pdb`, matching the `pdb_file` column of
`ProteinGym/reference_files/DMS_substitutions.csv`.

## Contents of this directory

The `*.pdb` files at the top level of `af2/` are the AF2 structures for the
generalization candidate proteins listed in
`analysis/results/generalization_candidate_proteins.csv`. They were moved here
from `ProteinGym_AF2_structures/` by looking up each candidate's `UniProt_ID`
against the `pdb_file` column of `DMS_substitutions.csv` (24 proteins).
