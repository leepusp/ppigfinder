# Guided input loading

The guided workflow supports multi-file input loading.

Use **Add input data** to select one or more files at the same time. ppigFinder detects each file by extension and routes it into the workflow:

- `.fa`, `.fasta`, `.fna`, `.ffn`, `.gb`, `.gbk`, `.genbank`, `.dna` → genome input;
- `.faa`, `.pep`, `.protein`, `.prot` → protein query input for BLAST;
- `.hmm` → HMM/domain profile input;
- `.ppigfinder.json` → ppigFinder project;
- `.json` → project/snapshot style input;
- AlphaFold 3 results are imported as folders.

The flow chips at the top of each module page are clickable. They can be used as shortcuts to move between Data, Genome, ORFs, Annotation, AlphaFold, DaVinci/HPC and Reports.
