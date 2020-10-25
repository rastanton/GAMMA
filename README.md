# GAMMA
GAMMA (Gene Annotation and Microbial Mutation Assessment) finds gene matches in microbial genomic data using protein coding (rather than nucleotide) identity, and then translates and annotates the match by providing the type (i.e., mutant, truncation, etc.) and a translated description (i.e., Y190S mutant, truncation at residue 110, etc.).

GAMMA runs in Python/2+ and requires the Biopython package (https://github.com/biopython) as well as Blat (https://genome.ucsc.edu/goldenPath/help/blatSpec.html), which has to be in your $PATH.

GAMMA requires a genome or assembly in fasta format and a multifasta database of the coding sequence of genes. The default output is a ".gamma" file that includes non-overlapping matches to any genes from your database found in your genome. It is run with this command:

> python GAMMA_Exe.py my_genome.fasta gene_db.fasta output_name [optional arguments]

There are five optional arguments:
  -a, --all           Returns all (including overlapping) gene matches
  -e, --extended      Returns all gene mutations, otherwise if there are more than 10 mutations present the count is given
  -f, --fasta         Writes out a multifasta file of the gene matches
  -g, -gff            Generates a general feature format (.gff) file of the output gene matches
  -i, --identity      The minimum nucleotide sequence identiy % used by the Blat search
