# GAMMA
**Introduction**

GAMMA (Gene Allele Mutation Microbial Assessment) is a command line tool that finds gene matches in microbial genomic data using protein coding (rather than nucleotide) identity, and then translates and annotates the match by providing the type (i.e., mutant, truncation, etc.) and a translated description (i.e., Y190S mutant, truncation at residue 110, etc.). Because microbial gene families often have multiple alleles and existing databases are rarely exhaustive, GAMMA is helpful in both identifying and explaining how unique alleles differ from their closest known matches.

**Quick Installation via Conda:**

GAMMA (and all the dependencies) can be installed via Conda:

https://bioconda.github.io/recipes/gamma/README.html

To create a new conda environment (called GAMMA) and install GAMMA into the environment:
```
conda create -n GAMMA gamma -y
```
Activate the Conda environment using the following command:
```
conda activate GAMMA
```
Run the GAMMA.py script to see usage:
```
GAMMA.py -h
```
Deactivate the Conda environment after you are finished (to protect the environment from clashing with other programs):
```
conda deactivate
```
If you download directly via Git, GAMMA requires Python 3+, the Biopython package (https://github.com/biopython), and Blat (https://genome.ucsc.edu/goldenPath/help/blatSpec.html), which has to be in your $PATH.

**Usage:**

The input for GAMMA is a genome or assembly in fasta format and a multifasta database of the coding sequence of genes. The default output is a ".gamma" file that includes non-overlapping matches to any genes from your database found in your genome. It is run with this command:
```
GAMMA.py my_genome.fasta gene_db.fasta output_name [optional arguments]
```
There are five optional arguments:
  -a,--all:            Returns all (including overlapping) gene matches
  
  -e, --extended:      Returns all gene mutations, otherwise if there are more than 10 mutations present the count is given
  
  -f, --fasta:         Writes out a multifasta file of the gene matches
  
  -g, -gff:            Generates a general feature format (.gff) file of the output gene matches
  
  -i, --identity:      The minimum nucleotide sequence identiy % used by the Blat search, input as an integer (i.e., "-i 95" for a 95% threshold), default is 90

**Output:**

The default output of GAMMA is a tab-delimited file with a “.gamma” extension with 14 columns:
1. Gene – The name of the closest matching gene (target) from the database
2. Contig – The name of the contig on which the match was found
3. Start – The start position of the sequence matching the gene on the contig
4. Stop – The end position of the sequence matching the gene on the contig 
5. Match_Type – The type of the gene match based on the translation of the sequence (i.e., the protein sequence). Can be native (for identical amino acid sequences to the target), mutant (for nonsynonymous mutations), truncation (for nonsense mutations), indels (for insertions/deletions), nonstop (for a missing stop codon), contig edge (for matches that are truncated at the start or stop of a contig), or a combination of multiple types (i.e., indel truncation).
6. Description – A short description of the match type.
7. Codon_Changes – The count of the codon changes in the sequence versus the closest match from the datbase.
8. BP_Changes - The count of the basepair changes in the sequence versus the closest match from the datbase.
9. Codon_Percent – The percent (expressed as a decimal value) of the codon similarity between the query and match sequence. Gene matches with large insertions may show a negative value.
10. BP_Percent - The percent (expressed as a decimal value) of the basepair similarity between the query and match sequence. Gene matches with large insertions may show a negative value.
11. Percent_Length - The percent (expressed as a decimal value) of the length of the target covered by the matching sequence, maximum of 1.
12. Match_Length – The length (in basepairs) of the matching sequence.
13. Target_Length - The length (in basepairs) of the target sequence.
14. Strand – The sense of the strand (+ or -) on which the match is found.

Additional outputs in the .gff format and a fasta of the gene matches (in the positive sense) can be generated using the -g and -f options, respectively.

Here’s a sample gamma output generated from running GAMMA on a drug resistant *Klebsiella pneumoniae* (Accession: SAMN11054834) using the ResFinder AR gene database (https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) downloaded on 05-20-2020:

Gene | Contig | Start | Stop | Match_Type | Description | Codon_Changes | BP_Changes | Codon_Percent | BP_Percent | Percent_Length | Match_Length | Target_Length | Strand |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| blaSHV-11_1_X98101 | DHQP1701672_chromosome | 2783310 | 2784171 | Native | No coding mutations | 0 | 8 | 1 | 0.9907 | 1 | 861 | 861 | + |
| oqxA_1_EU370913 | DHQP1701672_chromosome | 1175518 | 1176694 | Native | No coding mutations | 0 | 10 | 1 | 0.9915 | 1 | 1176 | 1176 | - |
| oqxB_1_EU370913 | DHQP1701672_chromosome | 1172342 | 1175495 | Mutant | G148N,G540S,D749E, | 3 | 39 | 0.9971 | 0.9876 | 1 | 3153 | 3153 | - |
| fosA_6_ACZD01000244 | DHQP1701672_chromosome | 4658498 | 4658918 | Native | No coding mutations | 0 | 12 | 1 | 0.9714 | 1 | 420 | 420 | - |
| blaTEM-1A_1_HM749966 | pDHQP1701672_amr_plasmid | 24988 | 25849 | Native | No coding mutations | 0 | 1 | 1 | 0.9988 | 1 | 861 | 861 | - |
| blaOXA-9_1_KQ089875 | pDHQP1701672_amr_plasmid | 26548 | 27373 | Truncation | truncation at codon 112 (of 275 codons),1 coding mutations | 1 | 1 | 0.9964 | 0.9988 | 1 | 825 | 825 | - |
| blaKPC-2_1_AY034847 | pDHQP1701672_amr_plasmid | 37034 | 37916 | Native | No coding mutations | 0 | 0 | 1 | 1 | 1 | 882 | 882 | - |
