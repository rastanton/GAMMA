# GAMMA
GAMMA (Gene Allele and Mutation Assessment) finds gene matches in microbial genomic data using protein coding (rather than nucleotide) identity, and then translates and annotates the match by providing the type (i.e., mutant, truncation, etc.) and a translated description (i.e., Y190S mutant, truncation at residue 110, etc.). Because microbial gene families often have multiple alleles and existing databases are rarely exhaustive, GAMMA is helpful in both identifying and explaining how unique alleles differ from their closest known matches.

GAMMA runs in Python 3+ and requires the Biopython package (https://github.com/biopython) as well as Blat (https://genome.ucsc.edu/goldenPath/help/blatSpec.html), which has to be in your $PATH.

The input for GAMMA is a genome or assembly in fasta format and a multifasta database of the coding sequence of genes. The default output is a ".gamma" file that includes non-overlapping matches to any genes from your database found in your genome. It is run with this command:

> GAMMA.py my_genome.fasta gene_db.fasta output_name [optional arguments]

There are five optional arguments:

  -a,--all:            Returns all (including overlapping) gene matches
  
  -e, --extended:      Returns all gene mutations, otherwise if there are more than 10 mutations present the count is given
  
  -f, --fasta:         Writes out a multifasta file of the gene matches
  
  -g, -gff:            Generates a general feature format (.gff) file of the output gene matches
  
  -i, --identity:      The minimum nucleotide sequence identiy % used by the Blat search, input as an integer (i.e., "-i 95" for a 95% threshold), default is 90
  
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

Here’s a sample gamma output generated from running GAMMA on a PacBio assembly (AR-0361, https://www.ncbi.nlm.nih.gov/assembly/GCF_002968495.1) using the ResFinder AR gene database (https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) downloaded on 05-20-2020:

| Gene | Contig | Start | Stop | Match_Type | Description | Codon_Changes | BP_Changes | Codon_Percent | BP_Percent | Percent_Length | Match_Length | Target_Length | Strand |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| fosA_6_ACZD01000244 | AR-0361_chromosome | 4930620 | 4931040 | Mutant | I91V,D138E, | 2 | 18 | 0.9857 | 0.9571 | 1 | 420 | 420 | + |
| blaSHV-11_1_X98101 | AR-0361_chromosome | 1548139 | 1549000 | Native | No coding mutations | 0 | 2 | 1 | 0.9977 | 1 | 861 | 861 | - |
| blaTEM-1A_1_HM749966 | AR-0361_plasmid_2 | 27269 | 28130 | Native | No coding mutations | 0 | 1 | 1 | 0.9988 | 1 | 861 | 861 | - |
| blaOXA-9_1_KQ089875 | AR-0361_plasmid_2 | 28829 | 29654 | Truncation | truncation at codon 112 (of 275 codons),1 coding mutations | 1 | 1 | 0.9964 | 0.9988 | 1 | 825 | 825 | - |
| blaKPC-2_1_AY034847 | AR-0361_plasmid_2 | 39315 | 40197 | Native | No coding mutations | 0 | 0 | 1 | 1 | 1 | 882 | 882 | - |
| aph(4)-Ia_1_V01499 | AR-0361_plasmid_6 | 37042 | 38068 | Native | No coding mutations | 0 | 0 | 1 | 1 | 1 | 1026 | 1026 | + |
| aac(3)-IV_1_DQ241380 | AR-0361_plasmid_6 | 36037 | 36814 | Mutant | W5L, | 1 | 1 | 0.9961 | 0.9987 | 1 | 777 | 777 | + |
| cmlA1_1_M64556 | AR-0361_plasmid_6 | 28540 | 29800 | Mutant | G412E, | 1 | 1 | 0.9976 | 0.9992 | 1 | 1260 | 1260 | + |
| aadA2b_1_D43625 | AR-0361_plasmid_6 | 27499 | 28279 | Mutant | G202D, | 1 | 1 | 0.9962 | 0.9987 | 1 | 780 | 780 | + |
| blaSHV-12_1_KF976405 | AR-0361_plasmid_7 | 54 | 915 | Native | No coding mutations | 0 | 0 | 1 | 1 | 1 | 861 | 861 | + |
