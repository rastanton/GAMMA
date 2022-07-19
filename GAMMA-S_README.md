**GAMMA-S:**

GAMMA-S (Gene Allele Mutation Microbial Assessment-Sequence) finds best matches from a gene database without translating them--so it will find the best match by nucleotides, rather by the translated protein sequence. However, it can perform protein-protein sequence matching as well, which requires two protein fastas as the input.

The usage is the same as with GAMMA:
```
GAMMA-S.py my_genome.fasta gene_db.fasta output_name [optional arguments]
```
There are five optional arguments: 
  
  -a,--all:            Returns all (including overlapping) gene matches
  
  -e, --extended:      Returns all gene mutations, otherwise if there are more than 10 mutations present the count is given
 
  -p, --protein:       For protein-protein comparisons, requires two protein sequence fastas as input
  
  -m, --minimum:       The minimum length percent match for output, input as an integer (i.e., "-m 50" for a 50% minimum match length to be reported), not active if the -a/--all option is used default is 20
  
  -i, --identity:      The minimum nucleotide sequence identity % used by the Blat search, input as an integer (i.e., "-i 95" for a 95% threshold), default is 90
  
The output of GAMMA-S is a tab-delimited file with a “.gamma” extension with 17 columns:
1. Gene – The name of the closest matching gene (target) from the database. If there are ambiguous gene matches (i.e., multiple target matches with the same number of basepair changes and transversions), the gene match will be appended with a "‡".
2. Contig – The name of the contig on which the match was found.
3. Start – The start position of the sequence matching the gene on the contig.
4. Stop – The end position of the sequence matching the gene on the contig.
5. Match_Type – The type of the gene match based on the sequence. Can be native (for identical sequences to the target), mutant, indels (for insertions/deletions), contig edge (for matches that are cut off at the start or stop of a contig), or a combination of multiple types (i.e., indel (contig edge)).
6. Description – A short description of the match type.
7. Mismatches – The count of nucleotide/protein sequence substitution mutations.
8. Transversions - The count of basepair changes that are transversions (i.e., purine to pyrimidine or vice versa, such as an A -> C or a T -> G)
9. Insertions – The count of the insertions in the matching sequence.
10. Insertion_BP – The count of the total bases/residues in the insertions in the matching sequence.
11. Deletions – The count of the deletions in the matching sequence.
12. Deletions_BP – The count of the total bases/residues in the insertions in the matching sequence.
13. Unweighted_Match_Percent – The percent match of overlapping sequences only (i.e., does not include sequences on the contig edges, insertions, or deletions).
14. Match_Percent – The percent match of the sequence, subtracing out sequences missing from contig edges, insertions, or deletions. Because insertions are subtracted out, this can lead to cases with negative values, if the insertion is larger than the gene itself.
15. Percent_Length – The percent (expressed as a decimal value) of the length of the target covered by the matching sequence, maximum of 1.
16. Target_Length – The length (in basepairs) of the target sequence.
17. Strand – The sense of the strand (+ or -) on which the match is found.
