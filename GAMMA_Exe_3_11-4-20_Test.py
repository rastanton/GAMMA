import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import MutableSeq
#from Bio.Alphabet import IUPAC
from Bio import SeqIO
import subprocess
from operator import itemgetter
from decimal import *
getcontext().prec = 4
import math
import argparse

##Written by Richard Stanton (njr5@cdc.gov)
##Requires Python/2+ and blat
##Usage: $ python GAMMA_Exe.py my_scaffolds.fasta gene_db.fasta My_output

def PSL_Type(PSL_Line):
    """Takes in a line from a PSL and returns its type"""
    List1 = PSL_Line.split('\t')
    Match_Length = int(List1[12]) - int(List1[11])
    Region = Genome_Region(PSL_Line)
    if In_Contig(Region, int(List1[10])) == False and Indel_Sum(PSL_Line) != 0:
        Type = 'Indel/Contig Edge'
    elif In_Contig(Region, int(List1[10])) == False and Indel_Sum(PSL_Line) == 0:
        Type = 'Contig Edge'
##    elif List1[0] == List1[14] and Match_Length == int(List1[14]):
##        Type = 'Native'
    elif  Match_Length == int(List1[14]):
        Type = 'Mutant'
    elif Indel_Sum(PSL_Line) != 0 and In_Contig(Region, int(List1[10])) == True:
        Type = 'Indel'
    else:
        Type = 'Mutant'
    return Type

def Indel_Sum(PSL_Line):
    """Returns the sum of the indels (+1 for insertions, -1 for deletions)"""
    Count = 0
    Indels = Indel_Base_Info(PSL_Line)
    for indels in Indels:
        if indels[0] == 'Insertion':
            Count = Count + indels[1]
        elif indels[0] == 'Deletion':
            Count = Count - indels[1]
    return Count

def Genome_Region(PSL_Line):
    """Determines the Genome start and stop positions of the match"""
    List1 = PSL_Line.split('\t')
    if List1[8] == '-':
        Start = int(List1[10]) - int(List1[12])
        End = int(List1[10]) - int(List1[11])
    else:
        Start = int(List1[11])
        End = int(List1[12])
    if List1[15] != '0':
        Start = Start - int(List1[15])
    if List1[16] != List1[14]:
        End = End + (int(List1[14]) - int(List1[16]))
    Sum = Indel_Sum(PSL_Line)
    if Sum < 0:
        End = End + abs(Sum)
    Output = [Start, End]
    return Output

def Indel_Match_Length(PSL_Line):
    """Determines the length of the match to the target gene"""
    List1 = PSL_Line.split('\t')
    Sum = Indel_Sum(PSL_Line)
    Target_Length = Match_Length_Maker(PSL_Line)
    Length = Target_Length + Sum
    return Length

def In_Contig(Start_Stop_List, Contig_Length):
    """Tells if a start and stop list (like from Genome_Region) is within a contig"""
    if Start_Stop_List[0] >= 0 and  Start_Stop_List[1] <= Contig_Length:
        return True
    else:
        return False
                
def Genome_Region_Extractor(PSL_Line, genome):
    """Returns the genome region of the full length match"""
    List1 = PSL_Line.split('\t')
    Positions = Genome_Region(PSL_Line)
    if List1[8] == '-':
        gene = genome[List1[9]]
        gene = gene.reverse_complement()
        gene = gene[Positions[0]:Positions[1]]
    else:
        gene = genome[List1[9]]
        gene = gene[Positions[0]:Positions[1]]
    return gene

def Indel_Typer(PSL_Line, genome_gene, gene):
    """Determines the type of Indel"""
    gene_pro = str(gene.seq.translate())
    genome_pro = str(genome_gene.seq.translate())
    for bases in genome_pro[0:-1]:
        if bases == '*':
            Type = 'Indel Truncation'
            return Type
    if genome_pro[-1] != '*' and gene_pro[-1] == '*':
        Type = 'Indel Nonstop'
        return Type
    else:
        Type = 'Indel'
        return Type

def N_Counter(input_gene):
    """Determines if gene has Ns"""
    Count = 0
    sequence = str(input_gene.seq)
    sequence = sequence.upper()
    for characters in sequence:
        if characters == 'N':
            Count = Count + 1
    return Count

def Truncation_Location(genome_gene):
    genome_pro = str(genome_gene.seq.translate())
    for positions in range(len(genome_pro)):
        if genome_pro[positions] == '*':
            return positions

def Mutant_Typer(PSL_Line, genome_gene, gene):
    """Determines the type of Indel"""
    gene_pro = str(gene.seq.translate())
    genome_pro = str(genome_gene.seq.translate())
    for bases in genome_pro[0:-1]:
        if bases == '*':
            Type = 'Truncation'
            return Type
    if genome_pro[-1] != '*' and gene_pro[-1] == '*':
        Type = 'Nonstop'
        return Type
    elif gene_pro == genome_pro:
        Type = 'Native'
        return Type
    else:
        Type = 'Mutant'
        return Type

def Is_Partial(PSL_Line):
    """If partial match (<90% length) returns True"""
    List1 = PSL_Line.split('\t')
    Match_Length = int(List1[12]) - int(List1[11])
    if Match_Length / float(List1[14]) < 0.9:
        return True
    else:
        return False

def Indel_Base_Info(PSL_Line):
    """Makes a list of Indel lengths and types from a PSL line"""
    Start_Stops_Blocks = Match_Start_Stop_Finder(PSL_Line)
    List1 = PSL_Line.split('\t')
    Lengths = Start_Stops_Blocks[2]
    Genome_Starts = Start_Stops_Blocks[3]
    Gene_Starts = Start_Stops_Blocks[4]
    Output = []
    for entries in range(1, len(Gene_Starts)):
        Gene_Difference = int(Gene_Starts[entries]) - (int(Gene_Starts[entries - 1]) + int(Lengths[entries - 1]))
        Genome_Difference = int(Genome_Starts[entries]) - (int(Genome_Starts[entries - 1]) + int(Lengths[entries - 1]))
        Difference = Gene_Difference - Genome_Difference
        Position = (int(Gene_Starts[entries - 1]) + int(Lengths[entries - 1]))
        if Difference > 0:
            Type = 'Deletion'
            Block = [Type, Difference, Position]
            Output.append(Block)
        elif Difference < 0:
            Type = 'Insertion'
            Difference = Difference * -1
            Block = [Type, Difference, Position]
            Output.append(Block)
    return Output

def Indel_Base_Output(PSL_Line):
    """Makes a readable output from Indel Info"""
    Info_List = Indel_Base_Info(PSL_Line)
    Output = ''
    for entries in Info_List:
        Info = str(entries[1]) + ' bp ' + entries[0] + ' at ' + str(entries[2] + 1) + ','
        Output = Output + Info
    return Output

def Indel_BP_Count(PSL_Line, genome_gene, gene):
    """Makes a count of mutants and indels"""
    genome_gene = str(genome_gene.seq)
    gene = str(gene.seq)
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info(PSL_Line)
    Start = 0
    Count = 0
    Offset = 0
    for entries in Base_Info:
        Block_Difference = Mutant_Count(genome_gene[Start + Offset:entries[2] + Offset], gene[Start:entries[2]])
        Count = Count + Block_Difference
        if entries[0] == 'Deletion':
            Offset = Offset - entries[1]
            Count = Count + entries[1]
            Start = entries[2] + Offset
        elif entries[0] == 'Insertion':
            Offset = Offset + entries[1]
            Count = Count + entries[1]
            Start = entries[2] + Offset
    Block_Difference = Mutant_Count(genome_gene[entries[2] + Offset:], gene[entries[2]:])
    Count = Count + Block_Difference
    return Count

def Indel_Transversion_Count(PSL_Line, genome_gene, gene):
    """Makes a count of transversions"""
    genome_gene = str(genome_gene.seq)
    gene = str(gene.seq)
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info(PSL_Line)
    Start = 0
    Count = 0
    Offset = 0
    for entries in Base_Info:
        Block_Difference = Transversion_Count(genome_gene[Start + Offset:entries[2] + Offset], gene[Start:entries[2]])
        Count = Count + Block_Difference
        if entries[0] == 'Deletion':
            Offset = Offset - entries[1]
            Start = entries[2] + Offset
        elif entries[0] == 'Insertion':
            Offset = Offset + entries[1]
            Start = entries[2] + Offset
    Block_Difference = Transversion_Count(genome_gene[entries[2] + Offset:], gene[entries[2]:])
    Count = Count + Block_Difference
    return Count

def Indel_Codon_Count(PSL_Line, genome_gene, gene):
    """Makes a count of mutants and indels"""
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info(PSL_Line)
    genome_pro  = str(genome_gene.seq.translate())
    gene_pro = str(gene.seq.translate())
    Start = 0
    Count = 0
    Offset = 0
    Codon = 0
    Codon_Difference = 0
    for entries in Base_Info:
        if entries[0] == 'Deletion':
            Offset = Offset - entries[1]
        elif entries[0] == 'Insertion':
            Offset = Offset + entries[1]
        if Offset != 0 and Offset % 3 == 0:
            Codon = math.ceil(entries[2] / float(3)) + Codon
            Codon = int(Codon)
            Block_Difference = Mutant_Count(genome_pro[Start + Codon_Difference:Codon + Codon_Difference], gene_pro[Start:Codon])
            Codon_Difference = Offset // 3
            Count = Count + Block_Difference + abs(Codon_Difference)
            Start = Codon
    Block_Difference = Mutant_Count(genome_pro[Start + Codon_Difference:], gene_pro[Start:])
    Count = Count + Block_Difference
    return Count

def Indel_Mutant_Count(PSL_Line, genome_gene, gene):
    """Makes a count of mutants in an indel"""
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info(PSL_Line)
    genome_pro  = str(genome_gene.seq.translate())
    gene_pro = str(gene.seq.translate())
    Start = 0
    Count = 0
    Offset = 0
    Codon = 0
    Codon_Difference = 0
    for entries in Base_Info:
        if entries[0] == 'Deletion':
            Offset = Offset - entries[1]
        elif entries[0] == 'Insertion':
            Offset = Offset + entries[1]
        if Offset != 0 and Offset % 3 == 0:
            Codon = math.ceil(entries[2] / float(3)) + Codon
            Codon = int(Codon)
            Block_Difference = Mutant_Count(genome_pro[Start + Codon_Difference:Codon + Codon_Difference], gene_pro[Start:Codon])
            Codon_Difference = Offset // 3
            Count = Count + Block_Difference#+ abs(Codon_Difference)
            Start = Codon
    Block_Difference = Mutant_Count(genome_pro[Start + Codon_Difference:], gene_pro[Start:])
    Count = Count + Block_Difference
    return Count

def Indel_Codon_Info(PSL_Line, genome_gene, gene):
    """Returns the mutation and indel info"""
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info(PSL_Line)
    genome_pro  = str(genome_gene.seq.translate())
    gene_pro = str(gene.seq.translate())
    Start = 0
    Offset = 0
    Codon = 0
    Codon_Difference = 0
    Info = ''
    for entries in Base_Info:
        if entries[0] == 'Deletion':
            Offset = Offset - entries[1]
        elif entries[0] == 'Insertion':
            Offset = Offset + entries[1]
        if Offset != 0 and Offset % 3 == 0:
            Codon = math.ceil(entries[2] / float(3)) + Codon
            Codon = int(Codon)
            Mutants = Mutant_Info_Offset(genome_pro[Start + Codon_Difference:Codon + Codon_Difference], gene_pro[Start:Codon], Start)
            Info = Info + Mutants
            Codon_Difference = Offset // 3
            Start = Codon
    Mutants = Mutant_Info_Offset(genome_pro[Start + Codon_Difference:], gene_pro[Start:], Start)
    Info = Info + Mutants
    Indel_Info = Indel_Base_Output(PSL_Line)
    Info = Indel_Info + Info
    return Info

def Frameshifted(PSL_Line):
    """Determines if Indels have caused a frameshift"""
    Count = Indel_Sum(PSL_Line)
    if Count % 3 != 0:
        return True
    else:
        return False

def Transversion_Count(mutant_gene, native_gene):
    """Takes in a mutant gene and a native gene and returns the # of mutations"""
    Count = 0
    mutant_gene = mutant_gene.upper()
    native_gene = native_gene.upper()
    Native_Length = len(native_gene)
    Mutant_Length = len(mutant_gene)
    Length = min([len(native_gene), len(mutant_gene)])
    for characters in range(Length):
        if native_gene[characters] != mutant_gene[characters]:
            if native_gene[characters] == 'G' or native_gene[characters] == 'A':
                if mutant_gene[characters] == 'C' or mutant_gene[characters] == 'T':
                    Count = Count + 1
            elif native_gene[characters] == 'C' or native_gene[characters] == 'T':
                if mutant_gene[characters] == 'G' or mutant_gene[characters] == 'A':
                    Count = Count + 1
    return Count
                                                                  
def Mutant_Count(mutant_gene, native_gene):
    """Takes in a mutant gene and a native gene and returns the # of mutations"""
    Count = 0
    mutant_gene = mutant_gene.upper()
    native_gene = native_gene.upper()
    Native_Length = len(native_gene)
    Mutant_Length = len(mutant_gene)
    Difference = max([len(native_gene), len(mutant_gene)]) - min([len(native_gene), len(mutant_gene)])
    Length = min([len(native_gene), len(mutant_gene)])
    for characters in range(Length):
        if native_gene[characters] != mutant_gene[characters]:
            Count = Count + 1
    Count = Count + Difference
    return Count

def Mutant_Info(mutant_gene, native_gene):
    mutant_gene = mutant_gene.upper()
    native_gene = native_gene.upper()
    Output = ''
    for characters in range(len(native_gene)):
        if native_gene[characters] != mutant_gene[characters]:
            Output = Output + native_gene[characters] + str(characters + 1) + mutant_gene[characters] + ','
    return Output

def Mutant_Info_Offset(mutant_gene, native_gene, offset):
    """Same as Mutant_Info but provides an offset value to match positions"""
    mutant_gene = mutant_gene.upper()
    native_gene = native_gene.upper()
    Output = ''
    for characters in range(len(native_gene)):
        if characters >= len(mutant_gene):
            Output = Output + native_gene[characters] + str(characters + 1 + offset) + ','
        elif native_gene[characters] != mutant_gene[characters]:
            Output = Output + native_gene[characters] + str(characters + 1 + offset) + mutant_gene[characters] + ','
    if len(Output) == 0:
        Output = ''
    return(Output)
    
def Indel_Line(PSL_Line, genome_gene, gene, verbose):
    """Makes a GAMA Line for an Indel"""
    Type = Indel_Typer(PSL_Line, genome_gene, gene)
    List1 = PSL_Line.split('\t')
    Description = Indel_Codon_Info(PSL_Line, genome_gene, gene)
    Codon_Changes = Indel_Codon_Count(PSL_Line, genome_gene, gene)
    Codon_Mutants = Indel_Mutant_Count(PSL_Line, genome_gene, gene)
    Sum = Indel_Sum(PSL_Line)
    Codon_Sum = Sum // 3
    Coding_Length = int(List1[14]) // 3
    if Type == 'Indel Truncation' and verbose == False and Codon_Changes - Codon_Sum > 10:
        Location = Truncation_Location(genome_gene)
        Location = str(Location + 1)
        Info = Indel_Base_Output(PSL_Line)
        Description = Info + 'truncation at codon ' + Location + ' (of ' + str(Coding_Length) + ' codons),' + str(Codon_Mutants) + ' coding mutations'
    elif (Type == 'Indel Truncation' and verbose == True) or (Type == 'Indel Truncation' and Codon_Changes - Codon_Sum >= 10):
        Location = Truncation_Location(genome_gene)
        Location = str(Location + 1)
        Info = Indel_Codon_Info(PSL_Line, genome_gene, gene)
        Description = Info + 'truncation at codon ' + Location + ' (of ' + str(Coding_Length) + ' codons),'
    elif Codon_Changes - Codon_Sum < 10:
        Description = Indel_Codon_Info(PSL_Line, genome_gene, gene)
    elif verbose==True:
        Description = Indel_Codon_Info(PSL_Line, genome_gene, gene)
    else:
        Info = Indel_Base_Output(PSL_Line)
        Description = Info + str(Codon_Mutants) + ' coding mutations'
    BP_Changes = Indel_BP_Count(PSL_Line, genome_gene, gene)
    Transversions = Indel_Transversion_Count(PSL_Line, genome_gene, gene)
    Percent_Codons = str(Decimal(Coding_Length - Codon_Changes) / Decimal(Coding_Length))
    Percent_Bases = str(Decimal(int(List1[14]) - BP_Changes) / Decimal(int(List1[14])))
    Match_Length = Match_Length_Maker(PSL_Line)
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Start = str(Blocks[0][0])
    Stop = str(Blocks[0][1])
    if List1[8] == '-':
        Start = str(int(List1[10]) - Blocks[0][1])
        Stop = str(int(List1[10]) - Blocks[0][0])
    Percent_Length = str(Decimal(Match_Length) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(Codon_Changes) + '\t' + str(BP_Changes) + '\t' + Percent_Codons + '\t' + Percent_Bases + '\t' + Percent_Length + '\t' + str(Match_Length) + '\t' + List1[14] + '\t' + str(Transversions)
    return Out

def Indel_Edge_Line(PSL_Line, genome_gene, gene, verbose):
    """Makes a GAMA Line for an Indel"""
    Type = Indel_Typer(PSL_Line, genome_gene, gene)
    List1 = PSL_Line.split('\t')
    Description = Indel_Codon_Info(PSL_Line, genome_gene, gene)
    Codon_Changes = Indel_Codon_Count(PSL_Line, genome_gene, gene)
    Codon_Missing = Edge_Codon_Missing(PSL_Line)
    Codon_Total = Codon_Changes + Codon_Missing
    Codon_Mutants = Indel_Mutant_Count(PSL_Line, genome_gene, gene)
    BP_Changes = Indel_BP_Count(PSL_Line, genome_gene, gene)
    BP_Missing = Edge_BP_Missing(PSL_Line)
    BP_Total = BP_Changes + BP_Missing
    Sum = Indel_Sum(PSL_Line)
    Codon_Sum = Sum // 3
    Coding_Length = int(List1[14]) // 3
    if Type == 'Indel Truncation' and Codon_Mutants > 10 and verbose != True:
        Location = Truncation_Location(genome_gene)
        Location = str(Location + 1)
        Info = Indel_Base_Output(PSL_Line)
        Description = Info + 'truncation at codon ' + Location + ' (of ' + str(Coding_Length) + ' codons),' + str(Codon_Mutants) + ' coding mutations'
    elif (Type == 'Indel Truncation' and Codon_Mutants < 10) or (Type == 'Indel Truncation' and verbose == True):
        Location = Truncation_Location(genome_gene)
        Location = str(Location + 1)
        Info = Indel_Base_Output(PSL_Line)
        Codon_Info = Indel_Codon_Info(PSL_Line, genome_gene, gene)
        Description = Info + 'truncation at codon ' + Location + ' (of ' + str(Coding_Length) + ' codons),' + Codon_Info + ' mutations'
    elif Codon_Mutants < 10 :
        Description = Indel_Codon_Info(PSL_Line, genome_gene, gene)
    elif verbose==True:
        Description = Indel_Codon_Info(PSL_Line, genome_gene, gene)
    else:
        Info = Indel_Base_Output(PSL_Line)
        Description = Info +  str(Codon_Mutants) + ' coding mutations'
    Type = Type + ' (contig edge)'
    Description = Description + ' for ' + str(int(List1[15]) + 1) + '-' + str(int(List1[16])) + ' of ' + List1[14] + ' bp'
    Transversions = Indel_Transversion_Count(PSL_Line, genome_gene, gene)
    Percent_Codons = str(Decimal(Coding_Length - Codon_Total) / Decimal(Coding_Length))
    Percent_Bases = str(Decimal(int(List1[14]) - BP_Total) / Decimal(int(List1[14])))
    Match_Length = Match_Length_Maker(PSL_Line)
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Start = str(Blocks[0][0])
    Stop = str(Blocks[0][1])
    if List1[8] == '-':
        Start = str(int(List1[10]) - Blocks[0][1])
        Stop = str(int(List1[10]) - Blocks[0][0])
    Percent_Length = str(Decimal(Match_Length) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(Codon_Total) + '\t' + str(BP_Total) + '\t' + Percent_Codons + '\t' + Percent_Bases + '\t' + Percent_Length + '\t' + str(Match_Length) + '\t' + List1[14] + '\t' + str(Transversions)
    return Out

def Mutant_Line(PSL_Line, genome_gene, gene, verbose):
    """Makes a GAMA Line for an Mutant"""
    Type = Mutant_Typer(PSL_Line, genome_gene, gene)
    List1 = PSL_Line.split('\t')
    genome_pro = str(genome_gene.seq.translate())
    gene_pro = str(gene.seq.translate())
    Description = Mutant_Info(genome_pro, gene_pro)
    Codon_Changes = Mutant_Count(genome_pro, gene_pro)
    if Description == '':
        Description = 'No coding mutations'
    elif Type == 'Truncation' and verbose==False:
        Location = Truncation_Location(genome_gene)
        Location = str(Location + 1)
        Description = 'truncation at codon ' + Location + ' (of ' + str(len(gene_pro)) + ' codons),' + str(Codon_Changes) + ' coding mutations'
    elif Type == 'Truncation' and verbose==True:
        Location = Truncation_Location(genome_gene)
        Location = str(Location + 1)
        Description = 'truncation at codon ' + Location + ' (of ' + str(len(gene_pro)) + ' codons),' + Description
    elif Codon_Changes > 10 and verbose==False:
        Description = str(Codon_Changes) + ' coding mutations'
    if Is_Partial(PSL_Line) == True:
        Type = Type + ' (partial match)'
        Description = List1[1] + ' SNPs in ' + str(int(List1[15]) + 1) + '-' + str(int(List1[16])) + ',' + Description 
    Coding_Length = int(List1[14]) // 3
    BP_Changes = Mutant_Count(str(genome_gene.seq), str(gene.seq))
    Transversions = Transversion_Count(str(genome_gene.seq), str(gene.seq))
    Percent_Codons = str(Decimal(Coding_Length - Codon_Changes) / Decimal(Coding_Length))
    Percent_Bases = str(Decimal(int(List1[14]) - BP_Changes) / Decimal(int(List1[14])))
    Match_Length = List1[14]
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Start = str(Blocks[0][0])
    Stop = str(Blocks[0][1])
    if List1[8] == '-':
        Start = str(int(List1[10]) - Blocks[0][1])
        Stop = str(int(List1[10]) - Blocks[0][0])
    Percent_Length = str(Decimal(Match_Length) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(Codon_Changes) + '\t' + str(BP_Changes) + '\t' + Percent_Codons + '\t' + Percent_Bases + '\t' + Percent_Length + '\t' + str(Match_Length) + '\t' + List1[14] + '\t' + str(Transversions)
    return Out

def Edge_Codon_Total(PSL_Line, genome_gene, gene):
    """Counts codon differences from edge matches"""
    Positions = Match_Start_Stop_Finder(PSL_Line)
    Gene_Start = Positions[1][0]
    Offset = int(Gene_Start) % 3
    if Offset != 0:
        Offset = 3 - Offset
    genome_pro = str(genome_gene[Offset:].seq.translate())
    gene_pro = str(gene[Offset:].seq.translate())
    Count = Mutant_Count(genome_pro, gene_pro)
    Missing = Edge_Codon_Missing(PSL_Line)
    Count = Count + Missing
    return Count

def Edge_Codon_Count(PSL_Line, genome_gene, gene):
    """Counts codon differences from edge matches"""
    Positions = Match_Start_Stop_Finder(PSL_Line)
    Gene_Start = Positions[1][0]
    Offset = int(Gene_Start) % 3
    if Offset != 0:
        Offset = 3 - Offset
    genome_pro = str(genome_gene[Offset:].seq.translate())
    gene_pro = str(gene[Offset:].seq.translate())
    Count = Mutant_Count(genome_pro, gene_pro)
##    Missing = Edge_Codon_Missing(PSL_Line)
##    Count = Count + Missing
    return Count

def Edge_Codon_Description(PSL_Line, genome_gene, gene):
    """Describes codon differences from edge matches"""
    Positions = Match_Start_Stop_Finder(PSL_Line)
    Gene_Start = Positions[1][0]
    Offset = int(Gene_Start) % 3
    if Offset != 0:
        Offset = 3 - Offset
    genome_pro = str(genome_gene[Offset:].seq.translate())
    gene_pro = str(gene[Offset:].seq.translate())
    Offset = int(math.ceil(float(Gene_Start) / 3.0))
    Info = Mutant_Info_Offset(genome_pro, gene_pro, Offset)
##    Missing = Edge_Codon_Missing(PSL_Line)
##    Count = Count + Missing
    return Info

def Edge_BP_Count(PSL_Line, genome_gene, gene):
    """Counts bp differences from edge matches"""
    Count = Mutant_Count(genome_gene, gene)
    return Count

def Edge_Transversion_Count(PSL_Line, genome_gene, gene):
    """Counts transversion differences from edge matches"""
    Count = Transversion_Count(genome_gene, gene)
    return Count

def Edge_BP_Total(PSL_Line, genome_gene, gene):
    """Counts bp differences from edge matches"""
    Count = Mutant_Count(genome_gene, gene)
    Missing = Edge_BP_Missing(PSL_Line)
    Count = Count + Missing
    return Count

def Edge_BP_Missing(PSL_Line):
    List1 = PSL_Line.split('\t')
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Total = Blocks[1][1] - Blocks[1][0]
    Missing = int(List1[14]) - Total
    return Missing

def Edge_Codon_Missing(PSL_Line):
    List1 = PSL_Line.split('\t')
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Codons = int(List1[14]) // 3
    Front = math.ceil(Blocks[1][0] / float(3))
    Back = Blocks[1][1] // 3
    Missing = Codons - (Back - int(Front))
    return Missing

def Edge_Line(PSL_Line, genome_gene, gene, verbose):
    """Makes a GAMA line from edge matches w/o indels"""
    Type = "Contig Edge"
    List1 = PSL_Line.split('\t')
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Codon_Changes = Edge_Codon_Total(PSL_Line, genome_gene, gene)
    Codon_Count = Edge_Codon_Count(PSL_Line, genome_gene, gene)
    Codon_Info = Edge_Codon_Description(PSL_Line, genome_gene, gene)
    Coding_Length = int(List1[14]) // 3
    BP_Changes = Edge_BP_Total(PSL_Line, genome_gene, gene)
    BP_Count = Edge_BP_Count(PSL_Line, genome_gene, gene)
    Transversions = Edge_Transversion_Count(PSL_Line, genome_gene, gene)
    Truncation = Indel_Edge_Truncation(Codon_Info)
    Total_Codons = str(len(gene) // 3)
    if Truncation != False:
        Codon_Info = Codon_Info + 'truncation at codon ' + Truncation + ' (of ' + Total_Codons + ' codons)'
        Type = "Contig Edge (truncation)"
    if Codon_Count > 10 and verbose != True and Truncation == False:
        Description = str(BP_Count) + ' SNPs in ' + str(Blocks[1][0] + 1) + '-' + str(Blocks[1][1]) + ',' + str(Codon_Count) + ' coding mutations'
    elif Codon_Count > 10 and verbose != True and Truncation != False:
        Description = str(BP_Count) + ' SNPs in ' + str(Blocks[1][0] + 1) + '-' + str(Blocks[1][1]) + ',' + str(Codon_Count) + ' coding mutations,truncation at codon ' + Truncation + ' (of ' + Total_Codons + ' codons)' 
    else:
        Description = str(BP_Count) + ' SNPs in ' + str(Blocks[1][0] + 1) + '-' + str(Blocks[1][1]) + ',' + Codon_Info
    Percent_Codons = str(Decimal(Coding_Length - Codon_Changes) / Decimal(Coding_Length))
    Percent_Bases = str(Decimal(int(List1[14]) - BP_Changes) / Decimal(int(List1[14])))
    Match_Length = Match_Length_Maker(PSL_Line)
    Start = str(Blocks[0][0])
    Stop = str(Blocks[0][1])
    if List1[8] == '-':
        Start = str(int(List1[10]) - Blocks[0][1])
        Stop = str(int(List1[10]) - Blocks[0][0])
    Percent_Length = str(Decimal(Match_Length) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(Codon_Changes) + '\t' + str(BP_Changes) + '\t' + Percent_Codons + '\t' + Percent_Bases + '\t' + Percent_Length + '\t' + str(Match_Length) + '\t' + List1[14] + '\t' + str(Transversions)
    return Out

def Indel_Edge_Truncation(input_info):
    if ('*,' in input_info) == True:
        Line = input_info.split(',')
        for mutants in Line:
            if mutants[-1] == '*':
                Truncation = mutants[1:-1]
                break
    else:
        Truncation = False
    return Truncation

def Match_Start_Stop_Finder(PSL_Line):
    """Finds the start and stop for the contig and gene"""
    List1 = PSL_Line.split('\t')
    Output = []
    Block_Lengths = List1[18].split(',')[0:-1]
    Blocks = []
    for lengths in Block_Lengths:
        Blocks.append(int(lengths))
    Block_Lengths = Blocks
    Blocks = []
    Gene_Starts = List1[-1].split(',')[0:-1]
    for lengths in Gene_Starts:
        Blocks.append(int(lengths))
    Gene_Starts = Blocks
    Blocks = []
    Genome_Starts = List1[-2].split(',')[0:-1]
    for lengths in Genome_Starts:
        Blocks.append(int(lengths))
    Genome_Starts = Blocks
    if List1[8] == '-':
        Genome_Start = int(List1[10]) - int(List1[12])
        Genome_End = int(List1[10]) - int(List1[11])
    else:
        Genome_Start = int(List1[11])
        Genome_End = int(List1[12])
    Genome_Length = int(List1[10])
    Gene_Start = int(List1[15])
    Gene_End = int(List1[16])
    Gene_Length = int(List1[14])
##    Sum = Indel_Sum(PSL_Line)
##    if Sum < 0:
##        Genome_End = Genome_End + abs(Sum)
    Gene_Starts.append(Gene_End)
    Genome_Starts.append(Genome_End)
    if Gene_Start > 0:
        Delta = Gene_Start
        if Genome_Start - Delta < 0:
            Delta = Genome_Start
            Genome_Start = 0
            Gene_Start = Gene_Start - Delta
            Block_Lengths[0] = Block_Lengths[0] + Delta
            Gene_Starts[0] = Gene_Start
            Genome_Starts[0] = Genome_Start
        else:
            Gene_Start = 0
            Genome_Start = Genome_Start - Delta
            Block_Lengths[0] = Block_Lengths[0] + Delta
            Gene_Starts[0] = Gene_Start
            Genome_Starts[0] = Genome_Start
    if Gene_End < Gene_Length:
        Delta = Gene_Length - Gene_End
        if Genome_End + Delta > Genome_Length:
            Delta = Genome_Length - Genome_End
            Genome_End = Genome_Length
            Gene_End = Gene_End + Delta
            Block_Lengths[-1] = Block_Lengths[-1] + Delta
            Gene_Starts[-1] = Gene_End
            Genome_Starts[-1] = Genome_End
        else:
            Genome_End = Genome_End + Delta
            Gene_End = Gene_Length
            Block_Lengths[-1] = Block_Lengths[-1] + Delta
            Gene_Starts[-1] = Gene_End
            Genome_Starts[-1] = Genome_End         
    Genome_Start_Stop = [Genome_Start, Genome_End]
    Output.append(Genome_Start_Stop)
    Gene_Start_Stop = [Gene_Start, Gene_End]
    Output.append(Gene_Start_Stop)
    Output.append(Block_Lengths)
    Output.append(Genome_Starts)
    Output.append(Gene_Starts)
    return Output

def Match_Length_Maker(PSL_Line):
    """Finds the length of the bp between the start and stop of a gene match on a contig"""
    Blocks = Match_Start_Stop_Finder(PSL_Line)
    Length = int(Blocks[0][1]) - int(Blocks[0][0])
    return Length

def GAMA_Line_Maker(PSL, genome_fasta, genes_fasta, verbose):
    """Makes a list of potential GAMA lines from a PSL file matching Genes to a Genome"""
    f = open(PSL, 'r')
    Lines = []
    String1 = f.readline()
    while String1 != '':
        Lines.append(String1)
        String1 = f.readline()
    f.close()
    Genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))
    Genes = SeqIO.to_dict(SeqIO.parse(genes_fasta, 'fasta'))
    Output = []
    for line in Lines:
        List1 = line.split('\t')
##        print(List1[13])
        Positions = Match_Start_Stop_Finder(line)
        gene = Genes[List1[13]][Positions[1][0]:Positions[1][1]]
        genome = Genome[List1[9]]
        if List1[8] == '-':
            genome = genome.reverse_complement()
        genome = genome[Positions[0][0]:Positions[0][1]]
        Ns = N_Counter(genome)
        if int(List1[0]) / float(List1[14]) >= 0.5:
            Type = PSL_Type(line)
            if Type == 'Mutant':
                Out = Mutant_Line(line, genome, gene, verbose)
            elif Type == 'Indel':
                Out = Indel_Line(line, genome, gene, verbose)
            elif Type == 'Contig Edge':
                Out = Edge_Line(line, genome, gene, verbose)
            elif Type == 'Indel/Contig Edge':
                Out = Indel_Edge_Line(line, genome, gene, verbose)
            Out = Out + '\t' + List1[8]
            if Ns > 0:
                Out_List = Out.split('\t')
                if Out_List[5][-1] == ',':
                    Out_List[5] = Out_List[5] + str(Ns) + ' Ns'
                else:
                    Out_List[5] = Out_List[5] + ',' + str(Ns) + ' Ns'
                Out = '\t'.join(Out_List)
##            print(Out)
            Output.append(Out)
    return Output
                
def Contig_Overlaps(input_list):
    """Takes in a GAMA list and makes a list of lists based on overlaps"""
    Out_List = []
    for lines in input_list:
        List1 = lines.split('\t')
        Add = 1
        for items in Out_List:
            if items[0].split('\t')[1] == List1[1]:
                items.append(lines)
                Add = 0
        if Add == 1:
            New_Contig = [lines]
            Out_List.append(New_Contig)
    return Out_List

def Internal(a,b):
    return (a[0] >= b[0] and a[1] <= b[1]) or (b[0] >= a[0] and b[1] <= a[1])

def Overlap(a,b):
    return (a[0] >= b[0] and a[1] <= b[1]) or (b[0] >= a[0] and b[1] <= a[1]) or (a[0] < b[1] and a[1] > b[0]) or (b[0] < a[1] and b[1] > a[0])

def Overlap_Fraction(a,b):
    if Overlap(a,b) == False:
        Fraction = 0
    Length_a = a[1] - a[0]
    Length_b = b[1] - b[0]
    Small_Distance = min(Length_a, Length_b)
    Distance = min(a[1], b[1]) - max(a[0], b[0])
    Fraction = float(Distance) / float(Small_Distance)
    return Fraction

def Frame_Finder(GAMA_Line):
    """Determines the frame of a GAMA_line"""
    List1 = GAMA_Line.split('\t')
    if List1[-1] == '-':
        Start = int(List1[3])
    else:
        Start = int(List1[2])
    Frame = Start % 3
    return Frame

def Frame_Overlap(GAMA_Line_1, GAMA_Line_2):
    """Determines if genes are coded on same frame"""
    Frame_1 = Frame_Finder(GAMA_Line_1)
    Frame_2 = Frame_Finder(GAMA_Line_1)
    if Frame_1 == Frame_2:
        return True
    else:
        return False

##def Edge_Codon_Changes(Edge_GAMA_line):
##    List1 = Edge_GAMA_line.split('\t')
##    Description = List1[5]
##    Codon_Info =  Description.split(',')[1]
##    Codons = int(Codon_Info.split(' ')[0])
##    return Codons

def Edge_Codon_Changes(Edge_GAMA_line):
    List1 = Edge_GAMA_line.split('\t')
    Description = List1[5]
    if ("coding" in Description):
        Codon_Info =  Description.split(',')[-1]
        Codons = int(Codon_Info.split(' ')[0])
    else:
        Codon_Info = Description.split(',')[1:]
        Codons = len(Codon_Info)
    return Codons

def Edge_BP_Changes(Edge_GAMA_line):
    List1 = Edge_GAMA_line.split('\t')
    Description = List1[5]
    BP_Info =  Description.split(',')[0]
    BP = int(BP_Info.split(' ')[0])
    return BP

def Edge_Line_Comp(Edge_1, Edge_2):
    """Determines if Edge_1 (1) is a better match than Edge_2 (0)"""
    Add = 1
    List1 = Edge_1.split('\t')
    List2 = Edge_2.split('\t')
    if Edge_Codon_Changes(Edge_1) >  Edge_Codon_Changes(Edge_2):
        Add = 0
    elif Edge_Codon_Changes(Edge_1) ==  Edge_Codon_Changes(Edge_2):
        if Edge_BP_Changes(Edge_1) > Edge_BP_Changes(Edge_2):
            Add = 0
        elif Edge_BP_Changes(Edge_1) == Edge_BP_Changes(Edge_2):
            if int(List1[13]) > int(List2[13]):
                Add = 0
            elif int(List1[13]) == int(List2[13]):
                Order_List = [Edge_1, Edge_2]
                Order_List.sort()
                if Order_List[0] != Edge_1:
                    Add = 0
    return Add

def Best_List(input_contig_list):
    """Finds the best matches for a set of matches to the same contig"""
    Output_List = []
    for items in input_contig_list:
        Add = 1
        List1 = items.split('\t')
        for items_2 in input_contig_list:
            List2 = items_2.split('\t')
            if items == items_2:
                continue
            elif Overlap_Fraction([int(List1[2]), int(List1[3])], [int(List2[2]), int(List2[3])]) > 0.5:
##                if ('Truncation' in List1[4]) == True and ('Truncation' in List2[4]) == False:
##                    Add = 0
                if List1[4] == 'Contig Edge' and List2[4] == 'Contig Edge':
                    Add = Edge_Line_Comp(items, items_2)
                    if Add == 0:
                        break
                elif float(List2[8]) > float(List1[8]):
                    Add = 0
                    break
                elif float(List2[8]) == float(List1[8]):
                    if float(List2[9]) > float(List1[9]):
                        Add = 0
                        break
                    elif float(List2[9]) == float(List1[9]):
                        if float(List2[10]) > float(List1[10]):
                            Add = 0
                            break
                        elif float(List2[10]) == float(List1[10]):
                            if int(List2[11]) > int(List1[11]):
                                Add = 0
                                break
                            elif int(List2[11]) == int(List1[11]):
                                if int(List2[13]) < int(List1[13]):
                                    Add = 0
                                    break
                                elif int(List2[13]) == int(List1[13]):
                                    Order_List = [items, items_2]
                                    Order_List.sort()
                                    if Order_List[0] != items:
                                        Add = 0
                                        break
        if Add == 1 and float(List1[10]) > 0.5:
            Out = '\t'.join(List1[0:13])
            Out = Out + '\t' + List1[-1]
            Output_List.append(Out)
    return(Output_List)
                        
def GAMA_List(PSL, genome_fasta, genes_fasta, verbose):
    Lines = GAMA_Line_Maker(PSL, genome_fasta, genes_fasta, verbose)
    Contig_Lines = Contig_Overlaps(Lines)
    Final_List = []
    for contigs in Contig_Lines:
        Best = Best_List(contigs)
        for lines in Best:
            Final_List.append(lines)
    return Final_List

def GAMA_Output(PSL, genome_fasta, genes_fasta, Out_File, verbose):
    List1 = GAMA_List(PSL, genome_fasta, genes_fasta, verbose)
    Output = open(Out_File, 'w')
    Output.write('Gene\tContig\tStart\tStop\tMatch_Type\tDescription\tCodon_Changes\tBP_Changes\tCodon_Percent\tBP_Percent\tPercent_Length\tMatch_Length\tTarget_Length\tStrand\n')
    for lines in List1:
        Output.write(lines + '\n')
    Output.close()

def GAMA_All(PSL, genome_fasta, genes_fasta, Out_File, verbose):
    List1= GAMA_Line_Maker(PSL, genome_fasta, genes_fasta, verbose)
    Output = open(Out_File, 'w')
    Output.write('Gene\tContig\tStart\tStop\tMatch_Type\tDescription\tCodon_Changes\tBP_Changes\tCodon_Percent\tBP_Percent\tPercent_Length\tMatch_Length\tTarget_Length\tStrand\tTransversions\n')
    for lines in List1:
        Output.write(lines + '\n')
    Output.close()

def GAMMA_All_Out(fasta, gene_db, output, identity, verbose):
    subprocess.call('blat' + ' ' + gene_db + ' '  + fasta + ' -noHead -minIdentity=' + str(identity) + ' ' + output + '.psl', shell=True)
    GAMA_All(output +'.psl', fasta, gene_db, output + '.gamma', verbose)

def GAMA_Fasta_Writer(input_fasta, input_GAMA, out_file):
    """Takes in a GAMA and a fasta file and makes a new fasta with just the GAMA matches"""
    Output = open(out_file, 'w')
    f = open(input_GAMA, 'r')
    Genome = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    Name = input_fasta.split('/')[-1]
    String1 = f.readline()
    String1 = f.readline()
    while String1 != '':
        List1 = String1.split('\t')
        gene = Genome[List1[1]][int(List1[2]):int(List1[3])]
        if ('-' in List1[-1]) == True:
            gene = gene.reverse_complement()
        gene.id = List1[0] + '_' + Name
        gene.description = List1[1]
        SeqIO.write(gene, Output, 'fasta')
        String1 = f.readline()
    f.close()
    Output.close()

def Reading_Frame(GAMMA_line):
    Frame = 0
    Info = GAMMA_line.split('\t')[4]
    if ('dge' in Info) == True or ('partial' in Info) == True:
        Info = GAMMA_line.split('\t')[5]
        Info = Info.split(',')[0]
        Info = Info.split(' ')[-1]
        Position = Info.split('-')[0]
        Position = int(Position)
        Position = Position - 1
        Frame = Position % 3
    return Frame

def GFF_Output(input_GAMA, output_file):
    Out = open(output_file, 'w')
    GAMA_list = []
    f = open(input_GAMA, 'r')
    String1 = f.readline()
    String1 = f.readline()
    while String1 != '':
        GAMA_list.append(String1)
        String1 = f.readline()
    f.close()
    for lines in GAMA_list:
        Frame = Reading_Frame(lines)
        List1 = lines.split('\t')
        Attributes = List1[0] + '; ' + List1[5] + '; ' + List1[7] + ' bp changes'
        Out_Line = List1[1] + ' ' + List1[4] + '\tGAMMA\tgene\t' + str(int(List1[2]) + 1) + '\t' + str(int(List1[3]) + 1) + '\t.\t' + List1[-1][0] + '\t' + str(Frame) + '\t' + Attributes + '\n' 
        Out.write(Out_Line)
    Out.close()

def GAMMA_Out(fasta, gene_db, output, identity, verbose):
    subprocess.call('blat' + ' ' + gene_db + ' '  + fasta + ' -noHead -minIdentity=' + str(identity) + ' ' + output + '.psl', shell=True)
    GAMA_Output(output +'.psl', fasta, gene_db, output + '.gamma', verbose)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_fasta", type=str, help="input fasta")
parser.add_argument("-d", "--database", type=str, help="input database")
parser.add_argument("-o", "--output", type=str, help="output name")
parser.add_argument("-a", "--all", action="store_true",
                    help="include all gene matches, even overlaps")
parser.add_argument("-e", "--extended", action="store_true",
                    help="writes out all protein mutations")
parser.add_argument("-f", "--fasta", action="store_true",
                    help="write fasta of gene matches")
parser.add_argument("-g", "--gff", action="store_true",
                    help="write gene matches as gff file")
parser.add_argument("-p", "--percent_identity", default=90, type=int,
                    help="minimum nucleotide identity %, default = 90%")
args = parser.parse_args()
if args.extended:
    Verbose = True
else:
    Verbose = False
if args.all:
    GAMMA_All_Out(args.input_fasta, args.database, args.output,  args.percent_identity, Verbose)
else:
    GAMMA_Out(args.input_fasta, args.database, args.output,  args.percent_identity, Verbose)
if args.fasta:
    GAMA_Fasta_Writer(args.input_fasta, args.output + ".gamma", args.output  + ".fasta")
if args.gff:
    GFF_Output(args.output + ".gamma", args.output + ".gff")
