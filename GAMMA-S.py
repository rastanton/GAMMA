#!/usr/bin/env python3

import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
from operator import itemgetter
from decimal import *
getcontext().prec = 4
import math
import argparse

##Written by Richard Stanton (njr5@cdc.gov)
##Requires Python/3+ and blat
##Usage: $ GAMMA-S.py my_scaffolds.fasta gene_db.fasta My_output

def PSL_Type(PSL_Line):
    """Takes in a line from a PSL and returns its type"""
    List1 = PSL_Line.split('\t')
    Match_Length = int(List1[12]) - int(List1[11])
    Region = Genome_Region(PSL_Line)
    if In_Contig(PSL_Line) == False and Indel_Truth(PSL_Line) == True:
        Type = 'Indel/Contig Edge'
    elif In_Contig(PSL_Line) == False and Indel_Truth(PSL_Line) == False:
        Type = 'Contig Edge'
##    elif List1[0] == List1[14] and Match_Length == int(List1[14]):
##        Type = 'Native'
    elif int(List1[16]) - int(List1[15]) < int(List1[14]) and Indel_Truth(PSL_Line) == True and In_Contig(PSL_Line) == True:
        Type = 'Indel/Partial'
    elif int(List1[16]) - int(List1[15]) < int(List1[14]):
        Type = 'Partial'
    elif Indel_Truth(PSL_Line) == True and In_Contig(PSL_Line) == True:
        Type = 'Indel'
    elif  Match_Length == int(List1[14]):
        Type = 'Match'
    else:
        Type = 'Match'
    return Type

def In_Contig(PSL_Line):
    """Tells if a start and stop list (like from Genome_Region) is within a contig"""
    List1 = PSL_Line.split('\t')
    if int(List1[15]) > int(List1[11]) or (int(List1[16]) < int(List1[14]) and int(List1[10]) == int(List1[12])):
        return False
    else:
        return True

def Genome_Region(PSL_Line):
    List1 = PSL_Line.split('\t')
    Out = [int(List1[11]), int(List1[12])]
    return Out

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

def Indel_Truth(PSL_Line):
    List1 = PSL_Line.split('\t')
    if (List1[4] != '0' or List1[6] != '0') and List1[5] != List1[7]:
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

def Indel_Base_Info_Offset(PSL_Line):
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
    if List1[15] != 0:
        for entry in Output:
            entry[2] = entry[2] - int(List1[15])
    return Output

def Indel_String(indel_list):
    Out = ''
    for entry in indel_list:
        Out = Out + str(entry[1]) + ' bp ' + str(entry[0]) + ' at ' + str(entry[2]) + ','
    return Out[0:-1]

def Indel_BP_Count(PSL_Line, genome_gene, gene):
    """Makes a count of mutants and indels"""
    genome_gene = str(genome_gene.seq)
    gene = str(gene.seq)
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info_Offset(PSL_Line)
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

def Transversion_Count_Mutant_List(mutant_list):
    Count = 0
    if mutant_list != [''] and len(mutant_list) > 0:
        for entry in mutant_list:
            if entry[0] == 'G' or entry[0] == 'A':
                if entry[-1] == 'C' or entry[-1] == 'T':
                    Count = Count + 1
            elif entry[0] == 'C' or entry[0] == 'T':
                if entry[-1] == 'G' or entry[-1] == 'A':
                    Count = Count + 1
    return Count

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

def Mutant_Info(mutant_gene, native_gene):
    mutant_gene = mutant_gene.upper()
    native_gene = native_gene.upper()
    Output = ''
    for characters in range(len(native_gene)):
        if native_gene[characters] != mutant_gene[characters]:
            Output = Output + native_gene[characters] + str(characters + 1) + mutant_gene[characters] + ','
    return Output[0:-1]

def Indel_BP_List(PSL_Line, genome_gene, gene):
    """Makes a list of mutants from indels"""
    String1 = Indel_BP_String(PSL_Line, genome_gene, gene)
    if String1 == '':
        Out = []
    else:
        Out = String1.split(',')
    return Out

def Indel_Totals(indel_list):
    Insertions = 0
    Insertion_Bases = 0
    Deletions = 0
    Deletion_Bases = 0
    for entry in indel_list:
        if entry[0] == 'Insertion':
            Insertions += 1
            Insertion_Bases += int(entry[1])
        else:
            Deletions += 1
            Deletion_Bases += int(entry[1])
    Out = str(Insertions) + '\t' + str(Insertion_Bases) + '\t' + str(Deletions) + '\t' + str(Deletion_Bases)
    return Out

def Indel_BP_String(PSL_Line, genome_gene, gene):
    """Makes a list of mutants from indels"""
    genome_gene = str(genome_gene.seq)
    gene = str(gene.seq)
    List1 = PSL_Line.split('\t')
    Length = len(gene)
    Base_Info = Indel_Base_Info_Offset(PSL_Line)
    Positions = Match_Start_Stop_Finder(PSL_Line)
    Start_Offset = Positions[1][0]
    Start = 0
    Mutant_All = ''
    Offset = 0
    for entries in Base_Info:
        Mutants = Mutant_Info_Offset(genome_gene[Start + Offset:entries[2] + Offset], gene[Start:entries[2]], Start_Offset)
        Mutant_All = Mutant_All + Mutants
        if entries[0] == 'Deletion':
            Offset = Offset - entries[1]
            Start = entries[2]
            Start_Offset = Start + Offset
        elif entries[0] == 'Insertion':
            Offset = Offset + entries[1]
            Start = entries[2]
            Start_Offset = Start + Offset
    Mutants = Mutant_Info_Offset(genome_gene[entries[2] + Offset:], gene[entries[2]:], Start_Offset)
    Mutant_All = Mutant_All + Mutants
    if len(Mutant_All) > 0 and Mutant_All[-1] == ',':
        Mutant_All = Mutant_All[0:-1]
    return Mutant_All

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
    Gene_Starts.append(Gene_End)
    Genome_Starts.append(Genome_End)       
    Genome_Start_Stop = [Genome_Start, Genome_End]
    Output.append(Genome_Start_Stop)
    Gene_Start_Stop = [Gene_Start, Gene_End]
    Output.append(Gene_Start_Stop)
    Output.append(Block_Lengths)
    Output.append(Genome_Starts)
    Output.append(Gene_Starts)
    return Output

    
def Match_Line(PSL_line, genome_gene, gene, verbose):
    List1 = PSL_line.split('\t')
    BP_Changes = Mutant_Count(genome_gene, gene)
    Type = 'Mutant'
    Description = Mutant_Info(genome_gene, gene)
    if BP_Changes == 0:
        Type = 'Native'
        Description = 'Exact match'
    elif verbose == True:
        Description = Mutant_Info(genome_gene, gene)
    elif BP_Changes > 10 and verbose == False:
        Description = str(BP_Changes) + ' mutations'
    Start = List1[11]
    Stop = List1[12]
    Transversions = Transversion_Count(genome_gene, gene)
    Percent_Bases = str(Decimal(int(List1[14]) - BP_Changes) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(BP_Changes) + '\t' + str(Transversions) + '\t0\t0\t0\t0\t' + Percent_Bases + '\t' + Percent_Bases + '\t1\t' + List1[14] + '\t' + List1[8]
    return Out

def Partial_Line(PSL_line, genome_gene, gene, verbose):
    List1 = PSL_line.split('\t')
    BP_Changes = Mutant_Count(genome_gene, gene)
    if BP_Changes == 0:
        Type = 'Partial'
        Description = 'No mutations for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
    elif BP_Changes <= 10 or verbose == True:
        Description = Mutant_Info_Offset(genome_gene, gene, int(List1[15]))
        Description = Description + ' for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Partial (Mutant)'
    elif BP_Changes > 10 and verbose == False:
        Description = str(BP_Changes) + ' mutations for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Partial (Mutant)'
    Start = List1[11]
    Stop = List1[12]
    Length = int(List1[12]) - int(List1[11])
    Transversions = Transversion_Count(genome_gene, gene)
    Percent_Bases = str(Decimal(int(List1[0])) / Decimal(int(List1[0]) + int(List1[1])))
    Percent_Match = str(Decimal(Length - BP_Changes) / Decimal(int(List1[14])))
    Percent_Length = str(Decimal(Length) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(BP_Changes) + '\t' + str(Transversions) + '\t0\t0\t0\t0\t'  + Percent_Bases + '\t' + Percent_Match + '\t' + Percent_Length + '\t' + List1[14] + '\t' + List1[8]
    return Out

def Edge_Line(PSL_line, genome_gene, gene, verbose):
    List1 = PSL_line.split('\t')
    BP_Changes = Mutant_Count(genome_gene, gene)
    if BP_Changes == 0:
        Type = 'Contig Edge'
        Description = 'No mutations for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
    elif BP_Changes <= 10 or verbose == True:
        Description = Mutant_Info_Offset(genome_gene, gene, int(List1[15]))
        Description = Description + ' for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Contig Edge (Mutant)'
    elif BP_Changes > 10 and verbose == False:
        Description = str(BP_Changes) + ' mutations for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Contig Edge (Mutant)'
    Start = List1[11]
    Stop = List1[12]
    Length = int(List1[12]) - int(List1[11])
    Transversions = Transversion_Count(genome_gene, gene)
    Percent_Bases = str(Decimal(int(List1[0])) / Decimal(int(List1[0]) + int(List1[1])))
    Percent_Match = str(Decimal(Length - BP_Changes) / Decimal(int(List1[14])))
    Percent_Length = str(Decimal(Length) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(BP_Changes) + '\t' + str(Transversions) + '\t0\t0\t0\t0\t'  + Percent_Bases + '\t' + Percent_Match + '\t' + Percent_Length + '\t' + List1[14] + '\t' + List1[8]
    return Out

def Indel_Line(PSL_line, genome_gene, gene, verbose):
    List1 = PSL_line.split('\t')
    BP_Changes = Indel_BP_Count(PSL_line, genome_gene, gene)
    Mutant_List = Indel_BP_List(PSL_line, genome_gene, gene)
    Indel_Info = Indel_Base_Info(PSL_line)
    Indel_Info_String = Indel_String(Indel_Info)
    if len(Mutant_List) == 0:
        Type = 'Indel'
        Description = 'No mutations, ' + Indel_Info_String
    elif len(Mutant_List) <= 10 or verbose == True:
        Description = Indel_BP_String(PSL_line, genome_gene, gene)
        Description = Description + ', ' + Indel_Info_String
        Type = 'Indel'
    elif len(Mutant_List)> 10 and verbose == False:
        Description = str(len(Mutant_List)) + ' mutations, ' + Indel_Info_String
        Type = 'Indel'
    Start = List1[11]
    Stop = List1[12]
    Length = int(List1[16]) - int(List1[15])
    Transversions = Transversion_Count_Mutant_List(Mutant_List)
    Percent_Bases = str(Decimal(int(List1[0])) / Decimal(int(List1[0]) + int(List1[1])))
##    Percent_Match = str(Decimal(Length - BP_Changes) / Decimal(int(List1[14])))
    Percent_Length = str(Decimal(Length) / Decimal(int(List1[14])))
    Indels_All_List = Indel_Base_Info(PSL_line)
    Indels_All = Indel_Totals(Indels_All_List)
    Indels_All_List = Indels_All.split('\t')
    Percent_Match = str(Decimal(int(List1[0]) - int(List1[1]) - int(Indels_All_List[1]) - int(Indels_All_List[3])) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(len(Mutant_List)) + '\t' + str(Transversions) + '\t' + Indels_All + '\t' + Percent_Bases + '\t' + Percent_Match + '\t' + Percent_Length + '\t' + List1[14] + '\t' + List1[8]
    return Out

def Indel_Partial_Line(PSL_line, genome_gene, gene, verbose):
    List1 = PSL_line.split('\t')
    BP_Changes = Indel_BP_Count(PSL_line, genome_gene, gene)
    Mutant_List = Indel_BP_List(PSL_line, genome_gene, gene)
    Indel_Info = Indel_Base_Info(PSL_line)
    Indel_Info_String = Indel_String(Indel_Info)
    if len(Mutant_List) == 0:
        Type = 'Indel (Partial)'
        Description = 'No mutations, ' + Indel_Info_String + ', for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
    elif len(Mutant_List) <= 10 or verbose == True:
        Description = Indel_BP_String(PSL_line, genome_gene, gene)
        Description = Description + ', ' + Indel_Info_String + ', for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Indel (Partial)'
    elif len(Mutant_List)> 10 and verbose == False:
        Description = str(len(Mutant_List)) + ' mutations, ' + Indel_Info_String + ', for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Indel (Partial)'
    Start = List1[11]
    Stop = List1[12]
    Length = int(List1[16]) - int(List1[15])
    Transversions = Transversion_Count_Mutant_List(Mutant_List)
    Percent_Bases = str(Decimal(int(List1[0])) / Decimal(int(List1[0]) + int(List1[1])))
##    Percent_Match = str(Decimal(Length - BP_Changes) / Decimal(int(List1[14])))
    Percent_Length = str(Decimal(Length) / Decimal(int(List1[14])))
    Indels_All_List = Indel_Base_Info(PSL_line)
    Indels_All = Indel_Totals(Indels_All_List)
    Indels_All_List = Indels_All.split('\t')
    Percent_Match = str(Decimal(int(List1[0]) - int(List1[1]) - int(Indels_All_List[1]) - int(Indels_All_List[3])) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(len(Mutant_List)) + '\t' + str(Transversions) + '\t' + Indels_All + '\t' + Percent_Bases + '\t' + Percent_Match + '\t' + Percent_Length + '\t' + List1[14] + '\t' + List1[8]
    return Out

def Indel_Edge_Line(PSL_line, genome_gene, gene, verbose):
    List1 = PSL_line.split('\t')
    BP_Changes = Indel_BP_Count(PSL_line, genome_gene, gene)
    Mutant_List = Indel_BP_List(PSL_line, genome_gene, gene)
    Indel_Info = Indel_Base_Info(PSL_line)
    Indel_Info_String = Indel_String(Indel_Info)
    if len(Mutant_List) == 0:
        Type = 'Indel (Contig Edge)'
        Description = 'No mutations, ' + Indel_Info_String + ', for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
    elif len(Mutant_List) <= 10 or verbose == True:
        Description = Indel_BP_String(PSL_line, genome_gene, gene)
        Description = Description + ', ' + Indel_Info_String + ', for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Indel (Contig Edge)'
    elif len(Mutant_List)> 10 and verbose == False:
        Description = str(len(Mutant_List)) + ' mutations, ' + Indel_Info_String + ', for bp ' + List1[15] + '-' + List1[16] + ' of '  + List1[14]
        Type = 'Indel (Contig Edge)'
    Start = List1[11]
    Stop = List1[12]
    Length = int(List1[16]) - int(List1[15])
    Transversions = Transversion_Count_Mutant_List(Mutant_List)
    Percent_Bases = str(Decimal(int(List1[0])) / Decimal(int(List1[0]) + int(List1[1])))
##    Percent_Match = str(Decimal(Length - BP_Changes) / Decimal(int(List1[14])))
    Percent_Length = str(Decimal(Length) / Decimal(int(List1[14])))
    Indels_All_List = Indel_Base_Info(PSL_line)
    Indels_All = Indel_Totals(Indels_All_List)
    Indels_All_List = Indels_All.split('\t')
    Percent_Match = str(Decimal(int(List1[0]) - int(List1[1]) - int(Indels_All_List[1]) - int(Indels_All_List[3])) / Decimal(int(List1[14])))
    Out = List1[13] + '\t' + List1[9] + '\t' + Start + '\t' + Stop + '\t' + Type + '\t' + Description + '\t' + str(len(Mutant_List)) + '\t' + str(Transversions) + '\t' + Indels_All + '\t' + Percent_Bases + '\t' + Percent_Match + '\t' + Percent_Length + '\t' + List1[14] + '\t' + List1[8]
    return Out

def Line_Lister(PSL, genome_fasta, gene_fasta, verbose):
    """Makes a list of potential GAMA lines from a PSL file matching Genes to a Genome"""
    f = open(PSL, 'r')
    Lines = []
    for line in f:
        Lines.append(line)
    f.close()
    Genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))
    Genes = SeqIO.to_dict(SeqIO.parse(gene_fasta, 'fasta'))
    Output = []
    for line in Lines:
        List1 = line.split('\t')
        Type = PSL_Type(line)
        Positions = Match_Start_Stop_Finder(line)
        genome = Genome[List1[9]]
        if List1[8] == '-':
            genome = genome.reverse_complement()
        genome = genome[Positions[0][0]:Positions[0][1]]
        gene = Genes[List1[13]][Positions[1][0]:Positions[1][1]]
        if Type == 'Match':
            Out = Match_Line(line, genome, gene, verbose)
        elif Type == 'Indel':
            Out = Indel_Line(line, genome, gene, verbose)
        elif Type == 'Partial':
            Out = Partial_Line(line, genome, gene, verbose)
        elif Type == 'Contig Edge':
            Out = Edge_Line(line, genome, gene, verbose)
        elif Type == 'Indel/Partial':
            Out = Indel_Partial_Line(line, genome, gene, verbose)
        elif Type == 'Indel/Contig Edge':
            Out = Indel_Edge_Line(line, genome, gene, verbose)
##        print(Out)
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

def Best_List(input_contig_list, length_minimum):
    """Finds the best matches for a set of matches to the same contig"""
    Output_List = []
    for items in input_contig_list:
        Add = 1
        Star = 0
        List1 = items.split('\t')
        for items_2 in input_contig_list:
            List2 = items_2.split('\t')
            if items == items_2:
                continue
            elif Overlap_Fraction([int(List1[2]), int(List1[3])], [int(List2[2]), int(List2[3])]) > 0.5:
                if (float(List2[12]) * float(List2[14])) > (float(List1[12]) * float(List1[14])):
                    Add = 0
                    break
                elif (float(List2[12]) * float(List2[14])) == (float(List1[12]) * float(List1[14])):
                    if float(List2[12]) > float(List1[12]):
                        Add = 0
                        break
                    elif float(List2[12]) == float(List1[12]):
                        if float(List2[14]) > float(List1[14]):
                            Add = 0
                            break
                        elif float(List2[14]) == float(List1[14]):
                            if int(List2[6]) < int(List1[6]):
                                Add = 0
                                break
                            elif int(List2[6]) == int(List1[6]):
                                if int(List2[7]) < int(List1[7]):
                                    Add = 0
                                    break
                                elif int(List2[7]) == int(List1[7]):
                                    if (int(List2[9]) + int(List2[11])) < (int(List2[9]) + int(List2[11])):
                                        Add = 0
                                        break
                                    elif (int(List2[9]) + int(List2[11])) == (int(List2[9]) + int(List2[11])):
                                        Star = 1
                                        Order_List = [items, items_2]
                                        Order_List.sort()
                                        if Order_List[0] != items:
                                            Add = 0
                                            break
        if Add == 1 and Star == 0 and float(List1[14]) > (length_minimum / 100):
            Output_List.append(items)
        elif Add == 1 and Star == 1 and float(List1[14]) > (length_minimum / 100):
            List1[0] = List1[0] + '‡'
            Out = '\t'.join(List1)
            Output_List.append(Out)
    return(Output_List)

def GAMMA_S_Output(PSL, genome_fasta, gene_fasta, out_file, verbose, minimum):
    Lines = Line_Lister(PSL, genome_fasta, gene_fasta, verbose)
    Contigs = Contig_Overlaps(Lines)
    Out = open(out_file, 'w')
    Out.write('Gene\tContig\tStart\tStop\tMatch_Type\tDescription\tMismatches\tTransversions\tInsertions\tInsertion_BP\tDeletions\tDeletions_BP\tUnweighted_Match_Percent\tMatch_Percent\tLength_Percent\tTarget_Length\tStrand\n')
    for items in Contigs:
        Best = Best_List(items, minimum)
        for entry in Best:
            Out.write(entry + '\n')
    Out.close()

def GAMMA_S_Output_All(PSL, genome_fasta, gene_fasta, out_file, verbose):
    Lines = Line_Lister(PSL, genome_fasta, gene_fasta, verbose)
    Out = open(out_file, 'w')
    Out.write('Gene\tContig\tStart\tStop\tMatch_Type\tDescription\tMismatches\tTransversions\tInsertions\tInsertion_BP\tDeletions\tDeletions_BP\tUnweighted_Match_Percent\tMatch_Percent\tLength_Percent\tTarget_Length\tStrand\n')
    for entry in Lines:
        Out.write(entry + '\n')
    Out.close()

def GAMMA_S_Out(fasta, gene_db, output, identity, verbose, minimum):
    subprocess.call('blat' + ' ' + gene_db + ' '  + fasta + ' -noHead -minIdentity=' + str(identity) + ' ' + output + '.psl', shell=True)
    GAMMA_S_Output(output + '.psl', fasta, gene_db, output + '.gamma', verbose, minimum)

def GAMMA_S_Out_All(fasta, gene_db, output, identity, verbose):
    subprocess.call('blat' + ' ' + gene_db + ' '  + fasta + ' -noHead -minIdentity=' + str(identity) + ' ' + output + '.psl', shell=True)
    GAMMA_S_Output_All(output + '.psl', fasta, gene_db, output + '.gamma', verbose)

def GAMMA_P_Out(fasta, gene_db, output, identity, verbose, minimum):
    subprocess.call('blat' + ' ' + gene_db + ' '  + fasta + ' -noHead -prot -minIdentity=' + str(identity) + ' ' + output + '.psl', shell=True)
    GAMMA_S_Output(output + '.psl', fasta, gene_db, output + '.gamma', verbose, minimum)

def GAMMA_P_Out_All(fasta, gene_db, output, identity, verbose):
    subprocess.call('blat' + ' ' + gene_db + ' '  + fasta + ' -noHead -prot -minIdentity=' + str(identity) + ' ' + output + '.psl', shell=True)
    GAMMA_S_Output_All(output + '.psl', fasta, gene_db, output + '.gamma', verbose)

def main():
    parser = argparse.ArgumentParser(description="""This scripts makes sequence (nucleotide or protein) match calls from matches in an assembly using a sequence database""")
    parser.add_argument("input_fasta", type=str, help="input fasta")
    parser.add_argument("database", type=str, help="input database")
    parser.add_argument("output", type=str, help="output name")
    parser.add_argument("-a", "--all", action="store_true",
                        help="include all gene matches, even overlaps")
    parser.add_argument("-e", "--extended", action="store_true",
                        help="writes out all protein mutations")
    parser.add_argument("-p", "--protein", action="store_true",
                        help="for protein-protein comparisons")
    parser.add_argument("-m", "--minimum", default=20, type=int,
                        help="minimum length percent match for output (default = 20)")
    parser.add_argument("-i", "--percent_identity", default=90, type=int,
                        help="minimum nucleotide identity for blat search (default = 90)")
    args = parser.parse_args()
    if args.extended:
        Verbose = True
    else:
        Verbose = False
    if args.all and args.protein:
        GAMMA_S_Out_All(args.input_fasta, args.database, args.output,  args.percent_identity, Verbose)
    elif args.all:
        GAMMA_S_Out_All(args.input_fasta, args.database, args.output,  args.percent_identity, Verbose)
    elif args.protein:
        GAMMA_P_Out(args.input_fasta, args.database, args.output,  args.percent_identity, Verbose, args.minimum)
    else:
        GAMMA_S_Out(args.input_fasta, args.database, args.output,  args.percent_identity, Verbose, args.minimum)

if __name__ == "__main__":
    sys.exit(main())
