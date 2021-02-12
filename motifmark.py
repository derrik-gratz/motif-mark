#!/usr/bin/env python3

import argparse
import re
import cairo
import colorsys
from random import random
from math import ceil


class Gene:
    '''
    Object for each fasta sequence passed in. Stores exons and motifs identified in each sequence.
    '''
    def __init__(self):
        self.exons = []
        self.motifs = {}
        self.sequence = ""
    def add_exon(self, exons):
        self.exons = exons
    def add_motif(self, motif, windows):
        self.motifs[motif] = windows
    def append_sequence(self, seq):
        self.sequence += seq.strip()


def get_args():
    '''
    Takes user arguments at runtime. One argument for the input fasta file to scan, one for the output
    directory, and one for the file of motifs to search for. 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', "--file", type=str, help="The input sequences in fasta format")
    parser.add_argument('-m', "--motifs", type=str, help="Text file containing motifs, 1 per line")
    return parser.parse_args()


def get_motifs(motifs):
    '''
    Translates the provided motifs into Regex search patterns using IUPAC nucleotide codes
    '''
    iupac_dict = {
        "A":"[Aa]",
        "C":"[Cc]",
        "G":"[Gg]",
        "T":"[TtUu]",
        "U":"[UuTt]",
        "W":"[AaTtUu]",
        "S":"[CcGg]",
        "M":"[AaCc]",
        "K":"[GgTtUu]",
        "R":"[AaGg]",
        "Y":"[CcTtUu]",
        "B":"[CcGgTtUu]",
        "D":"[AaGgTtUu]",
        "H":"[AaCcTtUu]",
        "V":"[AaCcGg]",
        "N":"[AaCcGgTtUu]",
        "Z":"[]",
    }
    translated_motifs = {}
    with open(motifs, "r") as fh:
        for line in fh:
            motif = line.strip()
            translated_motif = ""
            for character in motif.upper():
                translated_motif += iupac_dict[character]
            translated_motifs[translated_motif] = motif
    #outputs the regex pattern for the motif and a key to turn the regex string back into the original motif
    return translated_motifs


def exon_finder(sequence):
    '''
    Finds exons in fasta sequences by looking for capitalized nucleotides. Returns the start and
    stop positions of capitialized stretches within the sequence.
    '''
    iterator = re.finditer("[A-Z]+", sequence)
    # returns start and stop position, 0 indexed
    exons = []
    for match in iterator:
        exons.append(match.span())
    # exons as list of tuples with start stop position
    return exons


def motif_finder(sequence, motif):
    '''
    Finds regions in a sequence that match a motif. Return the start and stop positions 
    of each matching region. 
    '''
    #Number of nucleotides in motif
    motif_len = len(re.findall("\[[a-zA-Z]*\]", motif))
    #number of nucleotides in query sequence
    seq_len = len(sequence)
    #used for windows larger than the motif sequence
    still_in_motif = False
    motif_count = 0
    #store start and stop locations for each motif
    motif_locations = {}
    for i in range(0, seq_len-motif_len+1):
        sliding_window = sequence[i:i+motif_len]
        if re.fullmatch(motif, sliding_window):
            #sliding window in fasta matches the motif
            if not still_in_motif:
                #first instance of motif
                #adjusting positions for 0 base counting
                motif_start = i + 1
                motif_end = i + 1 + motif_len
                still_in_motif = True
                motif_count += 1
            elif still_in_motif:
                #consecutive encounters of this motif
                #motif window extends beyond length of motif, update the end position
                motif_end = i + 1 + motif_len
        elif still_in_motif:
            #end of motif window, store bounds
            motif_locations[motif_count] = [motif_start, motif_end]
            still_in_motif = False
    #return just the positions, as only 1 motif is provided at a time
    return [(locations) for number, locations in motif_locations.items()]


def generate_pallete(motifs):
    '''
    Create list of randomly generated RBG values for each motif, to be used in graphical representation
    '''
    colors = {}
    for motif in motifs:
        r,g,b,a = random(), random(), random(), 0.5
        colors[motif] = [r,g,b,a]
    return colors


def drawing(gene_objects, output_name, pallete, sequence_file, motif_file, motifs):
    '''
    All functions for the pycairo graphic
    Takes in each of the gene objects (with exons and motifs stored), the name of the file to
    write to, the color pallete, the original motifs in IUPAC form, and the names of the 
    input files.
    Outputs a SVG. Each fasta sequence from the input file gets its own drawing. The drawings
    are depicted as the same size, but exon and motif bounds are scaled in size relative to the 
    length of the gene. Each gene drawing shows the exons as grey boxes and introns as straight dark lines.
    Motifs are represented by colored boxes spanning the window which matched the motif.
    A key is provided at the bottom to associate motifs with their color. 
    The original file names are written to the output for record keeping. 
    '''
    #determine how many genes need to be represented
    gene_count = len(gene_objects) + 1
    #height for each gene
    base_height = 100
    #dimensions of total output, scaled for number of genes
    width, height = 600, base_height * (gene_count +1)
    #name of output
    filename = output_name + ".svg"
    #leaving whitespace for ease of viewing
    left_indent = width * .1
    #width of gene in visual, relative to width of whole image
    gene_width = width * .8
    with cairo.SVGSurface(filename, width, height) as surface:
        #create visual object
        context = cairo.Context(surface)
        #which gene were on, used for scaling height coordinates
        gene_count = 1
        for gene in gene_objects:
            #make a straight line for each gene
            context.set_source_rgb(0,0,0)
            x1, y1 = width * .1, base_height * .5 + base_height * gene_count
            x2, y2 = width * .9, base_height * .5 + base_height * gene_count
            context.set_line_width(1)
            context.move_to(x1,y1)
            context.line_to(x2,y2)
            #add the name of each gene/fasta entry above the visual
            x3, y3 = width * .09, base_height * .1 + base_height * gene_count
            context.move_to(x3, y3)
            gene_name = gene + "    {} BP".format(len(gene_objects[gene].sequence))
            context.show_text(gene_name)
            context.stroke()
            #find length of gene to scale motifs and exons accordingly
            gene_length = len(gene_objects[gene].sequence)
            #add the exon as a rectangle
            for exon in gene_objects[gene].exons:
                exon_start, exon_stop = exon
                #scale the start position for the graphic
                exon_x1 = (exon_start / gene_length) * gene_width + left_indent
                exon_y1 = base_height * 0.3 + base_height * gene_count
                exon_height = base_height * 0.4
                exon_width = (exon_stop - exon_start) / gene_length * gene_width
                #make a grey box that is wider than the intron lines
                context.rectangle(exon_x1, exon_y1, exon_width, exon_height)
                context.set_source_rgb(.6,.6,.6)
                context.fill()
            #now iterating through each motif found in this sequence
            for motif in gene_objects[gene].motifs:
                #grabbing BP coordinates for where the motifs matched the sequence (often more than one)
                for window in gene_objects[gene].motifs[motif]:
                    #output the motifs as boxes
                    re_x1 = (((window[0] / gene_length) * gene_width) + left_indent)
                    re_width = (window[1] - window[0]) / gene_length * gene_width
                    re_height = base_height * .2
                    context.rectangle(re_x1, base_height * .4 + base_height * gene_count, re_width, re_height)
                    #uses the same color for same motifs
                    r,g,b,a = pallete[motif]
                    context.set_source_rgba(r,g,b,a)
                    context.fill()
            gene_count += 1
        #adding figure title
        context.move_to(width * 0.1, base_height * 0.5)
        context.set_source_rgb(0,0,0)
        context.set_font_size(10.0)
        context.show_text("Identifying motifs described in the file '{}' found in the file '{}'".format(motif_file, sequence_file))
        context.move_to(width * 0.1, base_height * 0.7)
        context.set_source_rgb(0,0,0)
        context.set_font_size(7.0)
        context.show_text("Grey boxes represent exons, and colored boxes represent motifs. See motif legend below.")
        #adding the motif legend
        legend_x, legend_y = left_indent, base_height * 0.2 + base_height * gene_count
        legend_horiz_spacing = width * 0.2
        legend_vert_spacing = base_height * 0.1
        context.move_to(legend_x, legend_y - legend_vert_spacing)
        context.set_font_size(10)
        context.show_text("Motif legend")
        motif_count = 0
        line = 0
        for motif in motifs:
            #wrapping so 4 motifs appear in each line
            if motif_count % 4 == 0:
                line += 1
            x1, y1 = (legend_x + (motif_count % 4) * legend_horiz_spacing), (legend_y + legend_vert_spacing * line)
            context.move_to(x1, y1)
            #cute lil boxes
            context.rectangle(x1,y1,base_height*0.05, base_height * 0.05)
            r,g,b,a = pallete[motifs[motif]]
            context.set_source_rgba(r,g,b,a)
            context.fill()
            #write the original motif
            x1, y1 = x1 + base_height * 0.07, y1 + base_height * 0.07
            context.move_to(x1,y1)
            context.set_source_rgb(0,0,0)
            context.show_text(motifs[motif].upper())
            motif_count += 1


def main():
    #collect user parameters
    args = get_args()
    #get the regex form of motifs for searching 
    motifs = get_motifs(args.motifs)
    #initiallize dictionary of gene objects
    genes = {}
    with open(args.file, "r") as fh:
        #iterating through input fasta
        for line in fh:
            if line[0] == ">":
                #New fasta entry
                current_gene = line.strip()
                #define new object
                genes[current_gene] = Gene()
            else:
                #add the fasta sequence to the gene object
                genes[current_gene].append_sequence(line.strip())
    for gene in genes:
        #iterating through objects
        #find the exons in the gene's sequence
        exons = exon_finder(genes[gene].sequence)
        genes[gene].add_exon(exons)
        for regex_motif, nucleotide_motif in motifs.items():
            #loop through motifs
            motif_locations = motif_finder(genes[gene].sequence, regex_motif)
            if motif_locations:
                #if motif appears at least once in sequence, store it in object
                genes[gene].add_motif(nucleotide_motif, motif_locations)
    #making a color pallete for visual representation of motifs
    pallete = generate_pallete(motifs.values())
    #the name of the file to output
    output_name = args.file.split(".")[0]
    #passing in gene objects, color pallete, motifs, filenames, etc into drawing function for visual output
    drawing(genes, output_name, pallete, args.file, args.motifs, motifs)


main()
