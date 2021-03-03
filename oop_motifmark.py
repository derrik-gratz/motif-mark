#!/usr/bin/env python

import argparse
import cairo
from random import random
import re

# Global constants for image output
LEFT_MARGIN = 20
GENE_HEIGHT = 100
VERTICAL_BUFFER = GENE_HEIGHT / 5

class GeneGroup:
    ''' Object to hold other objects. One GeneGroup for each fasta sequence in the input.'''
    def __init__(self):
        self.gene = Gene()
        self.exon = Exon()
        self.header = FastaHeader()
        self.motifs = []
        self.number = 0


class FastaHeader():
    '''FASTA sequence header name. Contains method for writting visual output above each
    gene image.'''
    def __init__(self):
        self.name = ""
    def draw(self, context, gene_count):
        x = LEFT_MARGIN
        y = GENE_HEIGHT * gene_count
        context.set_source_rgb(0,0,0)
        context.move_to(x,y)
        context.set_font_size(10)
        context.show_text(self.name)


class Exon:
    '''Exon object stores location of exon in FASTA sequence. Can handle multiple exons.
    Also contains function for visualizing exons as dark boxes in output'''
    def __init__(self):
        self.locations = []
    def draw(self, context, gene_number):
        for location in self.locations:
            start, stop = location
            x = LEFT_MARGIN + start
            y = GENE_HEIGHT * gene_number + VERTICAL_BUFFER*2
            width = stop - start
            height = VERTICAL_BUFFER
            context.set_source_rgb(0,0,0)
            context.rectangle(x,y,width,height)
            context.fill()


class Gene:
    '''The gene sequence. Finds its width and scales the output image size accordingly.
    Also contains method to write gene visualization as a line.'''
    def __init__(self):
        self.sequence = ""
        self.width = 0
    def append_sequence(self, seq):
        self.sequence += seq
        self.width = len(self.sequence)
    def draw_gene(self, context, gene_number):
        x = LEFT_MARGIN
        y = GENE_HEIGHT * gene_number + (GENE_HEIGHT / 2)
        context.set_source_rgb(0,0,0)
        context.set_line_width(1)
        context.move_to(x,y)
        context.line_to(x + self.width, y)
        context.stroke()


class Motif:
    '''Object for each possible motif for each FASTA entry. If there are 4 motifs being
    searched for, each entry can have 4 motif objects. The motif object stores all 
    instances of that motif in that sequence. Visualizes motifs in the output as colored
    boxes.'''
    def __init__(self):
        self.sequence = ""
        self.locations = []
    def draw_motif(self, context, gene_count, pallete):
        for location in self.locations:
            start, stop = location
            width = stop - start
            x = LEFT_MARGIN + start
            y = GENE_HEIGHT * gene_count + (GENE_HEIGHT/2) - 20
            height = 40
            r,g,b,a = pallete[self.sequence]
            context.set_source_rgba(r,g,b,0.7)
            context.set_line_width(1)
            context.move_to(x,y)
            context.rectangle(x, y, width, height)
            context.fill()
            context.stroke()


def get_args():
    '''
    Takes user arguments at runtime. One argument for the input fasta file to scan, one for the output
    directory, and one for the file of motifs to search for. 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', "--file", type=str, help="The input sequences in fasta format", required=True)
    parser.add_argument('-m', "--motifs", type=str, help="Text file containing motifs, 1 per line", required=True)
    parser.add_argument('-c', "--color_randomizer", help="Optional flag to randomly generate colors. Needed for inputs of >10 motifs", action="store_true")
    return parser.parse_args()


def get_motifs(motifs):
    '''    Translates the provided motifs into Regex search patterns using IUPAC nucleotide codes'''
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
    motif_number = {}
    motif_count = 1
    with open(motifs, "r") as fh:
        for line in fh:
            motif = line.strip()
            motif_number[motif_count] = motif
            translated_motif = ""
            for character in motif.upper():
                translated_motif += iupac_dict[character]
            translated_motifs[motif] = translated_motif
            motif_count += 1
    #outputs the regex pattern for the motif and a key to turn the regex string back into the original motif
    return motif_number, translated_motifs


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


def generate_pallete(motifs, randomizer):
    '''
    Create list of pre-approved generated RBG values for each motif, to be used in graphical representation.
    If The user opts for random colors, they are randomly generated.
    '''
    pre_approved_colors = [[0.94, 0.57, 0.91, 1],
        [0.72, 0.89, 0.43, 1],
        [0.48, 0.8, 0.76, 1],
        [0.16, 0.54, 0.74, 1],
        [0.70, 0.34, 0.02, 1],
        [0.87, 0.50, 0.07, 1],
        [0.99, 0.72, 0.38, 1],
        [0.99, 0.87, 0.74, 1],
        [0.84, 0.85, 0.92, 1],
        [0.69, 0.67, 0.82, 1],
        [0.50, 0.45, 0.67, 1],
        [0.32, 0.15, 0.53,1],]
    colors = {}
    motif_count = 0
    for motif in motifs:
        if not randomizer:
            colors[motif] = pre_approved_colors[motif_count]
        else:
            r,g,b,a = random(), random(), random(), 0.7
            colors[motif] = [r,g,b,a]
        motif_count += 1 
    return colors


def create_genes(gene_fasta, motif_number, motif_regex):
    '''
    Iterates through the input FASTA. For each record, the objects are created.
    No output is generated at this time, just getting all information.
    '''
    #initiallize dictionary of gene objects
    gene_groups = {}
    gene_count = 0
    with open(gene_fasta, "r") as fh:
        #iterating through input fasta
        for line in fh:
            if line[0] == ">":
                #New fasta entry for a gene
                gene_count += 1
                #define new object
                gene_groups[gene_count] = GeneGroup()
                #add the gene name
                gene_groups[gene_count].header.name = line.strip()
                #which gene is this in the fasta
                gene_groups[gene_count].number = gene_count
            else:
                #add the fasta sequence to the gene object
                gene_groups[gene_count].gene.append_sequence(line.strip())
    gene_count = 0
    for a in gene_groups:
        #iterating through gene objects
        gene_group = gene_groups[a]
        #find the exons in the gene's sequence
        exons = exon_finder(gene_group.gene.sequence)
        #store exon locations in gene_group object
        gene_group.exon.locations = exons
        #loop through motifs
        for motif_count in motif_number:
            #ACTG motif
            standard_motif = motif_number[motif_count]
            #Regex search string
            regex_motif = motif_regex[standard_motif]
            #list of BP positions for this motif in this gene
            motif_locations = motif_finder(gene_group.gene.sequence, regex_motif)
            if motif_locations:
                #if motif appears at least once in sequence, store it in object
                motif = Motif()
                motif.sequence = standard_motif
                motif.locations = motif_locations
                #add the motif to a list for that gene_group
                gene_group.motifs.append(motif)
        gene_count += 1
    return gene_groups


def drawing(gene_objects, output_file, pallete, input_file, motif_number):
    '''
    With gene objects generated, now we iterate through the GeneGroup objects
    to visualize all the child objects in the groups. At the end, it calls
    add_legend to add the figure legend.
    '''
    #figure out how wide to make the output based on the input gene lengths
    output_width = 0
    for genes in gene_objects:
        if gene_objects[genes].gene.width > output_width:
            output_width = gene_objects[genes].gene.width
    with cairo.SVGSurface(output_file, output_width + LEFT_MARGIN, (len(gene_objects) +2) * GENE_HEIGHT) as surface:
        context = cairo.Context(surface)
        #loop through gene_groups
        gene_count = 0
        for count in gene_objects:
            gene_count += 1
            gene_group = gene_objects[count]
            gene_group.gene.draw_gene(context, count)
            gene_group.exon.draw(context, count)
            gene_group.header.draw(context, count)
            for motif in gene_group.motifs:
                motif.draw_motif(context, count, pallete)
        add_legend(context, motif_number, gene_count+1, pallete)


def add_legend(context, motif_number, gene_count, pallete):
    '''
    Messy, but only needs calling once. Adds a visual key at the bottom of the output
    to show what the colored motif boxes correspond to.
    '''
    legend_x, legend_y = LEFT_MARGIN, VERTICAL_BUFFER + GENE_HEIGHT * gene_count
    legend_horiz_spacing = 100
    legend_vert_spacing = GENE_HEIGHT * 0.1
    context.move_to(legend_x, legend_y - legend_vert_spacing)
    context.set_font_size(10)
    context.set_source_rgba(1,1,1,1)
    context.show_text("Motif legend")
    motif_count = 0
    line = 0
    for motif in motif_number:
        #wrapping so 4 motifs appear in each line
        if motif_count % 4 == 0:
            line += 1
        x1, y1 = (legend_x + (motif_count % 4) * legend_horiz_spacing), (legend_y + legend_vert_spacing * line)
        context.move_to(x1, y1)
        #cute lil boxes
        context.rectangle(x1,y1,GENE_HEIGHT*0.05, GENE_HEIGHT * 0.05)
        r,g,b,a = pallete[motif_number[motif]]
        context.set_source_rgba(r,g,b,a)
        context.fill()
        #write the original motif
        x1, y1 = x1 + GENE_HEIGHT * 0.07, y1 + GENE_HEIGHT * 0.07
        context.move_to(x1,y1)
        context.set_source_rgb(0,0,0)
        context.show_text(motif_number[motif].upper())
        motif_count += 1


def main():
    #collect user parameters
    args = get_args()
    #get the regex form of motifs for searching 
    motif_number, motif_regex = get_motifs(args.motifs)
    #making a color pallete for visual representation of motifs
    pallete = generate_pallete(motif_number.values(), args.color_randomizer)
    #create gene objects
    gene_objects = create_genes(args.file, motif_number, motif_regex)
    #the name of the file to output
    output_name = args.file.split(".")[0]
    output_name += ".svg"
    #passing in gene objects, color pallete, filenames, etc into drawing function for visual output
    drawing(gene_objects, output_name, pallete, args.file, motif_number)


main()