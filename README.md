# motif-mark

### Overview

A program to create visual representations for protein binding motifs in a gene. The user provides a fasta file of gene sequences and a text file of motifs to search for in the genes. The program outputs a .svg image with the location of motifs visualized for each gene. All genes are scaled to be visually represented at the same size, but relative sizes of motifs and exons are preserved. Only provided sequences are searched, reverse compliments are not considered. Motifs are represented with colored boxes spanning the full region of the gene that matches the motif, which may be longer than the motif sequence. No mismatches are allowed. Ambiguous bases (N) in the gene sequence are not treated as matches. Exons are represented by grey boxes, introns are black lines. Note that motif colors are randomly generated, and will change on reruns. This program was written for an assignment. 

### Runtime options

| Command | Input  |
|---------|--------|
| -f      | FASTA  |
| -m      | Motifs |

### Input requirements

Motifs must be encoded with [IUPAC nucleotides](https://www.bioinformatics.org/sms/iupac.html). Motifs should be provided in a .txt file with one motif per line. Motifs can be either RNA or DNA.
Gene sequences should be provided in FASTA format. Introns should be lowercase and exons uppercase characters. Gene sequences can be either DNA or RNA.

### Additional python packages required
cairo
