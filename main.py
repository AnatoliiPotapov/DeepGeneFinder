from Bio import SeqIO
from matplotlib import pyplot as plt
from BCBio import GFF
import numpy as np
from dna_features_viewer import GraphicFeature, GraphicRecord

# Reading file
def read_fasta(filename):
    seq = []
    with open(filename, "r") as handle:
        for s in SeqIO.parse(handle, "fasta"):
            seq.append(s)
    return seq

def read_gff(filename):
    gff = []
    with open(filename, "r") as in_handle:
        for rec in GFF.parse(in_handle):
            gff.append(rec)
    return gff

class Markup():
    def __init__(self, id, markup):
        self.id = id
        self.m = markup

class BacterialRefseq():
    def __init__(self, fasta_file, gff_file):
        self.seq = read_fasta(fasta_file)
        self.gff = read_gff(gff_file)
        self.markup = []

    def generate_markup(self):
        for seq in self.seq:
            print('Processing {}'.format(seq.id) )
            gff_record = next((x for x in self.gff if str(x.id) == str(seq.id)), None)
            if gff_record:
                # Filter 'gene' features
                genes = [x for x in gff_record.features if x.type=='gene']
                print('Found {} genes'.format(len(genes)))
                # Prepare numpy matrix for markup
                markup = np.zeros((6, len(seq.seq)))
                for gene in genes:
                    start, end, strand = gene.location.start, gene.location.end, gene.location.strand
                    if strand==1:
                        markup[start % 3, start:end] = 1
                    if strand==-1:
                        markup[start % 3 + 3, start:end] = 1
                self.markup.append(Markup(seq.id, markup))

    def visualize_markup(self, index=1):
        sequences = []
        for seq in self.seq:
            gff_record = next((x for x in self.gff if str(x.id) == str(seq.id)), None)
            if gff_record:
                # Filter 'gene' features
                genes = [x for x in gff_record.features if x.type=='gene']
                features = []
                for gene in genes:
                    start, end, strand = gene.location.start, gene.location.end, gene.location.strand
                    features.append(
                        GraphicFeature(start, end, strand, label=gene.qualifiers['Name'], color="#cffccc")
                    )
            sequences.append(features)
        output_file("test.html")
        record = GraphicRecord(sequence_length=1000, features=sequences[index])
        show(record.plot_with_bokeh(figure_width=5))

class TrainingDataGenerator():
    def __init__():
        pass

    def generate_dataset():
        pass



b = BacterialRefseq('data/Ecoli.fna', 'data/Ecoli.gff')
b.generate_markup()
b.visualize_markup()
print(b.markup)


class MLGeneFinder():
    def __init__():
        pass

    def train(self, seq, markup):
        pass

    def predict(self, seq):
        pass



for f in db.features_of_type('CDS', order_by='start'):
    print(f)
