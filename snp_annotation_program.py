
################################################### IMPORTS ###########################################################
import sqlite3, vcf, logging, gffutils, argparse, os, math
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from Bio import BiopythonWarning
################################################### PARSER ###########################################################
# parse arguments
parser = argparse.ArgumentParser(description='BCPy Assessment')
parser.add_argument('--vcfFile', required=True, help='Input VCF file')
parser.add_argument('--gffFile', required=True, help='Input GFF file')
parser.add_argument('--fastaFile', required=True, help='Input FASTA file')
parser.add_argument('--outputPrefix', required=True, help='Prefix for output files')
args = parser.parse_args()

################################################### LOGGER ###########################################################
# Check if the log file exists and delete it
log_file_path = f'{args.outputPrefix}.log'
if os.path.isfile(log_file_path):
    os.remove(log_file_path)

# Set up information LOGGER - file handler
infoLogger = logging.getLogger('infoLogger')
infoLogger.setLevel(logging.INFO)

# Create a file handler only if the log file does not exist
if not os.path.isfile(log_file_path):
    fh = logging.FileHandler(log_file_path)
    fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
    infoLogger.addHandler(fh)

    # Log that a new log file is being created
    infoLogger.info(f'Initialising program ...\n')
    infoLogger.info(f'Creating new log file {log_file_path} ...\n')

################################################### VCF File ###########################################################
infoLogger.info(f'Opening vcf file {args.vcfFile}\n')

try:
    vcfReader = vcf.Reader(filename=args.vcfFile)
except FileNotFoundError:
    infoLogger.error(f'File {args.vcfFile} cannot be opened for reading. Please check and try again')
    raise SystemExit(1)

# open the gff file 
infoLogger.info(f'Opening gff file {args.gffFile}\n')

################################################### GFF File ###########################################################
db = args.gffFile.replace('.gff', '.db')

# check if the database exists
if not os.path.isfile(db):
    infoLogger.info('Creating GFF database {}...\n'.format(db))
    try:
        db = gffutils.create_db(args.gffFile, dbfn=db, force=True, keep_order=True)
    except sqlite3.OperationalError:
        infoLogger.error('Cannot create database {}\n'.format(db))
        raise SystemExit(1)
    except ValueError:
        infoLogger.error('Cannot create dataabase {}\n'.format(db))
        raise SystemExit(1)
else:
    infoLogger.info('Connecting to existing database {}...\n'.format(db))
    try:
        db = gffutils.FeatureDB(db, keep_order=True)
    except ValueError:
        infoLogger.error('DB {} could not be read\n'.format(db))
        raise SystemExit(1)

################################################### TSV File ###########################################################
# open the output file
try:
    tsvOut = open(f'{args.outputPrefix}.tsv', 'wt')
    infoLogger.info(f'Opening TSV file {args.outputPrefix}.tsv\n')
except FileNotFoundError:
    infoLogger.error(f'File {args.outputPrefix} cannot be opened for writing. Please check and try again')
    raise SystemExit(1)

#writing the headers
tsvOut.write(f'Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein Location\tRef AA\tAlt AA\n')

################################################### CLASSES ###########################################################

class Sequence:
    def __init__(self, sequence):
        self._sequence = sequence

    @property
    def sequence(self):
        return self._sequence

    def protein_seq(self):
        return str(Seq(self.sequence).translate())

    def protein_pos(self, MmRNA_pos):
        return math.ceil(MmRNA_pos / 3)

    def reference_aa(self, protein_pos):
        codon_list = [self.sequence[i:i + 3] for i in range(0, len(self.sequence), 3)]
        codon_dict = {i + 1: codon for i, codon in enumerate(codon_list)}
        return codon_table[codon_dict[protein_pos]]

class FSequence(Sequence):
    def alt_aa(self, MmRNA_pos, alt_base):
        seq_dict = {index + 1: base for index, base in enumerate(self.sequence)}
        seq_dict[MmRNA_pos] = str(alt_base)
        alt_seq = ''.join(seq_dict.values())
        alt_codon_list = [alt_seq[i:i + 3] for i in range(0, len(alt_seq), 3)]
        alt_codon_dict = {i + 1: codon for i, codon in enumerate(alt_codon_list)}
        protein_pos = math.ceil(MmRNA_pos / 3)
        return codon_table[alt_codon_dict[protein_pos]]

class RSequence(Sequence):
    def alt_aa(self, MmRNA_pos, alt_base):
        #sequence dictionary (1 indexed)
        seq_dict = {index: base for index, base in enumerate(self.sequence, start=1)}
        # complement of the alt base
        alt_base_complement = str(Seq(str(alt_base)).complement())
        # changing the reference base to the alt base
        seq_dict[MmRNA_pos] = alt_base_complement
        # joining the dictionary values to form the alt sequence
        alt_seq = ''.join(seq_dict.values())
        #alternative codon list
        alt_codon_list = [alt_seq[i:i+3] for i in range(0, len(alt_seq), 3)]
        #alternative codon dictionary (1 indexed)
        alt_codon_dict = {i+1: codon for i, codon in enumerate(alt_codon_list)}
        protein_pos = math.ceil(MmRNA_pos / 3)
        return codon_table[alt_codon_dict[protein_pos]]


################################################### UTILS ###########################################################
#codon table
bases = "tcag".upper()
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
#supressing Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)
################################################### PROCESSING CODE ###########################################################

print('Program started...')

# less than QUAL 20 count
less_20_qual_count = 0
#psudogene count
psudogene_count = 0

#only consider high quality SNPs
for variant in vcfReader:
    if int(variant.QUAL) > 20:
        #feature is the SNP
        if len(list(db.region(seqid=variant.CHROM, start=variant.POS, end=variant.POS, featuretype='CDS')))>0:
            for feature in db.region(seqid=variant.CHROM, start=variant.POS, end=variant.POS, featuretype='CDS'):
                # find the CDS that the SNP is in
                for transcript in db.parents(feature.id, featuretype='mRNA'):
                    MmRNA_pos = 0 
                    seq = '' 
                    count = 0
                    # CDS on POSITIVE strand
                    if transcript.strand == '+':
                        for cds in db.children(transcript.id, featuretype='CDS', order_by='start'):

                            #get the sequence of each CDS and append it to the sequence container
                            seq = seq + cds.sequence(args.fastaFile, use_strand=True)

                            # coordinates of SNP relative to the Mature mRNA
                            if cds == feature:
                                cdsLength = int(variant.POS) - cds.start + 1
                                MmRNA_pos += cdsLength
                                # after the CDS where the SNP is break to exit the loop
                                break
                            else:
                                # for CDSs not coitaining the SNP, get the length and add to MmRNA_pos
                                cdsLength = cds.end - cds.start + 1
                                MmRNA_pos += cdsLength

                        #creating an instace of the FSequence class
                        cds_seq = FSequence(seq)
                        protein_seq = cds_seq.protein_seq()
                        
                        if protein_seq.startswith('M') and protein_seq.endswith('*') and protein_seq.count('*') == 1:                                 
                            #position in the protein where the SNP falls
                            protein_pos = cds_seq.protein_pos(MmRNA_pos=MmRNA_pos)
                            #reference AA
                            reference_aa = cds_seq.reference_aa(protein_pos=protein_pos)

                            #alternative AA
                            alt_aa = cds_seq.alt_aa(MmRNA_pos=MmRNA_pos, alt_base=variant.ALT[0])

                            if reference_aa != alt_aa:
                                variant_type = "Non-synonymous"
                            elif reference_aa == alt_aa:
                                variant_type = "Synonymous"
                                alt_aa = "NA"
                            
                            tsvOut.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF[0]}\t{variant.ALT[0]}\t{variant_type}\t{transcript.id}\t{protein_pos}\t{reference_aa}\t{alt_aa}\n')
                        else:
                            psudogene_count += 1
                            continue
                    # CDS on NEGATIVE strand
                    elif transcript.strand == '-':
                        for cds in db.children(transcript.id, featuretype='CDS', order_by='start', reverse=True):
                            
                            #get the sequence of each CDS and append it to the sequence container
                            seq = seq + cds.sequence(args.fastaFile, use_strand=True)

                            # coordinates of SNP relative to the Mature mRNA
                            if cds == feature:
                                cdsLength =  cds.end - int(variant.POS) + 1
                                MmRNA_pos += cdsLength
                                
                                # we don't want to consider any CDSs after this, so break to exit the loop
                                break
                            else:
                                # for CDSs not coitaining the SNP, get the length and add to MmRNA_pos
                                cdsLength = cds.end - cds.start + 1
                                MmRNA_pos += cdsLength

                        cds_seq = RSequence(seq)
                        protein_seq = cds_seq.protein_seq()
                        
                        if protein_seq.startswith('M') and protein_seq.endswith('*') and protein_seq.count('*') == 1:
                            #position in the protein where the SNP falls
                            protein_pos = cds_seq.protein_pos(MmRNA_pos=MmRNA_pos)
                            #reference AA
                            reference_aa = cds_seq.reference_aa(protein_pos=protein_pos)

                            #alternative AA
                            alt_aa = cds_seq.alt_aa(MmRNA_pos=MmRNA_pos, alt_base=variant.ALT[0])

                            if reference_aa != alt_aa:
                                variant_type = "Non-synonymous"
                            elif reference_aa == alt_aa:
                                variant_type = "Synonymous"
                                alt_aa = "NA"
                        
                            tsvOut.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF[0]}\t{variant.ALT[0]}\t{variant_type}\t{transcript.id}\t{protein_pos}\t{reference_aa}\t{alt_aa}\n')
                        else:
                            psudogene_count += 1
                            continue
        else:
            variant_type = "Non-coding"
            transcript_id = "NA"
            protein_pos = "NA"
            reference_aa = "NA"
            alt_aa = "NA"
            tsvOut.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF[0]}\t{variant.ALT[0]}\t{variant_type}\t{transcript_id}\t{protein_pos}\t{reference_aa}\t{alt_aa}\n')
    else:
        less_20_qual_count += 1

infoLogger.info(f'There are {less_20_qual_count} variants with quality < 20\n')
infoLogger.info(f'There are {psudogene_count} pseudogenes.\n')

################################################### GRAPH ###########################################################

# Read the TSV file into a DataFrame
try:
    tsvTable = pd.read_csv(f'{args.outputPrefix}.tsv', sep='\t')
except FileNotFoundError:
    infoLogger.error(f'File {args.outputPrefix} cannot be opened for reading. Please check and try again')
    raise SystemExit(1)

# Group by 'Type' and calculate the proportion
type_proportions = tsvTable['Type'].value_counts(normalize=True)
# Plotting the bar plot
type_proportions.plot(kind='bar', color=['yellowgreen', 'cornflowerblue', 'indianred'])
# labels and title
plt.xlabel('Variant Type')
plt.ylabel('Proportion')
plt.title('Proportion of Variants by Type')
# Save the plot as an image in the current directory
# Rotate x-axis labels for better readability
plt.xticks(rotation=0)  
plot_name = f'{args.outputPrefix}_type_proportions_plot.png'  # Name the plot
try:
    plt.savefig(plot_name)  # Save the plot as 'type_proportions_plot.png'
except FileNotFoundError:
    infoLogger.error(f'File {plot_name} cannot be saved. Please check and try again')
    raise SystemExit(1)

################################################### SAVING FILES ###########################################################
current_directory = os.getcwd()
infoLogger.info(f'The image {plot_name} was succesfully saved at: {current_directory}\n')
infoLogger.info(f'TSV file {args.outputPrefix}.tsv was sucessfully saved at: {current_directory}.\n')
infoLogger.info(f'TSV file {args.outputPrefix}.log was sucessfully saved at: {current_directory}.\n')
infoLogger.info(f'Program finished.\n')

print('Program finished successfully.')