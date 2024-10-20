# snp_annotation_program

This program produces annotations for several SNPs in coding regions. It also generates a log file with relevant information and an image with the proportions of non-coding, synonymous, and non-synonymous variants.

# How the Script Works

This program analyzes SNPs in a `Data.vcf.gz` type file, focusing on variants with a QUAL field higher than 20. The total count of SNPs that are QUAL ≤ 20 is reported to an `Output.log` type file.

It queries `Data.gff` for CDSs (Coding DNA Sequences). If an SNP is found to be located in a CDS, parent transcripts (mRNA) are called; otherwise, they are classified as "Non-coding."

The mature mRNA (MmRNA) sequence, composed only of CDSs, is reconstructed by querying all CDS for the previously called transcript and retrieving their sequences from a `Data.fasta` file. Thus, MmRNA sequences and SNP coordinates relative to MmRNA are retrieved.

The MmRNA sequence is translated to protein. If the translated protein is detected as a pseudogene, it is counted but not processed or annotated. A total count of the pseudogenes detected is reported in the `Output.log` file.

The alternative amino acid produced by the SNP in the protein is determined using the SNP location in the MmRNA. If the alternative amino acid is the same as the reference amino acid, the SNP is classified as "Synonymous"; otherwise, it is classified as "Non-synonymous."

Lastly, the annotations from each filtered SNP are written to a tab-separated `Output.tsv` file.

The generated `Output.tsv` file is read and used to produce an `Output.png` bar plot file illustrating the proportions of “Non-coding”, “Synonymous”, and “Non-synonymous” variants analyzed.

The locations of the `Output.log`, `Output.tsv`, and `Output.png` files are reported in the `Output.log` file.

# Usage

The program is implemented through a command line interface. Descriptions of the required input types are shown using the `-h` argument (Figure 1). The arguments are not positional, but all are required to run the program.


```bash
program.py [-h] --vcfFile VCFFILE --gffFile GFFFILE --fastaFile FASTAFILE --outputPrefix OUTPUTPREFIX

options:

- `-h`, `--help`  Show this help message and exit.

- `--vcfFile VCFFILE`  Input VCF file.

- `--gffFile GFFFILE`  Input GFF file.

- `--fastaFile FASTAFILE`  Input FASTA file.

- `--outputPrefix OUTPUTPREFIX`  Prefix for output files.
```

### Example
to correctly run the script, the following command line structure should be followed:
```bash
python program.py --vcfFile=Data.vcf.gz --gffFile=Data.gff    --fastaFile=Data.fasta --outputPrefix=1234567 
```

# API

## Sequence

### Class: `Sequence(sequence)`
This class creates a Sequence object, which is used to manipulate and extract information from a DNA sequence and avoid repetitive code in processing the files. Methods for this class and its subclasses use 1-indexed positions.

### Methods:
- `protein_seq()`  
  Returns a protein sequence based on the DNA sequence from the Sequence object.

- `protein_pos(MmRNA_pos)`  
  Returns the amino acid position in the protein sequence based on the SNP position relative to the MmRNA.

- `reference_aa(protein_pos)`  
  Returns the reference amino acid letter in the protein sequence based on an input position in the protein.

---

## FSequence

### Class: `FSequence(sequence)`
FSequence refers to Forward Sequence and represents sequences transcribed from a genomic sequence's forward strand (+). FSequence is a subclass that inherits from the Sequence class.

### Methods:
- `alt_aa(MmRNA_pos, alt_base)`  
  Returns the alternative amino acid letter in the protein sequence based on the input position in the MmRNA sequence and the alternative base to be overridden in the DNA sequence.

---

## RSequence

### Class: `RSequence(sequence)`
RSequence refers to Reverse Sequence and represents sequences transcribed from a genomic sequence's reverse strand (-). RSequence is a subclass that inherits from the Sequence class.

### Methods:
- `alt_aa(MmRNA_pos, alt_base)`  
  Same logic as in FSequence, but is used primarily for reverse strand encoded sequences because `alt_base` needs to be transformed to its complement base first.

# Input Files

## Data.vcf.gz
A compressed Variant Calling Format file. The file has the following fields (in order):
- **CHROM:** Chromosome name.
- **POS:** Position of the variant on the chromosome.
- **ID:** ID for the variant.
- **REF:** Reference base.
- **ALT:** Alternative base.
- **QUAL:** Quality score.
- **FILTER:** Filter status.
- **INFO:** Additional information about the variant.
- **FORMAT:** Format of the genotype information.
- **unknown:** Genotype information for a sample.

---

## Data.gff
A General Feature Format file. The file has the following fields (in order):
- **CHROM:** Chromosome name.
- **Source:** The source of the feature.
- **Feature Type:** The biological type of the feature (e.g., "CDS", "mRNA", "exon", etc.).
- **Start:** The starting position of the feature.
- **End:** The ending position of the feature.
- **Score:** A value representing the quality of the feature.
- **Strand:** The orientation of the feature on the sequence, “+” or “-”.
- **Frame:** If associated with CDS, the reading frame position.
- **Attributes:** Key-value pairs with additional information about the feature.

---

## Data.fasta
A FASTA file containing a single chromosome sequence.

---

## outputPrefix
An input string used to name all the output files.

---

# Output Files

## Output.tsv
A tab-separated file which has the following fields (in order):
- **Chrom:** Chromosome name.
- **Pos:** Position of the SNP on the chromosome.
- **Ref:** Reference base.
- **Alt:** Alternative base.
- **Transcript:** Transcript name/id.
- **Protein Location:** Position of the amino acid affected by the SNP.
- **Ref AA:** Reference amino acid.
- **Alt AA:** Alternative amino acid.

---

## Output.log
A log file which contains:
- Confirmation of the filenames given at the command line.
- The total count of SNPs which QUAL <= 20.
- The total count of SNPs which produce a pseudogene.
- The locations of the output files.

---

## Output.png
A bar plot file illustrating the proportions of “Non-coding”, “Synonymous”, and “Non-synonymous” variants analyzed.

# Issues and Observations

- There was a warning coming from the BioPython `Seq.translate()` method. Program performance was not affected, but it was hidden using the `warnings` and `BiopythonWarning` modules due to aesthetic preferences.

- Currently, pseudogenes are not addressed in any way other than counting them. They could be handled by obtaining a list and inputting their features into a TSV file.

- Some SNPs may cause a sequence change, producing a premature stop codon. These changes were overlooked because of how proteins were filtered to catch pseudogenes. It could be modified by changing `protein_seq.count('*')==1` to `protein_seq.count('*') >= 1`.

- Although a `Data.vcf.gz.tbi` file was provided, it was not used, and it is unknown if using this file as an input would produce the same results as using the `Data.vcf.gz` file.
