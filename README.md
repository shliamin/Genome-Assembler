

## SAM Files and Visualization of Mappings

For the list of mapping reads that we output in the last task, there is also a standard format in bioinformatics: The Sequence Alignment/Map Format, or SAM. It contains information about the reference sequence used in a mapping as well as the mapping reads. You can find the complete format specification with all possible additional information at the given link; we use only a minimal version. An example SAM file can be found at data/minimapping.sam.

The SAM format is line-based, with each line containing one or more tab-separated entries.

The SAM file begins with a header section, where each line starts with an @. The @SQ header contains information about the reference sequence used:

```text
@SQ	SN:GQ359764.1	LN:2445
```

The fields ```SN:``` and ```LN:``` indicate the name of the reference sequence (must be exactly the name from the reference FASTA file up to the first space) and its length.

Then, line by line, information about the mapped reads follows. Each line must contain at least 11 columns:

| Column | Field Name |N/A Value     | Description                                                     |
|--------|------------|--------------|-----------------------------------------------------------------|
| 1      | QNAME      | required     | Name of the read                                                |
| 2      | FLAG       | required     | Bitwise flag. 0 for "maps correctly"                            |
| 3      | RNAME      | *            | Name of the reference sequence (identical to SN: in @SQ header) |
| 4      | POS        | 0            | Position in the reference (first base is 1, not 0)              |
| 5      | MAPQ       | 255          | Mapping quality                                                 |
| 6      | CIGAR      | *            | CIGAR string of the alignment                                   |
| 7      | RNEXT      | *            | Reference where the next part of the read maps                  |
| 8      | PNEXT      | 0            | Position of the next part of the read                           |
| 9      | TLEN       | 0            | Insert length                                                   |
| 10     | SEQ        | *            | Read sequence                                                   |
| 11     | QUAL       | *            | Base quality of the reads (as in FASTQ file)                    |

For meaningful visualization, we provide information in the columns QNAME, FLAG, RNAME, POS, CIGAR, and SEQ.

The other columns are not important for us currently: The columns RNEXT and PNEXT are relevant for cases where the read is split into multiple parts during analysis and the different parts map to different locations on the reference or different references (e.g., to find structural variations). TLEN is relevant for paired-end sequencing. In QUAL, a mapper can indicate how confident they are in assigning the read to that exact location on the reference.

The relevant columns should be self-explanatory. The only exception is the CIGAR string: This is a condensed representation of how a read maps. It consists of a series of base numbers followed by status information. The status information can include, for example, M for "match" (base could be assigned to the reference sequence) or D (deletion - base is present in the reference sequence but missing in the read). 25M2D10M would thus mean: the first 25 bases fit, then 2 bases are deleted, and then 10 bases fit again. Since we do not consider deletions, the CIGAR string in our case is always <Read length>M.

The following line from the example SAM file:

```text
Read_95	0	GQ359764.1	10	255	50M	*	0	0	TCCATGGTGTATCCTGTTCCTGTTCCATGGCTGTATGGAGGATCTCCAGT	*
```

means: The read named "Read_95" (QNAME=Read_95) maps (FLAG=0) to the reference sequence "GQ359764.1" (RNAME=GQ359764.1) at position 10 (POS=10) without deletions or insertions (CIGAR=50M), and its sequence reads "TCCATGGTGTATCCTGTTCCTGTTCCATGGCTGTATGGAGGATCTCCAGT". The other fields are filled with their N/A values.

### Implementation of SAM Export

Implement a class ```SAMWriter``` with the following methods:
* ```__init__(self, mapping)```: Constructor, takes a ```Mapping```-object
* ```write_mapping(self, filename)```: Writes the mapping to the specified SAM file

### Visualization Using Tablet

Now map the file data/fluA_reads.fasta to data/fluA.fasta and save the result as fluA_mapping.sam.

Then download the program  [Tablet](https://ics.hutton.ac.uk/tablet/) and open (click the "Open Assembly" button in the top left) the files fluA_mapping.sam and the reference [data/fluA_reads.fasta](data/fluA.fasta). You should then receive a view like in this picture:

![t1](Bilder/Tablet1.png)

You see a schematic view of the entire reference with the mapped reads on it, below we see a detail view: First, the sequence translated into amino acids, then the nucleotide sequence of the reference, and then the individual reads.

If we hover the mouse pointer over a base, we get information about the read. The coordinate axis between the reference sequence and the reads also shows the position in red - in this case, it can be seen that the first base "T" from the read "Read_95" was mapped at position 10 of the reference sequence.

If we select the "Variants" option in the "Color Schemes" tab, all bases of the reads that match the reference sequence are grayed out:

![t2](Bilder/Tablet2.png)

You now recognize two red stripes at the top indicating differences from the reference sequence. The highlights in the overview are unfortunately not comprehensive. This means we must scroll through the entire mapping to find all variants. If we have selected the "Variants" scheme, the deviations will also be highlighted in red in the read view.

### Tablet Task

We can enter now in the format ```<Reference Base><Position><New Base>``` which four mutations we can recognize in the mapping (```T10A``` would mean, for example, that the base T stands in the reference at position 10, but according to the reads there is a mutation to A at this position):

```text
Mutation 1: G960C
Mutation 2: T1401G
Mutation 3: G1437T
Mutation 4: G1833T
```

## Antibiotic Resistances

One task where recognizing such mutations is particularly important is the treatment of bacterial infections. Bacteria can develop resistances to antibiotics - the administration of such antibiotics can then no longer contribute to healing. However, there are now many well-studied relationships between mutations in certain genes and the resulting antibiotic resistances. Accordingly, before treatment, a bacterium can be sequenced, and a treatment decision can be made based on the existing mutations.

A particularly prominent example is the bacterium Staphylococcus aureus, which quickly accumulates antibiotic resistances. Infections with multi-resistant S. aureus (MRSA) present a major challenge to medicine, as in the worst case, none of the available antibiotics may work against them (or only so-called "drugs of last resort" work - antibiotics that are held back for particularly severe cases, as bacteria have not yet been subjected to evolutionary pressure to develop resistances against these).

In this task, we will examine the samples of 4 individuals infected with S. aureus for mutations in the rpoB gene of the bacterium. The following two antibiotics are available for treatment - in brackets is the priority with which they should be used, if possible the antibiotic with the highest priority (the smallest number behind it) should be used:

* Daptomycin (1)
* Rifampicin (2)

The following three mutations are also known to convey resistances:

* C1862A: Resistance against Daptomycin
* T2858G: Resistance against Daptomycin
* C1402A: Resistance against Rifampicin

Map the read sequences of the 4 individuals ([data/patient1.fasta](data/patient1.fasta) - [data/patient4.fasta](data/patient4.fasta)) to the rpoB reference ([data/rpoB.fasta](data/rpoB.fasta)) and enter here which mutation(s) we could identify and which antibiotic we would recommend:

```text
Person 1 - Mutation(en): no mutations, Recommendation: Daptomycin 
Person 2 - Mutation(en): no mutations, Recommendation: Daptomycin 
Person 3 - Mutation(en): C1402A (many), Recommendation: Daptomycin 
Person 4 - Mutation(en): T2858G (one), Recommendation: Rifampicin 
```

Use a seed length of > 10 for the mapping. Please do not be confused by differences from the reference sequence that only occur in individual reads - this is a realistic dataset and the reads contain sequencing errors.

## Error Correction

As we may have noticed in the identification of antibiotic resistances, the real reads are affected by sequencing errors. These can complicate the analysis or even lead to misinterpretations of the data.

One way to correct these errors is the k-mer spectrum. It is assumed that because there is a large coverage with mostly correct reads, every sequenced k-mer should be represented several times in different reads. If a k-mer appears significantly less often than the other k-mers, it is probably not due to a mutation (which should be covered by several reads and whose k-mer should therefore appear several times) but due to a sequencing error.

Let's take an example of a short genome, which is sequenced error-free, and the 3-mer spectrum for it:

![sp1](Bilder/Spectrum1.png)

In this case, the genome was covered with 5 error-free reads, resulting in 4 3-mers with the frequencies 2, 4, 4, and 2.

However, if one of the reads contains an error (in this case, the 3rd base of Read 2 is erroneously read as C), the spectrum changes:

![sp2](Bilder/Spectrum2.png)


The error adds three new k-mers, each appearing only once.

Based on this information, a correction can be made: A threshold is defined, below which a k-mer is considered potentially faulty. For each k-mer that occurs less often than this threshold, the following steps are run through:

For each base X in the k-mer:
For each possible base Y (A, G, T, and C):
Generate a candidate k-mer by replacing the base X with the base Y
If the candidate k-mer also occurs in the dataset and at least as often as the threshold: Note it as a possible correction
If candidate k-mers were found: Replace the k-mer with the candidate k-mer that occurs most frequently in the dataset (in the case of two candidate k-mers with the same frequency, choose one at random)


For the mapping, the k-mers identified as correctable in all reads in which they occur must be replaced.

### Implementation

Implement the k-mer spectrum error correction as follows.

First, implement a class ```ReadPolisher``` with the following methods:
* ```__init__(self, kmerlen)```: Constructor, receives the k-mer length to be used
* ```add_read(self, readseq)```: Adds the provided read sequence to the k-mer spectrum
* ```get_replacements(self, minfreq)```: Calculates possible corrections for the k-mers that occur less often than ```minfreq``` in the k-mer spectrum and returns a corresponding dictionary. In it, keys are the correctable k-mers, values are the corrections (in the above example, a possible key-value pair might be "GCT":"GTT", which indicates that the k-mer "GCT" should be replaced with the k-mer "GTT")
  
Also, expand the ```Read``` class with the method ```replace_kmers(self, replacements)```, which receives the dictionary from ```get_replacements``` and replaces all k-mers occurring in the read that are listed as keys with their respective value. This is not the most efficient variant (more efficient would be to remember in ```ReadPolisher`` which reads contain which k-mers and then only make replacements there), but that would make the final task too long.

### Application

Use wer read correction to remap the read sequences of the 4 individuals ([data/patient1.fasta](data/patient1.fasta) - [data/patient4.fasta](data/patient4.fasta)) to the rpoB reference ([data/rpoB.fasta](data/rpoB.fasta)) using a seed length > 10 again. Do we see a difference? Which k-mer lengths and cutoffs seem sensible for the correction?

Enter here which mutation(s) we could identify and which antibiotic we would now recommend (use the parameters we find most sensible, but try at least once with [data/patient2.fasta](data/patient2.fasta) a k-mer length of 15 and a frequency cutoff of 3):

```text
Person 1 - Mutation(s): no mutations, Recommendation: Daptomycin 
Person 2 - Mutation(s): C1862A (many), Recommendation: Rifampicin 
Person 3 - Mutation(s): C1402A (many), Recommendation: Daptomycin 
Person 4 - Mutation(s): no mutations, Recommendation: Daptomycin 
```

Do we see a difference in the recommendations compared to those we gave without error correction? Briefly describe what the difference is and how it came about through error correction (no novel, 5-6 sentences are enough):

```text
The difference is that in the "corrected" version of the responses, more information is included. For example, we did not know before the corrected version of the code that patient number 2 had a resistance to the antibiotic Daptomycin, but now we do and can give more precise recommendations. Additionally, we probably gave patient number 4 an incorrect recommendation before the corrected version of the code, as in the corrected version of the code no mutations were found in this person.
```
