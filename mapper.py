class Sequence:
    def __init__(self, lines):
        self.name = lines[0].strip()[1:]
        self.bases = "".join([x.strip() for x in lines[1:]]).upper()

    def __str__(self):
        return self.name + ": " + self.bases[:20] + "..."

    def __repr__(self):
        return self.__str__()


class Read(Sequence):
    def get_seed(self, seedlength):
        return self.bases[:seedlength]

    def replace_kmers(self, replacements):
        for kmer, replacement in replacements.items():
            self.bases = self.bases.replace(kmer, replacement)


class Reference(Sequence):
    def __init__(self, lines):
        self.kmers = None
        super().__init__(lines)

    def calculate_kmers(self, kmersize):
        self.kmers = {}
        for pos in range(0, len(self.bases) - kmersize + 1):
            kmer = self.bases[pos:(pos + kmersize)]
            if kmer not in self.kmers:
                self.kmers[kmer] = []
            self.kmers[kmer] += [pos]

    def get_kmer_positions(self, kmer):
        if self.kmers is None or len(next(iter(self.kmers))) != len(kmer):
            self.calculate_kmers(len(kmer))
        if kmer not in self.kmers:
            return []
        return self.kmers[kmer]

    def count_mismatches(self, read, position):
        mismatches = 0
        for pos in range(position, position + len(read.bases)):
            if pos >= len(self.bases):
                break
            if read.bases[pos - position] != self.bases[pos]:
                mismatches += 1
        # Count every base of the read that goes out past the end of the reference as a mismatch
        mismatches += position + len(read.bases) - pos - 1
        return mismatches


class Mapping:
    def __init__(self, reference):
        self.reference = reference
        self.reads = {}

    def add_read(self, read, position):
        if position not in self.reads:
            self.reads[position] = []
        self.reads[position] += [read]

    def get_reads_at_position(self, position):
        if position not in self.reads:
            return []
        return self.reads[position]

    def __str__(self):
        res = ["Mapping to " + self.reference.name]
        for pos in self.reads:
            res += ["  " + str(len(self.reads[pos])) + " reads mapping at " + str(pos)]
        return "\n".join(res)


class SAMWriter:
    def __init__(self, mapping):
        self.mapping = mapping

    def write_mapping(self, filename):
        myfile = open(filename, "w")
        refname = self.mapping.reference.name.split(" ")[0]
        myfile.write("@SQ\tSN:" + refname + "\tLN:" + str(len(self.mapping.reference.bases)) + "\n")
        for pos in range(0, len(self.mapping.reference.bases)):
            for read in self.mapping.get_reads_at_position(pos):
                myfile.write("\t".join([read.name, "0", refname, str(pos + 1), "255",
                                        str(len(read.bases)) + "M", "*", "0", "0", read.bases, "*"]))
                myfile.write("\n")
        myfile.close()


class ReadPolisher:
    def __init__(self, kmerlen):
        self.kmer_length = kmerlen
        self.spectrum = {}

    def add_read(self, readseq):
        kmers = []
        k = self.kmer_length
        for i in range(len(readseq) - k + 1):
            kmer = readseq[i:i + k]
            kmers.append(kmer)
        self.create_spectrum(kmers)

    def create_spectrum(self, kmers):
        for kmer in kmers:
            if kmer in self.spectrum:
                self.spectrum[kmer] += 1
            else:
                self.spectrum[kmer] = 1

    def get_replacements(self, minfreq):
        corrections = {}
        for kmer, count in self.spectrum.items():
            if count < minfreq:
                candidate_kmers = []
                for i, base in enumerate(kmer):
                    for new_base in ['A', 'G', 'T', 'C']:
                        candidate_kmer = kmer[:i] + new_base + kmer[i + 1:]
                        if candidate_kmer in self.spectrum and self.spectrum[candidate_kmer] >= minfreq:
                            candidate_kmers.append(candidate_kmer)
                if candidate_kmers:
                    most_frequent_candidate = max(candidate_kmers, key=lambda x: self.spectrum[x])
                    corrections[kmer] = most_frequent_candidate
        return corrections


def read_fasta(fastafile, klassname):
    klass = globals()[klassname]
    f = open(fastafile, "r")
    readlines = []
    reads = []
    for line in f:
        if line[0] == '>' and len(readlines) != 0:
            reads += [klass(readlines)]
            readlines = []
        readlines += [line]
    reads += [klass(readlines)]
    f.close()
    return reads


def map_reads(reads, reference, kmersize, max_mismatches):
    mapping = Mapping(reference)
    reference.calculate_kmers(kmersize)
    for read in reads:
        seed = read.get_seed(kmersize)
        seed_positions = reference.get_kmer_positions(seed)
        for position in seed_positions:
            mismatches = reference.count_mismatches(read, position)
            if mismatches < max_mismatches:
                mapping.add_read(read, position)
    return mapping


def main():
    # ---------- Tablet-Aufgabe ---------
    # Mappen wir nun die Datei (z.B. data/fluA_reads.fasta auf data/fluA_reads.fasta )
    # und speichern wir das Ergebnis als z.B. fluA_mapping.sam.
    # reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
    # reference = read_fasta("data/fluA.fasta", Reference.__name__)[0]
    # mapping = map_reads(reads, reference, 8, 5)
    # print("Mapping reads: " + len(mapping.reads))
    # writer = SAMWriter(mapping)
    # writer.write_mapping("data/fluA_mapping.sam")

    # ---------- Antibiotika-Resistenzen ---------
    # Mappen wir die Read-Sequenzen der 4 Personen (1 bis 4):
    # reads = read_fasta("data/patient4.fasta", Read.__name__)
    # reference = read_fasta("data/rpoB.fasta", Reference.__name__)[0]
    # mapping = map_reads(reads, reference, 11, 5)
    # writer = SAMWriter(mapping)
    # writer.write_mapping("data/patient4.sam")

    # ---------- Anwendung ---------
    # Verwenden wir nun unsere ReadPolisher:
    reads = read_fasta("data/patient4.fasta", Read.__name__)
    reference = read_fasta("data/rpoB.fasta", Reference.__name__)[0]
    polisher = ReadPolisher(15)
    for read in reads:
        polisher.add_read(read.bases)
    replacements = polisher.get_replacements(3)
    nrep = 0
    for read in reads:
        nrep += 1
        if nrep % 1000 == 0:
            print(str(nrep) + "/" + str(len(reads)))
        read.replace_kmers(replacements)
    mapping = map_reads(reads, reference, 15, 3)
    writer = SAMWriter(mapping)
    writer.write_mapping("data/mapping_p4_corrected.sam")


if __name__ == "__main__":
    main()
