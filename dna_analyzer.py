"""
Enhanced DNA Sequence Analyzer
Analyzes DNA sequences for nucleotide composition, GC content,
transcription, translation, and motif searching.
Supports FASTA file input.
"""

CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met (Start)",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


def read_fasta(file_path):
    """Read a FASTA file and return the DNA sequence as a string."""
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):  # skip header lines
                sequence += line.strip()
    return sequence


def validate_sequence(seq):
    seq = seq.upper().strip()
    invalid = set(seq) - {"A", "T", "G", "C"}
    if invalid:
        return False, f"Invalid bases found: {', '.join(sorted(invalid))}"
    if len(seq) == 0:
        return False, "Sequence cannot be empty."
    return True, seq


def count_nucleotides(seq):
    return {base: seq.count(base) for base in "ATGC"}


def gc_content(seq):
    return ((seq.count("G") + seq.count("C")) / len(seq)) * 100


def at_content(seq):
    return ((seq.count("A") + seq.count("T")) / len(seq)) * 100


def dna_to_rna(seq):
    return seq.replace("T", "U")


def reverse_complement(seq):
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement[base] for base in reversed(seq))


def translate_rna(rna):
    codons = [rna[i:i+3] for i in range(0, len(rna) - 2, 3)]
    return [CODON_TABLE.get(codon, "???") for codon in codons]


def find_motif(seq, motif):
    positions = []
    start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)
        start = pos + 1
    return positions


def find_orfs(seq):
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    for i in range(len(seq) - 2):
        if seq[i:i+3] == "ATG":
            for j in range(i + 3, len(seq) - 2, 3):
                codon = seq[j:j+3]
                if codon in stop_codons:
                    orfs.append((i + 1, j + 3, seq[i:j+3]))
                    break
    return orfs


def melting_temperature(seq):
    a, t, g, c = seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")
    if len(seq) < 14:
        return 2 * (a + t) + 4 * (g + c)
    return 64.9 + 41 * (g + c - 16.4) / len(seq)


def display_bar(label, value, total, bar_width=20):
    filled = int((value / total) * bar_width)
    bar = "█" * filled + "░" * (bar_width - filled)
    print(f"  {label}  {bar}  {int(value):>4} ({value/total*100:5.1f}%)")


def print_section(title):
    print(f"\n{'─' * 52}")
    print(f"  {title}")
    print(f"{'─' * 52}")


def analyze(dna_sequence):
    print("\n" + "═" * 54)
    print("         🧬  DNA SEQUENCE ANALYZER  🧬")
    print("═" * 54)

    valid, result = validate_sequence(dna_sequence)
    if not valid:
        print(f"\n  ❌ Error: {result}")
        return
    seq = result

    print_section("SEQUENCE INFO")
    print(f"  Sequence : {seq[:60]}{'...' if len(seq) > 60 else ''}")
    print(f"  Length   : {len(seq)} bp")

    print_section("NUCLEOTIDE COMPOSITION")
    counts = count_nucleotides(seq)
    for base in "ATGC":
        display_bar(base, counts[base], len(seq))

    print_section("CONTENT METRICS")
    print(f"  GC Content : {gc_content(seq):.2f}%")
    print(f"  AT Content : {at_content(seq):.2f}%")
    print(f"  Est. Tm    : {melting_temperature(seq):.1f} °C  (Wallace rule)")

    print_section("TRANSCRIPTION  (DNA → RNA)")
    rna = dna_to_rna(seq)
    print(f"  RNA : {rna[:60]}{'...' if len(rna) > 60 else ''}")

    print_section("TRANSLATION  (RNA → Protein)")
    print(f"  Amino Acids : {' - '.join(translate_rna(rna))}")

    print_section("REVERSE COMPLEMENT")
    rev = reverse_complement(seq)
    print(f"  5' {seq[:50]}{'...' if len(seq) > 50 else ''} 3'")
    print(f"  3' {rev[:50]}{'...' if len(rev) > 50 else ''} 5'")

    print_section("OPEN READING FRAMES (ORFs)")
    orfs = find_orfs(seq)
    if orfs:
        for start, end, orf_seq in orfs:
            print(f"  Position {start}–{end}: {orf_seq[:40]}{'...' if len(orf_seq) > 40 else ''}  ({len(orf_seq)} bp)")
    else:
        print("  No complete ORFs found.")

    print_section("MOTIF SEARCH")
    for motif in ["ATG", "GCT", "CGA"]:
        positions = find_motif(seq, motif)
        if positions:
            print(f"  {motif} found at position(s): {positions[:10]}{'...' if len(positions) > 10 else ''}")
        else:
            print(f"  {motif} : not found")

    print("\n" + "═" * 54 + "\n")


if __name__ == "__main__":
    dna_sequence = read_fasta("fshb_sequence.fasta")
    analyze(dna_sequence)
