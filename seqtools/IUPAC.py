
"""Also see in Biopython: Bio.Alphabet.IUPAC"""


gaps = '-.'
nucleotides = 'ACGT'
unknown = 'N'
ambiguous = {'R': set('AG'),
             'Y': set('CT'),
             'S': set('GC'),
             'W': set('AT'),
             'K': set('GT'),
             'M': set('AC'),
             'B': set('CGT'),
             'D': set('AGT'),
             'H': set('ACT'),
             'V': set('ACG')}

aa = 'ABCDEFGHIKLMNPQRSTVWYZ'
unknown_aa = 'X'
stop_codons = set(('TGA', 'TAG', 'TAA'))  # In standard alphabet
#aa = {
#      'A':        Ala        Alanine
#      'B':        Asx        Aspartic acid or Asparagine
#      'C':        Cys        Cysteine
#      'D':        Asp        Aspartic Acid
#      'E':        Glu        Glutamic Acid
#      'F':        Phe        Phenylalanine
#      'G':        Gly        Glycine
#      'H':        His        Histidine
#      'I':        Ile        Isoleucine
#      'K':        Lys        Lysine
#      'L':        Leu        Leucine
#      'M':        Met        Methionine
#      'N':        Asn        Asparagine
#      'P':        Pro        Proline
#      'Q':        Gln        Glutamine
#      'R':        Arg        Arginine
#      'S':        Ser        Serine
#      'T':        Thr        Threonine
#      'V':        Val        Valine
#      'W':        Trp        Tryptophan
#      'X':        Xaa        Any amino acid
#      'Y':        Tyr        Tyrosine
#      'Z':        Glx        Glutamine or Glutamic acid
#      }
