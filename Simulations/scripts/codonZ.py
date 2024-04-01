std_table = {
    'TTT':'F','TCT':'S','TAT':'Y','TGT':'C',
    'TTC':'F','TCC':'S','TAC':'Y','TGC':'C',
    'TTA':'L','TCA':'S','TAA':'X','TGA':'X',
    'TTG':'L','TCG':'S','TAG':'X','TGG':'W',

    'CTT':'L','CCT':'P','CAT':'H','CGT':'R',
    'CTC':'L','CCC':'P','CAC':'H','CGC':'R',
    'CTA':'L','CCA':'P','CAA':'Q','CGA':'R',
    'CTG':'L','CCG':'P','CAG':'Q','CGG':'R',

    'ATT':'I','ACT':'T','AAT':'N','AGT':'S',
    'ATC':'I','ACC':'T','AAC':'N','AGC':'S',
    'ATA':'I','ACA':'T','AAA':'K','AGA':'R',
    'ATG':'M','ACG':'T','AAG':'K','AGG':'R',

    'GTT':'V','GCT':'A','GAT':'D','GGT':'G',
    'GTC':'V','GCC':'A','GAC':'D','GGC':'G',
    'GTA':'V','GCA':'A','GAA':'E','GGA':'G',
    'GTG':'V','GCG':'A','GAG':'E','GGG':'G',
}


# codon list in PAML mlc file
codonListY = [
    'TTT','TTC','TTA','TTG',
    'TCT','TCC','TCA','TCG',
    'TAT','TAC','TAA','TAG',
    'TGT','TGC','TGA','TGG',

    'CTT','CTC','CTA','CTG',
    'CCT','CCC','CCA','CCG',
    'CAT','CAC','CAA','CAG',
    'CGT','CGC','CGA','CGG',

    'ATT','ATC','ATA','ATG',
    'ACT','ACC','ACA','ACG',
    'AAT','AAC','AAA','AAG',
    'AGT','AGC','AGA','AGG',

    'GTT','GTC','GTA','GTG',
    'GCT','GCC','GCA','GCG',
    'GAT','GAC','GAA','GAG',
    'GGT','GGC','GGA','GGG',
]

codonListYwoTcodon = [
    'TTT','TTC','TTA','TTG',
    'TCT','TCC','TCA','TCG',
    'TAT','TAC',
    'TGT','TGC','TGG',

    'CTT','CTC','CTA','CTG',
    'CCT','CCC','CCA','CCG',
    'CAT','CAC','CAA','CAG',
    'CGT','CGC','CGA','CGG',

    'ATT','ATC','ATA','ATG',
    'ACT','ACC','ACA','ACG',
    'AAT','AAC','AAA','AAG',
    'AGT','AGC','AGA','AGG',

    'GTT','GTC','GTA','GTG',
    'GCT','GCC','GCA','GCG',
    'GAT','GAC','GAA','GAG',
    'GGT','GGC','GGA','GGG',
]

std_index_table = {
    0: 'TTT',  1: 'TTC',  2: 'TTA',  3: 'TTG',
    4: 'TCT',  5: 'TCC',  6: 'TCA',  7: 'TCG',
    8: 'TAT',  9: 'TAC',
    10: 'TGT', 11: 'TGC', 12: 'TGA', 13: 'TGG',
    14: 'CTT', 15: 'CTC', 16: 'CTA', 17: 'CTG',
    18: 'CCT', 19: 'CCC', 20: 'CCA', 21: 'CCG',
    22: 'CAT', 23: 'CAC', 24: 'CAA', 25: 'CAG',
    26: 'CGT', 27: 'CGC', 28: 'CGA', 29: 'CGG',
    30: 'ATT', 31: 'ATC', 32: 'ATA', 33: 'ATG',
    34: 'ACT', 35: 'ACC', 36: 'ACA', 37: 'ACG',
    38: 'AAT', 39: 'AAC', 40: 'AAA', 41: 'AAG',
    42: 'AGT', 43: 'AGC', 44: 'AGA', 45: 'AGG',
    46: 'GTT', 47: 'GTC', 48: 'GTA', 49: 'GTG',
    50: 'GCT', 51: 'GCC', 52: 'GCA', 53: 'GCG',
    54: 'GAT', 55: 'GAC', 56: 'GAA', 57: 'GAG',
    58: 'GGT', 59: 'GGC', 60: 'GGA', 61: 'GGG'
}

transversion = ['AC','CA','AT','TA','GC','CG','CT','TC']
transition = ['AG','GA','CT','TC']

AA_std_table = {
    'A': 0, 'R': 1, 'N': 2, 'D': 3,
    'C': 4, 'Q': 5, 'E': 6, 'G': 7,
    'H': 8, 'I': 9, 'L': 10, 'K': 11,
    'M': 12, 'F': 13, 'P': 14, 'S': 15,
    'T': 16, 'W': 17, 'Y': 18, 'V': 19,
    'X': 20, 'B': 21, 'Z': 22, 'O': 23,
    'U': 24
}

def reverse_AA_table():
    AA_reverse_table = {}
    for key, value in AA_std_table.items():
        AA_reverse_table[value] = key
    return AA_reverse_table