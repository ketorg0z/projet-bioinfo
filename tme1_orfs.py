import random
import re

CODEGENETIQUE = {
    "TTT": "F", "TTC": "F","TTA": "L","TTG": "L","TCT": "S","TCC": "S","TCA": "S","TCG": "S","TAT": "Y","TAC": "Y",
    "TAA": "*","TAG": "*","TGT": "C","TGC": "C","TGA": "*","TGG": "W","CTT": "L","CTC": "L","CTA": "L","CTG": "L",
    "CCT": "P","CCC": "P","CCA": "P","CCG": "P","CAT": "H","CAC": "H","CAA": "Q","CAG": "Q","CGT": "R","CGC": "R",
    "CGA": "R","CGG": "R","ATT": "I","ATC": "I","ATA": "I","ATG": "M","ACT": "T","ACC": "T","ACA": "T","ACG": "T",
    "AAT": "N","AAC": "N","AAA": "K","AAG": "K","AGT": "S","AGC": "S","AGA": "R","AGG": "R","GTT": "V","GTC": "V",
    "GTA": "V","GTG": "V","GCT": "A","GCC": "A","GCA": "A","GCG": "A","GAT": "D","GAC": "D","GAA": "E","GAG": "E",
    "GGT": "G","GGC": "G","GGA": "G","GGG": "G"
}

#Question 1
def remplace_non_identifies(seq: str):
  """
  Remplace les nucléotides non identifiés par une des possibilités de façon aléatoire.
  entrée seq: sequence ayant peut-etre des nucleotides non identifiés
  sortie    : sequence (nettoyé) sans nucléotides non identifiés

  """
  options = {
      "R": ["G", "A"],
      "Y": ["T", "C"],
      "K": ["G", "T"],
      "M": ["A", "C"],
      "S": ["G", "C"],
      "W": ["A", "T"],
      "B": ["G", "T", "C"],
      "D": ["G", "A", "T"],
      "H": ["A", "C", "T"],
      "V": ["G", "C", "A"],
      "N": ["A", "G", "C", "T"],
      "X": ["A", "G", "C", "T"],
      # avec cette solution, le temps d'exécution pourrait augmenter par rapport aux autres solutions possibles,
      # mais elle prend moins de lignes :)
      "A": ["A"],
      "T": ["T"],
      "G": ["G"],
      "C": ["C"]
  }

  return "".join([random.choice(options[s]) for s in seq])

# Question 2
def listecodon(seq:str):
    """
    Renvoie une liste de codons pour une séquence passée en paramètre.
    entrée seq : sequence de nucléotides
    sortie     : list de codons de la sequence d'entrée

    Si la longueur de la séquence n'est pas un multiple de 3 elle ne tiendra pas
    compte des 1 ou 2 nucléotides restant à la fin.

    >>> listecodon('AAACCC')
    ['AAA', 'CCC']
    >>> listecodon('AAACC')
    ['AAA']
    >>> listecodon('AAAC')
    ['AAA']
    """
    return [] if len(seq) < 3 else [seq[:3]] + listecodon(seq[3:])


# Question 3
def reversecompl(seq:str):
    """Renvoie le brin complémentaire d’une séquence.
    entrée seq : sequence de nucléotides (brin sens)
    sortie     : sequence de nucléotides (brin complementaire)
    >>> reversecompl('AACGTGGCA')
    'TGCCACGTT'
    """
    compl = {'A': 'T', 'C': 'G', 'G': 'C', 'T':'A'}
    return "".join(reversed([compl[i] for i in seq]))


def trouver_positions_orfs(codons:list):
    """Retourne les positions de cadre ouverts de la liste des codons.
    entrée codons: liste de codons
    sortie       : liste de positions [start] et [stop] de cadre ouverts

    >>> trouver_positions_orfs(['ATG', 'AAA', 'ATG', 'ATT', 'TAG', 'ATG', 'ACC', 'ATT', 'ACC', 'ACC', 'ACC', 'ATC', 'ACC', 'ATT', 'ACC', 'ACA', 'GGT', 'AAC', 'TGA', 'GGT', 'GCG', 'GGC'])
    ([0, 5], [4, 18])
    """
    starts = [i for i, codon in enumerate(codons) if codon in {'ATG', 'GTG', 'TTG'}] #positions de tous les starts dans codons
    stops =  [i for i, codon in enumerate(codons) if codon in {'TAA', 'TAG', 'TGA'}] #positions de tous les stops  dans codons

    orf_starts = [] #renvoie la liste de positions de codons start
    orf_stops = []  #renvoie la liste de positions de codons stop, les deux liste ont la meme taille et sont de pair start/stop

    lastStart = None

    for i in range(len(codons)):
      if i in stops:
        if lastStart != None:
          # ORF detected
          orf_starts.append(lastStart)
          orf_stops.append(i)
          lastStart = None
      elif i in starts and lastStart == None:
        lastStart = i

    return (orf_starts, orf_stops)


# Question 4
def liste_orfs_sens(seq:str):
    """
    Liste tous les cadres ouverts de lectures du brin sens.
    entrée : sequence de nucléotides
    sortie : liste contenant tous les cadres ouverts de lectures
    >>> sorted(liste_orfs_sens('AAATGATGTAATAGTGTTTTGATTAGGGCAT'))
    ['ATGATGTAA', 'GTGTTTTGA']
    """
    orfs = []
    for i in range(3):
      codons = listecodon(seq[i:])
      pos_orfs = trouver_positions_orfs(codons)
      for j in range(len(pos_orfs[0])):
        orfs.append(seq[
            i + pos_orfs[0][j] * 3 :
            i + (pos_orfs[1][j] + 1) * 3]
        )

    return orfs

def __liste_orfs_sens(seq):
    """Retourne la liste des cadres ouverts de lecture, sens 5' vers 3'.
    entrée : sequence de nucléotides
    sortie : liste contenant tous les cadres ouverts de lectures
    """
    return [
        match.group() for match in re.finditer(
            '((ATG|GTG|TTG)((?!TAA|TAG|TGA)([ACTG][ACTG][ACTG]))*(TAA|TAG|TGA))',
            seq)
    ]

def __liste_orfs(seq):
    """
    Retourne la liste de tous les cadres ouverts de lectures.
    entrée : sequence de nucléotides
    sortie : liste des cadres ouverts de lecture

    >>> sorted(__liste_orfs_sens('AAATGATGTAATAGTGTTTTGATTAGGGCAT'))
    ['ATGATGTAA', 'GTGTTTTGA']
    """
    liste = __liste_orfs_sens(seq)
    liste.extend(__liste_orfs_sens(reversecompl(seq)))
    return liste