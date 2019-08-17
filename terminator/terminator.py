from math import ceil

def ComplementarySeq(seq):
    """makes & returns a complementary sequence to the entered sequence"""
    
    nucleobases = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complementary_seq = ""
    for base in seq:
        complementary_seq += nucleobases[base]
        
    return complementary_seq


def GCRichRegionsSearch(seq):
    """searches guanine- and cytosine-rich regions in the entered sequence"""
    
    GC_regions = []
    GC_percentage = lambda s: (s.count("C") + s.count("G"))/len(s)  
    for N in range(18, 6, -1):              # a stem has a length of 6-18 bp
        for i in range(len(seq)-N+1):
            reg = seq[i:(i+N)]
            if GC_percentage(reg) < 0.8:    # most terminators contains 80% GC in their stems
                continue
            GC_regions.append([reg, (i, i+N-1)])
    
    return GC_regions


def Match_n_Find(GC_regions):
    """finds & matches with each other dyadically symmetrical sequences (whose base
    pair sequences are inverted repeats) among GC-rich regions"""
    
    palindromes = []
    length = lambda x: x[1][1]-x[1][0]
    complement = lambda x: [ComplementarySeq(x[0]), x[1]]
    complGC_regions = list(map(complement, GC_regions))
    for i in range(len(GC_regions)-1):
        for j in range(i+1, len(complGC_regions)):
            if length(GC_regions[i]) != length(complGC_regions[j]):
                continue
            for k in range(len(GC_regions[i][0])):
                if GC_regions[i][0][k] != complGC_regions[j][0][-(k+1)]:
                    break
            else:
                palindromes.append(GC_regions[i] + GC_regions[j])
    return palindromes


def PoliTScore(seq):
    """calculates the score for the given sequence"""
    
    x0 = 1
    score = 0
    for i in range(len(seq)):
        x0 = 0.9*x0 if seq[i] == "T" else 0.6*x0
        score -= x0
    return score


def PoliTSearch(seq):
    """searches poli-T regions (thymine-only) in the entered sequences"""
    
    poliT_regions = []
    min_score = -4.8
    for N in range(8, 5, -1):       # poli-T regions has an average length of 6-8 nt
        for i in range(len(seq)-N+1):
            reg = seq[i:(i+N)]
            if not reg.startswith("TTT"): # this kind of sequence can't be a poli-T region
                continue
            if PoliTScore(reg) > min_score:
                continue
            poliT_regions.append([reg, (i, i+N-1)])
        min_score += 0.3
    
    for reg in poliT_regions:       # deletes duplicates which are contained in longer sequences
        for reg_duplicate in reversed(poliT_regions):
            if reg[1][0] == reg_duplicate[1][0] and reg[1][1] != reg_duplicate[1][1]:
                del poliT_regions[poliT_regions.index(reg_duplicate)]
            
    return poliT_regions


def FoldAll(stem, poliT_regions):
    """fold pieces of sequence into a whole Rho-independent terminator (if it's possible)"""
    
    terminators = []
    no_interruption = lambda x, y: y[1][0]-x[3][1] == 1
    for s in stem:
        for t in poliT_regions:
            if no_interruption(s, t):
                terminators.append(s + t)
    return terminators

def GibbsFreeEnergy(t):
    """calculates Gibbs free energy ΔG [kcal/mol] of a terminator
    based on base pairing, base stacking and length of a loop
    
        --- the lower Gibbs energy, the more stable is a terminator --- """
    
    G = 0
    b_pairing = {"A": 0.57, "T": 0.57, "C": -0.11, "G": -0.11}
    b_stacking = {"AA": -1.2, "AT": -1.6, "TA": -1.6, "AG": -2.1,
                  "AC": -2.1, "CA": -2.1, "GA": -2.1, "CC": -4.8,
                  "CG": -3.0, "GC": -4.3, "TT": 1.0, "TC": -2.1,
                  "TG": 0.0, "CT": -2.1, "GG": 0.0, "GT": -0.3}   # Physical Chemistry of Nucleid Acid (authors: Bloomfield, Crothers, Tinocco)
    len_loop = {3: 8.0, 4: 5.0, 5: 5.0, 6: 4.0, 7: 4.5, 8: 5.0}
    loop_length = lambda x: (x[3][0]-1)-x[1][1]
    
    G += len_loop[loop_length(t)]
    for i in range(len(t[0])-1):
        G += b_pairing[t[0][i]] + b_stacking[t[0][i:(i+2)]]
    G += b_pairing[t[0][-1]]
    
    return G

def ShowTerminators(seq, terminators):
    """shows graphically how each terminator from the given list looks like"""
    
    ts = seq.replace("T", "U")
    l = lambda x: (x[3][0]-1)-x[1][1]
    for n, t in enumerate(terminators):
        print("TERMINATOR #{:03d}\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯".format(n+1))
        print("* loop length: {} nt\n* stem length: {} bp\n* energy: {} kcal/mol".format(l(t), len(t[0]), t[6]))
        print("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯")
        if l(t) == 3 or l(t) == 5 or l(t) == 7:
            i = t[1][1] + ceil(l(t)/2)
            print("{:>4}".format(ts[i]))
            if l(t) == 5 or l(t) == 7:
                print("{}{:>6}".format(ts[i-1], ts[i+1])) if l(t) == 5 else print("{}{:>6}\n{}{:>6}".format(ts[i-1], ts[i+1], ts[i-2], ts[i+2]))
        else:
            i = int(t[1][1] + l(t)/2)
            print("{:>3} {}".format(ts[i], ts[i+1]))
            if l(t) == 6 or l(t) == 8:
                print("{}{:>6}".format(ts[i-1], ts[i+1])) if l(t) == 6 else print("{}{:>6}\n{}{:>6}".format(ts[i-1], ts[i+1], ts[i-2], ts[i+2]))
        print("{:>2}{:>4}".format(ts[t[1][1]+1], ts[t[3][0]-1]))
    
        j, k = t[1][1], t[3][0]
        while j >= t[1][0] and k <= t[3][1]:
            s = "≡" if ts[j] == "C" or ts[j] == "G" else "="
            print("{:>3}{}{}".format(ts[j], s, ts[k])) if k < t[3][1] else (print("{:>3}{}{}{}".format(ts[j], s, ts[k], ts[t[5][0]:])) if t[5][1] >= len(ts)-15 else print("{:>3}{}{}{}".format(ts[j], s, ts[k], ts[t[5][0]:t[5][1]+16])))
            j, k = j-1, k+1
        s = " " * (len(t[4])+1)
        print("{:>3}{}{}\n".format(j, s, t[5][1]))
    return
        
    
def TerminatorSearch(seq):
    """searches Rho-independent transcription terminators in the entered sequence
    seq_cod - DNA coding sequence (5' -> 3'; +)
    strand - a type of strand where a potential terminator is searched"""
    
    seq_cod = seq.upper()
    GC_regions = GCRichRegionsSearch(seq_cod)
    palindromes = Match_n_Find(GC_regions)
    if palindromes == []:
        return "No terminators found"
    loop_length = lambda x: (x[3][0]-1)-x[1][1] >= 3 and (x[3][0]-1)-x[1][1] <= 8
    stem = list(filter(loop_length, palindromes)) # a loop in a terminator has an average length of 3-8 nt
    if stem == []:
        return "No terminators found"
    poliT_regions = PoliTSearch(seq_cod)
    terminators = FoldAll(stem, poliT_regions)
    if terminators == []:
        return "No terminators found"
    for t in terminators:
        G = GibbsFreeEnergy(t)
        t.append(G)
    ShowTerminators(seq_cod, terminators)
    return terminators

if __name__ == '__main__':
    DNA = "AACAAACCGAGCCCGCCTAATGAGCGGGCTTTTTTTTTGACAATGAAAAAACGACAAAGCAGCGCGGATTATCGCGCTGCTTTTTTATCCCTGT"
    TerminatorSearch(DNA)

#                                       OUTPUT
########################################################################################
########################################################################################
# TERMINATOR #001
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
# * loop length: 6 nt
# * stem length: 8 bp
# * energy: -19.299999999999997 kcal/mol
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
#  U U
#A     U
# G   U
#  G≡C
#  C≡G
#  G≡C
#  C≡G
#  G≡C
#  A=U
#  C≡G
#  G≡CUUUUUUAUCCCUGU
# 57         87

# TERMINATOR #002
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
# * loop length: 7 nt
# * stem length: 7 bp
# * energy: -18.89 kcal/mol
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
#   A
#A     U
#U     G
# C   A
#  C≡G
#  G≡C
#  C≡G
#  C≡G
#  C≡G
#  G≡C
#  A=UUUUUUUUUGACAAUGAAAAAACG
#  8         37

# TERMINATOR #003
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
# * loop length: 8 nt
# * stem length: 7 bp
# * energy: -15.19 kcal/mol
#¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
#  U U
#A     U
#G     A
# G   C
#  C≡G
#  G≡C
#  C≡G
#  G≡C
#  A=U
#  C≡G
#  G≡CUUUUUUAUCCCUGU
# 57         87
