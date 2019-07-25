import math as m

def Count(seq, consensus):
    """counts occurences of nucleobases at the specific positions of both sequences"""
    
    l = len(consensus)
    base_count = {'A': [0]*l, 'C': [0]*l, 'G': [0]*l, 'T': [0]*l }
    for s in [seq, consensus]:
        for i in range(l):
            for base in base_count.keys():
                if s[i] == "N": # if, according to consensus, a nucleobase at the i-th position 
                    base_count[base][i] = 2 # doesn't matter, then we don't want to change
                    # a score no matter what nucleobase is at the i-th position of a sequence seq
                if s[i] == base:
                    base_count[base][i] += 1
    return base_count


def Profile(seq, consensus):
    """makes a probability matrix for two given sequences: seq and consensus"""
    
    profile = Count(seq, consensus)
    for base in profile.keys():
        for i in range(len(consensus)):
            profile[base][i] /= 2 # max 2 nucleobases at the i-th position
    return profile


def Score(seq, consensus, pos_d, pos_m = [-1]):
    """calculates a score of the given seq relatively to the consensus
    (uses entropy H as a measure of alignment)
    seq - a query sequence
    consensus - an order of most frequent nucleotides at specific positions in a promoter
    pos_d - a list of positions where one of the nucleobases are dominant over three others
    for the given consensus
    pos_m - a list of positions where one of the nucleobases is moderately significant relatively
    to three others"""
    
    p = Profile(seq, consensus)
    H = 0
    for i in range(len(seq)):
        if i in pos_d:
            H += -(4*p[seq[i]][i]*m.log2(p[seq[i]][i]))
        elif i in pos_m:
            H += -(1.25*p[seq[i]][i]*m.log2(p[seq[i]][i]))
        else:
            H += -(p[seq[i]][i]*m.log2(p[seq[i]][i]))
    return H


def PotentialBoxSites(sequence, box, pos_d, pos_m = [-1], indexes = [-1]):
    """looks for potential box sites in a given sequence"""
    
    PBS = []
    n = len(sequence)
    l = len(box)
    if indexes != [-1]:
        for i, j in enumerate(indexes):
            if j[0]-3 < 0:
                continue
            s = Score(sequence[(j[0]-3):j[1]+1], box, pos_d, pos_m)
            if s <= 2:
                PBS.append([sequence[(j[0]-3):j[1]+1], (j[0]-3, j[1])])
        return PBS
    for i in range(n-l+1):
        s = Score(sequence[i:(i+l)], box, pos_d, pos_m)
        if s <= 2:
            PBS.append([sequence[i:i+l], (i, i+l-1), s])
    return PBS


def PotentialSigma70Promoter(P35BS, P10BS):
    """compares distances between one of the sequences in potential -35 box sites (P35BS)
    and all of the sequence in potential TATA-box sites (P10BS), then joins them into pairs
    (one promoter)"""
    
    PS70P = []
    
    for TTGbox in P35BS:
        a = TTGbox[1][1]
        for TATAbox in P10BS:
            b = TATAbox[1][0]
            if b-a <= 20 and b-a >= 15:
                PS70P.append([TTGbox[0], TTGbox[1], TATAbox[0], TATAbox[1]])
    return PS70P


def ExtractTGNTATASites(P10BS, promoters):
    """extract indexes of sequences from the list of potential TATA-box sites (P10BS)
    which weren't joined with any -35 box site"""
    
    last = 0
    while True:
        if promoters == []:
            break
        if last < len(P10BS):
            for i in range(last, len(P10BS)):
                for j in range(len(promoters)):
                    if P10BS[i][0] in promoters[j] and P10BS[i][1] in promoters[j]:
                        del P10BS[i]
                        last = i
                        break
                if i == len(P10BS)-1:
                    last = i+1
                    break
        if last == len(P10BS):
            break
    
    indexes = []
    for i in range(len(P10BS)):
        indexes.append(P10BS[i][1])
    return indexes


def ProkaryoticPromoterSearch(sequence):
    """looks for a gene prokaryotic promoter in a given sequence
    sequence - a region of gene where a promoter is assumed to be"""
    
    sequence = sequence.upper()
    box = ["TTGACA", "TATAAT", "TGNTATAAT"]
    P35BS = PotentialBoxSites(sequence, box[0], [0, 1, 5])
    # according to logo TTGACA, positions: 1st, 2nd and 6th are the most important
    P10BS = PotentialBoxSites(sequence, box[1], [0, 1, 5], [3, 4])
    # according to logo TATAAT, positions: 1st, 2nd and 6th are the most important, while
    # positions 4th and 5th are moderately significant
    promoters = PotentialSigma70Promoter(P35BS, P10BS)
    indexes = ExtractTGNTATASites(P10BS, promoters)
    PTGN10BS = PotentialBoxSites(sequence, box[2], [3, 4, 8], [6, 7], indexes)
    promoters += PTGN10BS
    print(promoters)
    return promoters


if __name__ == '__main__':
    test_seq = "ATGCTGTCATTGACACCAGTGCCTGTGAACGTTATAATGGTC"
    ProkaryoticPromoterSearch(test_seq)
