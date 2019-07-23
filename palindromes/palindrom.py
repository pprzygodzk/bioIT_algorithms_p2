def ComplementarySeq(seq):
    """makes & returns a complementary sequence to the entered sequence
    seq - sequence"""
    
    nucleobases = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complementary_seq = ""
    for base in seq:
        complementary_seq += nucleobases[base]
    
    return complementary_seq


def LongPalindromeSearch(seq1, seq2):
    """searches long palindrome sequences"""
    
    long_palindrome = []
    for N in range(20, 13, -1):
        for i in range(len(seq1)-N+1):
            pal1 = seq1[i:(i+N)]
            pal2 = seq2[i:(i+N)]
            ans = 1
            for j in range(N):
                if pal1[j] != pal2[-(j+1)]:
                    ans = 0
                    break
            if ans == 1:
                long_palindrome.append(pal1)
                indexes = (seq1.index(pal1), seq1.index(pal1)+N-1) # index at the start & the end of palindrome
                long_palindrome.append(indexes)
                break
        if ans == 1: # one iteration is enough cause it could search shorter palindromes which are contained in the longer ones
            break
        
    return long_palindrome


def ShortPalindromeSearch(seq1, seq2, initial_seq):
    """searches short palindrome sequences"""
    
    short_palindrome = []
    for N in range(8, 3, -1):
        for i in range(len(seq1)-N+1):
            pal1 = seq1[i:(i+N)]
            pal2 = seq2[i:(i+N)]
            ans = 1
            for j in range(N):
                if pal1[j] != pal2[-(j+1)]:
                    ans = 0
                    break
            if ans == 1:
                if pal1 in short_palindrome: # if found palindrome is already in the list but a new one is at further positions
                    if short_palindrome.count(pal1) == 1:
                        tmp = short_palindrome[short_palindrome.index(pal1)+1]
                    if short_palindrome.count(pal1) >= 2:
                        tmp = last_occur
                    short_palindrome.append(pal1)
                    indexes = (initial_seq.index(pal1, tmp[0]+1), initial_seq.index(pal1, tmp[0]+1)+N-1)
                    short_palindrome.append(indexes)
                    last_occur = indexes
                else: # if searched palindrome is not found before
                    short_palindrome.append(pal1)
                    indexes = (initial_seq.index(pal1), initial_seq.index(pal1)+N-1)
                    short_palindrome.append(indexes)
    
    return short_palindrome

def n_PalindromeSearch(seq1, seq2, n):
    """searches palindromic sequences of length n"""
    
    n_palindrome = []
    for i in range(len(seq1)-n+1):
        pal1 = seq1[i:(i+n)]
        pal2 = seq2[i:(i+n)]
        ans = 1
        for j in range(n):
            if pal1[j] != pal2[-(j+1)]:
                ans = 0
                break
        if ans == 1:
            if pal1 in n_palindrome: # if found palindrome is already in the list but a new one is at further positions
                if n_palindrome.count(pal1) == 1:
                    tmp = n_palindrome[n_palindrome.index(pal1)+1]
                if n_palindrome.count(pal1) >= 2:
                    tmp = last_occur
                n_palindrome.append(pal1)
                indexes = (seq1.index(pal1, tmp[0]+1), seq1.index(pal1, tmp[0]+1)+n-1)
                n_palindrome.append(indexes)
                last_occur = indexes
            else: # if searched palindrome is not found before
                n_palindrome.append(pal1)
                indexes = (seq1.index(pal1), seq1.index(pal1)+n-1)
                n_palindrome.append(indexes)
    
    return n_palindrome

def DeleteDuplicates(palindrome_list):
    """deletes duplicate palindrome (contained in longer ones)
    from the entered list of palindromes"""
    
    for i in range(0, len(palindrome_list), 2):
        for j in range(i+2, len(palindrome_list), 2):
            try:
                if str(palindrome_list[j]) in str(palindrome_list[i]) and len(palindrome_list[j]) < len(palindrome_list[i]):
                    palindrome_list.remove(palindrome_list[j+1])
                    palindrome_list.remove(palindrome_list[j])
                    j=j-2
            except IndexError:
                break
    
    return palindrome_list


def PalindromeSearch(seq, n = 0):
    """searches the longest of long and short palindrome sequences in the entered sequence
    seq_cod - sense (coding) strand sequence (+)
    seq_templ - antisense (template) strand sequence (-)"""
    
    assert seq != "", "The entered sequence cannot be empty!"
    
    seq_cod = seq.upper()
    seq_templ = ComplementarySeq(seq_cod)

    if n:
        n_palindrome = n_PalindromeSearch(seq_cod, seq_templ, n)
        n_palindrome = DeleteDuplicates(n_palindrome)
        
        return n_palindrome
    else:
        long_palindrome = LongPalindromeSearch(seq_cod, seq_templ)

        if long_palindrome != []:
            new_cod = seq_cod[0:long_palindrome[1][0]] + "B" + seq_cod[(long_palindrome[1][1]+1):]
            new_templ = seq_templ[0:long_palindrome[1][0]] + "R" + seq_templ[(long_palindrome[1][1]+1):]
        else:
            new_cod = seq_cod
            new_templ = seq_templ
        
        short_palindrome = ShortPalindromeSearch(new_cod, new_templ, seq_cod)
        short_palindrome = DeleteDuplicates(short_palindrome)
    
        return long_palindrome, short_palindrome

if __name__ == '__main__':
    test_seq = "TAATATACTGGACCTCGGATCCGAGGTGTGTAAACATGTGTGCAACCGGTT"
    test_long, test_short = PalindromeSearch(test_seq) # the long palindromic sequence is ACCTCGGATCCGAGGT in the test_seq
    print("Long p.:", test_long)
    print("Short p.:", test_short)
    
    test_n = PalindromeSearch(test_seq, 4)
    print("\nN-long p.:", test_n)
