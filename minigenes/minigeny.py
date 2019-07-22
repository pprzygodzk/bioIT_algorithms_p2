import re

def PotentialMGSearch(seq):
    """wyszukuje potencjalne minigeny majace dlugosc 2 aminokwasow
    w obrebie sekwencji antytoksyny
    seq - sekwencja kasety toksyna-antytoksyna
    at_start - miejsce startu ramki odczytu antytoksyny
    t_start - miejsce startu ramki odczytu toksyny"""
    
    assert seq != "", "sekwencja nie moze byc pusta"
    
    p = re.compile("atg|gtg|ttg")
    m = p.findall(seq)
    if m:
        i = 1
        for match in p.finditer(seq):
            tmp = match.span()
            if i:
                at_start = tmp[0]
                i = 0
            t_start = tmp[0]
    
    p = re.compile("ATG[ACGT]{3}TAA|ATG[ACGT]{3}TAG|ATG[ACGT]{3}TGA")
    m = p.findall(seq, at_start+3, t_start)
    i_list = [] # lista indeksow dopasowan
    m_list = [] # lista dopasowan
    if m:
        i = 0
        for match in p.finditer(seq, at_start+3, t_start):
            print(str(match.span()) + ' ' + str(match.group()) + '\tmatch nr ' + str(i) + '\n')
            i += 1
            i_list.append(match.span())
            m_list.append(match.group())
            
    if m_list == []:
        print("Nie znaleziono dopasowan")
        
    for indexes in i_list:
        if (t_start-indexes[0]) % 3 == 0:
            print(m_list[i_list.index(indexes)] + " jest zgodne z ramka odczytu toksyny")
        elif(indexes[0]-at_start) % 3 == 0:
            print(m_list[i_list.index(indexes)] + " jest zgodne z ramka odczytu antytoksyny")
        else:
            print(m_list[i_list.index(indexes)] + " - sekwencja ta nie jest minigenem")
            # nie jest zgodna z zadna ramka odczytu

    
if __name__ == '__main__':
    print("**********System TA u Enterococcus faecium**********")
    axe_txe = """AATTGTTTTATAGAAATAAATAAGGGGTGAAAGGAatgGAAGCAGTAGCTTATTCAAATTTCCGCCAAAATTTACGTAGTTATATGAAACAAGTTAATGAGGATGCTGAAACACTTATTGTAACAAGTAAAGATGTAGAAGATACAGTTGTTGTATTATCAAAAAGAGATTATGATTCTATGCAAGAAACGTTGAGAACACTTTCTAATAATTACGTCATGGAAAAAATTCGTCGAGGAGATGAACAATTCTCCAAAGGTGCATTTAAAACACATGACTTAATCGAGGTTGAATCTGatg"""
    PotentialMGSearch(axe_txe)
    
    print("\n\n**********System TA u Escherichia coli CFT073**********")
    Ecoli_CFT073 = """atgCGTACAATTAGCTACAGCGAAGCGCGTCAGAATTTGTCGGCAACAATGATGAAAGCCGTTGAAGATCATGCCCCGATCCTCATTACTCGTCAGAATGGAGAGGCTTGTGTTCTGATGTCACTCGAAGAATACAACTCGCTGGAAGAGACGGCTTATCTACTGCGTTCCCCCGCTAACGCCCGGAGATTGATGGACTCAATCGATAGCCTGAAATCAGGCAAAGGAACGGAAAAGGACATTATTGAgtgA"""
    PotentialMGSearch(Ecoli_CFT073)
    
    print("\n\n**********System TA u Lactobacillus rhamnosus Lc 705**********")
    Lacto_705 = """atgGAAGCAACGAATTATAGTGATTTCCGCCGCAACCTTAAGCATTATATGAGTCAAGTCAACGAAGACGCCGAACCGCTACTGGTTACCGCTAAAGATGATGATGACAATGTGGTGGTTATGAGCAAGCACGATTTTGACGCCATCGAAGAAACCCTGTATTTACTCAGCAATCCCAAGCTGATGGCCAAAATCAAACGTGGTGATGCCCAAATTGCCGCTGGAAAGGCTAAACAGCACGAGTTGTTAACGGACTTCGATCatgATTAA"""
    PotentialMGSearch(Lacto_705)