import palindrom as p
from datetime import datetime

try:
    file = open("palindromes_results.txt", 'a')
except IOError:
    print("Error: File cannot be open!")
else:
    file.write(str(datetime.now()) + '\n')
    file.write("****************************************************************************\n")
    name = input("Enter the name of tested gene(s) >>")
    file.write(name + '\n')
    
    seq_cod = input("Enter the sense strand's sequence (+) >>")
    file.write('5\' ' + seq_cod + ' 3\'\n')
    seq_templ = p.ComplementarySeq(seq_cod)
    file.write('3\' ' + seq_templ + ' 5\'\n\n')
    
    print("If you want to search palindromic sequences of specific length, please enter the length\n(number of nucleotides). If you don't have any specific length and want to search all palindromic sequences, enter \"0\".")
    n = input("Enter the length of palindromic sequence >>")
    n = int(n)
    
    if n != 0:
        n_palndrm = p.PalindromeSearch(seq_cod, n)
        file.write(str(n) + "-long palindromes: " + str(n_palndrm) + '\n\n\n')
    else:
        l_palndrm, sh_palndrm = p.PalindromeSearch(seq_cod)
        file.write("Long palindromes: " + str(l_palndrm) + '\n')
        file.write("Short palindromes: " + str(sh_palndrm) + '\n\n\n')
    
    while True:
        print("Do you want to test the next gene?")
        ans = input("Enter y or Y (if you do) and n or N (if not) >>")
        if ans == "y" or ans == "Y":
            name = input("Enter the name of tested gene(s) >>")
            file.write("\n\n" + name + '\n')
    
            seq_cod = input("Enter the sense strand's sequence (+) >>")
            file.write('5\' ' + seq_cod + ' 3\'\n')
            seq_templ = p.ComplementarySeq(seq_cod)
            file.write('3\' ' + seq_templ + ' 5\'\n\n')
    
            print("If you want to search palindromic sequences of specific length, please enter the length\n(number of nucleotides). If you don't have any specific length and want to search all palindromic sequences, enter \"0\".")
            n = input("Enter the length of palindromic sequence >>")
            n = int(n)
            
            if n != 0:
                n_palndrm = p.PalindromeSearch(seq_cod, n)
                file.write(str(n) + "-long palindromes: " + str(n_palndrm) + '\n\n\n')
            else:
                l_palndrm, sh_palndrm = p.PalindromeSearch(seq_cod)
                file.write("Long palindromes: " + str(l_palndrm) + '\n')
                file.write("Short palindromes: " + str(sh_palndrm) + '\n\n\n')
        if ans == "n" or ans == "N":
            break

finally:     
    file.close()
    