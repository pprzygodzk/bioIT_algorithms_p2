import promotor as p
from datetime import datetime

try:
    file = open("promoter_results.txt", 'a')
except IOError:
    print("Error: File cannot be open!")
else:
    file.write(str(datetime.now()) + '\n')
    file.write("****************************************************************************")
    print("PROKARYOTIC PROMOTER SEARCH")
    print("*********************************************************************************")
    name = input("Enter the name of a potential promoter region >>")
    file.write("\n****************************************************************************\nNAME:\t" + name + '\n')
    print("*********************************************************************************")
    print("Before you enter your sequence, remember that this program requires sequence running in the direction 5'->3'.")
    print("\nIf your gene is in the strand (-) and you have a sequence from the strand (+), make her reversed and complementary.")
    seq = input("Enter the sequence >>")
    file.write("SEQUENCE:\t" + seq)
    print("*********************************************************************************")
    print("FOUND PROMOTERS:")
    promoters = p.ProkaryoticPromoterSearch(seq)
    file.write("\nFOUND PROMOTERS:\t" + str(promoters))
    while True:
        print("*********************************************************************************")
        print("Do you want to test the next gene?")
        ans = input("Enter y or Y (if you do) and n or N (if not) >>")
        if ans == "y" or ans == "Y":
            print("*********************************************************************************\n")
            name = input("Enter the name of a potential promoter region >>")
            file.write('\n****************************************************************************\nNAME:\t' + name + '\n')
            print("*********************************************************************************\n")
            print("\nBefore you enter your sequence, remember that this program requires sequence running in the direction 5'->3'.")
            print("\nIf your gene is in the strand (-) and you have a sequence from the strand (+), make her reversed and complementary.")
            seq = input("\nEnter the sequence >>")
            file.write("SEQUENCE:\t" + seq)
            print("*********************************************************************************\n")
            print("\nFOUND PROMOTERS:")
            promoters = p.ProkaryoticPromoterSearch(seq)
            file.write("\nFOUND PROMOTERS:\t" + str(promoters))
            print("*********************************************************************************\n")
        if ans == "n" or ans == "N":
            break
finally:
    file.write("\n****************************************************************************\n\n\n")
    file.close()