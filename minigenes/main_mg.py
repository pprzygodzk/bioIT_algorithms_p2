import minigeny as mg

print("Przed podaniem sekwencji należy zaznaczyć kodony start dla genów\nantytoksyny i toksyny MAŁYMI LITERAMI!")
seq = input("Podaj sekwencję genu antytoksyny >>")
mg.PotentialMGSearch(seq)
i = 1
while i:
    print("Czy chcesz zbadać kolejną sekwencję?")
    ans = input("Wpisz y lub Y, jesli tak, albo n lub N, jesli nie >>")
    if ans == "y" or ans == "Y":
        seq = input("Podaj sekwencję genu antytoksyny >>")
        mg.PotentialMGSearch(seq)
    if ans == "n" or ans == "N":
        break
print("Dziekuje za skorzystanie z programu :)")