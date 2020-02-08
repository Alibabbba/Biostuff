
import random as rn

f=open("text.txt", "r")
g=open("text2.txt", "r")
h=open("text_newlines.txt","r")
text = str(f.read())
text2 = str(g.read())
text_newlines = [line.rstrip() for line in h.readlines()]
text_list = text_newlines[0].split(" ")
text_lists = []
#print(text_newlines[0].split(" "))
#for i in text_newlines:
#    text_lists.append(text_newlines[i]
#print(text_newlines)

def Count(text,pattern):
    count = 0
    for i in range (len(text) - len(pattern) + 1):
        if text[i : i + len(pattern)] == pattern:
            count += 1
    return count

#print(Count(text,text2))

def FrequentWords(text, lenght):
    frequentpatternts = []
    count = []
    for i in range (len(text) - lenght + 1):
        Pattern = text[i:i + lenght]
        count.append(Count(text,Pattern))
        print(count)
    maxCount = max(count)
    print(maxCount)
    for j in range (len(text) - lenght + 1):
        if count[j] == maxCount:
            frequentpatternts.append(text[j:j+lenght])
    return list(dict.fromkeys(frequentpatternts))
   
#print(FrequentWords(text, 4))

def DNACount(Text, Pattern):
        count = 0
        i = 0
        for i in range (len(Text) - len(Pattern)):
            if Text[i : i + len(Pattern)] == Pattern:
                count += 1
        return count

#print(DNACount(text, text2))
        
def SymbolToNumber(Symbol):
    if Symbol == "A":
        return 0
    elif Symbol == "C":
        return 1
    elif Symbol == "G":
        return 2
    elif Symbol == "T":
        return 3
    else:
        return "BLAND"
    
#print(SymbolToNumber("C"))
    
def NumberToSymbol(Number):
    if Number == 0:
        return "A"
    elif Number == 1:
        return "C"
    elif Number == 2:
        return "G"
    elif Number == 3:
        return "T"
    else:
        return "Bland"
    
def PatternToNumber(Pattern):
    x = len(Pattern)
    if Pattern == "":
        return 0
    symbol = Pattern[x-1]
    Prefix = Pattern[0: x - 1]
    return 4 * PatternToNumber(Prefix) + SymbolToNumber(symbol)
        
#print(PatternToNumber("ACGTAAGTTGCGGCAACT"))
    
def NumberToPattern(index, k):
    if k == 1:
        return NumberToSymbol(index)
    prefixindex = int(index/4)
    r = index % 4
    symbol = NumberToSymbol(r)
    PrefixPattern = NumberToPattern(prefixindex, k - 1)
    return PrefixPattern + symbol
        
#print(NumberToPattern(5461,9))

def Skew(text):
    skew = 0
    Skew = []
    Skew.append(skew)
    for i in text:
        if i == "G":
            skew +=1
            Skew.append(skew)
        elif i == "C":
            skew -= 1
            Skew.append(skew)
        else:
            skew = skew
            Skew.append(skew)
    return Skew

#print(Skew(text))

def Skew_min(text):
    skew = Skew(text)
    _min = []
    low = 0
    for i in skew:
        if i < low:
            low = i
    for j in range (len(skew)):
        if skew[j] == low:
            _min.append(j)
    return _min

#print(Skew_min(text))

def Skew_max(text):
    skew = Skew(text)
    _max = []
    high = 0
    for i in skew:
        if i > high:
            high = i
    for j in range (len(skew)):
        if skew[j] == high:
            _max.append(j)
    return _max

#print(Skew_max(text))

def Hammering(text,text2):
    if len(text) == len(text2):
        count = 0
        for i in range (len(text)):
            if text[i] != text2[i]:
                count += 1
        return count
    else:
        return "WTF diffrent lenght"

#print(Hammering(text,text2))
        
def All_Max_Hammering(text,pattern,max_mismatch):
    count_starts = []
    for i in range (len(text) - len(pattern) + 1):
        if Hammering(text[i:i+len(pattern)],pattern) <= max_mismatch:
            count_starts.append(str(i))
    return count_starts

#print(All_Max_Hammering(text,text2,2))
    
def ApproximatePatternCount(text,pattern, max_mismatch):
    a = All_Max_Hammering(text,pattern, max_mismatch)
    count = 0
    for i in a:
        count += 1
    return  count

#print(ApproximatePatternCount(text,text2,1))


def ImmitateNaighbors(pattern):
    aray= []
    DNA = ['A','C','G','T']
    for letter in pattern:
        aray.append(letter)
    neighborhood = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for j in DNA:
            arayB = aray.copy()
            if j != symbol:
                arayB[i] = j
            else:   
                continue
    return neighborhood

#print(ImmitateNaighbors(text2))

def Neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A','C','G','T'}
    Neighborhood = set()
    SuffixNeighbors = Neighbors(pattern[1:],d)
    for i in SuffixNeighbors:
        if Hammering(pattern[1:],i) < d:
            for j in ['A','C','G','T']:
                Neighborhood.add(j+i)
        else:
            Neighborhood.add(pattern[0] + i)
            
    return Neighborhood

#print(Neighbors("ATTTGGC",3))
    
def MotifEnumeration(Dna, k_mer, mismatches):
    results = set() 
    beta = []
    for part_Dna in Dna:
        alfa = set()
        beta.append(alfa)
        for pattern_start in range(len(part_Dna) - k_mer + 1):
            pattern = part_Dna[pattern_start:pattern_start + k_mer]
            pat = Neighbors(pattern, mismatches)
            alfa.update(pat)
    results = beta[0]
    for s in beta[1:]:
        results.intersection_update(s)
    return  results

#print(MotifEnumeration(text_newlines, 5, 1))


def DistanceBetweenMotifAndStrings(Pattern, Dna_strings):
    k = len(Pattern)
    distance = 0
    for string in Dna_strings:
#        print(string,Dna_strings, Pattern)
        HammingDistance = 1000000
        for k_mer_start in range( len(string) - k + 1) :
            k_mer = string[k_mer_start:k_mer_start + len(Pattern)]
            if HammingDistance > Hammering(Pattern, k_mer):
                HammingDistance = Hammering(Pattern, k_mer)
#                print(Pattern,k_mer,"\n")
        distance += HammingDistance
    return distance

#a = ['CGCCCCTC'] #'CCCTCTCG', 'TTCAGTAA', 'CTCTCGGG', 'GGGGGTGT']
#print(DistanceBetweenMotifAndStrings("CCCTCTCG", text_newlines))

def MedianString(Dna, k_mer):
    distance = 1000000
    for i in range(4 ** k_mer - 1):
        Pattern = NumberToPattern(i , k_mer)
        currentdistance = DistanceBetweenMotifAndStrings(Pattern, Dna)
        if distance > currentdistance:
            distance = currentdistance
            Median = Pattern
    return Median

#print(MedianString(text_newlines,7))

def ProbableKMer(Dna,k_mer,matrix):
    results = [[],[]]
    kil = ""
    for i in range(len(Dna) - k_mer + 1):
        kmer = Dna[i : i + k_mer]
        PossibleKmer = 1
        for symbol in range(k_mer):
            if kmer[symbol] == "A":
                PossibleKmer *= matrix[0][symbol]
            elif kmer[symbol] == "C":
                PossibleKmer *= matrix[1][symbol]
            elif kmer[symbol] == "G":
                PossibleKmer *= matrix[2][symbol]
            elif kmer[symbol] == "T":
                PossibleKmer *= matrix[3][symbol]
        results[0].append(kmer)
        results[1].append(PossibleKmer)
    for j in range(len(results[1])):
        if results[1][j] == max(results[1]):
            kil = results[0][j]
    return kil


def Profile_Dna(dna_list):
    results = [[],[],[],[]]
    for i in range (4):
        for j in range (len(dna_list[0])):
            results[i].append(0) 
    for dna in dna_list:
        for index, base in enumerate(dna):
            if base == 'A':
                results[0][index] += 1
            elif base == 'C':
                results[1][index] += 1
            elif base == 'G':
                results[2][index] += 1
            elif base == 'T':
                results[3][index] += 1
    for index, nucleotide in enumerate(results):
        for probability in range (len(nucleotide)):
            results[index][probability] /= len(dna_list)
    return results

# print(Profile_Dna(text_newlines))

def Random_starts(dna_list, k_mer):
    random_start = []
    for starts in range (len(dna_list)):
        a = rn.randrange(len(dna_list[0]) - k_mer)
        random_start.append(dna_list[0][a:a+k_mer])
    return random_start

    start = random_start
    print(start)   
    print(DistanceBetweenMotifAndStrings(start,dna_list))

    for max_iteration in range(1):
        results = [[],[],[],[]]
        posibility_matrix = []
        for i in range (4):
            for j in range (k_mer):
                results[i].append(0) 
        for dna in start:
            for index, base in enumerate(dna):
#                print(dna,index,base)
                if base == 'A':
                    results[0][index] += 1
                elif base == 'C':
                    results[1][index] += 1
                elif base == 'G':
                    results[2][index] += 1
                elif base == 'T':
                    results[3][index] += 1
#        print(results)
        for index, nucleotide in enumerate(results):
            for probability in range (len(nucleotide)):
                results[index][probability] /= len(dna_list)
                if results[index][probability] == 0.0:
                   results[index][probability] = 0.01 
        for mer in dna_list:
            posibility_matrix_mer = []
            for i in range (len(mer)- k_mer + 1):
                print('\n')
                posibility_matrix_kmer, posibility_matrix_posibility = "", 1
                for number,nucleotide in enumerate(mer[i:i+k_mer]):
                    print(posibility_matrix_posibility)
                    if nucleotide == 'A':
                        posibility_matrix_kmer += nucleotide
                        posibility_matrix_posibility *= results[0][number]
                    elif nucleotide == 'C':
                        posibility_matrix_kmer += nucleotide
                        posibility_matrix_posibility *= results[1][number]
                    elif nucleotide == 'G':
                        posibility_matrix_kmer += nucleotide
                        posibility_matrix_posibility *= results[2][number]
                    elif nucleotide == 'T':
                        posibility_matrix_kmer += nucleotide
                        posibility_matrix_posibility *= results[3][number]
                print(posibility_matrix_posibility, posibility_matrix_kmer, '\n' ,results)
                posibility_matrix_posibility = round(posibility_matrix_posibility,4)
                posibility_matrix_mer.append([posibility_matrix_kmer, posibility_matrix_posibility])
                print(posibility_matrix_mer)
            posibility_matrix.append(posibility_matrix_mer)
#        print(posibility_matrix)
    #            print(nucleotide,number)     
    #            print(mer[i:i+k_mer],mer)
    #            print(i)
        results_motifs = []
        results_motifs_new = []
        for dna_line in posibility_matrix:
            most_probable = 0
            most_probable_number = 0
            for inden,k_mer_new in enumerate(dna_line):
                if k_mer_new[1] == 0:
                    k_mer_new[1] = 0.0000001
                if k_mer_new[1] > most_probable:
                    most_probable = k_mer_new[1]
                    most_probable_number = inden
#                print(k_mer_new, inden)
            results_motifs.append(dna_line[most_probable_number])
#        print(results_motifs)
        for i in results_motifs:
            results_motifs_new.append(i[0])
        
#        print(start == results_motifs_new)
#        distance = DistanceBetweenMotifAndStrings(results_motifs_new[0],dna_list)
        distance = 0
        for pattern_number ,pattern in enumerate(results_motifs_new):
            distance = DistanceBetweenMotifAndStrings(pattern,dna_list[pattern_number])
#            print(distance, pattern, dna_list[pattern_number])
#        print("\n")
        distance = 0

        start = results_motifs_new
#        print(distance, results_motifs_new[0], "\n" ,dna_list)
#        print(results_motifs_new)
#        print(mer)
#    print(posibility_matrix_mer,"\n")
    print(dna_list)
    return None
#print(Random_starts(text_newlines, kamer))
#print(text_newlines)

f.close()
g.close()
h.close()