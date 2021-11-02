import os
import re

total_reads = 0
threshold = 0.7

first_line = 1
seq_list = []
seq_dict = {}

wordsize = 8
otucount= 1

path = '/Users/rashondaogletree/PycharmProjects/pythonProject1/sample.txt'

trimlength = input("Enter trimlength: ")
print("you entered " + trimlength)

def is_file_empty(path):
    # check if file exist and is it empty
    return os.path.exists(path) and os.stat(path).st_size == 0

is_empty = is_file_empty(path)

if is_empty:
    print('File is empty')


def word_list(sequence, wordsize):
    wordlist = []

    for i in range(0, len(sequence) - (int(wordsize) + 1)):
        word = sequence[i:int(wordsize)]
        wordlist.append(word)

        return wordlist

def make_otu(sequence, abundance, wordlist):
    otuname = "OTU"
    otu_dict = {}
    total_reads=0

    global otucount

    otucount += 1

    #otu_dict[otuname]= otu_dict.setdefault(sequence, "seedSeq")
    otu_dict[otuname]= otu_dict.setdefault("seedSeq", sequence)
    otu_dict[otuname]= otu_dict.setdefault("totalCount", 0)

    for word in wordlist:
        otu_dict["totalCount"]= otu_dict.setdefault("totalCount") +abundance
        otu_dict["word"] = otu_dict.setdefault("word", word)

        print("OTU NAME:" + otuname + "-" + str(otucount) + "," + "WORD:" + word + "," + "ABUNDANCE:" + str(abundance)) #WordFile

    total_reads +=1
    print(otu_dict["seedSeq"]) #SeedSequencesFile

def score_otu(sequence, wordList):
    bestotu = ""
    sumscore_for_otu = 0
    max_score = -999
    total_word_count= 0
    otuname= ""

    otu_dict = {}

    total_word_count = otu_dict.setdefault("totalCount", 0) +1
    for word in wordList:
        freq= otu_dict.setdefault(word, 0) +1
        curr_score_for_word = (freq / total_word_count)
        sumscore_for_otu += curr_score_for_word

    sedSeq= otu_dict.setdefault("seedSeq", sequence)
    sumscore_for_otu *= (len(sedSeq) / len(sequence))

    if (sumscore_for_otu >= max_score):
        max_score = sumscore_for_otu
        bestotu = otuname

        print(seq_dict["read"] + otuname + bestotu)

    return max_score, bestotu

def update_otu(sequence, abundance, wordList):
    otuname = ""

    for word in wordList:
        seq_dict[otuname] = seq_dict.setdefault("totalCount", 0) + 1
        seq_dict["word"] = seq_dict.setdefault(word, 0) +abundance

regex = r"^>([A-Z_0-9|a-z.]+)\s([A-Z]+)"
for line in open(path, 'r'):
    a1 = re.search(regex, line)
    if a1:
        header = a1.group(1)
        sequence = a1.group(2)
        if first_line == 1:
            length = len(sequence)
            if (length > int(trimlength)):
                seq = sequence[0:int(trimlength)]
                sequence = seq
                seq_list.append(sequence)

                for sample in seq_list:
                    seq_dict[sample] = seq_dict.setdefault(sample, 0) + 1
        seq_dict["read"]= seq_dict.setdefault("read", header)

    else:
        sequence = ""
        first_line = 0

for key in sorted(seq_dict):
    seq_len = len(key)
    if (seq_len > 10):
        # print ("%s: %s" % (key, seq_dict[key]))
        if key in seq_dict:
            abundance = seq_dict[key]
        else:
            abundance = 0
            print(key + "not find in data dictionary")

        currSeqW_List = word_list(key, wordsize)
        # print(currSeqW_List)

        if (otucount == 1):
            make_otu(key, abundance, currSeqW_List)
        else:
            temp_best_score, temp_best_otu = score_otu(key, currSeqW_List)
            if temp_best_score >= threshold:
                update_otu(key, abundance,currSeqW_List)
            else:
                make_otu(key, abundance, currSeqW_List)



