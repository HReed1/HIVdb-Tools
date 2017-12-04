#HIV_DB: EXTRACT

import csv, os.path

def remove_duplicates(List):
    output = []
    seen = set()
    for value in List:
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output

#Getting the data
#Takes in a list and a filename/path
def get_list(ls, fl):
    with open(fl, 'r') as f:
        reader = csv.reader(f)
        for i in reader:
            ls = list(reader)
    return ls

def extract_mut(ID,d,region):
    #Making region header
    for i in range(len(ID)):
        with open("Variants_{0}.txt".format(ID[i]), 'a') as f:
            f.write("\n" +region + "\n")
    #Extracting mutations
    for i in range(len(d)): #Tracks current sample per ID
        mut = ''
        for j in range(len(ID)): #Tracks current ID
            if ((ID[j] == d[i][0]) and (d[i][1] == region)): #if ID is the same and Gene hasn't changed...
                mut += (str(d[i][4]) + str(d[i][5]) + str(d[i][6]) +"\n") #...We'll add it to the gene list
                break #No need to iterate the whole list again...so we stop, and move on to the next sample.
        if (mut != ''):         #If the last ID had any mutations....
            with open("Variants_{0}.txt".format(ID[j]), 'a') as f:
                f.write(mut)    #We add those mutaions to the files.
#Main
user = input("Enter file with path: \n")
while(not os.path.exists(user)):
    user = input("File not found! Enter file with path or Ctrl-z to exit: ")

d = []
d = get_list(d,user)
ID = [i[0] for i in d]
ID = remove_duplicates(ID)

mut_RT = 'RT'
mut_PR = 'PR'
mut_IN = 'IN'

extract_mut(ID,d,mut_RT)
extract_mut(ID,d,mut_PR)
extract_mut(ID,d,mut_IN)
