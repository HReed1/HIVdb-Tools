import re, os.path, sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
def outRightMatched(DRM, line):
    #################################
    #        Match Codons           #
    #################################
    A = ['A','R','W','M','D','V' 'H']
    G = ['G','R','S','K','D','V','B']
    C = ['C','Y','S','M','V','H','B']
    T = ['T','Y','W','K','D','H','B']
    mismatch = 0
    DRM = re.findall('.',(''.join(DRM)))
    line = re.findall('.',(''.join(line)))
    Range = len(DRM)
    
    ################################
    #       Special Exceptions     #
    ################################
    if(len(DRM) != len(line)):
        return False
    if(DRM == ['C','C','C','A','T','T','A','G','C'] and line == ['C','C','C','A','T','T','A','G','T']):
        return True
    if(DRM == ['T','T','T','T','A','A','A','T','T'] and (line == ['T','T','T','T','A','A','A','T','C'] or line == ['C','T','T','T','A','A','G','T','T'] or line == ['C','T','T','T','A','A','A','T','T'] or line == ['T','T','T','T','A','A','A','T','C',])):
        return True
    if(DRM == ['T','T','T','T','T','A','G','A','T'] and (line == ['T','T','T','C','T','A','G','A','T'] or line == ['T','T','T','T','T','G','G','A','T'])):
        return True

    #### Add up mismatches ####
    for i in range(Range):
        if(DRM[i] == 'T'):
            if (line[i] not in T):
                mismatch += 1
        elif(DRM[i] == 'A'):
            if (line[i] not in A):
                mismatch += 1
        elif(DRM[i] == 'G'):
            if (line[i] not in G):
                mismatch += 1
        elif(DRM[i] == 'C'):
            if (line[i] not in C):
                mismatch += 1

    if(mismatch == 0):
        return True
    else:
        return False

def matched(PDRM, Pline):
    #Matches the front of DRM  at the AA level 
    change = False
    DRM = ''
    line = ''
    PDRM = ''.join(PDRM)
    Pline = ''.join(Pline)
    DRM = list(PDRM)
    line = list(Pline)
    if('nnn' in Pline):
        return False
    PDRM = (Seq(PDRM)).translate()
    Pline = (Seq(Pline)).translate()
    if not(Pline==PDRM):
        for i in range(len(Pline)):
            if(Pline[i] == 'X'):
                if(i ==0):
                    line[i:(i+3)] = DRM[i:(i+3)]
                    change = True
                else:
                    line[i*3:((i*3)+3)] = DRM[i*3:((i*3)+3)]
                    change = True
        if(change==True):
            line = ''.join(line)
            DRM = ''.join(DRM)
            line = (Seq(line)).translate()
            DRM = (Seq(DRM)).translate()
        return DRM == line
    else:
        return(Pline == PDRM)

def endMatched(PDRM, Pline):
    #Matches the back of DRM at the AA level
    change = False
    DRM = ''
    line = ''
    PDRM = ''.join(PDRM)
    PDRM = PDRM[::-1]
    Pline = ''.join(Pline)
    Pline = Pline[::-1]
    DRM = list(PDRM)
    line = list(Pline)
    if('nnn' in Pline):
        return False
    PDRM = (Seq(PDRM)).translate()
    Pline = (Seq(Pline)).translate()
    if not(Pline==PDRM):
        for i in range(len(Pline)):
            if(Pline[i] == 'X'):
                if(i ==0):
                    line[i:(i+3)] = DRM[i:(i+3)]
                    change = True
                else:
                    line[i*3:((i*3)+3)] = DRM[i*3:((i*3)+3)]
                    change = True
        if(change==True):
            line = ''.join(line)
            DRM = ''.join(DRM)
            line = (Seq(line)).translate()
            DRM = (Seq(DRM)).translate()
        return DRM == line
    else:
        return(Pline == PDRM)

def DRMwrite (string, ID, IDcount, line):
    print('Writing Sample ' + str(ID[IDcount]) + ' to ' + string)
    with open(string, 'a') as f:
        f.write(">" + ID[IDcount] + "\n" + line)

def find(DRM,line):
    slen = 2
    endSlen = 3
    stop = False
    frameShift = False
    backShift = False
    outRightMatch = False
    outRightEndMatch = False
    outShift = False
    outBackShift = False
    truthSeeker = [False,False, -1, -1, -1, -1]
    length = len(''.join(DRM))
    orgLength = length
    lineLength = len(''.join(line))
    DRM = ''.join(DRM)
    line = ''.join(line)
    endDRM = DRM[::-1]
    endLine = line[::-1]
    DRM = re.findall('...',DRM)
    temp = DRM
    line = re.findall('...',line)
    temp3 = line
    endDRM = re.findall('...',endDRM)
    endLine = re.findall('...',endLine)
    temp2 = endLine
    pop = 0
    snap = 0
   # print('\nStart finder\nDRM: ' + str(DRM) + '\nline: ' + str(line) + '\n')
    while(stop == False):
        match = 0
        matchPos = []
        endMatchPos = []
        endMatch = 0
   #     print('\nStart DRM: ' + str(DRM[0:3]) + ' Line: ' + str(line[0:3]))
    #    print('End DRM: ' + str(endDRM[0:3]) + " Line: " + str(endLine[0:3]) + "\n")
        if(frameShift == True and outShift == False):
    #        print('FrameShift: OUT, Invoked!')
            line = list(''.join(line))
            line.pop(0)
            pop += 1
            lineLength -= 1
            line = re.findall('...',(''.join(line)))
            frameShift = False
        if(frameShift == True and outShift == True):
    #        print('FrameShift Invoked! ')
            DRM = list(''.join(DRM))
            DRM.pop(0)
            pop +=1
            DRM = re.findall('...',(''.join(DRM)))
            length -= 1
            frameShift = False
        if(outShift == False):
  #          print('Section: OutRightMatch')
            for q in range(lineLength):
                if(outRightMatched(DRM[0:(slen+1)], line[q:q+(slen+1)])):
                    matchPos= []
                    match += 1
                    if(q==0):
   #                     print('OurRight Match! at POS: ' + str(q) + ':' + str(line[q:(q+(slen+1))]) + 'match: ' + str(match) + '\n')
                        matchPos.append(q)
                        truthSeeker[4] = (q)
                    else:
   #                     print('OurRight Match! at POS: ' + str((q*3)-pop) + ':' + str(line[q:(q+(slen+1))]) + 'match: ' + str(match) + '\n')
                        matchPos.append((q*3)+pop)
                        truthSeeker[4] = (q*3)+pop
                    outRightMatch = True
                elif(pop == 3):
                    lineLength += pop
                    pop = 0
                    line = temp3
                    outShift = True
        if(outRightMatch == False and outShift == True):
 #           print('Section: matched')
            for k in range(length):
                if(matched(DRM[k:(k+slen)],line[0:slen])):
                    matchPos = []
                    match += 1
                    if(k ==0):
    #                    print('PostMatch at POS: '+ str(k) + ':' + str(DRM[k:(k+slen)]) + 'match: ' + str(match)+ "\n")
                        matchPos.append(k)
                    else:
    #                    print('PostMatch at POS: '+ str((k*3)-pop+1) + ':' + str(DRM[k:(k+slen)]) + 'match: ' + str(match)+ "\n")
                       # matchPos.append((k*3)-pop+1)
                        matchPos.append((k*3)+pop+1)
        if(backShift == True):
    #        print('BackShift Invoked! ')
            endLine = list(''.join(endLine))
            endLine.pop(0)
            snap +=1
            endLine = re.findall('...',(''.join(endLine)))
            length -= 1
            backShift = False
        if(outBackShift == False):
  #          print('Section: OutRightBackMatched')
            for e in range(lineLength):
                if(outRightMatched(endDRM[0:endSlen],endLine[e:(e+endSlen)])):
                    endMatchPos = []
                    endMatch += 1
                    if(e == 0):
     #                   print('\n\nOutright EndMatch: at POS ' +str(e+len(''.join(temp2))-(snap+pop-2)) + ' : ' + str(endLine[e:(e+endSlen)]) + 'endMatch: ' + str(endMatch) + "\n\n")
                        endMatchPos.append(e+len(''.join(temp2))-(snap+pop-2))
                        truthSeeker[5] = (e+len(''.join(temp2))-(snap+pop-2))
                    else:
      #                  print('\n\nOutright EndMatch: at POS ' +str(-1*((e*3) - len(''.join(temp2)))) + ' : ' + str(endLine[e:(e+endSlen)]) + 'endMatch: ' + str(endMatch) + "\n\n")
                        endMatchPos.append(-1*((e*3) - len(''.join(temp2))))
                        truthSeeker[5] = (-1*((e*3) - len(''.join(temp2))))
                    outRightEndMatch = True
                elif(snap == 3):
                    lineLength += snap
                    snap = 0
                    endLine = temp2
                    outBackShift = True
        if(outRightEndMatch == False and outBackShift == True):
    #        print('Section: Back Match')
            for j in range(length):
                    if(endMatched(endDRM[0:endSlen],endLine[j:(j+endSlen)])):
                        endMatchPos = []
                        endMatch += 1
     #                   print('EndMatch at POS: '+ str(j) + ':' + str(endLine[0:endSlen]) + 'EndMatch: ' + str(endMatch)+ "\n")
    
                        if(j == 0):
       #                     print('\n\nEndMatch: at POS ' +str(j+len(''.join(endLine))) + ' : ' + str(endLine[j:(j+endSlen)]) + 'endMatch: ' + str(endMatch) + "\n\n")
                            endMatchPos.append(j+len(''.join(endLine)))
                        else:
      #                      print('\n\nEndMatch: at POS ' +str(-1*((j*3) - len(''.join(endLine)))) + ' : ' + str(endLine[j:(j+endSlen)]) + 'endMatch: ' + str(endMatch) + "\n\n")
                            endMatchPos.append(-1*((j*3) - len(''.join(endLine))))
#        print('Section: Checkup')
        if(match == 1 and endMatch == 1):
            stop = True
            if(endMatchPos[0] < matchPos[0]):
       #         print('\n\n---NO MATCH---\n\n')
                truthSeeker[0] = False
                truthSeeker[1] = False
            elif((endMatchPos[0]-matchPos[0]) > orgLength):
        #        print('\n\n***NO MATCH***\n\n')
       #         print('endMatch: ' + str(endMatchPos[0]) + '\nmatchPos' + str(matchPos[0]))
                truthSeeker[0] = False
                truthSeeker[1] = False
            else:
  #              print('\n\n----MATCH!!!---\n\n')
                truthSeeker[0] = True
                truthSeeker[1] = True		
        if(endMatch == 0 and snap == 3):
            stop = True
#            print("\n\n+++NO MATCH+++\n\n")
            truthSeeker[1] = False
            if(match == 1):
                truthSeeker[2] = matchPos[0]
        elif(endMatch == 0 and snap < 3):
            backShift = True
        if(match == 0 and pop == 3):
            stop = True
 #           print("\n\n...NO MATCH...\n\n")
            truthSeeker[0] = False
            if(endMatch == 1):
                truthSeeker[3] = endMatchPos[0]
        elif(match == 0 and pop < 3):
            frameShift = True
        if(match > 1):
            slen +=1
        if(endMatch > 1):
            endSlen += 1
    if(truthSeeker[0] == True and truthSeeker[1] == True):
        truthSeeker[2] = matchPos[0]
        truthSeeker[3] = endMatchPos[0]
      #  print('endMatch: ' + str(endMatchPos[0]) + '\nmatchPos' + str(matchPos[0]))
  
    return truthSeeker
##################TRUTHSEEKER###########################
#[0] = Front Match, Boolean                            #
#[1] = Back Match, Boolean                             #
#[2] = Front Match EndPos, int, -1 == no match         #
#[3] = Back Match EndPos, int, -1 == no match          #
#[4] = OutRightFront Match EndPos, int, -1 == no match #
#[5] = OutRightEnd Match EndPos, int, -1 == no match   #
########################################################

#MAIN
PR=['CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT']

RT=['CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTA']

IN=['TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT']

fast = input("Enter the cons filename\n")
samples = 'sampleList.txt'

ID = []
with open(samples, 'r') as f:
    ID = [line.rstrip('\n') for line in f]
isMatched = False
IDcount = 0
temp = ''
temp2 = ''
PR = ''.join(PR)
RT = ''.join(RT)
IN = ''.join(IN)
with open(fast, 'r') as f:
    for line in f:
        if(isMatched == True):   
            print('\n---------------' + ID[IDcount] + '---------------\n')
            RTforward = False
            PRmatched = find(PR, list(line))
            print('PRmatched: ' + str(PRmatched))
            RTmatched = find(RT, list(line))
            print('RTmatched: ' + str(RTmatched))
            INmatched = find(IN, list(line))
            print('INmatched: ' + str(INmatched))
            temp = line
            
            #If PR is truly matched...
            if(PRmatched[0] == True and PRmatched[1] == True and (PRmatched[3]-PRmatched[2]) > 150):
                if(RTmatched[4] != -1):
                #Checks for and corrects PR/RT overlap
                    PRmatched[3] = (PRmatched[3]+(RTmatched[4]-PRmatched[3]))
#                print('PR COORD ' + str(PRmatched[2]) + ':' + str(PRmatched[3]))
                line = line[PRmatched[2]:PRmatched[3]]
                DRMwrite('AllPR.fasta', ID, IDcount, line)
                #print('\nRTforward ' + str(RTforward))
                line = temp[(PRmatched[3]):] 
     
                #If PR matched & RT OutRightEndMatch & no IN StartMatch. Assume IN starts after RT.   
                if(RTmatched[5] != -1 and INmatched[0] == False):
                    line = temp[PRmatched[3]:RTmatched[5]]
                    DRMwrite('AllRT.fasta',ID,IDcount,line)
                    line = temp[RTmatched[5]:]
                    DRMwrite('AllIN.fasta',ID,IDcount,line)
                    line = temp[PRmatched[3]:]
                    DRMwrite('AllRT_IN.fasta',ID,IDcount,line)
                    line = temp[:RTmatched[5]]
                    DRMwrite('AllPR_RT.fasta',ID,IDcount,line)
                    line = temp
                    DRMwrite('AllReigons.fasta',ID,IDcount,line)
                    RTforward = True
                    #print('First Loop Easy Match')

                #Regain Line for More Searches. 
                RTmatched = find(RT, list(line)) 
                INmatched = find(IN,list(line)) 
                
                #If PR and IN is are  matched...
                if(INmatched[0] == True and INmatched[1] == True):
                    temp2 = line #temp 2 saves the PRend-to-end line
                    line = temp2[INmatched[2]:] #Only capturing IN
                    DRMwrite('AllIN.fasta', ID, IDcount, line)
                    line = temp2[:INmatched[2]]
                    #RT is Between PR and IN 
  #                  print('\nRT COORD: ' + str(PRmatched[3]) + ':' + str(INmatched[2]) + '\n')
                    DRMwrite('AllRT.fasta', ID, IDcount, line)
                    line = temp[:(INmatched[2]+PRmatched[3])]
                    DRMwrite('AllPR_RT.fasta', ID, IDcount, line)
                    line = temp[PRmatched[3]:]
                    DRMwrite('AllRT_IN.fasta', ID, IDcount, line)
                    line = temp #Regain original full line
                    DRMwrite('AllReigons.fasta', ID, IDcount, line)

                #If PR but not IN is matched...The rest is RT
                elif(RTforward == False):
                    line = temp[(PRmatched[3]):]
   #                 print('\nRT COORD: ' + str(PRmatched[3]) + ": (To the End)") 
                    DRMwrite('AllRT.fasta',ID,IDcount,line)
                    line = temp
                    DRMwrite('AllPR_RT.fasta',ID,IDcount,line)

            #If PR no match & RT OutRightEndMatch....The rest is IN
            elif(RTmatched[5] != -1 and INmatched[0] == False):
                    line = temp[:RTmatched[5]]
                    DRMwrite('AllRT.fasta',ID,IDcount,line)
                    line = temp[RTmatched[5]:]
                    DRMwrite('AllIN.fasta',ID,IDcount,line)
                    line = temp
                    DRMwrite('AllRT_IN.fasta',ID,IDcount,line)
                    RTforward = True

            #IF PR no match & IN matches.... 
            elif(INmatched[0] == True and INmatched[1] == True and (INmatched[3] -INmatched[2] > 800)):
                line = temp
                line = line[INmatched[2]:]
                DRMwrite('AllIN.fasta',ID,IDcount,line)

                #IF RT matches
                if(RTmatched[0] == True and RTmatched[1] == True):
                    line = temp
                    line = temp[:INmatched[2]]
    #                print('\nRT COORD: 0:' + str(INmatched[2]) + '\n')
                    DRMwrite('AllRT.fasta',ID,IDcount,line)
                    line = temp
                    DRMwrite('AllRT_IN.fasta',ID,IDcount,line)

            #IF PR and IN have not Matched....the whole thing is RT
            elif(RTforward == False):
                DRMwrite('AllRT.fasta',ID,IDcount,line)
            #Pushes forward the next sample or ends the program.
            isMatched = False 

        if(isMatched == False):
            it = re.finditer(r"^>(.+)", line)
            for match in it:
                for i in range(len(ID)):
                    if ID[i] == match.group(1):
                        isMatched = True
                        IDcount = i
                        break
        
