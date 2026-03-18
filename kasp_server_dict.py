import re
import math
import sys
import json
from Bio.SeqUtils import MeltingTemp as mt

def callSequence(chr, pos, radius):
    startPos, realPos = 0, 0
    
    fIdx = open('/storage1/sihoo219/web_source/Reference/Gmax_275_v2.0.idx', 'rt')
    for line in fIdx:
        item = line.strip().split('\t')
        if int(item[0][-2:]) == int(chr[-2:]):
            startPos = item[1]
    realPos = int(startPos) + int(pos) - 1
    fIdx.close()

    fSeq = open('/storage1/sihoo219/web_source/Reference/Gmax_275_v2.0.seq', 'rt')
    fSeq.seek(realPos - radius)
    sequence = fSeq.read(radius * 2  + 1)
    fSeq.close()

    return sequence, realPos

def checkTmValue(seq):
    tmValue = round(mt.Tm_NN(seq), 3)

    if 55 <= tmValue <= 65: return True, tmValue
    else:   return False, 0

def checkRepeatSequence(seq):
    # monoP = re.compile(r'([ATGC])\1{2,}')
    diP = re.compile(r'(AT|TA|CG|GC|AG|GA|CT|TC|AC|CA|GT|TG)\1{3,}')

    # monoRepeat = monoP.search(seq)
    diRepeat = diP.search(seq)

    # if not monoRepeat and not diRepeat: return True
    if not diRepeat: return True
    else: return False

def checkGCcontent(seq):
    GCcontent = round(((seq.count('G') + seq.count('C')) * 100 / len(seq)), 3)
    if 40 <= GCcontent <= 60:   return True, GCcontent
    else:   return False, 0

def convertComplementary(seq):
    compDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    convertedSeq = []
    
    for i in range(0, len(seq)):
        convertedSeq.append(compDict[seq[i]])

    return "".join(convertedSeq)


if __name__ == "__main__":
    selectedDatas_json = sys.argv[1]
    selectedDatas = json.loads(selectedDatas_json)      # selectedDatas = [ {"chr":"chr01", "pos":"12345", "marker":"chr01_12345"}, {}, {}, ... ] -> dict in array

    allPrimerList = []
    primerDict = {}
    
    for data in selectedDatas:
        # data = { "chr":"chr01", "pos":"12345", "marker":"chr01_12345" } -> dict
        marker_info = data['marker'].split("_")
        marker_chr, marker_pos = marker_info[0], marker_info[1]
        
        seq, realSnpPos = callSequence(marker_chr, marker_pos, 200)
        snpIdx = len(seq) // 2

        primerMin, primerMax = 18, 25
        ampliconMin, ampliconMax = 100, 150

        # Forward primer
        for lenF in range(primerMin, primerMax):
            primerF = seq[(snpIdx - lenF):snpIdx]
            tmBoolF, tmValF = checkTmValue(primerF)
            gcBoolF, gcValF = checkGCcontent(primerF) 
            if tmBoolF and gcBoolF and checkRepeatSequence(primerF):

                #Reverse primer
                for ampSize in range(ampliconMin, ampliconMax):
                    for lenR in range(primerMin, primerMax):
                        primerDict = {}
                        primerR = seq[snpIdx+ampSize+1:(snpIdx+ampSize+1)+lenR]
                        tmBoolR, tmValR = checkTmValue(primerR)
                        gcBoolR, gcValR = checkGCcontent(primerR)
                        tmDiff = round(abs(tmValF - tmValR), 3)
                        if tmBoolR and checkRepeatSequence(primerR) and gcBoolR and tmDiff <= 5:
                            # 1) template
                            primerDict['chr'] = data['chr']
                            primerDict['pos'] = data['pos']
                            primerDict['ref'] = data['ref']
                            primerDict['alt'] = data['alt']
                            # 2) forward primer
                            primerDict['fStart'] = (realSnpPos + 1) - (lenF - 1)
                            primerDict['fEnd'] = (realSnpPos + 1)
                            primerDict['fSize'] = lenF
                            primerDict['fTm'] = tmValF
                            primerDict['fGC'] = gcValF
                            primerDict['fSeq'] = primerF
                            # 3) reverse primer
                            primerDict['rStart'] = realSnpPos + ampSize
                            primerDict['rEnd'] = realSnpPos + ampSize + (lenR-1)
                            primerDict['rSize'] = lenR
                            primerDict['rTm'] = tmValR
                            primerDict['rGC'] = gcValR
                            primerDict['rSeq'] = convertComplementary(primerR[::-1])
                            # 4) relation of between 2 primers
                            primerDict['tmDiff'] = tmDiff
                            primerDict['ampSize'] = ampSize

                            allPrimerList.append(primerDict)

    result = allPrimerList
    
    print(json.dumps(result))  # JSON 문자열로 출력