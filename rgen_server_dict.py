import re
import subprocess
from Bio.SeqUtils import MeltingTemp as mt
from multiprocessing import Pool

import pickle
import argparse
import numpy as np
import sys, math, json
from functools import reduce

IUPAC_code = {
    "A": "A", "T": "T", "G": "G", "C": "C",
    "N": "[ATGC]", "R": "[AG]", "Y": "[CT]", "V":"[AGC]",
    "W": "[AT]", "S": "[GC]", "M": "[AC]", "K": "[GT]"
}

def callSequence(chr, pos, radius):
    startPos, realPos = 0, 0
    
    fIdx = open('/storage1/sihoo219/web_source/Reference/Gmax_275_v2.0.idx', 'rt')
    for line in fIdx:
        item = line.strip().split('\t')
        if item[0] == chr:
            startPos = item[1]
    realPos = int(startPos) + (int(pos) - 1)
    fIdx.close()

    fSeq = open('/storage1/sihoo219/web_source/Reference/Gmax_275_v2.0.seq', 'rt')
    fSeq.seek(realPos - radius)
    sequence = fSeq.read(radius * 2  + 1)
    fSeq.close()

    return sequence

def checkRepeatSequence(seq):
    diP = re.compile(r'(AT|TA|CG|GC|AG|GA|CT|TC|AC|CA|GT|TG)\1{3,}')

    diRepeat = diP.search(seq)

    if not diRepeat: return True
    else: return False

def calculateGCcontent(seq):
    GCcontent = (seq.count('G') + seq.count('C')) * 100 / len(seq)
    return GCcontent

def convertComplementary(seq):
    compDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    convertedSeq = []
    
    for i in range(0, len(seq)):
        convertedSeq.append(compDict[seq[i]])

    return "".join(convertedSeq)

def match_pam(seq, pam):
    regex = "".join(IUPAC_code[p] for p in pam)
    if re.match(regex, seq):
        return True
    return False

def init_result_dict(strand):
    return {
        "marker_chr": "", 
        "marker_pos": "",
        "marker_ref": "",
        "marker_alt": "",
        "spacer_chr": "", 
        "spacer_start_pos": "", 
        "spacer_end_pos": "",
        "spacer_size": "",
        "spacer_strand": strand,
        "spacer": "",
        "spacer_gc": "",
        "spacer_cleavage_site": "",
        "protospacer": "", 
        "protospacer_cleavage_site": "",
        "protospacer_gc": "",
        "pam_seq": "",
        "context": "",
        "mismatch_0": 0,
        "mismatch_1": 0,
        "mismatch_2": 0,
        "off_score": -1,
        "on_score": -1,
        "off_candidates": [],
    }

def init_result_dict_BE(strand):
    return {
        "marker_chr": "", 
        "marker_pos": "",
        "marker_ref": "",
        "marker_alt": "",
        "marker_relative_pos": "",
        "spacer_chr": "", 
        "spacer_start_pos": "", 
        "spacer_end_pos": "",
        "spacer_size": "",
        "spacer_strand": strand,
        "spacer": "",
        "spacer_gc": "",
        "spacer_window": "",
        "spacer_window_edited": "",
        "protospacer": "", 
        "protospacer_window": "",
        "protospacer_window_edited": "",
        "protospacer_gc": "",
        "pam_seq": "",
        "context": "",
        "mismatch_0": 0,
        "mismatch_1": 0,
        "mismatch_2": 0,
        "off_score": -1,
        "on_score": -1,
        "off_candidates": [],
    }

def handle_Cas(marker_dict, spacerLen, opt2):
    results, pamList = [], []
    cas9_dict = {
        "SpCas9":["NRG"], 
        "VQR-Cas9":["NGA"], 
        "VRER-Cas9":["NGCG"], 
        "EQR-Cas9":["NGAG"], 
        "xCas9":["GAW"], 
        "SpCas9-NG":["NG"],
        "StCas9":["NNAGAAW"], 
        "NmCas9":["NNNNGMTT"],
        "SaCas9":["NNGRRT"],
    }
    cas12_dict = {
        "SpCas12a":["TTTV"], 
        "AsCas12a/LbCas12a":["TTTN"], 
        "FnCas12a":["TTN"]
    }
    if opt2 in cas9_dict.keys():
        pamList = cas9_dict[opt2]
        results = design_Cas9(marker_dict, spacerLen, cas9_dict[opt2])
    elif opt2 in cas12_dict.keys():
        pamList = cas12_dict[opt2]
        results = design_Cas12(marker_dict, spacerLen, cas12_dict[opt2])

    return results, pamList

def handle_ABE(marker_dict, spacerLen, opt2):
    isCorrectKey = 0
    results = []

    cas9_ABE_dict = {
        "ABE(7.10)/ABEmax/SECURE": {"window_start": 4, "window_end": 7, "pams": ["NGG"]},
        "VQR/VRQR": {"window_start": 4, "window_end": 7, "pams": ["NGA"]},
        "VRER": {"window_start": 4, "window_end": 7, "pams": ["NGCG"]} ,
        "NG-ABE": {"window_start": 4, "window_end": 7, "pams": ["NG"]} ,
        "xCas9(3.7)": {"window_start": 4, "window_end": 8, "pams": ["GAW"]} ,
        "Spy-macABEmax": {"window_start": 4, "window_end": 7, "pams": ["NAA"]},
        "Sa-ABE": {"window_start": 4, "window_end": 14, "pams": ["NNGRRT"]},
        "SaKKH": {"window_start": 4, "window_end": 14, "pams": ["NNNRRT"]},
        "ScCas9": {"window_start": 4, "window_end": 7, "pams": ["NNG"]},
        "SauriABEmax": {"window_start": 6, "window_end": 14, "pams": ["NNGG"]},
        "ABE8e": {"window_start": 4, "window_end": 8, "pams": ["NGG"]},
        "SaABE8e": {"window_start": 3, "window_end": 14, "pams": ["NNGRRT"]},
        "ABE8s": {"window_start": 3, "window_end": 9, "pams": ["NGG"]},
        "NG-ABE8s": {"window_start": 3, "window_end": 9, "pams": ["NG"]},
        "Sa-ABE8s": {"window_start": 5, "window_end": 14, "pams": ["NNGRRT"]},
        "CP1012/CP1028/CP1041/CP1249": {"window_start": 4, "window_end": 12, "pams": ["NGG"]}
    }
    cas12_ABE_dict = {
        "LbABE8e/enAsABE8e": {"window_start": 8, "window_end": 14, "pams": ["TTTV"]},
    }

    if opt2 in cas9_ABE_dict.keys():
        isCorrectKey = 1
        pamList = cas9_ABE_dict[opt2]["pams"]
        results = design_Cas9_BE(marker_dict, spacerLen, cas9_ABE_dict[opt2])
    elif opt2 in cas12_ABE_dict.keys():
        isCorrectKey = 1
        pamList = cas12_ABE_dict[opt2]["pams"]
        results = design_Cas12_BE(marker_dict, spacerLen, cas12_ABE_dict[opt2])

    if isCorrectKey == 1:
        for res in results:
            res["spacer_window_edited"] = res["spacer_window"].replace("A", "G")
            res["protospacer_window_edited"] = res["protospacer_window"].replace("T", "C")
    else:
        return [], []

    return results, pamList
    
def handle_CBE(marker_dict, spacerLen, opt2):
    isCorrectKey = 0
    results = []

    cas9_CBE_dict = {
        "BE3/BE4/HF/A3A/eA3A": {"window_start": 4, "window_end": 8, "pams": ["NGG"]},
        "YE1/YE2/YEE/EE/SECURE(R33A/K34A)": {"window_start": 5, "window_end": 6, "pams": ["NGG"]},
        "CP1012/CP1028": {"window_start": 4, "window_end": 11, "pams": ["NGG"]},
        "xCas9(3.7)": {"window_start": 4, "window_end": 8, "pams": ["GAW"]},
        "VQR": {"window_start": 4, "window_end": 8, "pams": ["NGA"]},
        "EQR": {"window_start": 4, "window_end": 8, "pams": ["NGAG"]},
        "VRER": {"window_start": 4, "window_end": 8, "pams": ["NGCG"]},
        "NG-CBE": {"window_start": 4, "window_end": 8, "pams": ["NG"]},
        "Spy-macBE4max": {"window_start": 4, "window_end": 8, "pams": ["NAA"]},
        "Sa-CBE": {"window_start": 3, "window_end": 12, "pams": ["NNGRRT"]},
        "SaKKH": {"window_start": 3, "window_end": 12, "pams": ["NNNRRT"]},
        "SauriBE4max": {"window_start": 6, "window_end": 9, "pams": ["NNGG"]},
        "ScCas9": {"window_start": 4, "window_end": 8, "pams": ["NNG"]},
        "Target-AID": {"window_start": 2, "window_end": 4, "pams": ["NGG"]},
        "Target-AID-NG": {"window_start": 2, "window_end": 4, "pams": ["NG"]},
        "SsAPOBEC3B": {"window_start": 2, "window_end": 15, "pams": ["NGG"]}
    }
    cas12_CBE_dict = {
        "dLbCpf1/enAsCas12a": {"window_start": 8, "window_end": 13, "pams": ["TTTV"]},
        "dCpf1": {"window_start": 10, "window_end": 12, "pams": ["TTTV"]},
    }
    
    if opt2 in cas9_CBE_dict.keys():
        isCorrectKey = 1
        pamList = cas9_CBE_dict[opt2]["pams"]
        results = design_Cas9_BE(marker_dict, spacerLen, cas9_CBE_dict[opt2])
    elif opt2 in cas12_CBE_dict.keys():
        isCorrectKey = 1
        pamList = cas12_CBE_dict[opt2]["pams"]
        results = design_Cas12_BE(marker_dict, spacerLen, cas12_CBE_dict[opt2])

    if isCorrectKey == 1:
        for res in results:
            res["spacer_window_edited"] = res["spacer_window"].replace("C", "T")
            res["protospacer_window_edited"] = res["protospacer_window"].replace("G", "A")
    else: 
        return [], []

    return results, pamList

def handle_GBE(marker_dict, spacerLen, opt2):
    isCorrectKey = 0
    results = []

    cas9_CBE_dict = {
        "BE3/BE4/HF/A3A/eA3A": {"window_start": 4, "window_end": 8, "pams": ["NGG"]},
        "YE1/YE2/YEE/EE/SECURE(R33A/K34A)": {"window_start": 5, "window_end": 6, "pams": ["NGG"]},
        "CP1012/CP1028": {"window_start": 4, "window_end": 11, "pams": ["NGG"]},
        "xCas9(3.7)": {"window_start": 4, "window_end": 8, "pams": ["GAW"]} ,
        "VQR": {"window_start": 4, "window_end": 8, "pams": ["NGA"]},
        "EQR": {"window_start": 4, "window_end": 8, "pams": ["NGAG"]},
        "VRER": {"window_start": 4, "window_end": 8, "pams": ["NGCG"]},
        "NG-CBE": {"window_start": 4, "window_end": 8, "pams": ["NG"]},
        "Spy-macBE4max": {"window_start": 4, "window_end": 8, "pams": ["NAA"]},
        "Sa-CBE": {"window_start": 3, "window_end": 12, "pams": ["NNGRRT"]},
        "SaKKH": {"window_start": 3, "window_end": 12, "pams": ["NNNRRT"]},
        "SauriBE4max": {"window_start": 6, "window_end": 9, "pams": ["NNGG"]},
        "ScCas9": {"window_start": 4, "window_end": 8, "pams": ["NNG"]},
        "Target-AID": {"window_start": 2, "window_end": 4, "pams": ["NGG"]},
        "Target-AID-NG": {"window_start": 2, "window_end": 4, "pams": ["NG"]},
        "SsAPOBEC3B": {"window_start": 2, "window_end": 15, "pams": ["NGG"]}
    }
    cas12_CBE_dict = {
        "dLbCpf1/enAsCas12a": {"window_start": 8, "window_end": 13, "pams": ["TTTV"]},
        "dCpf1": {"window_start": 10, "window_end": 12, "pams": ["TTTV"]},
    }
    
    if opt2 in cas9_CBE_dict.keys():
        isCorrectKey = 1
        pamList = cas9_CBE_dict[opt2]["pams"]
        results = design_Cas9_BE(marker_dict, spacerLen, cas9_CBE_dict[opt2])
    elif opt2 in cas12_CBE_dict.keys():
        isCorrectKey = 1
        pamList = cas12_CBE_dict[opt2]["pams"]
        results = design_Cas12_BE(marker_dict, spacerLen, cas12_CBE_dict[opt2])

    if isCorrectKey == 1:
        for res in results:
            res["spacer_window_edited"] = res["spacer_window"].replace("C", "G")
            res["protospacer_window_edited"] = res["protospacer_window"].replace("G", "C")
    else:
        return [], []        

    return results, pamList

def design_Cas9(marker_dict, spacerLen, pams):
    results = []

    snp = marker_dict['marker']
    snpChr, snpPos = snp.split("_")
    fSeq = callSequence(snpChr, snpPos, 30)
    rSeq = convertComplementary(fSeq)
    snpIdx = len(fSeq) // 2

    for pam in pams:
        fPam = fSeq[(snpIdx + 3):(snpIdx + 3 + len(pam))]
        if match_pam(fPam, pam):
            res = init_result_dict("+")
            res.update({
                "marker_chr": snpChr,
                "marker_pos": snpPos,
                "marker_ref": marker_dict["ref"],
                "marker_alt": marker_dict["alt"], 
                "spacer_chr": snpChr,
                "spacer_start_pos": int(snpPos) - (spacerLen - 3),
                "spacer_end_pos": int(snpPos) + 2,
                "spacer_size": int(snpPos) - (int(snpPos) - (spacerLen - 3)) + 1,
                "spacer": fSeq[(snpIdx - (spacerLen - 3)):(snpIdx + 3)],
                "protospacer": rSeq[(snpIdx - (spacerLen - 3)):(snpIdx + 3)],
                "pam_seq": fPam,
                "spacer_cleavage_site": spacerLen - 2,
                "protospacer_cleavage_site": spacerLen - 2,
                "context": fSeq[(snpIdx + 3) - 24 : (snpIdx + 3) + 6],
            })
            results.append(res)
            
        rPam = rSeq[(snpIdx - 3):(snpIdx - 3 - len(pam)):-1]
        if match_pam(rPam, pam):
            res = init_result_dict("-")
            res.update({
                "marker_chr": snpChr,
                "marker_pos": snpPos,
                "marker_ref": marker_dict["ref"],
                "marker_alt": marker_dict["alt"], 
                "spacer_chr": snpChr,
                "spacer_start_pos": int(snpPos) + (spacerLen - 3),
                "spacer_end_pos": int(snpPos) - 2,
                "spacer_size": int(snpPos) + (spacerLen - 3) - (int(snpPos) - 2) + 1,
                "spacer": rSeq[(snpIdx + (spacerLen - 3)):(snpIdx - 3):-1],
                "protospacer": fSeq[(snpIdx + (spacerLen - 3)):(snpIdx - 3):-1],
                "pam_seq": rPam,
                "spacer_cleavage_site": spacerLen - 2,
                "protospacer_cleavage_site": spacerLen - 2,
                "context": rSeq[((snpIdx - 3) + 24) : ((snpIdx - 3) - 6) :-1],
            })
            results.append(res)

    return results

def design_Cas12(marker_dict, spacerLen, pams):
    results = []

    snp = marker_dict['marker']
    snpChr, snpPos = snp.split("_")
    fSeq = callSequence(snpChr, snpPos, 30)
    rSeq = convertComplementary(fSeq)
    snpIdx = len(fSeq) // 2

    for pam in pams:
        fPam = fSeq[(snpIdx - 18 - len(pam)) + 1:(snpIdx - 18) + 1]
        if match_pam(fPam, pam):
            res = init_result_dict("+")
            res.update({
                "marker_chr": snpChr,
                "marker_pos": snpPos,
                "marker_ref": marker_dict["ref"],
                "marker_alt": marker_dict["alt"],
                "spacer_chr": snpChr,
                "spacer_start_pos": int(snpPos) - 17,
                "spacer_end_pos": int(snpPos) + (spacerLen - 17 - 1),
                "spacer_size": int(snpPos) + (spacerLen - 17 - 1) - (int(snpPos) - 17) + 1,
                "spacer": fSeq[(snpIdx - 17):(snpIdx + spacerLen - 17)],
                "protospacer": rSeq[(snpIdx - 17):(snpIdx + spacerLen - 17)],
                "pam_seq": fPam,
                "spacer_cleavage_site": 18,
                "protospacer_cleavage_site": 23,
                "context": fSeq[(snpIdx - 18) - 8 + 1 : (snpIdx - 18) + 26 + 1]
            })
            results.append(res)

        # Reverse PAM
        rPam = rSeq[(snpIdx + 23 + len(pam)) - 1:(snpIdx + 23) - 1:-1]
        if match_pam(rPam, pam):
            res = init_result_dict("-")
            res.update({
                "marker_chr": snpChr,
                "marker_pos": snpPos,
                "marker_ref": marker_dict["ref"],
                "marker_alt": marker_dict["alt"],
                "spacer_chr": snpChr,
                "spacer_start_pos": int(snpPos) + (spacerLen - 1),
                "spacer_end_pos": int(snpPos),
                "spacer_size": int(snpPos) + (spacerLen - 1) - (int(snpPos)) + 1,
                "spacer": rSeq[(snpIdx + 23 - 1):(snpIdx - (spacerLen - 23) - 1):-1],
                "protospacer": fSeq[(snpIdx + 23 - 1):(snpIdx - (spacerLen - 23) - 1):-1],
                "pam_seq": rPam,
                "spacer_cleavage_site": 18,
                "protospacer_cleavage_site": 23,
                "context": rSeq[(snpIdx + 18) + 8 - 1 : (snpIdx + 18) - 26 - 1:-1]
            })
            results.append(res)

    return results

def design_Cas9_BE(marker_dict, spacerLen, infos):
    results = []

    snp = marker_dict['marker']
    snpChr, snpPos = snp.split("_")
    fSeq = callSequence(snpChr, snpPos, 30)
    rSeq = convertComplementary(fSeq)
    snpIdx = len(fSeq) // 2

    wStart, wEnd, pams = infos["window_start"], infos["window_end"], infos["pams"]

    for pam in pams:
        for wPos in range(wStart, wEnd + 1):
            # Forward PAM
            fPam = fSeq[snpIdx + (spacerLen - wPos) + 1 : snpIdx + (spacerLen - wPos) + 1 + len(pam)]
            if match_pam(fPam, pam):
                spacer = fSeq[snpIdx - wPos + 1 : snpIdx + (spacerLen - wPos) + 1]
                protospacer = rSeq[snpIdx - wPos + 1 : snpIdx + (spacerLen - wPos) + 1]
                context = fSeq[snpIdx - wPos + 1 - 4: snpIdx + (spacerLen - wPos) + 1 + 6]
                
                res = init_result_dict_BE("+")
                res.update({
                    "marker_chr": snpChr,
                    "marker_pos": snpPos,
                    "marker_ref": marker_dict["ref"],
                    "marker_alt": marker_dict["alt"], 
                    "spacer_chr": snpChr,
                    "marker_relative_pos": wPos,
                    "spacer_start_pos": int(snpPos) - wPos + 1,
                    "spacer_end_pos": int(snpPos) + (spacerLen - wPos),
                    "spacer_size": int(snpPos) + (spacerLen - wPos) - (int(snpPos) - wPos + 1) + 1,
                    "spacer": spacer,
                    "spacer_window": spacer[wStart - 1 : wEnd],
                    "protospacer": protospacer,
                    "protospacer_window": protospacer[wStart - 1 : wEnd],
                    "pam_seq": fPam,
                    "context": context
                })
                results.append(res)

            # Reverse PAM
            rPam = rSeq[snpIdx - (spacerLen - wPos) - 1 : snpIdx - (spacerLen - wPos) - 1 - len(pam) : -1]      # -17 ~ -19
            # print(f"{pam}, {rPam}")
            if match_pam(rPam, pam):
                spacer = rSeq[snpIdx + wPos - 1 : snpIdx - (spacerLen - wPos) - 1 : -1]                         # +3 ~ -16
                protospacer = fSeq[snpIdx + wPos - 1 : snpIdx - (spacerLen - wPos) - 1 : -1]
                context = rSeq[snpIdx + wPos - 1 + 4 : snpIdx - (spacerLen - wPos) - 1 - 6 : -1]
                
                res = init_result_dict_BE("-")
                res.update({
                    "marker_chr": snpChr,
                    "marker_pos": snpPos,
                    "marker_ref": marker_dict["ref"],
                    "marker_alt": marker_dict["alt"],
                    "spacer_chr": snpChr,
                    "marker_relative_pos": wPos,
                    "spacer_start_pos": int(snpPos) + wPos - 1,
                    "spacer_end_pos": int(snpPos) - (spacerLen - wPos),
                    "spacer_size": int(snpPos) + wPos - 1 - (int(snpPos) - (spacerLen - wPos)) + 1,
                    "spacer": spacer,
                    "spacer_window": spacer[wStart - 1 : wEnd],
                    "protospacer": protospacer,
                    "protospacer_window": protospacer[wStart - 1 : wEnd],
                    "pam_seq": rPam,
                    "context": context
                })
                results.append(res)

    return results

def design_Cas12_BE(marker_dict, spacerLen, infos):
    results = []

    snp = marker_dict['marker']
    snpChr, snpPos = snp.split("_")
    fSeq = callSequence(snpChr, snpPos, 30)
    rSeq = convertComplementary(fSeq)
    snpIdx = len(fSeq) // 2

    wStart, wEnd, pams = infos["window_start"], infos["window_end"], infos["pams"]

    for pam in pams:
        for wPos in range(wStart, wEnd + 1):
            # Forward PAM
            fPam = fSeq[snpIdx - wPos - len(pam) + 1 : snpIdx - wPos + 1]                                       # snpIdx - 7 ~ snpIdx - 4              23 ~ 26
            if match_pam(fPam, pam):
                spacer = fSeq[snpIdx - wPos + 1 : snpIdx + (spacerLen - wPos) + 1]
                protospacer = rSeq[snpIdx - wPos + 1 : snpIdx + (spacerLen - wPos) + 1]
                context = fSeq[snpIdx - wPos + 1 - 8 : snpIdx + (spacerLen - wPos) + 1 + 3]
                
                res = init_result_dict_BE("+")
                res.update({
                    "marker_chr": snpChr,
                    "marker_pos": snpPos,
                    "marker_ref": marker_dict["ref"],
                    "marker_alt": marker_dict["alt"],
                    "spacer_chr": snpChr,
                    "marker_relative_pos": wPos,
                    "spacer_start_pos": int(snpPos) - wPos + 1,
                    "spacer_end_pos": int(snpPos) + (spacerLen - wPos),
                    "spacer_size": int(snpPos) + (spacerLen - wPos) - (int(snpPos) - wPos + 1) + 1,
                    "spacer": spacer,
                    "spacer_window": spacer[wStart - 1 : wEnd],
                    "protospacer": protospacer,
                    "protospacer_window": protospacer[wStart - 1 : wEnd],
                    "pam_seq": fPam,
                    "context": context
                })
                results.append(res)

            # Reverse PAM
            rPam = rSeq[snpIdx + wPos + len(pam) - 1 : snpIdx + wPos - 1 : -1]                                   # snpIdx + 4 + 4 - 1 (+7) ~ snpIdx + 4 - 1 (+3)    7 6 5 4
            if match_pam(rPam, pam):
                spacer = rSeq[snpIdx + wPos - 1 : snpIdx - (spacerLen - wPos) - 1 : -1]                        
                protospacer = fSeq[snpIdx + wPos - 1 + 8 : snpIdx - (spacerLen - wPos) - 1 - 3 : -1]
                context = 0
                
                res = init_result_dict_BE("-")
                res.update({
                    "marker_chr": snpChr,
                    "marker_pos": snpPos,
                    "marker_ref": marker_dict["ref"],
                    "marker_alt": marker_dict["alt"],
                    "spacer_chr": snpChr,
                    "marker_relative_pos": wPos,
                    "spacer_start_pos": int(snpPos) + wPos - 1,
                    "spacer_end_pos": int(snpPos) - (spacerLen - wPos),
                    "spacer_size": int(snpPos) + wPos - 1 - (int(snpPos) - (spacerLen - wPos)) + 1,
                    "spacer": spacer,
                    "spacer_window": spacer[wStart - 1 : wEnd],
                    "protospacer": protospacer,
                    "protospacer_window": protospacer[wStart - 1 : wEnd],
                    "pam_seq": rPam,
                    "context": context
                })
                results.append(res)

    return results

def spacer_off_search(results, spacer_len, pam):
    genome_path = "/storage1/sihoo219/web_source/Reference/gmax_noscaff.fa"
    input_path = "/storage1/sihoo219/web_source/RGEN/01.Inputs/input.txt"
    output_path = "-"
    mode = "C"
    max_mismatch = 2

    spacerIdx = 0
    spacers = []
    cmd = ["conda", "run", "-n", "bio_env", "cas-offinder", input_path, mode, output_path]
    cas12_pams = ["TTTN", "TTTV", "TTN"]

    for res in results:
        if pam in cas12_pams:
            spacers.append("N"*len(pam) + res["spacer"] + " " + str(max_mismatch) + " " + str(spacerIdx))
            content = "\n".join([
                genome_path,
                pam + ("N" * spacer_len),
                "\n".join(spacers)
            ])
        else:
            spacers.append(res["spacer"] + "N"*len(pam) + " " + str(max_mismatch) + " " + str(spacerIdx))
            content = "\n".join([
                    genome_path,
                    ("N" * spacer_len) + pam, 
                    "\n".join(spacers)
            ])
        spacerIdx += 1

    with open(input_path, 'w') as fInp:
        fInp.write(content)

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    for line in process.stdout:
        item = line.strip().split("\t")
        if item[-1].isdigit():
            match_res = results[int(item[-1])]
        else:
            break

        if pam in cas12_pams and match_res["spacer_strand"] == item[-3] == "+" and match_res["marker_chr"] == item[1] and str(match_res["spacer_start_pos"] - len(pam)) == item[2]:
            continue
        elif pam in cas12_pams and match_res["spacer_strand"] == item[-3] == "-" and match_res["marker_chr"] == item[1] and str(match_res["spacer_end_pos"]) == item[2]:
            continue
        elif pam not in cas12_pams and match_res["spacer_strand"] == item[-3] == "+" and match_res["marker_chr"] == item[1] and str(match_res["spacer_start_pos"]) == item[2]:
            continue
        elif pam not in cas12_pams and match_res["spacer_strand"] == item[-3] == "-" and match_res["marker_chr"] == item[1] and str(match_res["spacer_end_pos"] - len(pam)) == item[2]:
            continue
        else:
            # match_res["off_candidates"].append(item)
            match_res["off_candidates"].append(item[3].upper())
            if item[-2] == '0':
                match_res["mismatch_0"] += 1
            elif item[-2] == '1':
                match_res["mismatch_1"] += 1
            elif item[-2] == '2':
                match_res["mismatch_2"] += 1

    return results

def spacer_on_scoring(results, pam):
    module_path = "/storage1/sihoo219/web_source/RGEN/02.Scoring/"
    onScore = ""

    for res in results:
        if pam not in ["TTTN", "TTTV", "TTN"]:
            # Cas9 on-target scoring
            exe_path = module_path + "on_cas9/ontarget-score-calculator.py"
            onScore = subprocess.run(["python3", exe_path, res["context"]], capture_output=True, text=True)
        else:
            # Cas12 on-target scoring
            exe_path = module_path + "on_cas12/DeepCpf1_sihoo.py"
            python_bin = "/home/sihoo219/anaconda3/envs/py39_env/bin/python"
            onScore = subprocess.run([python_bin, exe_path, res["context"]], capture_output=True, text=True)

        if onScore.stdout.strip():
            res["on_score"] = round(float(onScore.stdout.strip()), 2)

    return results

def spacer_off_scoring(results, pam):
    module_path = "/storage1/sihoo219/web_source/RGEN/02.Scoring/"

    offScore = ""
    offScores = []

    for res in results:
        for cand in res["off_candidates"]:
            if pam not in ["TTTN", "TTTV", "TTN"] and len(pam) == 3:                        # Cas9, PAM = 3
                exe_path = module_path + "off_cas9/cfd-score-calculator.py"
                cand_ref = res["spacer"] + res["pam_seq"]
                offScore = subprocess.run(["python3", exe_path, "--wt", cand_ref, "--off", cand], capture_output=True, text=True)
            elif pam not in ["TTTN", "TTTV", "TTN"] and len(pam) != 3:                      # Cas9, PAM != 3
                exe_path = module_path + "off_cas9/cfd-score-calculator.py"
                cand_ref = res["spacer"] + res["pam_seq"][0:3]
                cand = cand[0:23]
                offScore = subprocess.run(["python3", exe_path, "--wt", cand_ref, "--off", cand], capture_output=True, text=True)

            if offScore.stdout.strip():
                offScores.append(round(float(offScore.stdout.strip()), 2))

        res["off_score"] = round(reduce(lambda x, y: x * y, offScores), 2)

    return results
    
if __name__ == "__main__":
    selectedParams_json = sys.argv[1]
    selectedParams = json.loads(selectedParams_json)
    # selectedParams = {'markers': [{m1}, {m2}, {m3}], 'length': 20, 'method': 'Crispr', 'nuclease': 'SpCas9', 'checkOnCalc': true, 'checkOffCalc': false;}

    allGuideList = []
    guideDict = {}

    spacerLen = selectedParams['length']
    method = selectedParams['method']
    nuclease = selectedParams['nuclease']
    checkOnCalc = selectedParams['checkOnCalc']
    checkOffCalc = selectedParams['checkOffCalc']

    for marker in selectedParams['markers']:
        if method == "Crispr":
            allGuideList, pamList = handle_Cas(marker, spacerLen, nuclease)
        elif method == "BaseA":
            allGuideList, pamList = handle_ABE(marker, spacerLen, nuclease)
        elif method == "BaseC":
            allGuideList, pamList = handle_CBE(marker, spacerLen, nuclease)
        elif method == "BaseG":
            allGuideList, pamList = handle_GBE(marker, spacerLen, nuclease)

        if not allGuideList: continue

        if checkOnCalc and checkOffCalc:
            allGuideList = spacer_on_scoring(allGuideList, pamList[0])

            allGuideList = spacer_off_search(allGuideList, spacerLen, pamList[0])
            allGuideList = spacer_off_scoring(allGuideList, pamList[0])
        elif checkOnCalc:
            allGuideList = spacer_on_scoring(allGuideList, pamList[0])
        elif checkOffCalc:
            allGuideList = spacer_off_search(allGuideList, spacerLen, pamList[0])
            allGuideList = spacer_off_scoring(allGuideList, pamList[0])

        for guide in allGuideList:
            guide["marker_chr"] = int(guide["marker_chr"][2:])
            guide["spacer_gc"] = round(calculateGCcontent(guide["spacer"]), 2)
            guide["protospacer_gc"] = round(calculateGCcontent(guide["protospacer"]), 2)
        
    print(json.dumps(allGuideList))