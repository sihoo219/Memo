import sys, json
import random
import numpy as np
trait_dict = {
    'ARGININE': 'Seed aginine content',
    'CYSTEINE': 'Seed cysteine content',
    'FLOWERDATE': 'Plant flowering date',
    'HEIGHT': 'Plant height',
    'ISOLEUCINE': 'Seed Isoleucine content',
    'LEUCINE': 'Seed leucine content',
    'LINOLEIC': 'Seed linoleic content',
    'LINOLENIC': 'Seed linolenic content',
    'LYSINE': 'Seed lycine content',
    'LODGING': 'Plant lodging',
    'MATDATE': 'Plant maturity date',
    'METHIONINE': 'Seed methionine content',
    'OIL': 'Seed oil content',
    'OLEIC': 'Seed oleic acid content',
    'PALMITIC': 'Seed palmitic content',
    'PETUREIDE': 'Seed petiole ureide content',
    'PROTEIN': 'Seed protein content',
    'SEEDWEIGHT': 'Seed weight',
    'SHATEARLY': 'Early shattering',
    'SHATLATE': 'Late shattering',
    'STACHYOSE': 'Seed stachyose content',
    'STEARIC': 'Seed stearic content',
    'SUCROSE': 'Seed sucrose content',
    'THREONINE': 'Seed threonine content',
    'TRYPTOPHAN': 'Seed tryptophan content',
    'VALINE': 'Seed valine content',
    'YIELD': 'Plant yield',
    'BP': 'Resistance of bacterial pustule',
    'BRNSTEMROT': 'Resistance of brown stem rot',
    'FLWRCOLOR': 'Flower color',
    'MOTTLING': 'Resistance of seed mottling',
    'NEMATCYST_3': 'Resistance of cyst nematode race 3',
    'NEMATCYST_4': 'Resistance of cyst nematode race 4',
    'NEMATCYST_5': 'Resistance of cyst nematode race 5',
    'PHYTOROT_1': 'Resistance of phytophthora rot race 1',
    'PHYTOROT_17': 'Resistance of phytophthora rot race 17',
    'PHYTOROT_25': 'Resistance of phytophthora rot race 25',
    'PHYTOROT_3': 'Resistance of phytophthora rot race 3',
    'PHYTOROT_7': 'Resistance of phytophthora rot race 7',
    'PUBCOLOR': 'Plant pubescence color',
    'PUBDENSITY': 'Plant pubescence density',
    'SDS': 'Resistance to sudden death syndrome',
    'STEMTERM': 'Plant stem termination type',
    'HILUMCOLOR': 'Seed hilum color',
    'MATGROUP': 'Plant maturity group',
    'PODCOLOR': 'Plant pod color',
    'PUBFORM': 'Plant pubescence form',
    'SCOATCOLOR': 'Seed coat color',
    'SCOATLUST': 'Seed coat luster',
    'SEEDSHAPE1': 'Seed shape',
}

def compute_stats(genotype, v_list):
    if not v_list:
        return []

    q1 = np.percentile(v_list, 25)
    q3 = np.percentile(v_list, 75)
    median = np.median(v_list)
    iqr = q3 - q1
    lower_whisker = max(min(v_list), q1 - 1.5 * iqr)
    upper_whisker = min(max(v_list), q3 + 1.5 * iqr)

    return [lower_whisker, q1, median, q3, upper_whisker]
    
def handle_continuous(input_trait, input_marker):
    sampleDict = {}
    datas = []
    labels = []

    with open('/storage1/sihoo219/web_source/haplotype/traitMeta_noDup.txt', 'r') as fInp1:
        find = 0
        for line1 in fInp1: 
            item1 = line1.strip().split('\t')
            
            if item1[0] == input_trait:
                sampleDict[item1[1]] = {'value': float(item1[2]), 'idx': -1}
                find = 1
            
            if item1[0] != input_trait and find == 1:
                break

    with open('/storage1/sihoo219/web_source/haplotype/G2_totalSNP_nodup.hmp.txt', 'r') as fInp2:
            for line2 in fInp2:
                item2 = line2.strip().split('\t')

                if item2[0] == 'rs#':
                    for sample in sampleDict.keys():
                        if sample in item2:
                            sampleDict[sample]['idx'] = item2.index(sample)    
                    continue

                if item2[0] == input_marker.replace('Gm', 'Chr'):
                    ref, alt = item2[1].split('/')
                    homo1 = ref + ref
                    homo2 = alt + alt
                    hetero = ref + alt
                    labels = [homo1, hetero, homo2]

                    genotypeDict = {homo1: [], hetero: [], homo2: []}

                    for sample in sampleDict.keys():
                        genotype = item2[sampleDict[sample]['idx']]
                        
                        if genotype == homo1:
                            genotypeDict[homo1].append(sampleDict[sample]['value'])
                        elif genotype == homo2:
                            genotypeDict[homo2].append(sampleDict[sample]['value'])
                        elif genotype == (ref + alt) or genotype == (alt + ref):
                            genotypeDict[hetero].append(sampleDict[sample]['value'])
                    
                    break

    for genotype, value_list in genotypeDict.items():
        datas.append(compute_stats(genotype, value_list))

    boxplot_data = {
        "labels": labels,
        "datasets": [
            {
                "backgroundColor": "#90C0FA",
                "borderColor": "#595959",
                "label": trait_dict[input_trait],
                "data": datas
            }
        ]
    }

    return boxplot_data

def handle_binary(input_trait, input_marker):
    sampleDict = {}
    labels, datasets = [], []
    valueTypes = []
    valueCounts = {}
    basicColor = [200, 240, 255]

    with open('/storage1/sihoo219/web_source/haplotype/traitMeta_noDup.txt', 'r') as fInp1:
        find = 0
        for line1 in fInp1: 
            item1 = line1.strip().split('\t')
            
            if item1[0] == input_trait:
                newValue = int(item1[2])
                sampleDict[item1[1]] = {'value': newValue, 'idx': -1}
                find = 1
                if newValue not in valueTypes:
                    valueTypes.append(newValue)
                    color = f"rgb({basicColor[0]}, {basicColor[1]}, {basicColor[2]})"
                    dataset = {
                        'label': str(newValue),
                        'data': [0, 0, 0],
                        'backgroundColor': color
                    }
                    basicColor[0], basicColor[1] = max(0, basicColor[0] - 30), max(0, basicColor[1] - 30)
                    datasets.append(dataset)    # Label Sorting !!!!!!

            if item1[0] != input_trait and find == 1:
                break

    with open('/storage1/sihoo219/web_source/haplotype/G2_totalSNP_nodup.hmp.txt', 'r') as fInp2:
        for line2 in fInp2:
            item2 = line2.strip().split('\t')

            if item2[0] == 'rs#':
                for sample in sampleDict.keys():
                    if sample in item2:
                        sampleDict[sample]['idx'] = item2.index(sample)    
                continue

            if item2[0] == input_marker.replace('Gm', 'Chr'):
                ref, alt = item2[1].split('/')
                homo1, homo2, hetero1, hetero2 = ref + ref, alt + alt, ref + alt, alt + ref
                labels = [homo1, hetero1, homo2] 

                for v in valueTypes:
                    v = str(v)
                    valueCounts[v] = {homo1: 0, hetero1: 0, hetero2: 0, homo2: 0}      # {0: {AA: , AT:, TT: }, 1: {AA:, AT:, TT: }, ...}
                
                for sample, infos in sampleDict.items():
                    genotype = item2[infos['idx']]

                    if genotype not in [homo1, homo2, hetero1, hetero2]:
                        continue

                    valueCounts[str(infos['value'])][genotype] += 1

                break

    for d in datasets:
        d['data'][0] = valueCounts[d['label']][homo1]
        d['data'][1] = valueCounts[d['label']][hetero1] + valueCounts[d['label']][hetero2]
        d['data'][2] = valueCounts[d['label']][homo2]

    barplot_data = {
        'labels': labels,
        'datasets': datasets,
    }

    return barplot_data

# def handle_ordinary(input_trait, input_marker):
#     return 0


if __name__ == "__main__":
    selectedInfo_json = sys.argv[1]
    selectedInfo = json.loads(selectedInfo_json)

    input_trait = selectedInfo['trait']
    input_marker = selectedInfo['marker']
    input_type = selectedInfo['type']   # trait type = continuous, binary, ordinary

    if input_type == 'continuous':
        print(json.dumps(handle_continuous(input_trait, input_marker)))
    elif input_type == 'binary' or input_type == 'ordinary':
        print(json.dumps(handle_binary(input_trait, input_marker)))
    # elif input_type == 'ordinary':
    #     print(json.dumps(handle_ordinary(input_trait, input_marker)))
        