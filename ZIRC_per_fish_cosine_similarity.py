import csv
import numpy as np
from scipy import spatial

# Initiate empty lists and dictionaries

# Column names originally from Risa Askerooth's disease_columns_mutate_May17.csv
# once tranferred exclusively to hp_columns_binary_per_fish.csv

# Names of columns 1-22 (index 0-21) in order:
disease_list = [ 'PseudolomaNeurophilia',	'PleistophoraHyphessobryconis',	'MycobacteriumSpp',	\
'BacterialInfectionNonAcidFast',	'PseudocapillariaTomentosa',	'OtherHelminths',	\
'MyxidiumStreisingeri',	'FungalInfection',	'HepaticMegalocytosis',	\
'EggAssociatedInflammationOophoritis',	'Nephrocalcinosis',	'GillLesions',	\
'EdwardsiellaIctaluri',	'Seminoma',	'UltimobranchialAdenomaOrAdenocarcinoma',	\
'Chordoma',	'SpinalDeformityInHistologicSections',	'PiscinoodiniumPillulare',	\
'AerocystitisEtiologyUnknown',	'CoelomitisEtiologyUnknown',	'Splenomegaly',	'Bacilli' ]

def each_disease_cosin_into_diff_csv_file(disease_list):

    # dict to pair titles with index number
    title_dict = {0:'PseudolomaNeurophilia', 1:'PleistophoraHyphessobryconis', \
    2:'MycobacteriumSpp', 3:'BacterialInfectionNonAcidFast', 4:'PseudocapillariaTomentosa', \
    5:'OtherHelminths', 6:'MyxidiumStreisingeri', 7:'FungalInfection', 8:'HepaticMegalocytosis', \
    9:'EggAssociatedInflammationOophoritis', 10:'Nephrocalcinosis', 11:'GillLesions', \
    12:'EdwardsiellaIctaluri', 13:'Seminoma', 14:'UltimobranchialAdenomaOrAdenocarcinoma', \
    15:'Chordoma', 16:'SpinalDeformityInHistologicSections', 17:'PiscinoodiniumPillulare', \
    18:'AerocystitisEtiologyUnknown', 19:'CoelomitisEtiologyUnknown', 20: 'Splenomegaly', 21: 'Bacilli'}

    # will be populated by vector values (each value will be a list [read left to right]
    # containing the disease's column values from top to bottom)
    disease_dict = {'PseudolomaNeurophilia':[], 'PleistophoraHyphessobryconis':[], \
    'MycobacteriumSpp':[], 'BacterialInfectionNonAcidFast':[], 'PseudocapillariaTomentosa':[],'OtherHelminths':[], \
    'MyxidiumStreisingeri':[], 'FungalInfection':[], 'HepaticMegalocytosis':[], 'EggAssociatedInflammationOophoritis':[], \
    'Nephrocalcinosis':[], 'GillLesions':[], 'EdwardsiellaIctaluri':[], \
    'Seminoma':[], 'UltimobranchialAdenomaOrAdenocarcinoma':[], \
    'Chordoma':[], 'SpinalDeformityInHistologicSections':[], 'PiscinoodiniumPillulare':[],\
    'AerocystitisEtiologyUnknown':[], 'CoelomitisEtiologyUnknown':[], 'Splenomegaly':[], 'Bacilli':[]}

    # dict to pair disease_vectors with index number
    # will be poplated by vector values and used to populate disease_dict
    vector_dict = {0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],12:[],13:[],14:[],15:[],16:[],\
    17:[],18:[],19:[],20:[],21:[]}


    # will be populated by the result of a cosine similarity test between it and the UAC (disease of interest) vector
    cosin_sim_dict = {'PseudolomaNeurophilia':0, 'PleistophoraHyphessobryconis':0, \
    'MycobacteriumSpp':0, 'BacterialInfectionNonAcidFast':0, 'PseudocapillariaTomentosa':0,'OtherHelminths':0, \
    'MyxidiumStreisingeri':0, 'FungalInfection':0, 'HepaticMegalocytosis':0, 'EggAssociatedInflammationOophoritis':0, \
    'Nephrocalcinosis':0, 'GillLesions':0, 'EdwardsiellaIctaluri':0, \
    'Seminoma':0, 'UltimobranchialAdenomaOrAdenocarcinoma':0, \
    'Chordoma':0, 'SpinalDeformityInHistologicSections':0, 'PiscinoodiniumPillulare':0, \
    'AerocystitisEtiologyUnknown':0, 'CoelomitisEtiologyUnknown':0, 'Splenomegaly':0, 'Bacilli':0}

    # open and read the csv read_file
    read_file = csv.reader(open("hpcolumns_binary_per_fish.csv", "r", encoding="utf8"), delimiter=",",)
    data_in_list = list(read_file)
    print(len(data_in_list))

    # remove first row (titles) leaving only 11020 rows
    # each row corresponds to a specific caseID of binary data
    data_no_header = data_in_list[1:]
    print(len(data_no_header))

    single_row = data_no_header[0]

    num_columns = len(single_row)

    # populate vector_dict by actually grabbing csv data and putting them in lists containing column contents
    def populate_vector_dict(data_no_header, vector_dict):
        for row in data_no_header:
            for column in range(0,num_columns):
                column_data_pt = int(row[column])
                vector_dict[column].append(column_data_pt)
        return


    populate_vector_dict(data_no_header, vector_dict)

    # populate disease_dict
    # says "for every key in title_dict, aka every index num from 0-22
    # take the vector at that key, and populate disease_dict with that vector"

    def populate_disease_dict(title_dict, disease_dict):
        for key in title_dict:
            disease_dict[title_dict[key]] = np.array(vector_dict[key])
        return

    populate_disease_dict(title_dict, disease_dict)


    # case and edge case tests
    # print(type(UAC_vec))
    # print(UAC_vec[-1])
    # print(disease_dict['FungalInfection'][-1])
    # print(disease_dict['MycobacteriumSpp'][-1])
    # print(disease_dict['CoelomitisEtiologyUnknown'][-1])


    def run_cosin_sim(disease_list, disease_dict, cosin_sim_dict, disease_vec):
        # run cosine similarity tests between each disease and UAC
        for disease in disease_list:
            pair_vec = disease_dict[disease]
            cosin_sim = round((1 - spatial.distance.cosine(disease_vec, pair_vec)),5)
            cosin_sim_dict[disease] = cosin_sim
        return

    heat_map_matrix = []
    
    for disease in disease_list:
        disease_vec = disease_dict[disease]

        run_cosin_sim(disease_list, disease_dict, cosin_sim_dict, disease_vec)

        # cosin_sim_list = sorted(cosin_sim_dict.items(), key=lambda x: x[1], reverse=True)
        # print(cosin_sim_list)

        def sort_cosin_sim_dict(cosin_sim_dict, disease):

            cosin_sim_values = cosin_sim_dict.values()
            cosin_sim_keys = cosin_sim_dict.keys()
            sorted_cosin_sim_values = sorted(cosin_sim_values)

            print(sorted_cosin_sim_values)
            sorted_cosin_sim_dict = {}

            for value in sorted_cosin_sim_values:
                for disease in cosin_sim_dict:
                    if cosin_sim_dict[disease] == value:
                        sorted_cosin_sim_dict[disease] = cosin_sim_dict[disease]

            print(sorted_cosin_sim_dict)

            sorted_cosin_sim_keys = sorted_cosin_sim_dict.keys()

            return sorted_cosin_sim_dict, sorted_cosin_sim_keys, sorted_cosin_sim_values, cosin_sim_keys, cosin_sim_values

        sorted_cosin_sim_dict, sorted_cosin_sim_keys, sorted_cosin_sim_values, cosin_sim_keys, cosin_sim_values = sort_cosin_sim_dict(cosin_sim_dict, disease)

        heat_map_matrix.append(list(cosin_sim_values))
        
        # write to csv

        # NOTE: these csv files have to be pre-created to be populated, but they will be populated automatically once so.

        write_file = open('{}_per_fish_cosin_sim_results.csv'.format(disease), 'w+', newline ='')

        with write_file:
            write = csv.writer(write_file)
            write.writerow(sorted_cosin_sim_values)
            write.writerow(sorted_cosin_sim_keys)
            write.writerow(cosin_sim_keys)
            write.writerow(cosin_sim_values)

    return

# run the whole thing
each_disease_cosin_into_diff_csv_file(disease_list)
