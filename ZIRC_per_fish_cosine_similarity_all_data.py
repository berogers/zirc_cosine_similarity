import csv
import numpy as np
from scipy import spatial
import random

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
'AerocystitisEtiologyUnknown',	'CoelomitisEtiologyUnknown',	'Splenomegaly',	'Bacilli', 'SwimBladder' ]

def each_disease_cosin_into_diff_csv_file(disease_list):

    # dict to pair titles with index number
    title_dict = {0:'PseudolomaNeurophilia', 1:'PleistophoraHyphessobryconis', \
    2:'MycobacteriumSpp', 3:'BacterialInfectionNonAcidFast', 4:'PseudocapillariaTomentosa', \
    5:'OtherHelminths', 6:'MyxidiumStreisingeri', 7:'FungalInfection', 8:'HepaticMegalocytosis', \
    9:'EggAssociatedInflammationOophoritis', 10:'Nephrocalcinosis', 11:'GillLesions', \
    12:'EdwardsiellaIctaluri', 13:'Seminoma', 14:'UltimobranchialAdenomaOrAdenocarcinoma', \
    15:'Chordoma', 16:'SpinalDeformityInHistologicSections', 17:'PiscinoodiniumPillulare', \
    18:'AerocystitisEtiologyUnknown', 19:'CoelomitisEtiologyUnknown', 20: 'Splenomegaly', 21: 'Bacilli', 22:'SwimBladder'}

    # will be populated by vector values (each value will be a list [read left to right]
    # containing the disease's column values from top to bottom)
    disease_dict = {'PseudolomaNeurophilia':[], 'PleistophoraHyphessobryconis':[], \
    'MycobacteriumSpp':[], 'BacterialInfectionNonAcidFast':[], 'PseudocapillariaTomentosa':[],'OtherHelminths':[], \
    'MyxidiumStreisingeri':[], 'FungalInfection':[], 'HepaticMegalocytosis':[], 'EggAssociatedInflammationOophoritis':[], \
    'Nephrocalcinosis':[], 'GillLesions':[], 'EdwardsiellaIctaluri':[], \
    'Seminoma':[], 'UltimobranchialAdenomaOrAdenocarcinoma':[], \
    'Chordoma':[], 'SpinalDeformityInHistologicSections':[], 'PiscinoodiniumPillulare':[],\
    'AerocystitisEtiologyUnknown':[], 'CoelomitisEtiologyUnknown':[], 'Splenomegaly':[], 'Bacilli':[], 'SwimBladder':[]}

    rand_disease_dict = {'PseudolomaNeurophilia':[], 'PleistophoraHyphessobryconis':[], \
    'MycobacteriumSpp':[], 'BacterialInfectionNonAcidFast':[], 'PseudocapillariaTomentosa':[],'OtherHelminths':[], \
    'MyxidiumStreisingeri':[], 'FungalInfection':[], 'HepaticMegalocytosis':[], 'EggAssociatedInflammationOophoritis':[], \
    'Nephrocalcinosis':[], 'GillLesions':[], 'EdwardsiellaIctaluri':[], \
    'Seminoma':[], 'UltimobranchialAdenomaOrAdenocarcinoma':[], \
    'Chordoma':[], 'SpinalDeformityInHistologicSections':[], 'PiscinoodiniumPillulare':[],\
    'AerocystitisEtiologyUnknown':[], 'CoelomitisEtiologyUnknown':[], 'Splenomegaly':[], 'Bacilli':[], 'SwimBladder':[]}

    # dict to pair disease_vectors with index number
    # will be poplated by vector values and used to populate disease_dict
    vector_dict = {0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],12:[],13:[],14:[],15:[],16:[],\
    17:[],18:[],19:[],20:[],21:[],22:[]}

    rand_vector_dict = {0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],12:[],13:[],14:[],15:[],16:[],\
    17:[],18:[],19:[],20:[],21:[],22:[]}


    # will be populated by the result of a cosine similarity test between it and the UAC (disease of interest) vector
    cosin_sim_dict = {'PseudolomaNeurophilia':0, 'PleistophoraHyphessobryconis':0, \
    'MycobacteriumSpp':0, 'BacterialInfectionNonAcidFast':0, 'PseudocapillariaTomentosa':0,'OtherHelminths':0, \
    'MyxidiumStreisingeri':0, 'FungalInfection':0, 'HepaticMegalocytosis':0, 'EggAssociatedInflammationOophoritis':0, \
    'Nephrocalcinosis':0, 'GillLesions':0, 'EdwardsiellaIctaluri':0, \
    'Seminoma':0, 'UltimobranchialAdenomaOrAdenocarcinoma':0, \
    'Chordoma':0, 'SpinalDeformityInHistologicSections':0, 'PiscinoodiniumPillulare':0, \
    'AerocystitisEtiologyUnknown':0, 'CoelomitisEtiologyUnknown':0, 'Splenomegaly':0, 'Bacilli':0, 'SwimBladder':0}

    rand_cosin_sim_dict = {'PseudolomaNeurophilia':0, 'PleistophoraHyphessobryconis':0, \
    'MycobacteriumSpp':0, 'BacterialInfectionNonAcidFast':0, 'PseudocapillariaTomentosa':0,'OtherHelminths':0, \
    'MyxidiumStreisingeri':0, 'FungalInfection':0, 'HepaticMegalocytosis':0, 'EggAssociatedInflammationOophoritis':0, \
    'Nephrocalcinosis':0, 'GillLesions':0, 'EdwardsiellaIctaluri':0, \
    'Seminoma':0, 'UltimobranchialAdenomaOrAdenocarcinoma':0, \
    'Chordoma':0, 'SpinalDeformityInHistologicSections':0, 'PiscinoodiniumPillulare':0, \
    'AerocystitisEtiologyUnknown':0, 'CoelomitisEtiologyUnknown':0, 'Splenomegaly':0, 'Bacilli':0, 'SwimBladder':0}

    # open and read the csv read_file
    read_file = csv.reader(open("hpcolumns_binary_per_fish_all_data.csv", "r", encoding="utf8"), delimiter=",",)
    data_in_list = list(read_file)
    print(len(data_in_list))

    # remove first row (title row) leaving only 34253 rows
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

        # print("vector_dict[0][0:10]")
        # print(vector_dict[0][0:10]

        return vector_dict


    def populate_rand_vector_dict(vector_dict, rand_vector_dict):
        # print("len(rand_vector_dict)")
        # print(len(rand_vector_dict))

        # replicate:
        # vector_dict_values= vector_dict.values()
        # below:

        def take_values_from_keys(vector_dict):

            vector_dict_values = []

            for key in vector_dict:
                value = vector_dict[key]
                vector_dict_values.append(value)

            return vector_dict_values


        vector_dict_values = take_values_from_keys(vector_dict)

        print('len(vector_dict_values)')
        print(len(vector_dict_values))

        for v in range(0,len(vector_dict_values)):
            # print("v {}".format(v))
            vector = vector_dict_values[v]
            # print("vector {} 0:3: {}".format(v,vector[0:3]))
            vec_sum = sum(vector)
            print("vec {}_sum {}".format(v,vec_sum))
            length = len(vector)
            # print('len(vector) {} is {}'.format(v, length))

            rand_indexes = []
            for i in range(0,vec_sum):
                rand_index = random.randint(0,length-1)
                rand_indexes.append(rand_index)

            # print("len(rand_indexes)")
            # print(len(rand_indexes))
            # print('rand_indexes[0:3]')
            # print(rand_indexes[0:3])

            rand_vector = np.zeros(length, dtype=int)

            # print('rand_vector.shape is {}'.format(rand_vector.shape))

            for ri in range(0,len(rand_indexes)):
                rand_index = rand_indexes[ri]
                rand_vector[rand_index] = 1

            print('sum rand_vector: {}'.format(sum(rand_vector)))

            rand_vector_dict[v] = rand_vector

        # print("len(rand_vector_dict) after pop")
        # print(len(rand_vector_dict))
        # print("rand_vector_dict[0][0:10]")
        # print(rand_vector_dict[0][0:10])

        return rand_vector_dict

    # print('len(vector_dict)')
    # print(len(vector_dict))
    # print('len(vector_dict[2])')
    # print(len(vector_dict[2]))
    # print('sum(vector_dict[2])')
    # print(sum(vector_dict[2]))

    print("vector_dict[0][0:10]")
    print(vector_dict[0][0:10])


    vector_dict = populate_vector_dict(data_no_header, vector_dict)

    print("vector_dict[0][0:10] after pop")
    print(vector_dict[0][0:10])


    print("rand_vector_dict[0][0:10]")
    print(rand_vector_dict[0][0:10])


    rand_vector_dict = populate_rand_vector_dict(vector_dict, rand_vector_dict)

    print("rand_vector_dict[0][0:10] after pop")
    print(rand_vector_dict[0][0:10])

    print('len(vector_dict)')
    print(len(vector_dict))
    print('len(vector_dict[2])')
    print(len(vector_dict[2]))
    print('sum(vector_dict[2])')
    print(sum(vector_dict[2]))

    print('len(rand_vector_dict)')
    print(len(rand_vector_dict))
    print('len(rand_vector_dict[2])')
    print(len(rand_vector_dict[2]))
    print('sum(rand_vector_dict[2])')
    print(sum(rand_vector_dict[2]))

    # populate disease_dict
    # says "for every key in title_dict, aka every index num from 0-22
    # take the vector at that key, and populate disease_dict with that vector"

    def populate_disease_dict(title_dict, disease_dict):
        for key in title_dict:
            disease_dict[title_dict[key]] = np.array(vector_dict[key])
        return

    def populate_rand_disease_dict(title_dict, rand_disease_dict):
        for key in title_dict:
            rand_disease_dict[title_dict[key]] = np.array(rand_vector_dict[key])
        return

    populate_disease_dict(title_dict, disease_dict)
    populate_rand_disease_dict(title_dict, rand_disease_dict)


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

    def run_rand_cosin_sim(disease_list, rand_disease_dict, rand_cosin_sim_dict, rand_disease_vec):
        # run cosine similarity tests between each disease and UAC
        for disease in disease_list:
            rand_pair_vec = rand_disease_dict[disease]
            rand_cosin_sim = round((1 - spatial.distance.cosine(rand_disease_vec, rand_pair_vec)),5)
            rand_cosin_sim_dict[disease] = rand_cosin_sim
        return

    def sort_cosin_sim_dict(cosin_sim_dict, disease):

        cosin_sim_values = cosin_sim_dict.values()
        cosin_sim_keys = cosin_sim_dict.keys()
        sorted_cosin_sim_values = sorted(cosin_sim_values)

        # print(sorted_cosin_sim_values)
        sorted_cosin_sim_dict = {}

        for value in sorted_cosin_sim_values:
            for disease in cosin_sim_dict:
                if cosin_sim_dict[disease] == value:
                    sorted_cosin_sim_dict[disease] = cosin_sim_dict[disease]

        # print(sorted_cosin_sim_dict)

        sorted_cosin_sim_keys = sorted_cosin_sim_dict.keys()

        return sorted_cosin_sim_dict, sorted_cosin_sim_keys, sorted_cosin_sim_values, cosin_sim_keys, cosin_sim_values

    def sort_rand_cosin_sim_dict(rand_cosin_sim_dict, disease):

        rand_cosin_sim_values = rand_cosin_sim_dict.values()
        rand_cosin_sim_keys = rand_cosin_sim_dict.keys()
        sorted_rand_cosin_sim_values = sorted(rand_cosin_sim_values)

        # print(sorted_rand_cosin_sim_values)
        sorted_rand_cosin_sim_dict = {}

        for value in sorted_rand_cosin_sim_values:
            for disease in rand_cosin_sim_dict:
                if rand_cosin_sim_dict[disease] == value:
                    sorted_rand_cosin_sim_dict[disease] = rand_cosin_sim_dict[disease]

        # print(sorted_rand_cosin_sim_dict)

        sorted_rand_cosin_sim_keys = sorted_rand_cosin_sim_dict.keys()

        return sorted_rand_cosin_sim_dict, sorted_rand_cosin_sim_keys, sorted_rand_cosin_sim_values, rand_cosin_sim_keys, rand_cosin_sim_values

    for disease in disease_list:
        disease_vec = disease_dict[disease]

        run_cosin_sim(disease_list, disease_dict, cosin_sim_dict, disease_vec)
        run_rand_cosin_sim(disease_list, rand_disease_dict, rand_cosin_sim_dict, disease_vec)

        # cosin_sim_list = sorted(cosin_sim_dict.items(), key=lambda x: x[1], reverse=True)
        # print(cosin_sim_list)

        sorted_cosin_sim_dict, sorted_cosin_sim_keys, sorted_cosin_sim_values, cosin_sim_keys, cosin_sim_values = sort_cosin_sim_dict(cosin_sim_dict, disease)
        sorted_rand_cosin_sim_dict, sorted_rand_cosin_sim_keys, sorted_rand_cosin_sim_values, rand_cosin_sim_keys, rand_cosin_sim_values = sort_rand_cosin_sim_dict(rand_cosin_sim_dict, disease)

        # # calculate difference between real cosine sim value and cosine sim value between random vectors (negative control)
        #
        # sorted_differences = []
        # unsorted_differences = []
        #
        # for i in range(0,len(sorted_cosine_sim_values)):
        #     value = sorted_cosin_sim_values[i]
        #     rand_value = sorted_rand_cosin_sim_values[i]
        #     diff = value - sorted_rand-cosine_sim_values[


        # write to csv

        create_empty_csv = open('{}_per_fish_all_data_cosin_sim_results+randomvecs_onlypairisrandom.csv'.format(disease), 'w')

        write_file = open('{}_per_fish_all_data_cosin_sim_results+randomvecs_onlypairisrandom.csv'.format(disease), 'w+', newline ='')

        with write_file:
            write = csv.writer(write_file)
            write.writerow(sorted_cosin_sim_values)
            write.writerow(sorted_rand_cosin_sim_values)
            write.writerow(sorted_cosin_sim_keys)
            write.writerow(cosin_sim_keys)
            write.writerow(cosin_sim_values)
            write.writerow(rand_cosin_sim_values)


    return

# run the whole thing
each_disease_cosin_into_diff_csv_file(disease_list)

print('Finished')
