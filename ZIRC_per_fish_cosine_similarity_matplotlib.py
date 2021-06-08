import csv
import numpy as np
from scipy import spatial
import random
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Initiate empty lists and dictionaries

# Column names originally from Risa Askerooth's disease_columns_mutate_May17.csv
# once tranferred exclusively to hp_columns_binary_per_fish.csv

# Names of columns 1-22 (index 0-21) in order:
disease_list = [ 'PseudolomaNeurophilia',	'MycobacteriumSpp',	\
'BacterialInfectionNonAcidFast',	'PseudocapillariaTomentosa',	'OtherHelminths',	\
'MyxidiumStreisingeri',	'FungalInfection',	'HepaticMegalocytosis',	\
'EggAssociatedInflammationOophoritis',	'Nephrocalcinosis',	'GillLesions',	'Seminoma',	'UltimobranchialAdenomaOrAdenocarcinoma',	\
'Chordoma', 'SpinalDeformityInHistologicSections','AerocystitisEtiologyUnknown','CoelomitisEtiologyUnknown',	'Splenomegaly',	'Bacilli', "SwimBladder" ]

def each_disease_cosin_into_diff_csv_file(disease_list):

    # dict to pair titles with index number
    title_dict = {0:'PseudolomaNeurophilia', \
    1:'MycobacteriumSpp', 2:'BacterialInfectionNonAcidFast', 3:'PseudocapillariaTomentosa', \
    4:'OtherHelminths', 5:'MyxidiumStreisingeri', 6:'FungalInfection', 7:'HepaticMegalocytosis', \
    8:'EggAssociatedInflammationOophoritis', 9:'Nephrocalcinosis', 10:'GillLesions', 11:'Seminoma', 12:'UltimobranchialAdenomaOrAdenocarcinoma', \
    13:'Chordoma', 14:'SpinalDeformityInHistologicSections', 15:'AerocystitisEtiologyUnknown', 16:'CoelomitisEtiologyUnknown', 17: 'Splenomegaly', 18: 'Bacilli',19:'SwimBladder'}

    # will be populated by vector values (each value will be a list [read left to right]
    # containing the disease's column values from top to bottom)
    disease_dict = {'PseudolomaNeurophilia':[], \
    'MycobacteriumSpp':[], 'BacterialInfectionNonAcidFast':[], 'PseudocapillariaTomentosa':[],'OtherHelminths':[], \
    'MyxidiumStreisingeri':[], 'FungalInfection':[], 'HepaticMegalocytosis':[], 'EggAssociatedInflammationOophoritis':[], \
    'Nephrocalcinosis':[], 'GillLesions':[], \
    'Seminoma':[], 'UltimobranchialAdenomaOrAdenocarcinoma':[], \
    'Chordoma':[], 'SpinalDeformityInHistologicSections':[],'AerocystitisEtiologyUnknown':[], 'CoelomitisEtiologyUnknown':[], 'Splenomegaly':[], 'Bacilli':[], 'SwimBladder':[]}

    # dict to pair disease_vectors with index number
    # will be poplated by vector values and used to populate disease_dict
    vector_dict = {0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],12:[],13:[],14:[],15:[],16:[],\
    17:[],18:[],19:[]}


    # will be populated by the result of a cosine similarity test between it and the UAC (disease of interest) vector
    cosin_sim_dict = {'PseudolomaNeurophilia':0, \
    'MycobacteriumSpp':0, 'BacterialInfectionNonAcidFast':0, 'PseudocapillariaTomentosa':0,'OtherHelminths':0, \
    'MyxidiumStreisingeri':0, 'FungalInfection':0, 'HepaticMegalocytosis':0, 'EggAssociatedInflammationOophoritis':0, \
    'Nephrocalcinosis':0, 'GillLesions':0, \
    'Seminoma':0, 'UltimobranchialAdenomaOrAdenocarcinoma':0, \
    'Chordoma':0, 'SpinalDeformityInHistologicSections':0,'AerocystitisEtiologyUnknown':0, 'CoelomitisEtiologyUnknown':0, 'Splenomegaly':0, 'Bacilli':0, 'SwimBladder':0}

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
        
    disease_list = [ 'Pseudoloma Neurophilia',	'Mycobacterium Spp',	\
    'Bacterial Infection Non Acid Fast',	'Pseudocapillaria Tomentosa',	'Other Helminths',	\
    'Myxidium Streisingeri',	'Fungal Infection',	'Hepatic Megalocytosis',	\
    'Egg-Associated Inflammation Oophoritis',	'Nephrocalcinosis',	'Gill Lesions',	'Seminoma',	'Ultimobranchial Adenoma or Adenocarcinoma',	\
    'Chordoma',	'Spinal Deformity in Histologic Sections', "Aerocystitis Etiology Unknown",	'Coelomitis Etiology Unknown',	'Splenomegaly',	'Bacilli', 'Swim Bladder Infection' ]
       
    
    heat_map_matrix_np = np.array(heat_map_matrix)
    
    null_sim_score = 0.014
    
    # def subtract_null_sim_score(heat_map_matrix, null_sim_score):
    #     for row in heat_map_matrix:
    #         for column in row:
    #             column = column - null_sim_score
        
    #     return
                
    # subtract_null_sim_score(heat_map_matrix, null_sim_score)
    
    # heat_map_matrix = heat_map_matrix - null_sim_score
    
    for x in range(0,len(heat_map_matrix)):
        row = heat_map_matrix[x]
        print('row is')
        print(row)
        for y in range(0,len(row)):
            element = row[y]
            print('element is')
            print(element)
            subtracted = element - null_sim_score
            print('subtracted is ')
            print(subtracted)
            if subtracted < 0:
                subtracted += -1*(subtracted)
                
            if subtracted > 0.98:
                subtracted += null_sim_score
            
            heat_map_matrix[x][y] = subtracted
            
    print("heat_map_matrix[-2]")  
    print(heat_map_matrix[-2])

    # print('heat_map_matrix')
    # print(heat_map_matrix)
    # print('type(heat_map_matrix)')
    # print(type(heat_map_matrix))
    # print("heat_map_matrix_np")
    # print(heat_map_matrix_np)
    # print('type(heat_map_matrix_np)')
    # print(type(heat_map_matrix_np))
    # print('shape heat_map_matrix_np')
    # print(heat_map_matrix_np.shape)

    fig, ax = plt.subplots()
    im = ax.imshow(heat_map_matrix)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(disease_list)))
    ax.set_yticks(np.arange(len(disease_list)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(disease_list, fontsize=4)
    ax.set_yticklabels(disease_list,fontsize=4)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
              rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(disease_list)):
        for j in range(len(disease_list)):
            text = ax.text(j, i, round(heat_map_matrix[i][j],3),
                            ha="center", va="center", color="w", fontsize=2)

    ax.set_title("Cosine similarity between zebrafish condition pairs"+'\n'+
                 "2002-2021 per fish data",fontsize=6)
    fig.tight_layout()
    
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize=4)
    
    plt.savefig('cosine_sim_heatmap_2002-2021_per_fish_data.png',dpi=1200)
    plt.show()
    
    # sns.heatmap(heat_map_matrix, annot=True, linewidths=.5)

    
    return

# run the whole thing
each_disease_cosin_into_diff_csv_file(disease_list)
