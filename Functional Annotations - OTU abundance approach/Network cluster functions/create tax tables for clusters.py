create tax tables for clusters.py
# genera_otus = {}
#
# # write the genus and OTU ids to a dictionary
# with open("taxonomy table.csv", "r") as file:
#     for line in file:
#         line = line.strip("\n")
#         line = line.split(",")
#         # print(line)
#         if line[6] != "NA" and line[0] != "":
#             if line[6] not in genera_otus.keys():
#                 key = line[6]
#                 # initiate a list and put the value inside it at the same time.
#                 values = [line[0]]
#
#                 genera_otus[key] = values
#             else:
#                 genera_otus[line[6]].append(line[0])
#
# # remove duplicate OTUs from the dictionary
# for key in genera_otus:
#     genera_otus[key] = list(set(genera_otus[key]))
#
# # print(len(genera_otus))
#
# # file list
# file_list = ["cl1_75_95_Sparcc node data.csv", "cl2_55_75_Sparcc node data.csv",
#              "cl3_35_55_Sparcc node data.csv", "Large_Sparcc node data.csv", "Med_Small_Sparcc node data.csv"]
#
# otu_table_list = ["high_c_otu_table.tsv", "medium_c_otu_table.tsv", "low_c_otu_table.tsv", "large_otu_table.tsv",
#                   "medium_c_otu_table.tsv"]
# ############## -------- This code until line 105 works on generating cluster_function_differentOTUdoingFunction.
# for file, table in zip(file_list, otu_table_list):
#     # write the otu and cluster
#     cluster_otu = {}
#     with open(file, "r") as file_1:
#         for line in file_1:
#             line = line.strip("\n")
#             line = line.split(",")
#
#             if line[0] != '"index"':
#                 for key in genera_otus:
#                     otu_list = genera_otus[key]
#                     if key == line[1]:
#
#                         for otu in otu_list:
#
#                             if line[2] not in cluster_otu.keys():
#                                 cluster = line[2]
#                                 # initiate a list and put the value inside it at the same time.
#                                 otu_ids = [otu]
#                                 # print(values)
#                                 cluster_otu[cluster] = otu_ids
#                             else:
#                                 #
#                                 cluster_otu[line[2]].append(otu)
#
#     # # print the dictionary
#     # for key in cluster_otu:
#     #     print(key, cluster_otu[key])
#     #     for line in file_2:
#     #         line = line.strip("\n")
#     #         line = line.
#
#     # create the otu_taxonomy dictionary
#     otu_taxonomy = {}
#     with open("taxonomy table.csv", "r") as file_1:
#         for line in file_1:
#             line = line.strip("\n")
#             line = line.split(",")
#             # print(line)
#             # remove the column names
#             if line[1] != "Kingdom":
#                 if line[-1] == '':
#                     end = line.index('')
#                     print(end)
#                     otu_taxonomy[line[0]] = line[1:end]
#                 else:
#                     otu_taxonomy[line[0]] = line[1:]
#
#     # create the otu table for every cluster in every network
#     for cluster in cluster_otu:
#         file_name = file.split("_")[0]
#         with open(file_name + " cluster " + str(cluster) + " tax table.tsv", "w") as file_3:
#             for key in otu_taxonomy:
#                 if key in cluster_otu[cluster]:
#                     otu = str(key)
#
#                     values = otu_taxonomy[key]
#                     print(values)
#                     rank_length = len(values)
#                     print(rank_length)
#
#                     values = str(otu_taxonomy[key]).replace("\'", "")
#                     values = values.replace("[", "")
#                     values = values.replace("]", "")
#                     values = values.replace(" ", "")
#                     values = values.split(",")
#                     # print(values)
#                     if rank_length == 7:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#                         ranks1 = str(values[1])
#                         ranks2 = str(values[2])
#                         ranks3 = str(values[3])
#                         ranks4 = str(values[4])
#                         ranks5 = str(values[5])
#                         ranks6 = str(values[6])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
#                                           ";f_" + ranks4 + ";g_" + ranks5 + ";s_" + ranks6 + "\n")
#
#                     elif rank_length == 6:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#                         ranks1 = str(values[1])
#                         ranks2 = str(values[2])
#                         ranks3 = str(values[3])
#                         ranks4 = str(values[4])
#                         ranks5 = str(values[5])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
#                                           ";f_" + ranks4 + ";g_" + ranks5 + "\n")
#
#                     elif rank_length == 5:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#                         ranks1 = str(values[1])
#                         ranks2 = str(values[2])
#                         ranks3 = str(values[3])
#                         ranks4 = str(values[4])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
#                                           ";f_" + ranks4 + "\n")
#
#                     elif rank_length == 4:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#                         ranks1 = str(values[1])
#                         ranks2 = str(values[2])
#                         ranks3 = str(values[3])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
#                                           "\n")
#
#                     elif rank_length == 3:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#                         ranks1 = str(values[1])
#                         ranks2 = str(values[2])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + "\n")
#
#                     elif rank_length == 2:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#                         ranks1 = str(values[1])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + "\n")
#
#                     elif rank_length == 1:
#                         # format the necessary data
#                         ranks0 = str(values[0])
#
#                         file_3.writelines(otu + "\t" + "k_" + ranks0 + "\n")

otu_table_list = ["full cluster 1 otu_table.tsv", "full cluster 2 otu_table.tsv", "full cluster 3 otu_table.tsv",
                  "full cluster 4 otu_table.tsv", "full cluster 5 otu_table.tsv", "full cluster 6 otu_table.tsv",
                  "full cluster 7 otu_table.tsv", "full cluster 8 otu_table.tsv", "full cluster 9 otu_table.tsv",
                  "cl1 cluster 1 otu_table.tsv", "cl1 cluster 2 otu_table.tsv", "cl1 cluster 3 otu_table.tsv",
                  "cl1 cluster 4 otu_table.tsv", "cl2 cluster 1 otu_table.tsv", "cl2 cluster 2 otu_table.tsv",
                  "cl2 cluster 3 otu_table.tsv", "cl2 cluster 4 otu_table.tsv", "cl2 cluster 5 otu_table.tsv",
                  "cl2 cluster 6 otu_table.tsv", "cl2 cluster 7 otu_table.tsv", "cl3 cluster 1 otu_table.tsv",
                  "cl3 cluster 2 otu_table.tsv", "cl3 cluster 3 otu_table.tsv", "cl3 cluster 4 otu_table.tsv",
                  "Large cluster 1 otu_table.tsv", "Large cluster 2 otu_table.tsv", "Large cluster 3 otu_table.tsv",
                  "Large cluster 4 otu_table.tsv", "Large cluster 5 otu_table.tsv", "Large cluster 6 otu_table.tsv",
                  "Large cluster 7 otu_table.tsv", "Med cluster 1 otu_table.tsv", "Med cluster 2 otu_table.tsv",
                  "Med cluster 3 otu_table.tsv", "Med cluster 4 otu_table.tsv", "Med cluster 5 otu_table.tsv",
                  "Med cluster 6 otu_table.tsv", "Med cluster 7 otu_table.tsv", "Med cluster 8 otu_table.tsv"]

# create the otu_taxonomy dictionary
otu_taxonomy = {}
with open("taxonomy table.csv", "r") as file_1:
    for line in file_1:
        line = line.strip("\n")
        line = line.split(",")
        # print(line)
        # remove the column names
        if line[1] != "Kingdom":
            if line[-1] == '':
                end = line.index('')
                print(end)
                otu_taxonomy[line[0]] = line[1:end]
            else:
                otu_taxonomy[line[0]] = line[1:]

for file in otu_table_list:
    otu_data = []
    # create the otu table for every cluster in every network
    with open("input/" + file, "r") as file_4:
        for line in file_4:
            line = line.strip("\n")
            line = line.split("\t")

            if line[0] != "OTU_ID":
                otu_data.append(line[0])

    file_name = file.split("_")[0]
    file_name = file_name.split(" ")
    file_name = file_name[0] + " " + file_name[1] + " " + file_name[2]
    with open("output/" + file_name + " tax table.tsv", "w") as file_3:
        for key in otu_taxonomy:
            if key in otu_data:
                otu = str(key)

                values = otu_taxonomy[key]
                print(values)
                rank_length = len(values)
                print(rank_length)

                values = str(otu_taxonomy[key]).replace("\'", "")
                values = values.replace("[", "")
                values = values.replace("]", "")
                values = values.replace(" ", "")
                values = values.split(",")
                # print(values)
                if rank_length == 7:
                    # format the necessary data
                    ranks0 = str(values[0])
                    ranks1 = str(values[1])
                    ranks2 = str(values[2])
                    ranks3 = str(values[3])
                    ranks4 = str(values[4])
                    ranks5 = str(values[5])
                    ranks6 = str(values[6])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
                                      ";f_" + ranks4 + ";g_" + ranks5 + ";s_" + ranks6 + "\n")

                elif rank_length == 6:
                    # format the necessary data
                    ranks0 = str(values[0])
                    ranks1 = str(values[1])
                    ranks2 = str(values[2])
                    ranks3 = str(values[3])
                    ranks4 = str(values[4])
                    ranks5 = str(values[5])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
                                      ";f_" + ranks4 + ";g_" + ranks5 + "\n")

                elif rank_length == 5:
                    # format the necessary data
                    ranks0 = str(values[0])
                    ranks1 = str(values[1])
                    ranks2 = str(values[2])
                    ranks3 = str(values[3])
                    ranks4 = str(values[4])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
                                      ";f_" + ranks4 + "\n")

                elif rank_length == 4:
                    # format the necessary data
                    ranks0 = str(values[0])
                    ranks1 = str(values[1])
                    ranks2 = str(values[2])
                    ranks3 = str(values[3])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + ";o_" + ranks3 +
                                      "\n")

                elif rank_length == 3:
                    # format the necessary data
                    ranks0 = str(values[0])
                    ranks1 = str(values[1])
                    ranks2 = str(values[2])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + ";c_" + ranks2 + "\n")

                elif rank_length == 2:
                    # format the necessary data
                    ranks0 = str(values[0])
                    ranks1 = str(values[1])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + ";p_" + ranks1 + "\n")

                elif rank_length == 1:
                    # format the necessary data
                    ranks0 = str(values[0])

                    file_3.writelines(otu + "\t" + "k_" + ranks0 + "\n")
