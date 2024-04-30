create otu tables for clusters.py
genera_otus = {}

# write the genus and OTU ids to a dictionary
with open("taxonomy table.csv", "r") as file:
    for line in file:
        line = line.strip("\n")
        line = line.split(",")
        # print(line)
        if line[6] != "" and line[0] != "":
            if line[6] not in genera_otus.keys():
                key = line[6]
                # initiate a list and put the value inside it at the same time.
                values = [line[0]]

                genera_otus[key] = values
            else:
                genera_otus[line[6]].append(line[0])

# remove duplicate OTUs from the dictionary
for key in genera_otus:
    genera_otus[key] = list(set(genera_otus[key]))

# print(len(genera_otus))

# file list
file_list = ["full_network_Sparcc node data.csv","cl1_75_95_Sparcc node data.csv", "cl2_55_75_Sparcc node data.csv",
             "cl3_35_55_Sparcc node data.csv", "Large_Sparcc node data.csv", "Med_Small_Sparcc node data.csv"]

otu_table_list = ["full_otu_table.tsv","high_c_otu_table.tsv", "medium_c_otu_table.tsv", "low_c_otu_table.tsv", "large_otu_table.tsv",
                  "med_small_otu_table.tsv"]
############## -------- This code until line 105 works on generating cluster_function_differentOTUdoingFunction.
for file, table in zip(file_list, otu_table_list):
    # write the otu and cluster
    cluster_otu = {}
    with open(file, "r") as file_1:
        for line in file_1:
            line = line.strip("\n")
            line = line.split(",")

            if line[0] != '"index"':
                for key in genera_otus:
                    otu_list = genera_otus[key]
                    if key == line[1]:

                        for otu in otu_list:

                            if line[2] not in cluster_otu.keys():
                                cluster = line[2]
                                # initiate a list and put the value inside it at the same time.
                                otu_ids = [otu]
                                # print(values)
                                cluster_otu[cluster] = otu_ids
                            else:
                                #
                                cluster_otu[line[2]].append(otu)

    # # print the dictionary
    # for key in cluster_otu:
    #     print(key, cluster_otu[key])
    #     for line in file_2:
    #         line = line.strip("\n")
    #         line = line.

    # read the relevant otu table file and create the otu_abundance dictionary
    otu_abundance = {}
    column_names = []
    header = []
    with open(table, "r") as file_2:
        for line in file_2:

            line = line.strip("\n")

            line = line.split("\t")

            if line[0] != "OTU_ID":
                key = str(line[0]).replace("\'", "")
                key = key.replace("\"", "")

                values = str(line[1:]).replace("\'", "")
                values = values.replace("[", "")
                values = values.replace("]", "")
                otu_abundance[key] = values
            else:
                # global header
                header.append(line[1:])
                # print(line[1:-1])

    #remove duplicate otus from dictionary
    cluster_otu_unique = {}
    for cluster in cluster_otu:
        cluster_otu[cluster] = list(set(cluster_otu[cluster]))
    # cluster_otu_unique = {}
    # for i in range(len(cluster_otu)):  # Use a numerical index to iterate through the list
    #     cluster_otu[i] = list(set(cluster_otu[i]))

        # create the otu table for every cluster in every network
    for cluster in cluster_otu:
        file_name = file.split("_")[0]
        with open(file_name + " cluster " + str(cluster) + " otu_table.tsv", "w") as file_3:
            # remove list marks from column_names

            header = str(header).replace("[", "")
            header = header.replace("]", "")
            # print(header)
            header = header.replace("\"", "")
            header = header.replace("'", "")
            header = header.replace(",", "\t")
            print(header)
            file_3.write("OTU_ID" + "\t" + header + "\n")

            for otu in cluster_otu[cluster]:
                if otu in otu_abundance.keys():
                    file_3.writelines(otu + "\t" + str(otu_abundance[otu]).replace(",","\t") + "\n")
