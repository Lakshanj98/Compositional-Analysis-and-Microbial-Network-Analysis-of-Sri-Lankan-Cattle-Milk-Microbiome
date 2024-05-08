generate tax table for networks.py
otu_table_list = ["full net otu_table.tsv", "high c otu_table.tsv", "medium c otu_table.tsv",
                  "low c otu_table.tsv", "large net otu_table.tsv", "med small otu_table.tsv"]

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
    file_name = file_name[0] + " " + file_name[1]
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
