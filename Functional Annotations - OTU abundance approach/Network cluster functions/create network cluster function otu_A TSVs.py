create network cluster function otu_A TSVs.py
# file lists for networks

full_list = ["full cluster 1_faprotax.csv", "full cluster 2_faprotax.csv", "full cluster 3_faprotax.csv",
             "full cluster 4_faprotax.csv", "full cluster 5_faprotax.csv", "full cluster 6_faprotax.csv",
             "full cluster 7_faprotax.csv", "full cluster 8_faprotax.csv", "full cluster 9_faprotax.csv"]

cl1_list = ["cl1 cluster 1_faprotax.csv", "cl1 cluster 2_faprotax.csv", "cl1 cluster 3_faprotax.csv",
            "cl1 cluster 4_faprotax.csv"]

cl2_list = ["cl2 cluster 1_faprotax.csv", "cl2 cluster 2_faprotax.csv", "cl2 cluster 3_faprotax.csv",
            "cl2 cluster 4_faprotax.csv", "cl2 cluster 5_faprotax.csv", "cl2 cluster 6_faprotax.csv",
            "cl2 cluster 7_faprotax.csv"]

cl3_list =["cl3 cluster 1_faprotax.csv", "cl3 cluster 2_faprotax.csv", "cl3 cluster 3_faprotax.csv",
           "cl3 cluster 4_faprotax.csv"]

large_list = ["large cluster 1_faprotax.csv", "large cluster 2_faprotax.csv", "large cluster 3_faprotax.csv",
              "large cluster 4_faprotax.csv", "large cluster 5_faprotax.csv", "large cluster 6_faprotax.csv",
              "large cluster 7_faprotax.csv"]

med_list = ["med cluster 1_faprotax.csv", "med cluster 2_faprotax.csv", "med cluster 3_faprotax.csv",
            "med cluster 4_faprotax.csv", "med cluster 5_faprotax.csv", "med cluster 6_faprotax.csv",
            "med cluster 7_faprotax.csv", "med cluster 8_faprotax.csv"]

final_list = [full_list, cl1_list, cl2_list, cl3_list, large_list, med_list]

for each_list in final_list:
    cluster_function_abundance = {}

    results_file_name = ""
    for file in each_list:

        # prepare the results file name
        results_file_name = file.split("_")[0]
        results_file_name = results_file_name.split(" ")[0]

        with open("input/" + file, "r") as file_1:
            # get the cluster no
            cluster_no = file.split("_")[0]
            cluster_no = cluster_no.split(" ")[2]


            for line in file_1:
                line = line.strip("\n")
                line = line.split(",")
                if line[0] != "group":
                    if line[1] != "0":

                        cluster_function_abundance[cluster_no + line[0]] = [cluster_no,line[0], line[1]]

    # sort the dictionary on alphabetical order
    sorted_c_f_a = dict(sorted(cluster_function_abundance.items()))

    # create files for every network and write the data
    with open("output/" + results_file_name + " network cluster functions and otu abundance.tsv", "w") as file_2:
        file_2.write("Cluster\tFunction\tOTU_abundance\n")
        for key in sorted_c_f_a:

            data = sorted_c_f_a[key]

            file_2.writelines(data[0] + "\t" + data[1] + "\t" + data[2] + "\n")


