create network function otu_A TSVs.py
# file lists for networks

file_list = ["full net faprotax.csv", "high c faprotax.csv", "medium c faprotax.csv", "low c faprotax.csv",
             "large fs faprotax.csv", "med small fs faprotax.csv"]

for file in file_list:
    cluster_function_abundance = {}

    # prepare the results file name
    results_file_name = file.split(" ")
    results_file_name = results_file_name[0] + " " + results_file_name[1]

    with open("input/" + file, "r") as file_1:

        for line in file_1:
            line = line.strip("\n")
            line = line.split(",")
            if line[0] != "group":
                if line[1] != "0":
                    cluster_function_abundance[line[0]] = [line[0], line[1]]

    # sort the dictionary on alphabetical order
    sorted_c_f_a = dict(sorted(cluster_function_abundance.items()))

    # create files for every network and write the data
    with open("output/" + results_file_name + " network cluster functions and otu abundance.tsv", "w") as file_2:
        file_2.write("Function\tOTU_abundance\n")
        for key in sorted_c_f_a:
            data = sorted_c_f_a[key]

            file_2.writelines(data[0] + "\t" + data[1] + "\n")
