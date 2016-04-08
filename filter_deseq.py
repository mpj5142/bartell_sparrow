count_data=open("counts_all_norm_deseq2.txt", "rU")
filter_name=open("names_all_100.txt", "rU")
out_file=open("deseq2_all_100.tdms", "w")

count_lines=count_data.readlines()
filter_lines=filter_name.readlines()

for line in count_lines:
	line_split=line.split("\t")
	count_sample=line_split[0]
	for filter_sample in filter_lines:
		if count_sample+"\n"==filter_sample:
			out_file.write(line)
			break

