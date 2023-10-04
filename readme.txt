###
platform：
Windows 10 professional
Intel® Xeon® W-2145 CPU@ 3.70GHz
128 GB memory 

###
compile:
make
clean:
make clean
###

usage:
test.exe start_of_conditionAttribute end_of_conditionAttribute decisionAttributePosition datasets initsize alpha  foldsize >> result_path

eg:
test.exe 1 19 20 demo 2310 0.3 5 >> result\demo\output.txt

# we split the original datasets into six subsets and in each loop one subsets was removed. 
