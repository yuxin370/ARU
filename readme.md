## Background
In the big data era, some data, becoming meaningless or illegal over time and space, need to be deleted from historical knowledge. It is a challenging problem, called machine unlearning, to efficiently forget the information of those outdated data from historical models. Some unlearning techniques have been proposed in loss-well-defined classification models, such as SVM, Random Forest, and Federated learning model. Yet, it is under study to remove outdated data from learned feature selection in fuzzy rough philosophy. To narrow this gap, we propose a fuzzy rough unlearning model for feature selection. 

## Platform：
Windows 10 professional

Intel® Xeon® W-2145 CPU@ 3.70GHz

128 GB memory 

## Usage
- compile: make
- clean: make clean
- usage: test.exe start_of_conditionAttribute end_of_conditionAttribute decisionAttributePosition datasets initsize alpha  foldsize >> result_path

  eg: test.exe 1 19 20 demo 2310 0.3 5 >> result\demo\output.txt

## Data split
we split the original datasets into six subsets and in each loop one subsets was removed. 
