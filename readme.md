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

## Cite as
@article{TANG2024109102,
title = {Fuzzy rough unlearning model for feature selection},
journal = {International Journal of Approximate Reasoning},
volume = {165},
pages = {109102},
year = {2024},
issn = {0888-613X},
doi = {https://doi.org/10.1016/j.ijar.2023.109102},
url = {https://www.sciencedirect.com/science/article/pii/S0888613X23002335},
author = {Yuxin Tang and Suyun Zhao and Hong Chen and Cuiping Li and Junhai Zhai and Qiangjun Zhou},
keywords = {Fuzzy rough sets, Unlearning, Feature selection, Positive region}
}
