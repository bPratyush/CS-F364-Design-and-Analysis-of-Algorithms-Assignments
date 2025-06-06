
# CS F364 Design and Analysis of Algorithms Assignment II

The repository contains source code for the two algorithms implemented mentioned in the paper:
- Algorithm 1: Exact Algorithm
- Algorithm 4: CoreExact Algorithm

## Webpage Link

[DAA Assignment II Website](https://glittery-croissant-0480a0.netlify.app/)


## Dataset Preparation

All datasets are provided in ```.txt``` format. No extensive preprocessing has been performed on the given datasets. Instead, the source code dynamically processes the data as required.  

[Dataset Folder](https://drive.google.com/drive/folders/1PMXBEOQy44Y198SDR4oAbJyCjG57gqef?usp=sharing)

## Execution Instruction
- Ensure the source code and datasets are in the same directory.
- We compiled all source codes on ```clang++``` or ```g++``` compiler and using C++ version 17 or C++ version 14.
- Kindly clear the input buffer and allocated recursion stack space before proceeding with the next input or source code to ensure that a segmentation fault does not occur
- All executions were performed on a 2020 MacBook Pro running MacOS Sequoia 15.3.2 with 8GB of RAM and 8-core Apple ARM M2 CPU or with a MacBook Pro having a 10-core Apple M1 Pro CPU and 16 GB of RAM
- The codes written have been supported and verified for Linux Ubuntu distribution and MacOS Operating Systems

### Compilation
For compiling ```src.cpp``` stored in ```dir``` directory with C++ version 17 and using clang++ compiler, run the following command
 ```bash
clang++ -std=c++17 dir/src.cpp -o src
```
- For using ```g++``` compiler, replace ```clang++``` with ```g++```
- For changing C++ version, either change the value in ```-std=c++17``` or remove for default version execution present in the system

### Execution
For execution, run the following command:
```bash
./src t1.txt h
```
```t1.txt``` represents the dataset file in ```.txt``` format and ```h``` represents the clique size that the algorithm uses to find the densest subgraph

### Output File

The Exact Algorithm stores results in ```algo_1_result.txt``` file
The CoreExact Algorithm results are stored in ```algo-4-results.txt``` file


## Contributors

- ```Pratyush Bindal (2022A7PS0119H)```: Pratyush did a comprehensive walkthrough of the Exact Algorithm mentioned in the paper. He worked upon preparing the Netscience dataset and the README file. He also contributed in the final report and ensured coordination among the group members.

- ```Kalash Bhattad (2022A7PS0065H)```: Kalash was deeply involved in the extensive implementation of the Exact algorithm, meticulously resolving bugs and ensuring its correctness. He also carefully prepared the AS-733 and Yeast dataset on Exact Algorithm and contributed meaningfully to the preparation of the final report.

- ```RVS Aashrey Kumar (2022A7PS0160H)```: Aashrey has conducted comprehensive implementation of the CoreExact Algorithm as mentioned in the paper. He has also worked on curating relevant datasets (AS-733, NetScience and CA-HepTh), collecting results and benchmarking with the paper's results.

- ```Venkata Saketh Dakuri (2022A7PS0056H)```: Saketh played a key role in designing and preparing the website for presenting the source code and results. He also significantly contributed to refining and optimizing both algorithms, ensuring improved efficiency and clarity.

- ```Kavya Ganatra (2022A7PS0057H)```: Kavya made significant contributions to the project report and verified the code and results to ensure accuracy. He also conducted a theoretical review of the assigned research paper. His efforts strengthened the overall depth and quality of the work.
