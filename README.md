
# CS F364 Design and Analysis of Algorithms Assignment I

The repository contains source code for the three algorithms implemented for Maximal Clique Enumeration (MCE):
- Eppstein, Löffler & Strash (2010) algorithm
- Tomita, Tanata & Takahashi (2006) algorithm
- Chiba & Nishizeki (1985) algorithm

## Webpage Link

[Maximal Cliques Algorithms](https://serene-gnome-fc5c3f.netlify.app/)


## Dataset Preparation

All datasets are provided in ```.txt``` format. No extensive preprocessing has been performed on the given datasets. Instead, the source code dynamically processes the data as required.  

[Dataset Folder](https://drive.google.com/drive/folders/1PMXBEOQy44Y198SDR4oAbJyCjG57gqef?usp=sharing)

## Execution Instruction
- Ensure the source code and datasets are in the same directory.
- We compiled all source codes on ```clang++``` compiler using optimisation level 3 (```O3```) and C++ version 17 or C++ version 14.
- Kindly clear the input buffer and allocated recursion stack space before proceeding with the next input or source code to ensure that a segmentation fault does not occur
- All executions were performed on a 2020 MacBook Pro running MacOS Sequoia 15.3.2 with 8GB of RAM and 8-core Apple ARM M2 CPU or with a MacBook Pro having a 10-core Apple M1 Pro CPU and 16 GB of RAM
- The codes written have been supported and verified for Linux Ubuntu distribution and MacOS Operating Systems

### Compilation
For compiling ```src.cpp``` stored in ```dir``` directory, run the following command
 ```bash
clang++ -O3 -std=c++17 dir/src.cpp -o src
```

### Execution with Increased Stack Space
For source codes requiring large recursion stack space such as Chiba & Nishizeki (1985) algorithm and optimised Tomita, Tanata & Takahashi (2006) algorithm, execute the code with root privileges for increasing the default stack space of the machine by running the following command
```bash
sudo ./src t1.txt
```
After running the command, enter your login password when prompted.  

```t1.txt``` represents the dataset file in ```.txt``` format

### Execution
For all other algorithms execution, either continue execution with increased stack space or run the following command
```bash
./src t1.txt
```
```t1.txt``` represents the dataset file in ```.txt``` format

### Output File

The source codes for all algorithms store results in ```output.txt``` file


## Contributors

- ```Pratyush Bindal (2022A7PS0119H)```: Pratyush did a comprehensive walkthrough and worked on the source code implementation as well as optimisation of Eppstein, Löffler & Strash (2010) and Tomita, Tanata & Takahashi (2006) algorithms. He also contributed in optimising Chiba and Nishizeki (1985) algorithm and worked upon preparing the Readme file and ensured coordination among the group members.

- ```Kalash Bhattad (2022A7PS0065H)```: Kalash worked on the implementation of Chiba & Nishizeki (1985) algorithm, extensively contributed to the final project report, and wrote the introduction along with detailed explanations of two algorithms including the english description, source code implementation walkthrough, relevant results and in-depth analysis.

- ```RVS Aashrey Kumar (2022A7PS0160H)```: Aashrey did a comprehensive walkthrough of Tomita, Tanata & Takahashi (2006) algorithm and did a detailed comparative analysis of all the three algorithms. He contributed in implementation of Chiba and Nishizeki (1985) algorithm. He also listed key findings from maximal cliques in graph theory and applications of maximal cliques in real life scenarios.

- ```Venkata Saketh Dakuri (2022A7PS0056H)```: Saketh extensively worked on the website ensuring it to be fully functional and updated it with the source codes, visualisation and presentation of obtained results, background and theoretical ideas about the algorithms from the papers. He also did a comprehensive walkthrough of Eppstein, Löffler & Strash (2010) algorithm by preparing a report.

- ```Kavya Ganatra (2022A7PS0057H)```: Kavya worked extensively on the Chiba & Nishizeki (1985) algorithm, exploring and experimenting with various optimization techniques to improve its efficiency without altering its core execution logic. He prepared an in-depth explanation of the Chiba & Nishizeki (1985) algorithm, detailing its core concepts, code implementation, and the optimizations explored to enhance its efficiency and additionally played a key role in drafting and refining the final project report.
