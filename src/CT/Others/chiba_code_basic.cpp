#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <map>
#include <sys/resource.h>
#include <cstring>
#include <cstdlib>
using namespace std;
using namespace chrono;

int maxCliqueSize = 0;
int totalMaximalCliques = 0;
int maxVertex = 0;
int counter = 0; 
vector<vector<int> > adjacencyList; 
map<int, int> cliqueSizeDistribution;
vector<int> T, S;

void addEdge(int u, int v, vector<vector<int> > & adj) {
    if(u >= (int)adj.size() || v >= (int)adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    
    adj[u].push_back(v);
    adj[v].push_back(u);
}


void printClique(vector<int> &clique){
    int size = (int)clique.size();
    
    cliqueSizeDistribution[size]++;
    
    if (size > maxCliqueSize) {
        maxCliqueSize = size;
    }
    
    cout << counter++ << endl;
}


vector<int> setDiff(vector<int> &A, vector<int> &B) {
    vector<int> r;
    for (int x : A) {
        if (find(B.begin(), B.end(), x) == B.end()) {
            r.push_back(x);
        }
    }
    return r;
}


vector<int> setIntersection(vector<int> &A, vector<int> &B) {
    vector<int> r;
    if (A.size() < B.size()) {
        for (int x : A) {
            if (find(B.begin(), B.end(), x) != B.end()) {
                r.push_back(x);
            }
        }
    } else {
        for (int x : B) {
            if (find(A.begin(), A.end(), x) != A.end()) {
                r.push_back(x);
            }
        }
    }
    return r;
}

void UPDATE(int i, vector<int> & C) {
    
    if (!C.empty()) {
        int m = *max_element(C.begin(), C.end());
        if (i <= m) {
            UPDATE(m + 1, C);
            return;
        }
    }
    

    if(i == maxVertex + 1){
        
        printClique(C);
    }
    else {
        
        vector<int> diff = setDiff(C, adjacencyList[i]);
        
        vector<int> intersection = setIntersection(C, adjacencyList[i]);

        
        if(!diff.empty()){
            UPDATE(i + 1, C);
        }

        
        for(int x : intersection){
            vector<int> diff1 = setDiff(adjacencyList[x], C);
            for(int y : diff1){
                if(y != i){
                    T[y]++;
                }
            }
        }
        
        for(int x : diff){
            vector<int> diff1 = setDiff(adjacencyList[x], C);
            for(int y : diff1){
                S[y]++;
            }
        }

        
        bool flag = true;
        
        vector<int> diff1 = setDiff(adjacencyList[i], C);
        for(int y : diff1){
            if(y < i && T[y] == (int)intersection.size()){
                flag = false;
                break;
            }
        }

        
        vector<int> v(diff.begin(), diff.end());
        sort(v.begin(), v.end());
        int p = (int)v.size();

        
        for(int k = 0; k < p && flag; k++){
            int x = v[k];
            vector<int> diff2 = setDiff(adjacencyList[x], C);
            
            for(int y : diff2) {
                
                if(y >= i) continue;
                
                if(T[y] == (int)intersection.size()){
                    if(y >= x){
                        
                        S[y]--;
                    }
                    else {
                        
                        if((k == 0 || y >= v[k - 1]) && (S[y] + k == p)) {
                            flag = false;
                            break;
                        }
                    }
                }
            }
        }

       
        if(!intersection.empty()){
            
            for(int x : intersection){
                vector<int> diffX = setDiff(adjacencyList[x], C);
                for(int y : diffX){
                    if(y < i && y != i && T[y] == (int)intersection.size() && S[y] == 0){
                        
                        if(!v.empty() && v[p-1] < y){
                            flag = false;
                        }
                    }
                }
            }
        }
        else {
            
            if(!v.empty() && v[p-1] < i - 1){
                flag = false;
            }
        }

        
        for(int x : intersection){
            vector<int> diffX = setDiff(adjacencyList[x], C);
            for(int y : diffX){
                if(y != i){
                    T[y] = 0;
                }
            }
        }
        
        for(int x : diff){
            vector<int> diffX = setDiff(adjacencyList[x], C);
            for(int y : diffX){
                S[y] = 0;
            }
        }

        
        if(flag){
            
            vector<int> save = diff;
            
            vector<int> result = intersection;
            result.push_back(i);

            
            C = result;
            UPDATE(i + 1, C);

            
            C.erase(find(C.begin(), C.end(), i));
            C.insert(C.end(), save.begin(), save.end());
        }
    }
}

int main(int argc, char* argv[]){
    
    struct rlimit rl;
    rlim_t stack_size = 512 * 1024 * 1024;
    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result != 0) {
        cerr << "Error getting stack limit: " << strerror(errno) << endl;
        return 1;
    }
    if (rl.rlim_cur < stack_size) {
        rl.rlim_cur = stack_size;
        if (rl.rlim_max < rl.rlim_cur)
            rl.rlim_max = rl.rlim_cur;
        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0)
            cerr << "Error setting stack limit: " << strerror(errno) << endl;
        else
            cerr << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << endl;
    }
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [input file]" << "\n";
        return 1;
    }
    
    
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    if (!infile) {
        cerr << "Error: Unable to open input file " << argv[1] << "\n";
        return 1;
    }
    int u, v;
    string line;
    vector<pair<int,int> > edgeList;
    while(getline(infile, line)){
        if(line.empty() || line[0]=='#'){
            continue;
        }
        istringstream iss(line);
        if(iss >> u >> v){
            edgeList.push_back(std::make_pair(u,v));
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    
    
    adjacencyList.resize(maxVertex+1);
    for(auto &edge : edgeList){
        addEdge(edge.first, edge.second, adjacencyList);
    }
    
    
    int n = adjacencyList.size();
    vector<pair<int, int> > degreeIndex;
    degreeIndex.reserve(n);
    for (int i = 0; i < n; i++) {
        degreeIndex.push_back(std::make_pair((int)adjacencyList[i].size(), i));
    }
    
    sort(degreeIndex.begin(), degreeIndex.end());
    
    
    vector<int> oldToNew(n), newToOld(n);
    for (int newIndex = 0; newIndex < n; newIndex++) {
        int oldIndex = degreeIndex[newIndex].second;
        oldToNew[oldIndex] = newIndex;
        newToOld[newIndex] = oldIndex;
    }
    
   
    vector<vector<int> > newAdjacencyList(n);
    for (int oldU = 0; oldU < n; oldU++) {
        for (int oldV : adjacencyList[oldU]) {
            int newU = oldToNew[oldU];
            int newV = oldToNew[oldV];
            newAdjacencyList[newU].push_back(newV);
        }
    }
    
    for (int i = 0; i < n; i++) {
        sort(newAdjacencyList[i].begin(), newAdjacencyList[i].end());
    }
    adjacencyList = std::move(newAdjacencyList);
    maxVertex = n - 1;
   
    T.resize(maxVertex+1, 0);
    S.resize(maxVertex+1, 0);
    
    
    vector<int> C;
    C.push_back(0);
    
    auto startTime = high_resolution_clock::now();
    UPDATE(1, C);
    auto stopTime = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stopTime - startTime).count();
    
    
    outfile << "Largest size of a clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << counter << "\n";
    outfile << "Execution time (ms): " << duration << "\n";
    outfile << "Distribution of clique sizes:\n";
    for(auto &kv : cliqueSizeDistribution){
        outfile << "Size " << kv.first << ": " << kv.second << "\n";
    }
    outfile.close();
    
    cout << "\nFinal Results:\n";
    cout << "Largest size of a clique: " << maxCliqueSize << "\n";
    cout << "Total number of maximal cliques: " << counter << "\n";
    cout << "Execution time (ms): " << duration << "\n";
    cout << "Distribution of clique sizes:\n";
    for(auto &kv : cliqueSizeDistribution){
        cout << "Size " << kv.first << ": " << kv.second << "\n";
    }
    cout << endl;

    return 0;
}