#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <set>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <sys/resource.h>
using namespace std;
using namespace chrono;

vector<int> T_, S_;
unordered_map<int, int> mp; 
int numCliques = 0;
int maxCliqueSize = 0;
vector<unordered_set<int>> adj;
int n;
bool FLAG;
set<vector<int>> printedCliques;

int myMax(int a, int b) {
    return a > b ? a : b;
}

void printClique(const unordered_set<int>& C) {
    vector<int> v(C.begin(), C.end());
    maxCliqueSize = myMax(maxCliqueSize, (int)v.size());
    sort(v.begin(), v.end());
    if (printedCliques.count(v) == 0) {
        mp[v.size()]++;  
        printedCliques.insert(v);
        numCliques++;
    }
}

unordered_set<int> setdiff(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> r;
    for (int x : A) {
        if (!B.count(x)) {
            r.insert(x);
        }
    }
    return r;
}

unordered_set<int> setintersect(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> r;
    if (A.size() < B.size()) {
        for (int x : A) {
            if (B.count(x)) {
                r.insert(x);
            }
        }
    } else {
        for (int x : B) {
            if (A.count(x)) {
                r.insert(x);
            }
        }
    }
    return r;
}

void UPDATE(int i, unordered_set<int>& C) {
    if (i % 1000 == 0) { 
        cout << "[DEBUG] Processing node " << i << endl;
    }
    
    if (i == n) { 
        printClique(C);
    } else {
        unordered_set<int> diff = setdiff(C, adj[i]);
        if (!diff.empty()) { 
            UPDATE(i + 1, C);
        }
        unordered_set<int> cap = setintersect(C, adj[i]);
  
        for (int x : cap) { 
            for (int y : adj[x]) {
                if (y != i && !C.count(y)) {
                    T_[y]++;
                }
            }
        }
        for (int x : diff) {
            for (int y : adj[x]) {
                if (!C.count(y)) { 
                    S_[y]++; 
                }
            }
        }
     
        FLAG = true; 
        int cs = (int)cap.size();
        for (int y : adj[i]) {
            if (!C.count(y) && y < i && T_[y] == cs) { 
                FLAG = false; 
                break;
            }
        }
    
        vector<int> dv(diff.begin(), diff.end());
        sort(dv.begin(), dv.end());
    
        for (int k = 0; k < (int)dv.size(); k++) {
            int j_k = dv[k];
            bool dummyflag = false;
            for (int y : adj[j_k]) {
                if (y < i && !C.count(y) && T_[y] == cs) { 
                    if (y >= j_k) { 
                        S_[y]--;
                    } else {            
                       if (!dummyflag) {
                           dummyflag = true;
                           if (k > 0 && S_[y] + k == (int)dv.size() && y >= dv[k-1]) {   
                                FLAG = false;
                                break;
                           }
                       }       
                    }
                }
            }
        }
     
        int jp = 0;
        if (!dv.empty())
            jp = dv[dv.size() - 1];
        if (!cap.empty()) { 
            for (int y : dv) {
                if (y != i && !C.count(y) && y < i && (int)cap.size() == T_[y] && S_[y] == 0) {
                    if (jp < y) {
                        FLAG = false;
                        break;
                    }
                }
            }
        } else if (jp < i - 1) {
            FLAG = false;
        }
    
        for (int x : cap) {
            for (int y : adj[x]) {
                if (y != i && C.count(y) == 0) {
                    T_[y] = 0; 
                }
            }
        }
    
        for (int x : cap) {
            for (int y : adj[x]) {
                if (C.count(y) == 0) {
                    S_[y] = 0; 
                }
            }
        }
    
        if (FLAG) {
            unordered_set<int> save = diff;
            cap.insert(i);      
            C = cap;
            UPDATE(i + 1, C);
            C.erase(i);
            C.insert(save.begin(), save.end());
        }
    }
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [input file]" << endl;
        return 1;
    }
    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error: Unable to open input file " << argv[1] << endl;
        return 1;
    }
    
    int u, v;
    int maxVertex = -1;
    vector<pair<int, int>> edges;
    string line;
    while(getline(infile, line)) {
        if(line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        if(iss >> u >> v) {
            maxVertex = max(maxVertex, max(u, v));
            edges.push_back({u, v});
        }
    }
    infile.close();
    n = maxVertex + 1;
    adj.resize(n);
    for (auto &edge : edges) {
        u = edge.first;
        v = edge.second;
        if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
            adj[u].insert(v);
            adj[v].insert(u);
        }
    }
    
    T_.resize(n, 0);
    S_.resize(n, 0);
    unordered_set<int> C;
    C.insert(0);
    auto start_time = high_resolution_clock::now();
    UPDATE(1, C);
    auto stop_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop_time - start_time).count();
    
    ofstream outfile("output.txt");
    if (!outfile) {
        cerr << "Error: Unable to open output.txt for writing" << endl;
        return 1;
    }
    
    outfile << "Largest size of a clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << numCliques << "\n";
    outfile << "Execution time (ms): " << duration << "\n";
    outfile << "Distribution of clique sizes:\n";
    for (int i = 1; i <= maxCliqueSize; i++) {
        outfile << "Size " << i << ": " << mp[i] << "\n";
    }
    outfile.close();
    return 0;
}