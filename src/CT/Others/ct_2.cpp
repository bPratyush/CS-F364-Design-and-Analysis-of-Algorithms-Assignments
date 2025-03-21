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

// Global variables (removed printedCliques)
vector<int> T_, S_;
unordered_map<int, int> mp; 
int numCliques = 0;
int maxCliqueSize = 0;
vector<unordered_set<int>> adj;
int n;
bool FLAG;

int myMax(int a, int b) {
    return a > b ? a : b;
}

// Record clique C by counting it; does not store the clique itself.
void printClique(const unordered_set<int>& C) {
    vector<int> v(C.begin(), C.end());
    maxCliqueSize = myMax(maxCliqueSize, (int)v.size());
    sort(v.begin(), v.end());
    mp[v.size()]++;  
    numCliques++;
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
    // Enforce canonical ordering:
    if (!C.empty()) {
        int m = *max_element(C.begin(), C.end());
        if (i <= m) {
            UPDATE(m + 1, C);
            return;
        }
    }
    
    if (i == n) {
        printClique(C);
        return;
    }
    
    // diff = elements in C that are NOT neighbors of i.
    unordered_set<int> diff;
    for (int x : C) {
        if (!adj[i].count(x))
            diff.insert(x);
    }
    
    // Recurse if some vertices of C are not adjacent to i.
    if (!diff.empty())
        UPDATE(i + 1, C);
    
    // cap = elements in C that ARE neighbors of i.
    unordered_set<int> cap;
    for (int x : C) {
        if (adj[i].count(x))
            cap.insert(x);
    }
    
    // Update T_ and S_.
    for (int x : cap)
        for (int y : adj[x])
            if (y != i && !C.count(y))
                T_[y]++;
    for (int x : diff)
        for (int y : adj[x])
            if (!C.count(y))
                S_[y]++;
    
    FLAG = true;
    int cs = cap.size();
    for (int y : adj[i]) {
        if (!C.count(y) && y < i && T_[y] == cs) {
            FLAG = false;
            break;
        }
    }
    
        // Prepare sorted diff (dv) for lexicographic tests.
    vector<int> dv(diff.begin(), diff.end());
    sort(dv.begin(), dv.end());
    const int p_int = static_cast<int>(dv.size());
    for (int k = 0; k < p_int && FLAG; k++) {
        int j_k = dv[k];  // current vertex from sorted diff (C - N(i))
        for (int y : adj[j_k]) {
            if ((y < i) && !C.count(y) && (T_[y] == cs)) {
                if (y == j_k) {
                    S_[y]--;
                } else if (y < j_k) {
                    S_[y] = S_[j_k];
                    if ((S_[y] + k) == p_int && (y > (j_k - 1))) { 
                        FLAG = false;
                        break; // exit inner loop if duplicate extension detected.
                    }
                }
            }
        }
    }
    
    // Lexicographic test – Step 7: further remove non-maximal cliques.
    int jp = dv.empty() ? 0 : dv.back();
    if (!cap.empty()) {
        for (int y = 0; y < i && FLAG; y++) {
            if (!C.count(y) && (T_[y] == cs) && (S_[y] == 0)) {
                if (jp < y) {
                    FLAG = false;
                    break;
                }
            }
        }
    } else {
        if (jp < (i - 1))
            FLAG = false;
    }
    
    // Reset T_ and S_ for nodes in adjacency of cap.
    for (int x : cap)
        for (int y : adj[x])
            if ((y != i) && !C.count(y))
                T_[y] = 0;
    for (int x : cap)
        for (int y : adj[x])
            if (!C.count(y))
                S_[y] = 0;
    
    if (FLAG) {
        unordered_set<int> save = diff;
        cap.insert(i);
        C = cap;
        UPDATE(i + 1, C);
        C.erase(i);
        for (int val : save)
            C.insert(val);
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
        if (line.empty() || line[0] == '#')
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