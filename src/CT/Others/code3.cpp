#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace chrono;

vector<int> T_, S_;
set<vector<int>> printedCliques;
unordered_map<int, int> cliqueSizeCount;
int maxCliqueSize = 0;
int totalCliques = 0;

auto start_time = high_resolution_clock::now();

void printClique(const unordered_set<int> &C) {
    vector<int> v(C.begin(), C.end());
    sort(v.begin(), v.end());
    if (printedCliques.count(v) == 0) {
        printedCliques.insert(v);
        int size = v.size();
        maxCliqueSize = max(maxCliqueSize, size);
        totalCliques++;
        cliqueSizeCount[size]++;
        if (totalCliques % 50000 == 0) {
            auto elapsed = duration_cast<seconds>(high_resolution_clock::now() - start_time).count();
            cout << "[Checkpoint] Found " << totalCliques 
                 << " unique cliques so far. Time elapsed: " 
                 << elapsed << " sec\n" << flush;
        }
    }
}

unordered_set<int> setdiff(const unordered_set<int> &A, const unordered_set<int> &B) {
    unordered_set<int> r;
    for (int x : A) {
        if (!B.count(x)) {
            r.insert(x);
        }
    }
    return r;
}

unordered_set<int> setintersect(const unordered_set<int> &A, const unordered_set<int> &B) {
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

void UPDATE(int i, unordered_set<int> &C, int n, const vector<unordered_set<int>> &adj) {
    if (i == n) {
        printClique(C);
        return;
    }
    
    if (i % 500 == 0) {
        auto elapsed = duration_cast<seconds>(high_resolution_clock::now() - start_time).count();
        cout << "[Checkpoint] Recursion i=" << i << ", |C|=" << C.size() 
             << ", Time elapsed: " << elapsed << " sec\n" << flush;
    }
    
    vector<int> oldT = T_;
    vector<int> oldS = S_;
    
    unordered_set<int> diff = setdiff(C, adj[i]);
    if (!diff.empty()) {
        // Continue without adding vertex i.
        UPDATE(i + 1, C, n, adj);
        T_ = oldT;
        S_ = oldS;
    }
    
    // Reset T_ and S_ for non-members of C.
    for (int y = 0; y < n; y++) {
        if (!C.count(y)) {
            T_[y] = 0;
            S_[y] = 0;
        }
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
    
    bool FLAG = true;
    int cs = (int)cap.size();
    for (int y : adj[i]) {
        if (!C.count(y) && y < i && T_[y] == cs) {
            FLAG = false;
            break;
        }
    }
    
    vector<int> dv(diff.begin(), diff.end());
    sort(dv.begin(), dv.end());
    
    for (int j : dv) {
        if (!FLAG) break;
        for (int y : adj[j]) {
            if (y < i && !C.count(y) && T_[y] == cs) {
                if (adj[y].count(j)) {
                    S_[y]--;
                } else {
                    FLAG = false;
                    break;
                }
            }
        }
    }
    
    if (!cap.empty() && !dv.empty()) {
        int jb = dv.back();
        for (int y : cap) {
            if (y < i && (int)adj[i].size() == T_[y] && S_[y] == 0) {
                if (jb < y || jb < i - 1) {
                    FLAG = false;
                    break;
                }
            }
        }
    }
    
    if (FLAG) {
        // Add vertex i to the clique and recurse.
        unordered_set<int> oc = C;
        unordered_set<int> nc = cap;
        nc.insert(i);
        C = nc;
        UPDATE(i + 1, C, n, adj);
        C = oc;
    }
    
    T_ = oldT;
    S_ = oldS;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    ifstream infile("wiki-Vote.txt");
    if (!infile) {
        cerr << "Error: Unable to open wiki-Vote.txt\n";
        return 1;
    }
    cout << "[Checkpoint] Successfully opened wiki-Vote.txt\n" << flush;
    
    // Read edges without requiring a node count.
    vector<pair<int, int>> edges;
    int u, v;
    int maxNode = -1;
    int edgeCount = 0;
    while (infile >> u >> v) {
        edges.push_back({u, v});
        maxNode = max(maxNode, max(u, v));
        edgeCount++;
        if (edgeCount % 100000 == 0) {
            cout << "[Checkpoint] Read " << edgeCount << " edges so far...\n" << flush;
        }
    }
    infile.close();
    cout << "[Checkpoint] Finished reading " << edgeCount << " edges.\n" << flush;
    
    int n = maxNode + 1;
    cout << "[Checkpoint] Computed number of nodes: " << n << "\n" << flush;
    
    // Build graph as a vector of unordered_set for each node.
    vector<unordered_set<int>> adj(n);
    for (auto &e : edges) {
        int a = e.first, b = e.second;
        if (a != b && a >= 0 && a < n && b >= 0 && b < n) {
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }
    
    T_.resize(n, 0);
    S_.resize(n, 0);
    
    start_time = high_resolution_clock::now();
    // Iterate over all vertices as starting points.
    for (int v = 0; v < n; v++) {
        unordered_set<int> C;
        C.insert(v);
        UPDATE(v + 1, C, n, adj);
    }
    auto stop_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop_time - start_time).count();
    
    cout << "\nFinal Results:\n";
    cout << "Largest size of the clique: " << maxCliqueSize << "\n";
    cout << "Total number of maximal cliques: " << totalCliques << "\n";
    cout << "Execution time (ms): " << duration << "\n";
    cout << "Distribution of different size cliques:\n";
    for (int i = 1; i <= maxCliqueSize; i++) {
        cout << "Size " << i << ": " << cliqueSizeCount[i] << "\n";
    }
    
    return 0;
}

