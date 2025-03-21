#include <iostream>
#include <vector>
#include <bitset>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <fstream>
using namespace std;
using namespace chrono;

const int MAX_N = 10000; // Upper bound on number of nodes; adjust if needed

// Global arrays for bookkeeping
vector<int> T_, S_;
set<vector<int>> printedCliques;
unordered_map<int,int> cliqueSizeCount;
int maxCliqueSize = 0;
int totalCliques = 0;

auto start_time = high_resolution_clock::now();

// Utility: convert a bitset representing a clique into a sorted vector of indices.
vector<int> bitsetToVector(const bitset<MAX_N>& B, int n) {
    vector<int> v;
    for (int i = 0; i < n; i++) {
        if (B.test(i)) v.push_back(i);
    }
    return v;
}

// Print clique (stored as bitset) if it hasn't been printed before.
void printClique(const bitset<MAX_N>& C, int n) {
    vector<int> v = bitsetToVector(C, n);
    sort(v.begin(), v.end());
    if (printedCliques.count(v) == 0) {
        printedCliques.insert(v);
        int sz = v.size();
        maxCliqueSize = max(maxCliqueSize, sz);
        totalCliques++;
        cliqueSizeCount[sz]++;
    }
}

// Bitset set difference: returns A \ B
bitset<MAX_N> setDiff(const bitset<MAX_N>& A, const bitset<MAX_N>& B) {
    return A & ~B;
}

// Bitset set intersection: returns A & B
bitset<MAX_N> setIntersect(const bitset<MAX_N>& A, const bitset<MAX_N>& B) {
    return A & B;
}

// UPDATE: recursively enumerate maximal cliques using the Chiba‐style approach.
// C is the current clique (as a bitset); i is the current vertex index we are considering.
// adj is the graph, stored as a vector of bitset.
void UPDATE(int i, bitset<MAX_N>& C, int n, const vector<bitset<MAX_N>>& adj) {
    if (i == n) {
        printClique(C, n);
        return;
    }
    if (i % 1000 == 0) {
        auto elapsed = duration_cast<seconds>(high_resolution_clock::now() - start_time).count();
        cout << "[CP] Entering UPDATE: i=" << i << ", |C|=" << C.count() 
             << ", Time=" << elapsed << " sec" << "\n" << flush;
    }
    
    vector<int> oldT = T_;
    vector<int> oldS = S_;
    
    // Compute diff = (C \ adj[i]): vertices in C that are not adjacent to i.
    bitset<MAX_N> diff = setDiff(C, adj[i]);
    cout << "[CP] i=" << i << " diff.size()=" << diff.count() << "\n" << flush;
    if (diff.any()) {
        UPDATE(i+1, C, n, adj);
        T_ = oldT;
        S_ = oldS;
    }
    
    // For all vertices not in C, reset T_ and S_
    for (int y = 0; y < n; y++) {
        if (!C.test(y)) {
            T_[y] = 0;
            S_[y] = 0;
        }
    }
    
    // Compute cap = C ∩ adj[i]: vertices in C that are adjacent to i.
    bitset<MAX_N> cap = setIntersect(C, adj[i]);
    cout << "[CP] i=" << i << " cap.size()=" << cap.count() << "\n" << flush;
    
    // For each vertex x in cap, update T_ for neighbors not in C.
    for (int x = 0; x < n; x++) {
        if (cap.test(x)) {
            for (int y = 0; y < n; y++) {
                if (adj[x].test(y) && !C.test(y) && y != i) {
                    T_[y]++;
                }
            }
        }
    }
    // For each vertex x in diff, update S_ for neighbors not in C.
    for (int x = 0; x < n; x++) {
        if (diff.test(x)) {
            for (int y = 0; y < n; y++) {
                if (adj[x].test(y) && !C.test(y)) {
                    S_[y]++;
                }
            }
        }
    }
    
    bool FLAG = true;
    int cs = cap.count();
    for (int y = 0; y < n; y++) {
        if (!C.test(y) && y < i && adj[i].test(y) && T_[y] == cs) {
            FLAG = false;
            break;
        }
    }
    
    // Build dv = sorted vector of indices from diff.
    vector<int> dv = bitsetToVector(diff, n);
    sort(dv.begin(), dv.end());
    for (int j : dv) {
        if (!FLAG) break;
        for (int y = 0; y < n; y++) {
            if (y < i && !C.test(y) && T_[y] == cs) {
                if (adj[j].test(y))
                    S_[y]--;
                else {
                    FLAG = false;
                    break;
                }
            }
        }
    }
    
    if (!cap.none() && !dv.empty()) {
        int jb = dv.back();
        for (int y = 0; y < n; y++) {
            if (cap.test(y) && y < i && ((int)adj[i].count() == T_[y]) && S_[y] == 0) {
                if (jb < y || jb < i - 1) {
                    FLAG = false;
                    break;
                }
            }
        }
    }
    
    if (FLAG) {
        auto cp_elapsed = duration_cast<seconds>(high_resolution_clock::now()-start_time).count();
        cout << "[CP] i=" << i << " recursing: adding i, current |C|=" << C.count() 
             << ", Time=" << cp_elapsed << " sec" << "\n" << flush;
        bitset<MAX_N> oc = C;
        // Build new clique = (C ∩ adj[i]) U { i }
        bitset<MAX_N> nc = setIntersect(C, adj[i]);
        nc.set(i, true);
        C = nc;
        UPDATE(i+1, C, n, adj);
        cout << "[CP] i=" << i << " returning from recursion, restoring C" << "\n" << flush;
        C = oc;
    }
    
    T_ = oldT;
    S_ = oldS;
    
    if (i % 1000 == 0) {
        cout << "[CP] Leaving UPDATE: i=" << i << ", |C|=" << C.count() << "\n" << flush;
    }
}
 
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
 
    ifstream infile("wiki-Vote.txt");
    if (!infile) {
        cerr << "Error: Unable to open wiki-Vote.txt\n";
        return 1;
    }
    int n;
    infile >> n;
    vector<bitset<MAX_N>> adj(n);
    int u, v;
    while(infile >> u >> v) {
        if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
            adj[u].set(v, true);
            adj[v].set(u, true);
        }
    }
    infile.close();
 
    T_.resize(n, 0);
    S_.resize(n, 0);
 
    bitset<MAX_N> C;
    C.reset();
    C.set(0, true);
 
    start_time = high_resolution_clock::now();
    UPDATE(1, C, n, adj);
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