#pragma GCC optimize("Ofast")
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <chrono>
#include <algorithm>
#include <functional>
#include <sstream>
using namespace std;
using namespace chrono;
 
// Uncomment the next line to enable debug prints
//#define DEBUG

#include <sys/resource.h> // For increasing stack size on POSIX systems
 
// Global variables for membership testing and clique printing
vector<int> T_, S_;
set<vector<int>> printedCliques;
unordered_map<int,int> cliqueSizeCount;
int maxCliqueSize = 0;
int totalCliques = 0;
vector<bool> inC;  // fast membership indicator for clique C
 
auto start_time = high_resolution_clock::now();
 
// --- INLINE SET OPERATIONS ---
// Using two-pointer techniques (both A and B are assumed sorted)
inline vector<int> setDiff(const vector<int>& A, const vector<int>& B) {
    vector<int> r;
    r.reserve(A.size());
    size_t i = 0, j = 0;
    while(i < A.size() && j < B.size()){
        if(A[i] < B[j]) { r.push_back(A[i]); i++; }
        else if(B[j] < A[i]) { j++; }
        else { i++; j++; }
    }
    while(i < A.size()){
        r.push_back(A[i]);
        i++;
    }
    return r;
}
 
inline vector<int> setIntersect(const vector<int>& A, const vector<int>& B) {
    vector<int> r;
    r.reserve(min(A.size(), B.size()));
    size_t i = 0, j = 0;
    while(i < A.size() && j < B.size()){
        if(A[i] < B[j]) i++;
        else if(B[j] < A[i]) j++;
        else { r.push_back(A[i]); i++; j++; }
    }
    return r;
}
 
// --- Membership update helper ---
// Resets the global membership vector 'inC' to reflect the current clique C.
inline void updateMembership(const vector<int>& C, int n) {
    // For large n, instead of filling the whole vector each time, you might
    // consider using a version that only updates the changed indices.
    fill(inC.begin(), inC.end(), false);
    for (int v : C)
        inC[v] = true;
}
 
// --- Optimized printClique ---
inline void printClique(const vector<int>& C) {
    // We use the same structure; note that using a hash of C could be even faster.
    if(printedCliques.count(C) == 0){
        printedCliques.insert(C);
        int sz = C.size();
        maxCliqueSize = max(maxCliqueSize, sz);
        totalCliques++;
        cliqueSizeCount[sz]++;
    }
}
 
// --- The UPDATE procedure ---
void UPDATE(int i, vector<int>& C, int n, const vector<vector<int>>& adj) {
    if(i == n) {
        printClique(C);
        return;
    }
    
#ifdef DEBUG
    if(i % 1000 == 0) {
        auto elapsed = duration_cast<seconds>(high_resolution_clock::now() - start_time).count();
        cout << "[CP] Entering UPDATE: i=" << i << ", |C|=" << C.size() 
             << ", Time=" << elapsed << " sec\n" << flush;
    }
#endif
 
    vector<int> oldT = T_;
    vector<int> oldS = S_;
 
    vector<int> diff = setDiff(C, adj[i]);
#ifdef DEBUG
    cout << "[CP] i=" << i << " diff.size()=" << diff.size() << "\n" << flush;
#endif
    if(!diff.empty()){
        UPDATE(i+1, C, n, adj);
        T_ = oldT;
        S_ = oldS;
    }
 
    // Instead of using binary_search repeatedly, we use our inC flag.
    for(int y = 0; y < n; y++){
        if(!inC[y]){
            T_[y] = 0;
            S_[y] = 0;
        }
    }
 
    vector<int> cap = setIntersect(C, adj[i]);
#ifdef DEBUG
    cout << "[CP] i=" << i << " cap.size()=" << cap.size() << "\n" << flush;
#endif
 
    for (int x : cap) {
        // Iterate over neighbors of x and update T_
        for (int y : adj[x]) {
            if (y != i && !inC[y]) { // replaced binary_search(C.begin(),C.end(),y)
                T_[y]++;
            }
        }
    }
    for (int x : diff) {
        for (int y : adj[x]) {
            if (!inC[y]) { // replaced binary_search(C.begin(),C.end(),y)
                S_[y]++;
            }
        }
    }
 
    bool FLAG = true;
    int cs = (int)cap.size();
    for (int y : adj[i]) {
        if (!inC[y] && y < i && T_[y] == cs) {
            FLAG = false;
            break;
        }
    }
 
    vector<int> dv = diff;
    sort(dv.begin(), dv.end());
    for (int j : dv) {
        if (!FLAG) break;
        for (int y : adj[j]) {
            if (y < i && !inC[y] && T_[y] == cs) {
                // Here, we still use binary_search on adj[y] because we must
                // check membership in a sorted neighbor list.
                if (binary_search(adj[y].begin(), adj[y].end(), j)) {
                    S_[y]--;
                } else {
                    FLAG = false;
                    break;
                }
            }
        }
    }
 
    if(!cap.empty() && !dv.empty()){
        int jb = dv.back();
        for (int y : cap) {
            if(y < i && (int)adj[i].size() == T_[y] && S_[y] == 0){
                if(jb < y || jb < i - 1){
                    FLAG = false;
                    break;
                }
            }
        }
    }
 
    if(FLAG){
#ifdef DEBUG
        cout << "[CP] i=" << i << " recursing: adding i, current |C|=" << C.size() << "\n" << flush;
#endif
        vector<int> oc = C;
        // Form new clique: new clique = (C ∩ adj[i]) ∪ {i}
        vector<int> nc = cap;
        nc.push_back(i);
        sort(nc.begin(), nc.end());
        C = nc;
        updateMembership(C, n); // update membership flags for new C
        UPDATE(i+1, C, n, adj);
#ifdef DEBUG
        cout << "[CP] i=" << i << " returning from recursion, restoring C\n" << flush;
#endif
        C = oc;
        updateMembership(C, n); // restore membership flags
    }
 
    T_ = oldT;
    S_ = oldS;
 
#ifdef DEBUG
    if(i % 1000 == 0) {
        cout << "[CP] Leaving UPDATE: i=" << i << ", |C|=" << C.size() << "\n" << flush;
    }
#endif
}
 
// --- The main function ---
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    // Increase stack size (for POSIX systems)
    const rlim_t kStackSize = 256 * 1024 * 1024; // 256 MB
    struct rlimit rl;
    if(getrlimit(RLIMIT_STACK, &rl) == 0) {
        if(rl.rlim_cur < kStackSize) {
            rl.rlim_cur = kStackSize;
            if(setrlimit(RLIMIT_STACK, &rl) != 0) {
                cerr << "Error setting stack limit.\n";
            }
        }
    }
    
    ifstream infile("wiki-Vote.txt");
    if(!infile){
        cerr << "Error: Unable to open wiki-Vote.txt\n";
        return 1;
    }
    int n;
    infile >> n;
    vector<vector<int>> adj(n);
    int u, v;
    while(infile >> u >> v){
        if(u >= 0 && u < n && v >= 0 && v < n && u != v){
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    infile.close();
    
    // Sort each vertex's neighbor list (required for binary_search and set operations)
    for(auto &neighbors : adj){
        sort(neighbors.begin(), neighbors.end());
    }
    
    // Initialize global arrays (T_ and S_) and membership vector
    T_.resize(n, 0);
    S_.resize(n, 0);
    inC.assign(n, false);
    
    // Starting clique: using vertex 0 as the initial clique member.
    vector<int> C;
    C.push_back(0);
    updateMembership(C, n);
    
    start_time = high_resolution_clock::now();
    UPDATE(1, C, n, adj);
    auto stop_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop_time - start_time).count();
    
    cout << "\nFinal Results:\n";
    cout << "Largest size of the clique: " << maxCliqueSize << "\n";
    cout << "Total number of maximal cliques: " << totalCliques << "\n";
    cout << "Execution time (ms): " << duration << "\n";
    cout << "Distribution of different size cliques:\n";
    for (int i = 1; i <= maxCliqueSize; i++){
        cout << "Size " << i << ": " << cliqueSizeCount[i] << "\n";
    }
    
    return 0;
}