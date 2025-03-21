#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <set>
#include <unordered_map>
using namespace std;
using namespace chrono;

#include <sys/resource.h> // For increasing stack size on POSIX systems

vector<int> T_, S_;
set<vector<int>> printedCliques;
unordered_map<int,int> cliqueSizeCount;
int maxCliqueSize = 0;
int totalCliques = 0;
 
auto start_time = high_resolution_clock::now();
 
vector<int> setDiff(const vector<int>& A, const vector<int>& B) {
    vector<int> r;
    for (int x : A) {
        if (!binary_search(B.begin(), B.end(), x))
            r.push_back(x);
    }
    return r;
}
 
vector<int> setIntersect(const vector<int>& A, const vector<int>& B) {
    vector<int> r;
    int i = 0, j = 0;
    while(i < A.size() && j < B.size()){
        if(A[i] < B[j]) i++;
        else if(B[j] < A[i]) j++;
        else { r.push_back(A[i]); i++; j++; }
    }
    return r;
}
 
void printClique(const vector<int>& C) {
    vector<int> v = C;
    if(printedCliques.count(v) == 0){
        printedCliques.insert(v);
        int sz = v.size();
        maxCliqueSize = max(maxCliqueSize, sz);
        totalCliques++;
        cliqueSizeCount[sz]++;
    }
}
 
void UPDATE(int i, vector<int>& C, int n, const vector<vector<int>>& adj) {
    if(i == n) {
        printClique(C);
        return;
    }
    if(i % 1000 == 0) {
        auto elapsed = duration_cast<seconds>(high_resolution_clock::now()-start_time).count();
        cout << "[CP] Entering UPDATE: i=" << i << ", |C|=" << C.size() 
             << ", Time=" << elapsed << " sec\n" << flush;
    }
 
    vector<int> oldT = T_;
    vector<int> oldS = S_;
 
    vector<int> diff = setDiff(C, adj[i]);
    cout << "[CP] i=" << i << " diff.size()=" << diff.size() << "\n" << flush;
    if(!diff.empty()){
        UPDATE(i+1, C, n, adj);
        T_ = oldT;
        S_ = oldS;
    }
 
    for(int y = 0; y < n; y++){
        if(!binary_search(C.begin(), C.end(), y)){
            T_[y] = 0;
            S_[y] = 0;
        }
    }
 
    vector<int> cap = setIntersect(C, adj[i]);
    cout << "[CP] i=" << i << " cap.size()=" << cap.size() << "\n" << flush;
 
    for(int x : cap){
        for(int y : adj[x]){
            if(y != i && !binary_search(C.begin(), C.end(), y)){
                T_[y]++;
            }
        }
    }
    for(int x : diff){
        for(int y : adj[x]){
            if(!binary_search(C.begin(), C.end(), y)){
                S_[y]++;
            }
        }
    }
 
    bool FLAG = true;
    int cs = (int)cap.size();
    for(int y : adj[i]){
        if(!binary_search(C.begin(), C.end(), y) && y < i && T_[y] == cs){
            FLAG = false;
            break;
        }
    }
 
    vector<int> dv = diff;
    sort(dv.begin(), dv.end());
    for (int j : dv){
        if(!FLAG) break;
        for (int y : adj[j]){
            if(y < i && !binary_search(C.begin(), C.end(), y) && T_[y] == cs){
                if(binary_search(adj[y].begin(), adj[y].end(), j)){
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
        for (int y : cap){
            if(y < i && (int)adj[i].size() == T_[y] && S_[y] == 0){
                if(jb < y || jb < i - 1){
                    FLAG = false;
                    break;
                }
            }
        }
    }
 
    if(FLAG){
        cout << "[CP] i=" << i << " recursing: adding i, current |C|=" << C.size() << "\n" << flush;
        vector<int> oc = C;
        vector<int> nc = cap;
        nc.push_back(i);
        sort(nc.begin(), nc.end());
        C = nc;
        UPDATE(i+1, C, n, adj);
        cout << "[CP] i=" << i << " returning from recursion, restoring C\n" << flush;
        C = oc;
    }
 
    T_ = oldT;
    S_ = oldS;
 
    if(i % 1000 == 0) {
        cout << "[CP] Leaving UPDATE: i=" << i << ", |C|=" << C.size() << "\n" << flush;
    }
}
 
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
 
    for(auto &neighbors : adj){
        sort(neighbors.begin(), neighbors.end());
    }
 
    T_.resize(n, 0);
    S_.resize(n, 0);
 
    vector<int> C;
    C.push_back(0);
 
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