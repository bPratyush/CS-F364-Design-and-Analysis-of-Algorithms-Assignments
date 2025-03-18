#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <thread>
#include <mutex>
using namespace std;
using namespace std::chrono;

// --- Original Helper Functions (unchanged) ---

void addEdge(int u, int v, vector<unordered_set<int> >& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

unordered_set<int> setintersect(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for(int a : A) {
        if(B.find(a) != B.end())
            res.insert(a);
    }
    return res;
}

unordered_set<int> setdiff(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for(int a : A) {
        if(B.find(a) == B.end())
            res.insert(a);
    }
    return res;
}

// Global clique statistics
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
unordered_map<int, int> cliqueSizeDistribution;
mutex clique_mutex;  // Protects updates below

// Original BronKerboschPivot – modified to lock when reporting a clique.
void BronKerboschPivot(unordered_set<int> P, unordered_set<int> R, unordered_set<int> X, const vector<unordered_set<int> >& adj) {
    unordered_set<int> unionPX = P;
    unionPX.insert(X.begin(), X.end());
    if(unionPX.empty()){
        int cliqueSize = R.size();
        if(cliqueSize > 1) {
            lock_guard<mutex> lock(clique_mutex);
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    }
    int u = *unionPX.begin();
    unordered_set<int> diff = setdiff(P, adj[u]);
    for(int v : diff) {
        unordered_set<int> newP = setintersect(P, adj[v]);
        unordered_set<int> newX = setintersect(X, adj[v]);
        unordered_set<int> newR = R;
        newR.insert(v);
        BronKerboschPivot(newP, newR, newX, adj);
        P.erase(v);
        X.insert(v);
    }
}

// This function (degeneracyorder) is not used in the parallel version
// since we manually compute the degeneracy ordering below.
vector<int> degeneracyorder(const vector<unordered_set<int> >& adj){
    int n = adj.size();
    vector<bool> used(n, false);
    vector<int> degree(n, 0);
    for(int i = 0; i < n; i++)
        degree[i] = adj[i].size();
    vector<int> ordering;
    for(int k = 0; k < n; k++){
        int u = -1;
        int minDeg = INT_MAX;
        for(int i = 0; i < n; i++){
            if(!used[i] && degree[i] < minDeg){
                minDeg = degree[i];
                u = i;
            }
        }
        if(u == -1)
            break;
        used[u] = true;
        ordering.push_back(u);
        for(int w : adj[u]){
            if(!used[w])
                degree[w]--;
        }
    }
    reverse(ordering.begin(), ordering.end());
    return ordering;
}

// --- Main Function (Modified for Parallelism with limited threads) ---

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    auto start = high_resolution_clock::now();
    
    // Read input file and build graph.
    ifstream ip("wiki-Vote.txt");
    if(!ip){
        cerr << "Failed to open file." << endl;
        return 1;
    }
    
    unordered_set<int> nodes;
    int a, b;
    int maxVertex = 0;
    string line;
    while(getline(ip, line)){
        if(line.empty() || line[0] == '#')
            continue;
        istringstream ss(line);
        ss >> a >> b;
        if(!ss.fail()){
            maxVertex = max(maxVertex, max(a, b));
            nodes.insert(a);
            nodes.insert(b);
        }
    }
    ip.clear();
    ip.seekg(0, ios::beg);  // Reset stream to reread
    // Initialize graph vector and read edges using addEdge.
    vector<unordered_set<int>> adj(maxVertex + 1);
    while(getline(ip, line)){
        if(line.empty() || line[0] == '#')
            continue;
        istringstream ss(line);
        ss >> a >> b;
        if(!ss.fail()){
            addEdge(a, b, adj);
            nodes.insert(a);
            nodes.insert(b);
        }
    }
    ip.close();
    
    cout << "Adjacency list created." << endl;
    
    // --- Compute degeneracy ordering using the same manual method as before ---
    unordered_map<int, unordered_set<int>> tempadj;
    for (int i = 0; i < adj.size(); i++){
        tempadj[i] = adj[i];
    }
    
    unordered_map<int, int> deg;
    unordered_map<int, unordered_set<int>> deglist;
    for (int node : nodes) {
        int d = 0;
        if(node < adj.size())
            d = adj[node].size();
        deg[node] = d;
        deglist[d].insert(node);
    }
    cout << "Degree list created." << endl;
    
    int global_maxdeg = 0;
    for(auto &p : deg){
        if(p.second > global_maxdeg)
            global_maxdeg = p.second;
    }
    
    vector<int> degenlist;
    degenlist.reserve(nodes.size());
    int minDegree = 0;
    while(degenlist.size() < nodes.size()){
        while(minDegree <= global_maxdeg && deglist[minDegree].empty()){
            minDegree++;
        }
        if(minDegree > global_maxdeg)
            break; 
        int u = *deglist[minDegree].begin();
        deglist[minDegree].erase(u);
        degenlist.push_back(u);
        for(auto v : tempadj[u]){
            tempadj[v].erase(u);
            int curr = deg[v];
            deglist[curr].erase(v);
            deg[v]--;
            int newdeg = deg[v];
            deglist[newdeg].insert(v);
            if(newdeg < minDegree)
                minDegree = newdeg;
        }
    }
    cout << "Degeneracy ordering created." << endl;
    
    // Build a lookup for positions: pos[v] is the position of v in degenlist.
    unordered_map<int, int> pos;
    for(int i = 0; i < degenlist.size(); i++){
        pos[degenlist[i]] = i;
    }
    
    // --- Parallelize the outer loop of BronKerboschDegeneracy ---
    int totalNodes = degenlist.size();
    // Instead of creating one thread per small chunk, use a fixed number of threads.
    unsigned int nThreads = thread::hardware_concurrency();
    if(nThreads == 0) nThreads = 4;
    // Partition indices among threads.
    vector<thread> workers;
    int blockSize = totalNodes / nThreads;
    int remainder = totalNodes % nThreads;
    mutex print_mutex; // Protects concurrent printing to cout
    
    int start_idx = 0;
    for(unsigned int t = 0; t < nThreads; t++){
        int end_idx = start_idx + blockSize + (t < remainder ? 1 : 0);
        workers.emplace_back([start_idx, end_idx, &degenlist, &pos, &adj, &print_mutex](){
            for (int i = start_idx; i < end_idx; i++){
                { // Print in a critical section.
                    lock_guard<mutex> guard(print_mutex);
                    cout << "node " << i << " being processed" << endl;
                }
                int vi = degenlist[i];
                unordered_set<int> P, X, R;
                // For each neighbor of vi, use pos map.
                for(auto neighbor : adj[vi]){
                    if(pos.find(neighbor) != pos.end()){
                        if(pos[neighbor] > i)
                            P.insert(neighbor);
                        else if(pos[neighbor] < i)
                            X.insert(neighbor);
                    }
                }
                R.insert(vi);
                BronKerboschPivot(P, R, X, adj);
            }
        });
        start_idx = end_idx;
    }
    
    for(auto &t : workers){
        t.join();
    }
    
    // --- Output Results ---
    ofstream op("output.txt");
    if(!op){
        cerr << "Failed to open output file." << endl;
        return 1;
    }
    
    op << "\nNo. of maximal cliques: " << totalMaximalCliques << "\n";
    op << "\nClique size distribution:\n";
    for(auto &p : cliqueSizeDistribution){
        op << "Clique size: " << p.first << ", Frequency: " << p.second << "\n";
    }
    op << "Maximum clique size: " << maxCliqueSize << "\n";
    op.close();
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start);
    cout << "Number of nodes: " << nodes.size() << "\n";
    cout << "No. of elements in degeneracy ordering: " << degenlist.size() << "\n";
    cout << "Time taken: " << duration.count() << " milliseconds" << "\n";
    
    return 0;
}