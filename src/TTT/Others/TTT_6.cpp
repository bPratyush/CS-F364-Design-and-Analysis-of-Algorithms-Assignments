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
using namespace std;
using namespace std::chrono;

// Global variables for Tomita algorithm
vector<int> Q;    // current clique (recursion stack)
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
map<int, int> cliqueSizeDistribution;

// Add an edge into the graph (0-indexed vertices)
void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

// Optimized EXPAND using Tomita’s pivot strategy with iterator-based loops.
// Temporary containers reserve capacity to reduce reallocations.
void EXPAND(unordered_set<int> SUBG, unordered_set<int> CAND, vector<unordered_set<int>>& adj) {
    if(SUBG.empty()) {
        int cliqueSize = Q.size();
        if(cliqueSize > 1) {
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    } 
    // Pivot selection: choose u in SUBG that maximizes count(neighbors in CAND)
    int u = -1, maxCount = 0;
    for(auto it = SUBG.begin(); it != SUBG.end(); ++it) {
        int cnt = 0;
        for(auto jt = CAND.begin(); jt != CAND.end(); ++jt) {
            if(adj[*it].find(*jt) != adj[*it].end())
                ++cnt;
        }
        if(cnt > maxCount) {
            maxCount = cnt;
            u = *it;
        }
    }
    // Compute Extu = { v in CAND such that v is NOT adjacent to u }
    unordered_set<int> Extu;
    Extu.reserve(CAND.size());
    for(auto it = CAND.begin(); it != CAND.end(); ++it) {
        if(adj[u].find(*it) == adj[u].end())
            Extu.insert(*it);
    }
    
    // Process each vertex q in Extu
    while(!Extu.empty()){
        int q = *Extu.begin();
        Q.push_back(q);
        
        // Build SUBGq = { v in SUBG : v adjacent to q }
        unordered_set<int> SUBGq;
        SUBGq.reserve(SUBG.size());
        for(auto it = SUBG.begin(); it != SUBG.end(); ++it) {
            if(adj[q].find(*it) != adj[q].end())
                SUBGq.insert(*it);
        }
        
        // Build CANDq = { v in CAND : v adjacent to q }
        unordered_set<int> CANDq;
        CANDq.reserve(CAND.size());
        for(auto it = CAND.begin(); it != CAND.end(); ++it) {
            if(adj[q].find(*it) != adj[q].end())
                CANDq.insert(*it);
        }
        
        // Recurse using fully-qualified std::move to avoid redundant copying.
        EXPAND(std::move(SUBGq), std::move(CANDq), adj);
        
        // Remove q from CAND and from Extu.
        CAND.erase(q);
        Extu.erase(q);
        Q.pop_back();
        // Extu is now a subset of the (shrunk) CAND; no full recomputation is needed.
    }
}
    
// CLIQUES launches the Tomita algorithm over the full vertex set.
// Debug output is disabled when NDEBUG is defined.
void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    vector<int> ordering;
    for (int i = 0; i < V; i++){
        Vset.insert(i);
        ordering.push_back(i);
    }
#ifndef NDEBUG
    for(auto it = ordering.begin(); it != ordering.end(); ++it) {
        cout << "[DEBUG] Processing node " << *it 
             << " (loop index " << distance(ordering.begin(), it) << ")" << endl;
    }
#endif
    // Initially, both SUBG and CAND are the full vertex set.
    EXPAND(Vset, Vset, adj);
}
    
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " [input file]" << endl;
        return 1;
    }
    
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            edgeList.push_back({u, v});
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    
    // Build graph as an adjacency list using unordered_set.
    vector<unordered_set<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    
    auto start = high_resolution_clock::now();
    CLIQUES(adj, maxVertex + 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    outfile << "Largest size of the clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << "\n";
    outfile << "Execution time (ms): " << duration.count() << "\n";
    outfile << "Distribution of different size cliques:\n";
    for(auto it = cliqueSizeDistribution.begin(); it != cliqueSizeDistribution.end(); ++it) {
        outfile << "Size " << it->first << ": " << it->second << "\n";
    }
    outfile.close();
    return 0;
}