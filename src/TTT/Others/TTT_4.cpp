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
// Graph is represented as a vector of unordered_set<int>
void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

// Optimized EXPAND using Tomita’s pivot strategy.
// This version computes Extu once and then updates it incrementally.
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
    // Pivot selection: choose u in SUBG that maximizes the count of neighbors in CAND.
    int u = -1, maxCount = 0;
    for (int x : SUBG) {
        int cnt = 0;
        for (int y : CAND) {
            if(adj[x].count(y))
                ++cnt;
        }
        if(cnt > maxCount) {
            maxCount = cnt;
            u = x;
        }
    }
    // Compute Extu = { v in CAND such that v is NOT adjacent to u } only once.
    unordered_set<int> Extu;
    Extu.reserve(CAND.size());
    for (int v : CAND) {
        if(!adj[u].count(v))
            Extu.insert(v);
    }
    
    // Process each q in Extu.
    while(!Extu.empty()){
        int q = *Extu.begin();
        Q.push_back(q);
        
        // Compute SUBGq = { v in SUBG that are adjacent to q }.
        unordered_set<int> SUBGq;
        SUBGq.reserve(SUBG.size());
        for (int v : SUBG)
            if(adj[q].count(v))
                SUBGq.insert(v);
        
        // Compute CANDq = { v in CAND that are adjacent to q }.
        unordered_set<int> CANDq;
        CANDq.reserve(CAND.size());
        for (int v : CAND)
            if(adj[q].count(v))
                CANDq.insert(v);
        
        // Recurse with move-semantics to avoid copying.
        EXPAND(move(SUBGq), move(CANDq), adj);
        
        // Remove q from CAND and from Extu.
        CAND.erase(q);
        Extu.erase(q);
        Q.pop_back();
        // Since CAND shrank, Extu remains valid as it is a subset of CAND.
    }
}
    
// CLIQUES launches the Tomita algorithm over the full vertex set.
// It also prints debug messages similar to the ELS variant.
void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    vector<int> ordering;
    for(int i = 0; i < V; i++){
        Vset.insert(i);
        ordering.push_back(i);
    }
    // Debug output: print processing order.
    for (int i = 0; i < ordering.size(); i++){
        cout << "[DEBUG] Processing node " << ordering[i]
             << " (loop index " << i << ")" << endl;
    }
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
    ofstream outfile("output_1.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0]=='#')
            continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            edgeList.push_back({u, v});
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    
    // Build graph as an adjacency list using unordered_set for quick lookups.
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
    for(const auto &p : cliqueSizeDistribution)
        outfile << "Size " << p.first << ": " << p.second << "\n";
    outfile.close();
    return 0;
}