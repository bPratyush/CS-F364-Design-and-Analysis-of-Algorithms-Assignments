#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <map>
using namespace std;
using namespace std::chrono;
 
// Global variables to collect the clique metrics
int totalMaximalCliques = 0;
int maxCliqueSize = 0;
map<int, int> cliqueSizeDistribution;
 
// Function to add an edge to an undirected graph
void addEdge(int u, int v, vector<unordered_set<int> > &adj) {
    adj[u].insert(v);
    adj[v].insert(u);
}
 
// The Bron-Kerbosch algorithm with pivot selection
void bronKerbosch2(unordered_set<int> R, unordered_set<int> P, unordered_set<int> X, vector<unordered_set<int> > &adj) {
    if (P.empty() && X.empty()) {
        cout << "Maximal Clique: ";
        for (int v : R)
            cout << v << " ";
        cout << endl;
        // Update clique metrics
        totalMaximalCliques++;
        int cliqueSize = R.size();
        if (cliqueSize > maxCliqueSize)
            maxCliqueSize = cliqueSize;
        cliqueSizeDistribution[cliqueSize]++;
        return;
    }
    // Compute P ∪ X
    unordered_set<int> unionPX = P;
    unionPX.insert(X.begin(), X.end());
    
    // Choose pivot u in unionPX that maximizes |P ∩ N(u)|
    int pivot = -1;
    int maxNeighbors = -1;
    for (int u : unionPX) {
        int count = 0;
        for (int v : P) {
            if (adj[u].find(v) != adj[u].end())
                count++;
        }
        if (count > maxNeighbors) {
            maxNeighbors = count;
            pivot = u;
        }
    }
    
    // Compute vertices in P that are not neighbors of pivot.
    unordered_set<int> candi;
    for (int v : P) {
        if (adj[pivot].find(v) == adj[pivot].end())
            candi.insert(v);
    }
    for (int v : candi) {
        unordered_set<int> newR = R;
        newR.insert(v);
        unordered_set<int> newP, newX;
        for (int n : adj[v]) {
            if (P.find(n) != P.end())
                newP.insert(n);
            if (X.find(n) != X.end())
                newX.insert(n);
        }
        bronKerbosch2(newR, newP, newX, adj);
        P.erase(v);
        X.insert(v);
    }
}
 
// CLIQUES: Initializes the Bron-Kerbosch algorithm call 
void CLIQUES(vector<unordered_set<int> > &adj, int V) {
    unordered_set<int> R, X, P;
    for (int i = 0; i < V; i++) {
        P.insert(i);
    }
    bronKerbosch2(R, P, X, adj);
}
 
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        exit(1);
    }
    
    ifstream infile(argv[1]);
    if (!infile.is_open()) {
        cerr << "Error opening file " << argv[1] << endl;
        exit(1);
    }
    
    ofstream outfile("output_alt.txt");
    string line;
    int u, v;
    vector<pair<int,int>> edgeList;
    int maxVertex = 0;
    
    // Read the file line by line. Ignore all lines starting with '#'
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        if (iss >> u >> v) {
            edgeList.push_back({u, v});
            // Update maximum vertex id.
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    
    // Allocate adjacency list based on the maximum vertex found
    vector<unordered_set<int>> adj(maxVertex + 1);
    for (auto &edge : edgeList) {
        addEdge(edge.first, edge.second, adj);
    }
    
    // Debug: Print the full adjacency list with sorted neighbors for consistency
    for (int i = 0; i < (int)adj.size(); ++i) {
        cout << "Vertex " << i << ": ";
        if (adj[i].empty()) {
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for (int neighbor : neighbors) {
                cout << neighbor << " ";
            }
        }
        cout << endl;
    }
    
    auto start = high_resolution_clock::now();
    CLIQUES(adj, maxVertex + 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;
    for (const auto& pair : cliqueSizeDistribution)
        outfile << "Size " << pair.first << ": " << pair.second << endl;
    outfile.close();
    
    return 0;
}