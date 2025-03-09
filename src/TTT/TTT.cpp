#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <map> // for consistent ordering
using namespace std;
using namespace std::chrono;

// Tomita, Tanata & Takahashi (2006) algorithm for finding all maximal cliques in an undirected graph
vector<int> Q;
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
map<int, int> cliqueSizeDistribution;

void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    // Treat edge as undirected
    adj[u].insert(v);
    adj[v].insert(u);
}

void EXPAND(unordered_set<int> SUBG, unordered_set<int> CAND, vector<unordered_set<int>>& adj) {
    if(SUBG.empty()){
        int cliqueSize = Q.size();
        maxCliqueSize = max(maxCliqueSize, cliqueSize);
        totalMaximalCliques++;
        cliqueSizeDistribution[cliqueSize]++;
        return;
    } else {
        int u = -1;
        int maxCount = 0;
        for (int x : SUBG) {
            int cnt = 0;
            for (int y : CAND) {
                if(adj[x].find(y) != adj[x].end())
                    ++cnt;
            }
            if(cnt > maxCount){
                maxCount = cnt;
                u = x;
            }
        }
        unordered_set<int> Extu;
        for (int v : CAND) {
            if(adj[u].find(v) == adj[u].end())
                Extu.insert(v);
        }
        unordered_set<int> FINI;
        while(!Extu.empty()){
            int q = *Extu.begin();
            Q.push_back(q);
            unordered_set<int> SUBGq;
            for (int v : SUBG) {
                if(adj[q].find(v) != adj[q].end())
                    SUBGq.insert(v);
            }
            unordered_set<int> CANDq;
            for (int v : CAND) {
                if(adj[q].find(v) != adj[q].end())
                    CANDq.insert(v);
            }
            EXPAND(SUBGq, CANDq, adj);
            CAND.erase(q);
            FINI.insert(q);
            Q.pop_back();
            Extu.clear();
            for (int v : CAND) {
                if(adj[u].find(v) == adj[u].end())
                    Extu.insert(v);
            }
        }
    }
}

void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    for(int i = 0; i < V; i++)
        Vset.insert(i);
    EXPAND(Vset, Vset, adj);
}

int main(int argc, char* argv[]) {
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        exit(1);
    }
    
    ifstream infile(argv[1]);
    if(!infile.is_open()){
        cerr << "Error opening file " << argv[1] << endl;
        exit(1);
    }
    
    ofstream outfile("output.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    
    // Read file line by line and ignore lines beginning with '#'
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
    
    // Allocate adjacency list based on the maximum vertex found
    vector<unordered_set<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    
    // Debug: Print the full adjacency list with sorted neighbors for consistency
    for(int i = 0; i < adj.size(); ++i){
        cout << "Vertex " << i << ": ";
        if(adj[i].empty()){
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for(int neighbor : neighbors){
                cout << neighbor << " ";
            }
        }
        cout << endl;
        cout.flush();
    }
    
    auto start = high_resolution_clock::now();
    // Use maxVertex+1 as the vertex count
    CLIQUES(adj, maxVertex + 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;
    
    for(const auto& pair : cliqueSizeDistribution)
        outfile << "Size " << pair.first << ": " << pair.second << endl;
    outfile.close();
    
    return 0;
}