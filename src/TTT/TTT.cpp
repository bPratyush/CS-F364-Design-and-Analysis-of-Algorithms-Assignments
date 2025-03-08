#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
using namespace std;
using namespace std::chrono;

// Tomita, Tanata & Takahashi (2006) algorithm for finding all maximal cliques in an undirected graph
vector<int> Q;
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
unordered_map<int, int> cliqueSizeDistribution;

void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if (u >= adj.size() || v >= adj.size()) {
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

void EXPAND(unordered_set<int> SUBG, unordered_set<int> CAND, vector<unordered_set<int>>& adj) {
    if (SUBG.empty()) {
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
                if (adj[x].find(y) != adj[x].end()) ++cnt;
            }
            if (cnt > maxCount) {
                maxCount = cnt;
                u = x;
            }
        }
        unordered_set<int> Extu;
        for (int v : CAND) {
            if (adj[u].find(v) == adj[u].end()) Extu.insert(v);
        }
        unordered_set<int> FINI;
        while (!Extu.empty()) {
            int q = *Extu.begin();
            Q.push_back(q);
            unordered_set<int> SUBGq;
            for (int v : SUBG) {
                if (adj[q].find(v) != adj[q].end()) SUBGq.insert(v);
            }
            unordered_set<int> CANDq;
            for (int v : CAND) {
                if (adj[q].find(v) != adj[q].end()) CANDq.insert(v);
            }
            EXPAND(SUBGq, CANDq, adj);
            CAND.erase(q);
            FINI.insert(q);
            Q.pop_back();
            Extu.clear();
            for (int v : CAND) {
                if (adj[u].find(v) == adj[u].end()) Extu.insert(v);
            }
        }
    }
}

void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    for (int i = 0; i < V; i++) Vset.insert(i);
    EXPAND(Vset, Vset, adj);
}

int main(int argc, char* argv[]) {
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    while (getline(infile, line)) {
        if (line[0] != '#') break;
    }
    istringstream iss(line);
    int V, E;
    iss >> V >> E;
    vector<unordered_set<int>> adj(V);
    int u, v;
    while (infile >> u >> v) addEdge(u, v, adj);
    infile.close();

    // Debug: Print the adjacency list
    for (int i = 0; i < adj.size(); ++i) {
        cout << "Vertex " << i << ": ";
        for (int neighbor : adj[i]) {
            cout << neighbor << " ";
        }
        cout << endl;
    }

    auto start = high_resolution_clock::now();
    CLIQUES(adj, V);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;

    for (const auto& pair : cliqueSizeDistribution) outfile << "Size " << pair.first << ": " << pair.second << endl;
    outfile.close();
    return 0;
}