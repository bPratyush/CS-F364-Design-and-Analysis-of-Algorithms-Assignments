#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <sstream>

using namespace std;
using namespace std::chrono;

int n; // Number of vertices
vector< set<int> > G; // Adjacency list for the undirected graph

// Global statistics
long long total_maximal_cliques = 0;
int largest_clique_size = 0;
map<int, int> clique_distribution;

// Utility function to check if a clique is maximal
bool isMaximal(const set<int>& C) {
    for (int u = 1; u <= n; u++) {
        if (C.find(u) == C.end()) {
            bool isNeighborToAll = true;
            for (int v : C) {
                if (G[u].find(v) == G[u].end()) {
                    isNeighborToAll = false;
                    break;
                }
            }
            if (isNeighborToAll) return false; // Not maximal
        }
    }
    return true;
}

// Process a maximal clique
void processClique(const set<int>& C) {
    if (!isMaximal(C)) return; // Ensure only maximal cliques are counted
    int size = C.size();
    if (size < 2) return;
    total_maximal_cliques++;
    largest_clique_size = max(largest_clique_size, size);
    clique_distribution[size]++;
}

// Recursive function for finding maximal cliques
void findCliques(set<int> C, set<int> candidates) {
    if (candidates.empty()) {
        processClique(C);
        return;
    }

    while (!candidates.empty()) {
        int v = *candidates.begin();
        candidates.erase(v);

        set<int> newC = C;
        newC.insert(v);

        set<int> newCandidates;
        for (int u : candidates) {
            if (G[v].find(u) != G[v].end()) {
                newCandidates.insert(u);
            }
        }

        findCliques(newC, newCandidates);
    }
}

// Main function to start finding maximal cliques
void CLIQUE() {
    set<int> emptySet;
    set<int> allVertices;
    for (int i = 1; i <= n; i++) {
        if (!G[i].empty()) { // Only consider nodes with edges
            allVertices.insert(i);
        }
    }
    findCliques(emptySet, allVertices);
}

// Read input file and construct the undirected graph
void loadGraph(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }

    string line;
    set<int> nodes;

    // Read edges
    while (getline(file, line)) {
        if (line[0] == '#') continue; // Ignore comments
        stringstream ss(line);
        int u, v;
        if (ss >> u >> v) {
            G[u].insert(v);
            G[v].insert(u); // Convert to undirected
            nodes.insert(u);
            nodes.insert(v);
        }
    }

    file.close();

    n = nodes.size(); // Set number of vertices
    cout << "Graph loaded: " << n << " nodes." << endl;
}

// Main function
int main() {
    string filename = "/Users/bpratyush/Documents/GitHub/CS-F364-Design-and-Analysis-of-Algorithms-Assignments/src/CT/wiki-Vote.txt";
    
    // Allocate large graph adjacency list
    G.assign(40000, set<int>()); // Slightly larger than 7115 to handle index shifts

    loadGraph(filename);

    auto start = high_resolution_clock::now();
    CLIQUE();
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();

    cout << "Largest size of the clique: " << largest_clique_size << "\n";
    cout << "Total number of maximal cliques: " << total_maximal_cliques << "\n";
    cout << "Execution time (ms): " << duration << "\n";
    cout << "Distribution of different size cliques:\n";
    for (auto &entry : clique_distribution) {
        cout << "Size " << entry.first << ": " << entry.second << "\n";
    }

    return 0;
}