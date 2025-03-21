#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include <chrono>

using namespace std;
using namespace std::chrono;

int n; // Number of vertices
vector< set<int> > G; // Graph adjacency list

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
        allVertices.insert(i);
    }
    findCliques(emptySet, allVertices);
}

// Main function
int main() {
    n = 5;
    G.assign(n + 1, set<int>());

    // Graph input (edges)
    G[1] = {2, 3};
    G[2] = {1, 3,4};
    G[3] = {1,2,5};
    G[4] = {5};
    G[5] = {4};

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
