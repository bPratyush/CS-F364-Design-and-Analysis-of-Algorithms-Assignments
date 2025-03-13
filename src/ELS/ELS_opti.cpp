#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <omp.h> // OpenMP header
using namespace std;
using namespace std::chrono;

// Eppstein, Löffler & Strash (2010) algorithm for finding all maximal cliques in an undirected graph
void addEdge(int u, int v, vector<unordered_set<int> >& adj) {
    if (u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

unordered_set<int> setintersect(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for (int a : A) {
        if (B.find(a) != B.end())
            res.insert(a);
    }
    return res;
}

unordered_set<int> setdiff(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for (int a : A) {
        if (B.find(a) == B.end())
            res.insert(a);
    }
    return res;
}

int maxCliqueSize = 0;
int totalMaximalCliques = 0;
unordered_map<int, int> cliqueSizeDistribution;

// Recursive Bron–Kerbosch algorithm with pivoting.
void BronKerboschPivot(unordered_set<int> P, unordered_set<int> R, unordered_set<int> X,
                         const vector<unordered_set<int> >& adj) {
    unordered_set<int> unionPX = P;
    unionPX.insert(X.begin(), X.end());
    
    // If P ∪ X is empty, R is a maximal clique.
    if (unionPX.empty()){
        int cliqueSize = R.size();
        #pragma omp critical
        {
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    }
    
    // Choose a pivot from P ∪ X arbitrarily.
    int u = *unionPX.begin();
    // Compute P \ neighbors(u)
    unordered_set<int> diff = setdiff(P, adj[u]);
    for (int v : diff) {
        unordered_set<int> newP = setintersect(P, adj[v]);
        unordered_set<int> newX = setintersect(X, adj[v]);
        unordered_set<int> newR = R;
        newR.insert(v);
        BronKerboschPivot(newP, newR, newX, adj);
        P.erase(v);
        X.insert(v);
    }
}

/* 
  Degeneracy ordering computed by repeatedly removing a vertex with the smallest
  degree (and updating neighbors' degrees accordingly) until the graph is empty.
*/
vector<int> degeneracyorder(const vector<unordered_set<int> >& adj) {
    int n = adj.size();
    vector<bool> used(n, false);
    vector<int> degree(n, 0);
    for (int i = 0; i < n; i++)
        degree[i] = adj[i].size();
    vector<int> ordering;
    for (int k = 0; k < n; k++) {
        int u = -1;
        int minDeg = INT_MAX;
        for (int i = 0; i < n; i++) {
            if (!used[i] && degree[i] < minDeg) {
                minDeg = degree[i];
                u = i;
            }
        }
        if (u == -1)
            break;
        used[u] = true;
        ordering.push_back(u);
        for (int w : adj[u]) {
            if (!used[w])
                degree[w]--;
        }
    }
    reverse(ordering.begin(), ordering.end());
    return ordering;
}

// Parallelized outer loop using OpenMP in the degeneracy-based Bron–Kerbosch algorithm.
void BronKerboschDegeneracy(const vector<unordered_set<int> >& adj) {
    int n = adj.size();
    vector<int> ordering = degeneracyorder(adj);
    vector<int> pos(n, 0);
    for (int i = 0; i < n; i++)
        pos[ordering[i]] = i;
    
    // The outer loop is parallelized. Dynamic scheduling helps balance uneven workloads.
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        int vi = ordering[i];
        unordered_set<int> P, X, R;
        // Build P: neighbors of vi with ordering greater than vi.
        for (int w : adj[vi]) {
            if (pos[w] > pos[vi])
                P.insert(w);
        }
        // Build X: neighbors of vi with ordering less than vi.
        for (int w : adj[vi]) {
            if (pos[w] < pos[vi])
                X.insert(w);
        }
        // R is initialized with {vi}.
        R.insert(vi);
        BronKerboschPivot(P, R, X, adj);
    }
}

int main(int argc, char* argv[]) {
    if(argc < 2) {
        cerr << "Usage: " << argv[0] << " [input file]" << endl;
        return 1;
    }
    
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    
    // Read input file; skip lines starting with '#'
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        int u, v;
        if (iss >> u >> v) {
            edgeList.push_back({u, v});
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    
    // Allocate and construct the adjacency list.
    vector<unordered_set<int>> adj(maxVertex + 1);
    for (auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    
    // Debug: Print the full adjacency list with sorted neighbors (optional).
    for (int i = 0; i < adj.size(); ++i) {
        cout << "Vertex " << i << ": ";
        if (adj[i].empty())
            cout << "No neighbors";
        else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for (int neighbor : neighbors)
                cout << neighbor << " ";
        }
        cout << endl;
    }
    
    auto start = high_resolution_clock::now();
    BronKerboschDegeneracy(adj);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    // Output results.
    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;
    for (const auto& p : cliqueSizeDistribution)
        outfile << "Size " << p.first << ": " << p.second << endl;
    
    outfile.close();
    return 0;
}