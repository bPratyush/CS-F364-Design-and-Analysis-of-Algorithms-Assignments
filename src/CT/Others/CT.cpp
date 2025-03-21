#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <unordered_map>
using namespace std;
using namespace std::chrono;

// Chiba & Nishizeki (1985) algorithm for finding all maximal cliques in a graph
// We assume vertices are numbered from 1 to n (index 0 is unused).

// (Optional) Mapping of new vertex labels to original values.
// For now we simply print the vertex as is.
vector<int> newToOld;

void printclique(const unordered_set<int>& C) {
    cout << "Clique: ";
    for (int v : C)
        cout << v << " ";
    cout << endl;
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

// UPDATE procedure (recursive) for clique enumeration.
// i runs from 1 to n (inclusive)
void UPDATE(int i, unordered_set<int>& C, int n, const vector<unordered_set<int> >& adj, vector<int>& S, vector<int>& T) {
    if (i == n + 1) {
        printclique(C);
        int sizeC = C.size();
        maxCliqueSize = max(maxCliqueSize, sizeC);
        totalMaximalCliques++;
        cliqueSizeDistribution[sizeC]++;
        return;
    }
    // diff = C \ adj[i]: vertices in C that are not adjacent to i.
    unordered_set<int> diff = setdiff(C, adj[i]);
    if (!diff.empty()) {
        UPDATE(i + 1, C, n, adj, S, T);
    }
    // For each vertex y not in C and not equal to i, set T[y] = |adj[y] ∩ diff|
    for (int y = 1; y <= n; y++) {
        if (C.find(y) == C.end() && y != i) {
            unordered_set<int> inter = setintersect(adj[y], diff);
            T[y] = inter.size();
        }
    }
    // cap = C ∩ adj[i]: vertices in C adjacent to i.
    unordered_set<int> cap = setintersect(C, adj[i]);
    for (int x : cap) {
        unordered_set<int> Nx = setdiff(adj[x], C);
        if (Nx.find(i) != Nx.end())
            Nx.erase(i);
        for (int y : Nx) {
            unordered_set<int> inter = setintersect(adj[y], diff);
            T[y] += inter.size();
        }
    }
    for (int x : diff) {
        unordered_set<int> Nx = setdiff(adj[x], C);
        for (int y : Nx) {
            S[y]++;
        }
    }
    bool FLAG = true;
    int capSize = cap.size();
    for (int y : adj[i]) {
        if (C.find(y) == C.end() && y < i && T[y] == capSize) {
            FLAG = false;
            break;
        }
    }
    vector<int> diffVec(diff.begin(), diff.end());
    sort(diffVec.begin(), diffVec.end());
    for (size_t k = 0; k < diffVec.size(); k++) {
        int j = diffVec[k];
        unordered_set<int> Nj = setdiff(adj[j], C);
        for (int y : Nj) {
            if (y < i && T[y] == capSize) {
                if (adj[y].find(j) != adj[y].end()) {
                    S[y]--;
                } else if (k == 0 && S[y] > 0) {
                    FLAG = false;
                }
            }
        }
    }
    if (!cap.empty() && !diffVec.empty()) {
        int jp = diffVec.back();
        for (int y : cap) {
            if (y < i && T[y] == (int)adj[i].size() && S[y] == 0) {
                if (jp < y || jp < i - 1) {
                    FLAG = false;
                    break;
                }
            }
        }
    }
    // Reset counters in T for vertices in cap and S for vertices in diff.
    for (int x : cap) {
        for (int y : adj[x]) {
            if (C.find(y) == C.end() && y != i)
                T[y] = 0;
        }
    }
    for (int x : diff) {
        for (int y : adj[x]) {
            if (C.find(y) == C.end())
                S[y] = 0;
        }
    }
    if (FLAG) {
        unordered_set<int> SAVE = diff;
        unordered_set<int> newC = cap;
        newC.insert(i);
        unordered_set<int> oldC = C;
        C = newC;
        UPDATE(i + 1, C, n, adj, S, T);
        // Recover C: remove i and reinsert SAVE
        unordered_set<int> recovered;
        for (int x : C) {
            if (x != i)
                recovered.insert(x);
        }
        recovered.insert(SAVE.begin(), SAVE.end());
        C = recovered;
    }
}

void addEdge(int u, int v, vector<unordered_set<int> >& adj) {
    if (u >= adj.size() || v >= adj.size()) {
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    // Treat edge as undirected
    adj[u].insert(v);
    adj[v].insert(u);
}

int main(int argc, char* argv[]){
    if(argc < 2) {
        cerr << "Usage: " << argv[0] << " [input file]\n";
        return 1;
    }
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    // Read input (ignoring blank and comment lines).
    // Adjusted for 0-based input (we add 1 to each vertex).
    while(getline(infile, line)){
        if(line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            u++; v++; // convert from 0-based to 1-based indexing
            edgeList.push_back({u, v});
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    
    // Allocate adjacency list based on the maximum vertex found.
    // We allocate size maxVertex+1 and use indices 1..maxVertex (ignore index 0).
    vector<unordered_set<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    
    // Debug: Print the full adjacency list with sorted neighbors.
    for(int i = 1; i < adj.size(); ++i){
        cout << "Vertex " << i << ": ";
        if(adj[i].empty()){
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for (int neighbor : neighbors) {
                cout << neighbor << " ";
            }
        }
        cout << endl;
        cout.flush();
    }
    
    // Begin clique enumeration.
    auto start = high_resolution_clock::now();
    {
        int n = maxVertex;  // vertices 1..n
        unordered_set<int> C; // initial clique is empty
        vector<int> S(n + 1, 0), T(n + 1, 0);
        // Call the UPDATE procedure starting with i = 1.
        UPDATE(1, C, n, adj, S, T);
    }
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