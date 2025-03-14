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

// Tomita, Tanaka & Takahashi (2006) algorithm for finding all maximal cliques in an undirected graph
vector<int> Q;
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
map<int, int> cliqueSizeDistribution;
// Global vector to store all maximal cliques (each as a sorted vector of vertices).
vector<vector<int>> allCliques;

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
        // Save the current clique: sort Q and store a copy.
        vector<int> clique = Q;
        sort(clique.begin(), clique.end());
        allCliques.push_back(clique);
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

// Compute degeneracy ordering by repeatedly removing the vertex of minimum degree.
vector<int> degeneracyorder(const vector<unordered_set<int>>& adj) {
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
    // (Optional) reverse(ordering.begin(), ordering.end());
    return ordering;
}

// New function using degeneracy ordering to drive the search.
void CLIQUES_DEGEN(vector<unordered_set<int>>& adj) {
    vector<int> ordering = degeneracyorder(adj);
    int n = adj.size();
    // Build positions to compare ordering.
    vector<int> pos(n, 0);
    for (int i = 0; i < ordering.size(); i++) {
        pos[ordering[i]] = i;
    }
    // For each vertex v in degeneracy order, search for cliques where v is the smallest member.
    for (int i = 0; i < ordering.size(); i++) {
        int v = ordering[i];
        unordered_set<int> cand;
        for(auto w : adj[v]) {
            if(pos[w] > pos[v])
                cand.insert(w);
        }
        unordered_set<int> subg = cand;
        Q.clear();
        Q.push_back(v);
        EXPAND(subg, cand, adj);
    }
}

int main(int argc, char* argv[]) {
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
    
    // Allocate adjacency list based on the maximum vertex found.
    vector<unordered_set<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    
    // Debug: Print the full adjacency list with sorted neighbors.
    for (int i = 0; i < adj.size(); ++i) {
        cout << "Vertex " << i << ": ";
        if(adj[i].empty()){
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for (int neighbor : neighbors)
                cout << neighbor << " ";
        }
        cout << endl;
        cout.flush();
    }
    
    auto start = high_resolution_clock::now();
    // Use degeneracy ordering based CLIQUES_DEGEN
    CLIQUES_DEGEN(adj);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;
    for(const auto& pair : cliqueSizeDistribution)
        outfile << "Size " << pair.first << ": " << pair.second << endl;
    
    outfile << "\nMaximal Cliques:" << endl;
    for (const auto& clique : allCliques) {
        for (int v : clique)
            outfile << v << " ";
        outfile << "\n";
    }
    outfile.close();
    
    return 0;
}