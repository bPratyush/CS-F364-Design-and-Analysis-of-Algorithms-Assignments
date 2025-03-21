#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <functional>
using namespace std;

// This code reads the wiki-Vote dataset (ignoring comment lines),
// builds an undirected graph, computes a degeneracy ordering, and
// then counts triangles using Chiba and Nishizeki’s (1985) idea.

// New: Bron–Kerbosch algorithm (pivot version) to count maximal cliques.
void bronKerbosch(unordered_set<int> R, unordered_set<int> P, unordered_set<int> X, 
                    const unordered_map<int, unordered_set<int>> &graph, 
                    long long &maximalCliqueCount) {
    if (P.empty() && X.empty()) {
        maximalCliqueCount++;
        return;
    }
    int u = -1;
    unordered_set<int> unionPX;
    for (int v : P) unionPX.insert(v);
    for (int v : X) unionPX.insert(v);
    if (!unionPX.empty())
        u = *unionPX.begin();
    unordered_set<int> P_without_neighbors;
    if (u != -1) {
        for (int v : P) {
            if (graph.at(u).find(v) == graph.at(u).end())
                P_without_neighbors.insert(v);
        }
    } else {
        P_without_neighbors = P;
    }
    for (int v : P_without_neighbors) {
        unordered_set<int> R_new = R;
        R_new.insert(v);
        unordered_set<int> P_new, X_new;
        for (int w : P) {
            if (graph.at(v).find(w) != graph.at(v).end())
                P_new.insert(w);
        }
        for (int w : X) {
            if (graph.at(v).find(w) != graph.at(v).end())
                X_new.insert(w);
        }
        bronKerbosch(R_new, P_new, X_new, graph, maximalCliqueCount);
        P.erase(v);
        X.insert(v);
    }
}

int main() {
    string filename = "/Users/bpratyush/Documents/GitHub/CS-F364-Design-and-Analysis-of-Algorithms-Assignments/src/CT/wiki-Vote.txt";
    ifstream infile(filename);
    if (!infile) {
        cerr << "Cannot open file: " << filename << "\n";
        return 1;
    }
  
    // Build the undirected graph
    unordered_map<int, unordered_set<int>> graph;
    string line;
    int u, v;
    while(getline(infile, line)) {
        if (line.empty() || line[0]=='#')
            continue;
        istringstream iss(line);
        if (!(iss >> u >> v))
            continue;
        // Treat the input as undirected.
        graph[u].insert(v);
        graph[v].insert(u);
    }
  
    // Get the list of vertices.
    vector<int> vertices;
    for (const auto &p : graph)
        vertices.push_back(p.first);
  
    // Compute initial degrees.
    unordered_map<int, int> degree;
    for (const auto &p: graph)
        degree[p.first] = p.second.size();
  
    // Determine the maximum degree.
    int max_deg = 0;
    for (const auto &p: degree)
        max_deg = max(max_deg, p.second);
  
    // Set up buckets to support the degeneracy ordering
    vector<vector<int>> buckets(max_deg + 1);
    for (const auto &p: degree) {
        buckets[p.second].push_back(p.first);
    }
  
    vector<int> ordering;
    ordering.reserve(vertices.size());
    unordered_map<int, bool> removed;
    // Initially none is removed.
    for (const auto &p: graph)
        removed[p.first] = false;
  
    // Compute degeneracy order (smallest degree first)
    int n = vertices.size();
    for (int i = 0; i < n; i++) {
        int cur = -1;
        // Find a vertex that is not removed.
        for (int d = 0; d <= max_deg; d++) {
            while (!buckets[d].empty()) {
                int candidate = buckets[d].back();
                buckets[d].pop_back();
                if (!removed[candidate]) {
                    cur = candidate;
                    break;
                }
            }
            if (cur != -1)
                break;
        }
        if (cur == -1)
            break;
        removed[cur] = true;
        ordering.push_back(cur);
        // For each neighbor not removed, lower its degree.
        for (int nb : graph[cur]) {
            if (!removed[nb]) {
                int d_old = degree[nb];
                degree[nb] = d_old - 1;
                // Place nb into the bucket for degree d_old - 1.
                buckets[d_old - 1].push_back(nb);
            }
        }
    }
  
    // Map each vertex to its position in the ordering.
    unordered_map<int, int> pos;
    for (int i = 0; i < ordering.size(); i++) {
        pos[ordering[i]] = i;
    }
  
    // Count triangles; for each vertex u, consider only those neighbors v with pos[u] < pos[v]
    // then count intersections among these "forward" neighbors.
    long long triangleCount = 0;
    for (int u : ordering) {
        vector<int> Nu;
        for (int w : graph[u]) {
            if (pos[u] < pos[w])
                Nu.push_back(w);
        }
        sort(Nu.begin(), Nu.end(), [&](int a, int b) { return pos[a] < pos[b]; });
        for (size_t i = 0; i < Nu.size(); i++) {
            for (size_t j = i + 1; j < Nu.size(); j++) {
                int v = Nu[i], w = Nu[j];
                if (graph[v].find(w) != graph[v].end())
                    triangleCount++;
            }
        }
    }
  
    // --- New Part: Count maximal cliques using Bron–Kerbosch ---
    long long maximalCliqueCount = 0;
    {
        unordered_set<int> R, P, X;
        for (const auto &p : graph)
            P.insert(p.first);
        bronKerbosch(R, P, X, graph, maximalCliqueCount);
    }
  
    // --- New Part: Count total cliques of size ≥ 3 (unique cliques via ordering) ---
    vector<int> verticesSorted;
    for (const auto &p: graph)
        verticesSorted.push_back(p.first);
    sort(verticesSorted.begin(), verticesSorted.end());
    long long totalCliquesCount = 0;
    // Recursive lambda to extend cliques.
    function<void(vector<int>&, vector<int>&)> extendClique = [&](vector<int> &clique, vector<int> &cand) {

            totalCliquesCount++;
        for (size_t i = 0; i < cand.size(); i++) {
            int v = cand[i];
            vector<int> newCand;
            for (size_t j = i+1; j < cand.size(); j++) {
                int w = cand[j];
                if (graph[v].find(w) != graph[v].end())
                    newCand.push_back(w);
            }
            clique.push_back(v);
            extendClique(clique, newCand);
            clique.pop_back();
        }
    };
    for (int u : verticesSorted) {
        vector<int> clique = {u};
        vector<int> cand;
        for (int w : graph[u]) {
            if (w > u)
                cand.push_back(w);
        }
        extendClique(clique, cand);
    }
  
    cout << "Triangle Count: " << triangleCount << "\n";
    cout << "Maximal Cliques Count: " << maximalCliqueCount << "\n";
    cout << "Total Cliques (size >= 3) Count: " << totalCliquesCount << "\n";
    return 0;
}
