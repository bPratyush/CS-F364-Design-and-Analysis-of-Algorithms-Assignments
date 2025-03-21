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

long long maximalCliqueCount = 0;

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
  
    // Set up buckets to support the degeneracy ordering.
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
                buckets[d_old - 1].push_back(nb);
            }
        }
    }
  
    // Map each vertex to its position in the ordering.
    unordered_map<int, int> pos;
    for (int i = 0; i < ordering.size(); i++) {
        pos[ordering[i]] = i;
    }
  
    // Count triangles using forward neighbors (Chiba–Nishizeki idea)
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
   
        auto startTime = chrono::steady_clock::now();
        int largestCliqueSize = 0;  // track largest maximal clique size
        unordered_map<int, long long> cliqueDistribution;  // key: clique size, value: count
    
        function<void(vector<int>&, vector<int>&)> extendClique =
            [&](vector<int> &clique, vector<int> &cand) {
                if (cand.empty()){
                    // Found a maximal clique.
                    int sz = clique.size();
                    maximalCliqueCount++;
                    largestCliqueSize = max(largestCliqueSize, sz);
                    cliqueDistribution[sz]++;
                    return;
                }
                // Extend the clique with each candidate.
                for (size_t i = 0; i < cand.size(); i++) {
                    int v = cand[i];
                    vector<int> newCand;
                    // Only consider candidates after v that are also adjacent to v.
                    for (size_t j = i + 1; j < cand.size(); j++) {
                        int w = cand[j];
                        if (graph.at(v).find(w) != graph.at(v).end())
                            newCand.push_back(w);
                    }
                    clique.push_back(v);
                    extendClique(clique, newCand);
                    clique.pop_back();
                }
            };
    
        vector<int> verticesSorted;
        for (const auto &p: graph)
            verticesSorted.push_back(p.first);
        // Use degeneracy (pos) ordering for uniqueness: u is the smallest vertex
        sort(verticesSorted.begin(), verticesSorted.end(),
             [&](int a, int b) { return pos[a] < pos[b]; });
        for (int u : verticesSorted) {
            vector<int> clique = {u};
            vector<int> cand;
            // Only neighbors with a higher order than u.
            for (int w : graph[u])
                if (pos[u] < pos[w])
                    cand.push_back(w);
            extendClique(clique, cand);
        }
        auto endTime = chrono::steady_clock::now();
        auto execTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    
        cout << "Largest size of the clique: " << largestCliqueSize << "\n";
        cout << "Total number of maximal cliques: " << maximalCliqueCount << "\n";
        cout << "Execution time (ms): " << execTime << "\n";
        cout << "Distribution of different size cliques:\n";
        vector<int> sizes;
        for (const auto &p : cliqueDistribution)
            sizes.push_back(p.first);
        sort(sizes.begin(), sizes.end());
        for (int sz : sizes) {
            cout << "Size " << sz << ": " << cliqueDistribution[sz] << "\n";
        }
    }