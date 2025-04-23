#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <climits>

using namespace std;

// Class to represent a graph
class Graph {
private:
    int n;                              // Number of vertices
    vector<vector<int>> adjList;        // Adjacency list
    vector<vector<int>> cliques;        // List of h-cliques
    
public:
    // Constructor
    Graph(int vertices) {
        n = vertices;
        adjList.resize(n);
    }
    
    // Add an edge to the graph
    void addEdge(int u, int v) {
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }
    
    // Find all h-cliques in the graph (for simplicity, we'll implement for h=3, i.e., triangles)
    void findCliques() {
        // For triangles (h=3)
        for (int u = 0; u < n; u++) {
            for (size_t i = 0; i < adjList[u].size(); i++) {
                int v = adjList[u][i];
                if (v > u) { // Avoid counting the same triangle multiple times
                    for (size_t j = 0; j < adjList[v].size(); j++) {
                        int w = adjList[v][j];
                        if (w > v) { // Ensure we're not double counting
                            // Check if w is also connected to u
                            if (find(adjList[u].begin(), adjList[u].end(), w) != adjList[u].end()) {
                                cliques.push_back({u, v, w});
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Compute clique degree for each vertex
    vector<int> computeCliqueDegrees() {
        vector<int> cliqueDegrees(n, 0);
        for (const auto& clique : cliques) {
            for (int v : clique) {
                cliqueDegrees[v]++;
            }
        }
        return cliqueDegrees;
    }
    
    // Perform (k, Ψ)-core decomposition
    vector<int> kCliqueCore() {
        vector<int> coreNumbers(n, 0);
        vector<int> cliqueDegrees = computeCliqueDegrees();
        
        // Maximum possible core number
        int maxDegree = *max_element(cliqueDegrees.begin(), cliqueDegrees.end());
        
        // Bin sort vertices by their degrees
        vector<vector<int>> bins(maxDegree + 1);
        for (int v = 0; v < n; v++) {
            bins[cliqueDegrees[v]].push_back(v);
        }
        
        // Process vertices in non-decreasing order of degrees
        vector<bool> processed(n, false);
        int remaining = n;
        
        for (int k = 0; k <= maxDegree && remaining > 0; k++) {
            while (!bins[k].empty()) {
                int v = bins[k].back();
                bins[k].pop_back();
                
                if (processed[v]) continue;
                processed[v] = true;
                remaining--;
                coreNumbers[v] = k;
                
                // Update neighbors
                for (int u : adjList[v]) {
                    if (!processed[u] && cliqueDegrees[u] > k) {
                        // Find common triangles and decrement
                        for (const auto& clique : cliques) {
                            if (find(clique.begin(), clique.end(), v) != clique.end() && 
                                find(clique.begin(), clique.end(), u) != clique.end()) {
                                // Remove v from a clique containing u
                                bins[cliqueDegrees[u]].erase(
                                    remove(bins[cliqueDegrees[u]].begin(), bins[cliqueDegrees[u]].end(), u),
                                    bins[cliqueDegrees[u]].end()
                                );
                                cliqueDegrees[u]--;
                                bins[cliqueDegrees[u]].push_back(u);
                            }
                        }
                    }
                }
            }
        }
        
        return coreNumbers;
    }
    
    // Extract the (k, Ψ)-core for a given k
    Graph extractKCore(int k) {
        vector<int> coreNumbers = kCliqueCore();
        Graph kCore(n);
        
        for (int u = 0; u < n; u++) {
            if (coreNumbers[u] >= k) {
                for (int v : adjList[u]) {
                    if (coreNumbers[v] >= k && v > u) { // Avoid adding edges twice
                        kCore.addEdge(u, v);
                    }
                }
            }
        }
        
        return kCore;
    }
    
    // Build a flow network for the densest subgraph computation
    void buildFlowNetwork(double alpha, vector<vector<pair<int, double>>>& flowNetwork, int s, int t) {
        int nodes = n + 2; // n vertices + source + sink
        flowNetwork.resize(nodes);
        
        // Add edges from source to all vertices
        vector<int> cliqueDegrees = computeCliqueDegrees();
        for (int v = 0; v < n; v++) {
            flowNetwork[s].push_back({v, cliqueDegrees[v]});
        }
        
        // Add edges from vertices to sink
        for (int v = 0; v < n; v++) {
            flowNetwork[v].push_back({t, alpha * n});
        }
        
        // Add edges between vertices based on cliques
        for (const auto& clique : cliques) {
            for (size_t i = 0; i < clique.size(); i++) {
                for (size_t j = i + 1; j < clique.size(); j++) {
                    int u = clique[i];
                    int v = clique[j];
                    flowNetwork[u].push_back({v, 1.0});
                    flowNetwork[v].push_back({u, 1.0});
                }
            }
        }
    }
    
    // Find minimum s-t cut using Ford-Fulkerson algorithm
    pair<double, vector<int>> findMinCut(vector<vector<pair<int, double>>>& flowNetwork, int s, int t) {
        int nodes = flowNetwork.size();
        vector<vector<double>> residual(nodes, vector<double>(nodes, 0.0));
        
        // Initialize residual graph
        for (int u = 0; u < nodes; u++) {
            for (const auto& edge : flowNetwork[u]) {
                int v = edge.first;
                double capacity = edge.second;
                residual[u][v] += capacity;
            }
        }
        
        // Ford-Fulkerson algorithm
        double maxFlow = 0.0;
        vector<int> parent(nodes);
        
        while (true) {
            fill(parent.begin(), parent.end(), -1);
            queue<int> q;
            q.push(s);
            parent[s] = -2;
            
            while (!q.empty() && parent[t] == -1) {
                int u = q.front();
                q.pop();
                
                for (int v = 0; v < nodes; v++) {
                    if (parent[v] == -1 && residual[u][v] > 0) {
                        parent[v] = u;
                        q.push(v);
                    }
                }
            }
            
            if (parent[t] == -1) break; // No augmenting path
            
            double pathFlow = INT_MAX;
            for (int v = t; v != s; v = parent[v]) {
                int u = parent[v];
                pathFlow = min(pathFlow, residual[u][v]);
            }
            
            for (int v = t; v != s; v = parent[v]) {
                int u = parent[v];
                residual[u][v] -= pathFlow;
                residual[v][u] += pathFlow;
            }
            
            maxFlow += pathFlow;
        }
        
        // Find vertices in the source side of the cut
        vector<bool> visited(nodes, false);
        queue<int> q;
        q.push(s);
        visited[s] = true;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (int v = 0; v < nodes; v++) {
                if (!visited[v] && residual[u][v] > 0) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }
        
        // Collect vertices in the source side (excluding s)
        vector<int> sourceSet;
        for (int v = 0; v < n; v++) {
            if (visited[v]) {
                sourceSet.push_back(v);
            }
        }
        
        return {maxFlow, sourceSet};
    }
    
    // CoreExact algorithm implementation
    vector<int> coreExact() {
        // Step 1: Perform (k, Ψ)-core decomposition
        vector<int> coreNumbers = kCliqueCore();
        int kMax = *max_element(coreNumbers.begin(), coreNumbers.end());
        
        // Step 2: Locate the CDS in the (k'', Ψ)-core using Pruning1 and Pruning2
        // For simplicity, we'll use k'' = kMax
        int kDoublePrime = kMax;
        
        // Step 3-4: Initialize bounds
