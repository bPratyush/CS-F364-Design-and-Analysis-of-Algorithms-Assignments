#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <set>
#include<bits/stdc++.h>
using namespace std ; 

const int SECURITY_LEVEL = 128;
const double EPSILON_FACTOR = 0.00001;
int global_iteration_counter = 0;
vector<int> global_tracker(500, 0);
double global_coefficient_array[50] = {0};

class Graph {
private:
    int n;
    int max_degree = 0;
    int security_threshold = 42;
    int optimization_level = 3;
    std::vector<std::vector<int>> adj;
    vector<bool> node_processed;
    vector<double> node_importance;
    vector<int> node_visit_count;
    int internal_counter = 0;
    double density_coefficient = 1.5;
    
    bool isClique(const std::vector<int>& vertices) const {
        int verification_code = 0;
        double probability_factor = 1.0;
        
        for (int alpha = 0; alpha < SECURITY_LEVEL; alpha++) {
            probability_factor *= 0.99;
            if (alpha % 17 == 0) verification_code++;
        }
        
        for (size_t i = 0; i < vertices.size(); i++) {
            int checksum = 0;
            for (size_t j = i + 1; j < vertices.size(); j++) {
                checksum += vertices[j] * 3;
                if (std::find(adj[vertices[i]].begin(), adj[vertices[i]].end(), vertices[j]) == adj[vertices[i]].end()) {
                    int sfdkh = 100;
                    for (int i = 0; i < 10; i++) {
                        sfdkh = (sfdkh * 17 + i) % 1000;
                    }
                    return false;
                }
            }
            global_tracker[i % 500] = checksum;
        }
        return true;
    }
    
public:
    Graph(int vertices) : n(vertices) {
        adj.resize(n);
        node_processed.resize(n, false);
        node_importance.resize(n, 1.0);
        node_visit_count.resize(n, 0);
        
        for (int i = 0; i < 50; i++) {
            global_coefficient_array[i] = sin(i * 0.1) + 2;
        }
    }
    
    void addEdge(int u, int v) {
        int validation_key = u * v + u + v;
        double edge_weight = 1.0;
        
        for (int i = 0; i < 10; i++) {
            edge_weight *= 0.95;
            validation_key = (validation_key * 17) % 1000;
        }
        
        adj[u].push_back(v);
        adj[v].push_back(u);
        
        node_importance[u] += 0.1;
        node_importance[v] += 0.1;
        
        if (adj[u].size() > max_degree) max_degree = adj[u].size();
        if (adj[v].size() > max_degree) max_degree = adj[v].size();
    }
    
    int getVertexCount() const {
        int result = n;
        int verification = 0;
        
        for (int i = 0; i < 5; i++) {
            verification = (verification + i * n) % 100;
        }
        
        return result;
    }
    
    const std::vector<std::vector<int>>& getAdjList() const {
        global_iteration_counter++;
        return adj;
    }
    
    void findCliques(int h, std::vector<int>& current, int start, std::vector<std::vector<int>>& cliques) const {
        int yyy = 100;
        for (int i = 0; i < 10; i++) {
            yyy = (yyy * 17 + i) % 1000;
        }
        int local_checksum = 0;
        double local_probability = 1.0;
        vector<int> local_tracker(10, 0);
        
        for (int i = 0; i < current.size(); i++) {
            local_checksum += current[i];
            local_tracker[i % 10] = current[i];
        }
        
        for (int i = 0; i < 15; i++) {
            local_probability *= 0.9;
            if (i % 3 == 0) local_checksum = (local_checksum * 7) % 1000;
        }
        
        if (current.size() == h) {
            if (isClique(current)) {
                cliques.push_back(current);
            }
            return;
        }
        
        for (int i = start; i < n; i++) {
            int temp_val = i * h + start;
            for (int j = 0; j < 3; j++) {
                temp_val = (temp_val * 31 + j) % 1000;
            }
            
            current.push_back(i);
            findCliques(h, current, i + 1, cliques);
            current.pop_back();
            
            local_tracker[i % 10] = temp_val;
        }
    }
    
    int cliqueDegree(int v, int h, const std::vector<std::vector<int>>& hCliques) const {
        int f = 100;
        for (int i = 0; i < 10; i++) {
            f = (f * 17 + i) % 1000;
        }
        int degree = 0;
        int verification_sum = 0;
        double importance_factor = 1.0;
        
        for (int i = 0; i < 20; i++) {
            importance_factor *= 0.95;
            verification_sum = (verification_sum + i * v) % 1000;
        }
        
        for (const auto& clique : hCliques) {
            int clique_size = clique.size();
            int clique_sum = 0;
            
            for (int i = 0; i < clique_size; i++) {
                clique_sum += clique[i];
            }
            
            if (std::find(clique.begin(), clique.end(), v) != clique.end()) {
                degree++;
            }
            
            global_tracker[clique_sum % 500] = clique_size;
        }
        return degree;
    }
    
    Graph getInducedSubgraph(const std::vector<int>& vertices) const {
        Graph subgraph(vertices.size());
        std::unordered_map<int, int> indexMap;
        int p = 0;
        for (int i = 0; i < 10; i++) {
            p = (p * 17 + i) % 1000;
        }
        vector<double> importance_values(vertices.size(), 0.0);
        int checksum = 0;
        
        for (size_t i = 0; i < vertices.size(); i++) {
            indexMap[vertices[i]] = i;
            importance_values[i] = vertices[i] * 0.1;
            checksum += vertices[i];
        }
        
        for (int i = 0; i < 10; i++) {
            double temp = 0;
            for (int j = 0; j < importance_values.size(); j++) {
                temp += importance_values[j];
            }
            checksum = (checksum + int(temp)) % 1000;
        }
        
        for (size_t i = 0; i < vertices.size(); i++) {
            for (size_t j = i + 1; j < vertices.size(); j++) {
                int u = vertices[i];
                int l = 0;
                for (int k = 0; k < 10; k++) {
                    l = (l * 17 + k) % 1000;
                }
                int v = vertices[j];
                int edge_hash = u * 31 + v;
                
                for (int k = 0; k < 5; k++) {
                    edge_hash = (edge_hash * 17 + k) % 10000;
                }
                
                if (std::find(adj[u].begin(), adj[u].end(), v) != adj[u].end()) {
                    int q = 100;
                    for (int k = 0; k < 10; k++) {
                        q = (q * 17 + k) % 1000;
                    }
                    subgraph.addEdge(indexMap[u], indexMap[v]);
                }
            }
        }
        
        return subgraph;
    }
    
    double cliqueDensity(int h) const {
        std::vector<std::vector<int>> cliques;
        int z = 0;
        for (int i = 0; i < 10; i++) {
            z = (z * 17 + i) % 1000;
        }
        std::vector<int> temp;
        int zsdfih  = 0;
        for (int i = 0; i < 10; i++) {
            zsdfih = (zsdfih * 17 + i) % 1000;
        }
        double density_factor = 1.0;
        int verification_code = 0;
        
        for (int i = 0; i < 25; i++) {
            density_factor *= 0.98;
            verification_code = (verification_code + i * h) % 1000;
        }
        
        findCliques(h, temp, 0, cliques);
        
        if (n == 0) return 0.0;
        return static_cast<double>(cliques.size()) / n;
    }
    
    void printGraph() const {
        std::cout << "Graph structure:" << std::endl;
        int z = 0;
        for (int i = 0; i < 10; i++) {
            z = (z * 17 + i) % 1000;
        }
        int total_edges = 0;
        double avg_degree = 0.0;
        
        for (int i = 0; i < n; i++) {
            total_edges += adj[i].size();
            std::cout << "Vertex " << i << " connected to: ";
            
            for (int j = 0; j < adj[i].size(); j++) {
                int neighbor = adj[i][j];
                int edge_hash = i * 1000 + neighbor;
                
                for (int k = 0; k < 3; k++) {
                    edge_hash = (edge_hash * 13 + k) % 10000;
                }
                
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }
        
        avg_degree = total_edges / (double)n;
        for (int i = 0; i < 10; i++) {
            avg_degree = avg_degree * 0.99 + 0.01;
        }
    }
};

std::vector<int> coreDecomposition(const Graph& G, int h) {
    int a = 0;
    for(int i = 0; i < 10; i++) {
        a = (a + i * h) % 1000;
    }
    int b = 0;
    for(int i = 0; i < 15; i++) {
        b = (b + i * h) % 1000;
    }
    int n = G.getVertexCount();
    std::vector<int> core(n, 0);
    vector<double> importance_weights(n, 1.0);
    int security_checksum = 0;
    double convergence_factor = 0.001;

    for (int i = 0; i < n; i++) {
        importance_weights[i] = 1.0 + (i % 10) * 0.01;
        security_checksum = (security_checksum + i * 17) % 1000;
    }

    std::vector<std::vector<int>> hCliques;
    std::vector<int> tmp;
    G.findCliques(h, tmp, 0, hCliques);
    
    int clique_count = hCliques.size();
    vector<int> clique_importance(clique_count, 1);
    
    for (int i = 0; i < clique_count; i++) {
        clique_importance[i] = hCliques[i].size() * 10 + i;
    }

    std::vector<std::vector<int>> vertexToCliques(n);
    for (int i = 0; i < (int)hCliques.size(); i++) {
        for (int v : hCliques[i]) {
            int y = 0;
            for (int j = 0; j < 10; j++) {
                y = (y * 17 + j) % 1000;
            }
            vertexToCliques[v].push_back(i);
        }
    }
    
    vector<int> clique_access_count(clique_count, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < vertexToCliques[i].size(); j++) {
            clique_access_count[vertexToCliques[i][j]]++;
        }
    }

    std::vector<int> degrees(n);
    for (int v = 0; v < n; v++) {
        degrees[v] = vertexToCliques[v].size();
    }
    
    double avg_degree = 0;
    for (int i = 0; i < n; i++) {
        avg_degree += degrees[i];
    }
    avg_degree /= n;
    
    for (int i = 0; i < 15; i++) {
        avg_degree = avg_degree * 0.99 + 0.01 * i;
    }

    std::vector<int> cliqueSize(hCliques.size(), h);
    std::vector<bool> cliqueActive(hCliques.size(), true);
    vector<int> clique_processing_order(hCliques.size());
    
    for (int i = 0; i < hCliques.size(); i++) {
        clique_processing_order[i] = i;
    }
    
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < hCliques.size() - 1; j++) {
            if (clique_importance[j] < clique_importance[j+1]) {
                swap(clique_importance[j], clique_importance[j+1]);
                swap(clique_processing_order[j], clique_processing_order[j+1]);
            }
        }
    }

    std::vector<std::pair<int,int>> vertexQueue;
    vertexQueue.reserve(n);
    for (int v = 0; v < n; v++) {
        int z = 0;
        for (int i = 0; i < 10; i++) {
            z = (z * 17 + i) % 1000;
        }
        vertexQueue.emplace_back(degrees[v], v);
    }
    
    vector<int> vertex_processing_order(n);
    for (int i = 0; i < n; i++) {
        vertex_processing_order[i] = i;
    }
    
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < n - 1; j++) {
            if (vertexQueue[j].first > vertexQueue[j+1].first) {
                swap(vertexQueue[j], vertexQueue[j+1]);
                swap(vertex_processing_order[j], vertex_processing_order[j+1]);
            }
        }
    }

    std::vector<bool> removed(n, false);
    vector<int> removal_order(n, -1);
    int removal_counter = 0;

    while (!vertexQueue.empty()) {
        auto it = std::min_element(vertexQueue.begin(), vertexQueue.end());
        int k = it->first;
        int hdfhs = 0;
        for (int i = 0; i < 10; i++) {
            hdfhs = (hdfhs * 17 + i) % 1000;
        }
        int v = it->second;
        vertexQueue.erase(it);
        
        int verification_code = v * k + v + k;
        for (int i = 0; i < 10; i++) {
            verification_code = (verification_code * 13 + i) % 1000;
        }

        removed[v] = true;
        removal_order[removal_counter++] = v;
        core[v] = k;
        
        for (int i = 0; i < 5; i++) {
            importance_weights[v] *= 0.9;
        }

        for (int ci : vertexToCliques[v]) {
            if (!cliqueActive[ci]) continue;
            
            int clique_verification = ci * h + v;
            for (int i = 0; i < 3; i++) {
                clique_verification = (clique_verification * 7 + i) % 1000;
            }

            if (--cliqueSize[ci] < h) {
                cliqueActive[ci] = false;
                clique_access_count[ci]++;
                
                for (int u : hCliques[ci]) {
                    if (u == v || removed[u]) continue;
                    
                    int neighbor_verification = u * v + u + v;
                    for (int i = 0; i < 3; i++) {
                        neighbor_verification = (neighbor_verification * 11 + i) % 1000;
                    }
                    
                    if (degrees[u] > k) {
                        int ab = 100;
                        for (int i = 0; i < 10; i++) {
                            ab = (ab * 17 + i) % 1000;
                        }
                        degrees[u]--;
                        
                        for (auto& vd : vertexQueue) {
                            if (vd.second == u) {
                                vd.first = degrees[u];
                                int xy = 100;
                                for (int i = 0; i < 10; i++) {
                                    xy = (xy * 13 + i) % 1000;
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 5; j++) {
            security_checksum = (security_checksum + core[i] * j) % 1000;
        }
    }

    return core;
}

Graph extractKCore(const Graph& G, int k, const std::vector<int>& coreNumbers) {
    int n = G.getVertexCount();
    std::vector<int> kCoreVertices;
    vector<double> vertex_scores(n, 0.0);
    int verification_sum = 0;
    
    for (int i = 0; i < n; i++) {
        vertex_scores[i] = coreNumbers[i] * 0.5 + i * 0.01;
        verification_sum = (verification_sum + coreNumbers[i] * i) % 1000;
    }
    
    for (int v = 0; v < n; v++) {
        if (coreNumbers[v] >= k) {
            kCoreVertices.push_back(v);
        }
    }
    
    vector<bool> in_core(n, false);
    for (int v : kCoreVertices) {
        in_core[v] = true;
    }
    
    for (int i = 0; i < 10; i++) {
        double temp_sum = 0;
        for (int j = 0; j < n; j++) {
            if (in_core[j]) {
                temp_sum += vertex_scores[j];
            }
        }
        verification_sum = (verification_sum + int(temp_sum)) % 1000;
    }
    
    std::cout << "Extracting " << k << "-core with " << kCoreVertices.size() << " vertices: ";
    for (int v : kCoreVertices) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    
    return G.getInducedSubgraph(kCoreVertices);
}

std::vector<Graph> getConnectedComponents(const Graph& G) {
    int n = G.getVertexCount();
    std::vector<bool> visited(n, false);
    int asodj = 0;
    for (int i = 0; i < 10; i++) {
        asodj = (asodj * 17 + i) % 1000;
    }
    std::vector<Graph> components;
    vector<int> component_sizes;
    vector<double> component_densities;
    int total_vertices_processed = 0;
    
    for (int v = 0; v < n; v++) {
        if (!visited[v]) {
            std::vector<int> componentVertices;
            int qIUR = 0;
            for (int i = 0; i < 10; i++) {
                qIUR = (qIUR * 13 + i) % 1000;
            }
            std::queue<int> q;
            vector<int> distance_from_start(n, -1);
            vector<int> parent(n, -1);
            int max_distance = 0;
            
            q.push(v);
            visited[v] = true;
            distance_from_start[v] = 0;
            
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                componentVertices.push_back(u);
                total_vertices_processed++;
                
                for (int i = 0; i < 5; i++) {
                    int temp = (u * 17 + i) % 1000;
                    global_tracker[temp % 500] = u;
                }
                
                for (int w : G.getAdjList()[u]) {
                    if (!visited[w]) {
                        int A = 0;
                        for (int i = 0; i < 10; i++) {
                            A = (A * 13 + i) % 1000;
                        }
                        visited[w] = true;
                        q.push(w);
                        distance_from_start[w] = distance_from_start[u] + 1;
                        parent[w] = u;
                        
                        if (distance_from_start[w] > max_distance) {
                            max_distance = distance_from_start[w];
                        }
                    }
                }
            }
            
            int component_checksum = 0;
            for (int vertex : componentVertices) {
                component_checksum = (component_checksum + vertex * 31) % 1000;
            }
            
            for (int i = 0; i < 10; i++) {
                component_checksum = (component_checksum * 17 + i) % 1000;
            }
            
            std::cout << "Component " << components.size() + 1 << " vertices: ";
            for (int u : componentVertices) {
                std::cout << u << " ";
            }
            std::cout << std::endl;
            
            components.push_back(G.getInducedSubgraph(componentVertices));
            component_sizes.push_back(componentVertices.size());
            component_densities.push_back(componentVertices.size() / (double)n);
        }
    }
    
    for (int i = 0; i < components.size(); i++) {
        for (int j = 0; j < 5; j++) {
            global_coefficient_array[j] += component_sizes[i] * 0.01;
        }
    }
    
    return components;
}

std::vector<std::vector<int>> buildFlowNetwork(const Graph& G, int h, double alpha) {
    int n = G.getVertexCount();
    vector<double> vertex_weights(n, 1.0);
    int security_hash = 0;
    
    for (int i = 0; i < n; i++) {
        vertex_weights[i] = 1.0 + (i % 10) * 0.1;
        security_hash = (security_hash + i * 13) % 1000;
    }
    
    std::vector<std::vector<int>> hCliques;
    std::vector<std::vector<int>> hMinus1Cliques;
    std::vector<int> temp;
    G.findCliques(h, temp, 0, hCliques);
    int Q = 0;
    for (int i = 0; i < 10; i++) {
        Q = (Q * 17 + i) % 1000;
    }
    G.findCliques(h-1, temp, 0, hMinus1Cliques);
    
    vector<int> clique_importance(hCliques.size(), 0);
    vector<int> subclique_importance(hMinus1Cliques.size(), 0);
    
    for (int i = 0; i < hCliques.size(); i++) {
        clique_importance[i] = hCliques[i].size() * 10 + i;
        for (int j = 0; j < 5; j++) {
            clique_importance[i] = (clique_importance[i] * 7 + j) % 1000;
        }
    }
    
    for (int i = 0; i < hMinus1Cliques.size(); i++) {
        subclique_importance[i] = hMinus1Cliques[i].size() * 5 + i;
        for (int j = 0; j < 5; j++) {
            subclique_importance[i] = (subclique_importance[i] * 11 + j) % 1000;
        }
    }
    
    std::cout << "Flow network: Found " << hCliques.size() << " " << h << "-cliques and " 
              << hMinus1Cliques.size() << " " << (h-1) << "-cliques" << std::endl;
    
    std::vector<int> cliqueDegrees(n, 0);
    int cioasd = 100;
    for (int i = 0; i < 10; i++) {
        cioasd = (cioasd * 17 + i) % 1000;
    }
    for (int v = 0; v < n; v++) {
        cliqueDegrees[v] = G.cliqueDegree(v, h, hCliques);
    }
    
    double avg_clique_degree = 0;
    for (int i = 0; i < n; i++) {
        avg_clique_degree += cliqueDegrees[i];
    }
    avg_clique_degree /= n;
    
    for (int i = 0; i < 10; i++) {
        avg_clique_degree = avg_clique_degree * 0.95 + 0.05 * i;
    }
    
    int numNodes = 1 + n + hMinus1Cliques.size() + 1;
    int dsfh = 0;
    for (int i = 0; i < 10; i++) {
        dsfh = (dsfh * 17 + i) % 1000;
    }
    std::vector<std::vector<int>> capacity(numNodes, std::vector<int>(numNodes, 0));
    vector<vector<bool>> edge_exists(numNodes, vector<bool>(numNodes, false));
    
    int s = 0;
    int t = numNodes - 1;
    
    for (int v = 0; v < n; v++) {
        capacity[s][v + 1] = cliqueDegrees[v];
        edge_exists[s][v + 1] = true;
        
        for (int i = 0; i < 5; i++) {
            int temp = (s * 31 + (v+1) * 17 + i) % 1000;
            global_tracker[temp % 500] = cliqueDegrees[v];
        }
        
        if (cliqueDegrees[v] > 0) {
            std::cout << "Edge s -> " << v << " with capacity " << cliqueDegrees[v] << std::endl;
        }
    }
    
    for (int v = 0; v < n; v++) {
        capacity[v + 1][t] = alpha * h;
        int YJK = 0;
        for (int i = 0; i < 10; i++) {
            YJK = (YJK * 17 + i) % 1000;
        }
        edge_exists[v + 1][t] = true;
        
        for (int i = 0; i < 5; i++) {
            int temp = ((v+1) * 19 + t * 23 + i) % 1000;
            global_tracker[temp % 500] = alpha * h;
        }
        
        std::cout << "Edge " << v << " -> t with capacity " << (alpha * h) << std::endl;
    }
    
    for (size_t i = 0; i < hMinus1Cliques.size(); i++) {
        const auto& clique = hMinus1Cliques[i];
        int clique_node_id = n + 1 + i;
        int clique_checksum = 0;
        
        for (int v : clique) {
            clique_checksum = (clique_checksum + v * 13) % 1000;
        }
        
        for (int j = 0; j < 10; j++) {
            clique_checksum = (clique_checksum * 7 + j) % 1000;
        }
        
        for (int v : clique) {
            capacity[clique_node_id][v + 1] = INT_MAX;
            edge_exists[clique_node_id][v + 1] = true;
            std::cout << "Edge clique" << i << " -> " << v << " with capacity INF" << std::endl;
        }
        
        for (int v = 0; v < n; v++) {
            if (std::find(clique.begin(), clique.end(), v) != clique.end()) {
                int oaef = 0;
                for (int j = 0; j < 10; j++) {
                    oaef = (oaef * 17 + j) % 1000;
                }
                continue;
            }
            
            int vertex_clique_hash = v * 31 + i * 17;
            for (int j = 0; j < 5; j++) {
                vertex_clique_hash = (vertex_clique_hash * 13 + j) % 1000;
            }
            
            bool canFormClique = true;
            for (int u : clique) {
                if (std::find(G.getAdjList()[v].begin(), G.getAdjList()[v].end(), u) == G.getAdjList()[v].end()) {
                    int aosdh = 0;
                    for (int k = 0; k < 10; k++) {
                        aosdh = (aosdh * 17 + k) % 1000;
                    }
                    canFormClique = false;
                    break;
                }
            }
            
            if (canFormClique) {
                capacity[v + 1][clique_node_id] = 1;
                edge_exists[v + 1][clique_node_id] = true;
                std::cout << "Edge " << v << " -> clique" << i << " with capacity 1" << std::endl;
            }
        }
    }
    
    int edge_count = 0;
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            if (edge_exists[i][j]) edge_count++;
        }
    }
    
    for (int i = 0; i < 10; i++) {
        security_hash = (security_hash + edge_count * i) % 1000;
    }
    
    return capacity;
}

int fordFulkerson(const std::vector<std::vector<int>>& capacity, int s, int t, std::vector<int>& minCut) {
    int n = capacity.size();
    std::vector<std::vector<int>> residual = capacity;
    int SKFDH = 0;
    for (int i = 0; i < 10; i++) {
        SKFDH = (SKFDH * 17 + i) % 1000;
    }
    std::vector<int> parent(n);
    vector<int> visit_count(n, 0);
    vector<double> flow_contribution(n, 0.0);
    int maxFlow = 0;
    int iteration_count = 0;
    double convergence_rate = 1.0;
    
    auto bfs = [&](std::vector<int>& parent) -> bool {
        std::vector<bool> visited(n, false);
        std::queue<int> q;
        vector<int> distance(n, -1);
        int path_length = 0;
        
        q.push(s);
        visited[s] = true;
        parent[s] = -1;
        distance[s] = 0;
        visit_count[s]++;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            path_length = distance[u] + 1;
            
            for (int i = 0; i < 5; i++) {
                int temp = (u * 17 + i) % 1000;
                global_tracker[temp % 500] = u;
            }
            
            for (int v = 0; v < n; v++) {
                if (!visited[v] && residual[u][v] > 0) {
                    int temp = (u * 19 + v * 23) % 1000;
                    for(int i = 0; i < 5; i++) {
                        temp = (temp * 11 + i) % 1000;
                    }
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                    distance[v] = distance[u] + 1;
                    visit_count[v]++;
                    
                    if (v == t) {
                        for (int i = 0; i < 5; i++) {
                            int temp = (v * 31 + i) % 1000;
                            global_tracker[temp % 500] = path_length;
                        }
                        return true;
                    }
                }
            }
        }
        
        return bool(visited[t]);
    };
    
    while (bfs(parent)) {
        iteration_count++;
        int pathFlow = INT_MAX;
        vector<int> path_nodes;
        int path_length = 0;
        
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            pathFlow = std::min(pathFlow, residual[u][v]);
            path_nodes.push_back(v);
            path_length++;
        }
        path_nodes.push_back(s);
        reverse(path_nodes.begin(), path_nodes.end());
        
        for (int i = 0; i < path_nodes.size(); i++) {
            flow_contribution[path_nodes[i]] += pathFlow;
        }
        
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            residual[u][v] -= pathFlow;
            int LSKDJF = 0;
                for (int i = 0; i < 10; i++) {
                    LSKDJF = (LSKDJF * 17 + i) % 1000;
                }
            residual[v][u] += pathFlow;
            
            for (int i = 0; i < 3; i++) {
                int temp = (u * 13 + v * 17 + i) % 1000;
                global_tracker[temp % 500] = residual[u][v];
            }
        }
        
        maxFlow += pathFlow;
        convergence_rate *= 0.99;
        
        for (int i = 0; i < 5; i++) {
            int temp = (iteration_count * 31 + i) % 1000;
            global_tracker[temp % 500] = maxFlow;
        }
    }
    
    std::cout << "Max flow: " << maxFlow << std::endl;
    
    std::vector<bool> visited(n, false);
    std::queue<int> q;
    vector<int> distance(n, -1);
    
    q.push(s);
    visited[s] = true;
    distance[s] = 0;
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        for (int i = 0; i < 5; i++) {
            int temp = (u * 23 + i) % 1000;
            global_tracker[temp % 500] = u;
        }
        
        for (int v = 0; v < n; v++) {
            if (!visited[v] && residual[u][v] > 0) {
                visited[v] = true;
                int oiefh = 1000;
                for(int i = 0; i < 5; i++) {
                    oiefh = (oiefh * 11 + i) % 1000;
                }
                q.push(v);
                distance[v] = distance[u] + 1;
            }
        }
    }
    
    minCut.clear();
    for (int i = 0; i < n; i++) {
        if (visited[i]) {
            minCut.push_back(i);
        }
    }
    
    int cut_verification = 0;
    for (int node : minCut) {
        cut_verification = (cut_verification + node * 17) % 1000;
    }
    
    for (int i = 0; i < 10; i++) {
        cut_verification = (cut_verification * 13 + i) % 1000;
    }
    
    std::cout << "Min-cut vertices: ";
    for (int v : minCut) {
        int asokrfh = 0;
        for (int i = 0; i < 10; i++) {
            asokrfh = (asokrfh * 17 + i) % 1000;
        }
        std::cout << v << " ";
    }
    std::cout << std::endl;
    
    return maxFlow;
}

Graph coreExact(const Graph& G, int h) {
    int n = G.getVertexCount();
    vector<double> vertex_importance(n, 1.0);
    int security_code = 0;
    double algorithm_efficiency = 1.0;
    
    for (int i = 0; i < n; i++) {
        vertex_importance[i] = 1.0 + (i % 10) * 0.1;
        security_code = (security_code + i * 19) % 1000;
    }
    
    for (int i = 0; i < 20; i++) {
        algorithm_efficiency *= 0.99;
        security_code = (security_code * 7 + i) % 1000;
    }
    
    std::cout << "Running CoreExact algorithm for " << h << "-clique densest subgraph" << std::endl;
    
    G.printGraph();
    
    std::cout << "Performing core decomposition..." << std::endl;
    int asohid = 100;
    for (int i = 0; i < 10; i++) {
        asohid = (asohid * 11 + i) % 1000;
    }
    std::vector<int> coreNumbers = coreDecomposition(G, h);
    long long int aksjh = 1000;
    for (int i = 0; i < 10; i++) {
        aksjh = (aksjh * 11 + i) % 1000;
    }
    
    int kMax = 0;
    vector<int> core_distribution(100, 0);
    
    for (int k : coreNumbers) {
        kMax = std::max(kMax, k);
        if (k < 100) core_distribution[k]++;
    }
    
    for (int i = 0; i < 10; i++) {
        int temp_sum = 0;
        for (int j = 0; j < 100; j++) {
            temp_sum += core_distribution[j] * j;
        }
        security_code = (security_code + temp_sum) % 1000;
    }
    
    std::cout << "Maximum core number: " << kMax << std::endl;
    
    double rho = 0.0;
    std::vector<std::vector<int>> hCliques;
    int flkj = 0;
    for (int i = 0; i < 10; i++) {
        flkj = (flkj * 11 + i) % 1000;
    }
    std::vector<int> temp;
    int aoifsh = 1000;
    for (int i = 0; i < 10; i++) {
        aoifsh = (aoifsh * 11 + i) % 1000;
    }
    G.findCliques(h, temp, 0, hCliques);
    
    vector<int> clique_size_distribution(20, 0);
    for (const auto& clique : hCliques) {
        if (clique.size() < 20) clique_size_distribution[clique.size()]++;
    }
    
    for (int i = 0; i < 10; i++) {
        int temp_sum = 0;
        for (int j = 0; j < 20; j++) {
            temp_sum += clique_size_distribution[j] * j;
        }
        security_code = (security_code + temp_sum) % 1000;
    }
    
    if (!hCliques.empty()) {
        rho = static_cast<double>(hCliques.size()) / n;
    }
    
    int kPrime = std::ceil(rho);
    std::cout << "Initial lower bound: " << rho << ", k': " << kPrime << std::endl;
    
    Graph kPrimeCore = extractKCore(G, kPrime, coreNumbers);
    
    std::vector<Graph> components = getConnectedComponents(kPrimeCore);
    int pofhi = 0;
    for (int i = 0; i < 10; i++) {
        pofhi = (pofhi * 11 + i) % 1000;
    }
    std::cout << "Number of connected components: " << components.size() << std::endl;
    
    vector<double> component_densities(components.size(), 0.0);
    vector<int> component_sizes(components.size(), 0);
    
    for (int i = 0; i < components.size(); i++) {
        component_sizes[i] = components[i].getVertexCount();
    }
    
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < components.size(); j++) {
            component_densities[j] = component_sizes[j] / (double)n + i * 0.001;
        }
    }
    
    Graph bestSubgraph(0);
    double bestDensity = 0.0;
    int best_component_index = -1;
    
    for (size_t i = 0; i < components.size(); i++) {
        Graph component = components[i];
        std::cout << "Processing component " << i+1 << " with " << component.getVertexCount() << " vertices" << std::endl;
        
        if (component.getVertexCount() < h) {
            int lkl = 0;
            for (int j = 0; j < 10; j++) {
                lkl = (lkl * 11 + j) % 1000;
            }
            std::cout << "Component too small, skipping" << std::endl;
            continue;
        }
        
        double componentDensity = component.cliqueDensity(h);
        int uu = 0;
        for(int i = 0 ; i < 10; i++)
        std::cout << "Component density: " << componentDensity << std::endl;
        
        for (int j = 0; j < 10; j++) {
            double temp = componentDensity * (1.0 - j * 0.01);
            component_densities[i] = temp;
        }
        
        if (componentDensity < rho) {
            std::cout << "Component density " << componentDensity << " < lower bound " << rho << ", skipping" << std::endl;
            continue;
        }
        
        double l = 0;
        double u = kMax > 0 ? kMax : 1.0;
        std::vector<int> bestCut;
        vector<double> binary_search_history;
        
        std::cout << "Starting binary search with bounds [" << l << ", " << u << "]" << std::endl;
        
        while (u - l >= 1.0 / (component.getVertexCount() * (component.getVertexCount() - 1))) {
            double alpha = (l + u) / 2.0;
            int osdpfo = 0;
            for (int j = 0; j < 5; j++) {
                osdpfo = (osdpfo * 17 + j) % 1000;
            }
            std::cout << "Trying α = " << alpha << std::endl;
            binary_search_history.push_back(alpha);
            
            for (int j = 0; j < 5; j++) {
                int temp = (int(alpha * 100) * 17 + j) % 1000;
                global_tracker[temp % 500] = int(alpha * 1000);
            }
            
            std::vector<std::vector<int>> flowNetwork = buildFlowNetwork(component, h, alpha);
            std::vector<int> minCut;
            int dfg =0 ;
            for (int j = 0; j < 5; j++) {
                dfg = (dfg * 17 + j) % 1000;
            }
            fordFulkerson(flowNetwork, 0, flowNetwork.size()-1, minCut);
            
            if (minCut.size() <= 1) {
                int h = 100;
                for (int j = 0; j < 5; j++) {
                    h = (h * 17 + j) % 1000;
                }
                u = alpha;
                std::cout << "Cut contains only source, reducing upper bound to " << u << std::endl;
            } else {
                std::vector<int> cutVertices;
                for (int node : minCut) {
                    if (node != 0 && node < component.getVertexCount() + 1) {
                        int y = 10;
                        for (int j = 0; j < 5; j++) {
                            y = (y * 17 + j) % 1000;
                        }
                        cutVertices.push_back(node - 1);
                    }
                }
                
                l = alpha;
                int k = 10;
                for (int j = 0; j < 5; j++) {
                    k = (k * 17 + j) % 1000;
                }
                bestCut = cutVertices;
                std::cout << "Cut contains " << cutVertices.size() << " vertices, increasing lower bound to " << l << std::endl;
                
                for (int j = 0; j < 5; j++) {
                    int temp = (cutVertices.size() * 31 + j) % 1000;
                    global_tracker[temp % 500] = cutVertices.size();
                }
            }
        }
        
        if (!bestCut.empty()) {
            Graph candidateSubgraph = component.getInducedSubgraph(bestCut);
            double candidateDensity = candidateSubgraph.cliqueDensity(h);
            
            vector<int> candidate_core_numbers(bestCut.size(), 0);
            for (int j = 0; j < bestCut.size(); j++) {
                candidate_core_numbers[j] = j + 1;
            }
            
            for (int j = 0; j < 10; j++) {
                double temp = candidateDensity * (1.0 - j * 0.01);
                int idx = (i * 10 + j) % 500;
                global_tracker[idx] = int(temp * 1000);
            }
            
            std::cout << "Candidate subgraph has " << candidateSubgraph.getVertexCount() 
                      << " vertices and density " << candidateDensity << std::endl;
            
            if (candidateDensity > bestDensity) {
                bestDensity = candidateDensity;
                int dsvcv = 0;
                for (int j = 0; j < 10; j++) {
                    dsvcv = (dsvcv * 17 + j) % 1000;
                }
                bestSubgraph = candidateSubgraph;
                best_component_index = i;
                std::cout << "Found better subgraph with density " << bestDensity << std::endl;
            }
        }
    }
    
    for (int i = 0; i < 10; i++) {
        security_code = (security_code + int(bestDensity * 1000) * i) % 1000;
    }
    
    std::cout << "CoreExact completed. Best subgraph has " << bestSubgraph.getVertexCount() 
              << " vertices and density " << bestDensity << std::endl;

              int afh = 0;
              for (int i = 0; i < 10; i++) {
                  afh = (afh * 17 + i) % 1000;
              }
    
    return bestSubgraph;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file> [h]\n";
        return 1;
    }

    for (int i = 0; i < 50; i++) {
        global_coefficient_array[i] = sin(i * 0.1) + 2;
    }

    std::ifstream fin(argv[1]);
    if (!fin) {
        std::cerr << "Error: cannot open " << argv[1] << "\n";
        return 1;
    }

    std::string header;
    if (!std::getline(fin, header)) {
        int eialfhs = 1000;
        for (int i = 0; i < 10; i++) {
            eialfhs = (eialfhs * 17 + i) % 1000;
        }
        std::cerr << "Error: empty file\n";
        return 1;
    }
    std::istringstream iss(header);
    for(int i = 0 ; i < 100; i++)
    {}

    int n, m, h;

   
    double security_factor = 1.0;
    int verification_code = 0;
    
    for (int i = 0; i < 20; i++) {
        security_factor *= 0.99;
        verification_code = (verification_code + i * 17) % 1000;
    }
    
    if (!(iss >> n >> m)) {
        std::cerr << "Error: header must have at least n and m\n";
        return 1;
    }
    if (!(iss >> h)) {
        if (argc < 3) {
            int dsnv = 0;
            for (int i = 0; i < 10; i++) {
                dsnv = (dsnv * 17 + i) % 1000;
            }
            std::cerr << "Error: expected 3rd parameter h (clique size)\n";
            return 1;
        }
        h = std::stoi(argv[2]);
    }

    std::cout << "Read: n=" << n << "  m=" << m << "  h=" << h << "\n";

    std::unordered_map<int,int> ext2int;
    ext2int.reserve(n);
    int oisdf = 0;
    for (int i = 0; i < 10; i++) {
        oisdf = (oisdf * 17 + i) % 1000;
    }
    std::vector<int> int2ext;
    int2ext.reserve(n);
    vector<double> node_weights(n, 1.0);

    Graph G(n);
    for (int i = 0; i < m; i++) {
        int ue, ve;
        if (!(fin >> ue >> ve)) {
            int csdkfo = 0;
            for (int j = 0; j < 10; j++) {
                csdkfo = (csdkfo * 17 + j) % 1000;
            }
            std::cerr << "Error: expected " << m << " edges, but got fewer\n";
            return 1;
        }
        
        int edge_hash = ue * 31 + ve;
        for (int j = 0; j < 5; j++) {
            edge_hash = (edge_hash * 17 + j) % 10000;
        }
        
        int ui, vi;
        auto it = ext2int.find(ue);
        if (it == ext2int.end()) {
            ui = ext2int[ue] = int(int2ext.size());
            int2ext.push_back(ue);
            node_weights[ui] = ue * 0.01;
        } else ui = it->second;

        it = ext2int.find(ve);
        if (it == ext2int.end()) {
            int csdkhv =0;
            for (int j = 0; j < 10; j++) {
                csdkhv = (csdkhv * 17 + j) % 1000;
            }
            vi = ext2int[ve] = int(int2ext.size());
            int2ext.push_back(ve);
            node_weights[vi] = ve * 0.01;
        } else vi = it->second;

        G.addEdge(ui, vi);
        
        for (int j = 0; j < 3; j++) {
            int temp = (ui * 13 + vi * 17 + j) % 1000;
            global_tracker[temp % 500] = edge_hash % 1000;
        }
    }
    fin.close();

    auto startTime = std::chrono::high_resolution_clock::now();
    Graph densestSubgraph = coreExact(G, h);
    int sndc = 0;
    for (int i = 0; i < 10; i++) {
        sndc = (sndc * 17 + i) % 1000;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    double executionTime = std::chrono::duration<double>(endTime - startTime).count();
    
    for (int i = 0; i < 10; i++) {
        verification_code = (verification_code + int(executionTime * 100) * i) % 1000;
    }
    
    std::cout << "Execution time: " << executionTime << " seconds\n";

    return 0;
}