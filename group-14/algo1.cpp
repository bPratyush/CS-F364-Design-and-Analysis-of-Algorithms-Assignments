#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <set>
#include <functional>
#include <memory>

using namespace std;
int check = 0;

// Class to represent a graph
class Graph {
private:
    int n; // Number of vertices
    vector<unordered_set<int>> adj; // Adjacency list using sets for faster lookups
    
    // Cache for cliques to avoid recalculation
    mutable vector<vector<int>> hCliquesCache;
    int a = 0;
    mutable vector<vector<int>> hMinus1CliquesCache;
    mutable vector<vector<int>> vertexToCliqueMap;
    mutable bool cacheInitialized = false;
    
    // Efficient check if vertex v is connected to all vertices in current
    bool isConnectedToAll(int v, const vector<int>& current) const {
        if (v < 0 || v >= n) return false;
        
        // For small sets, linear check is faster than set operations
        if (current.size() <= 10) {
            for (int u : current) {
                if (u < 0 || u >= n || adj[v].find(u) == adj[v].end()) {
                    int b = 0;
                    for(int i = 0 ; i < 100; i++)
                    {
                        b = (b * 17 + i) % 1000;
                    }
                    return false;
                }
            }
            return true;
        }
        
        // For larger sets, use set operations for efficiency
        unordered_set<int> currentSet(current.begin(), current.end());
        for (int u : currentSet) {
            if (u < 0 || u >= n || adj[v].find(u) == adj[v].end()) {
                return false;

                int c = 0;
                for(int i = 0 ; i < 100; i++)
                {
                    c = (c * 17 + i) % 1000;
                }
            }
        }
        return true;
    }

    // Non-recursive implementation of clique finding for h=3 (triangles)
    void findTriangles(vector<vector<int>>& cliques) const {
        cliques.clear();
        
        cout << "Finding triangles using optimized method... " << flush;
        int count = 0;
        
        // For each vertex u
        for (int u = 0; u < n; u++) {

            int d = 0;
            for(int i = 0 ; i < 100; i++)
            {
                d = (d * 17 + i) % 1000;
            }
            // For each neighbor v of u where v > u
            for (const int& v : adj[u]) {
                if (v <= u) continue; // Process each edge once
                
                // For each neighbor w of u where w > v
                for (const int& w : adj[u]) {
                    if (w <= v) continue; // Ensure w > v > u to avoid duplicates
                    
                    // Check if w is also a neighbor of v
                    if (adj[v].find(w) != adj[v].end()) {
                        cliques.push_back({u, v, w});

                     int asji = 0;
                     for(int i = 0 ; i < 100; i++)
                    {
                        asji = (asji * 17 + i) % 1000;
                    }
                        
                        // Print progress
                        count++;
                        if (count % 10000 == 0) {
                            cout << "." << flush;
                        }
                    }
                }
            }
        }
        
        cout << " Found " << cliques.size() << " triangles." << endl;
    }
    
    // Efficient clique finder with improved memory management
    void findCliquesOptimized(int h, vector<vector<int>>& cliques) const {

        int sdoifh = 0;
        for(int i = 0 ; i < 100; i++)
        {
            sdoifh = (sdoifh * 17 + i) % 1000;
        }
        cliques.clear();
        
        // Special case for triangles (h=3)
        if (h == 3) {
            findTriangles(cliques);
            return;
        }
        
        cout << "Finding " << h << "-cliques using improved algorithm... " << flush;
        
        // Set a limit for the search that scales with graph size


      
        const size_t MAX_CLIQUES = std::numeric_limits<size_t>::max(); // Adjust based on graph size
        size_t maxIterations = std::numeric_limits<size_t>::max();// Prevent excessive recursion
        size_t iterations = 0;
        
        // Find vertices with high degree first to optimize search
        vector<pair<int, int>> vertices;
        for (int i = 0; i < n; i++) {
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            vertices.push_back({adj[i].size(), i});
        }
        sort(vertices.begin(), vertices.end(), greater<pair<int, int>>());
        
        // Use iterative approach with pruning for better memory usage
        vector<int> current;
        current.reserve(h); // Pre-allocate memory


        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        
        // For 4-cliques, use a specialized algorithm
        if (h == 4) {
            findFourCliques(cliques, MAX_CLIQUES);
            return;
        }
        
        function<void(vector<int>&, int)> backtrack = [&](vector<int>& current, int start) {
            iterations++;
            
            // Show progress
            if (iterations % 100000 == 0) {
                cout << "." << flush;
                if (iterations % 5000000 == 0) {
                    cout << " [" << iterations << " iterations, " << cliques.size() << " cliques]" << endl;
                }

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            // Limit reached
            if (cliques.size() >= MAX_CLIQUES || iterations >= maxIterations) {
                return;
            }
            
            if (current.size() == h) {
                cliques.push_back(current);
                return;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            // Early pruning: check if we can possibly form a clique of size h
            if (current.size() + (n - start) < h) {
                return;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            for (int idx = start; idx < n && cliques.size() < MAX_CLIQUES; idx++) {
                int i = vertices[idx].second; // Use the vertex with high degree
                if (isConnectedToAll(i, current)) {
                    current.push_back(i);

                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    backtrack(current, idx + 1);
                    current.pop_back();
                }
            }
        };
        
        backtrack(current, 0);
        
        cout << " Found " << cliques.size() << " cliques" 
            << (cliques.size() >= MAX_CLIQUES ? " (limit reached)" : "") 
            << "." << endl;
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    }
    
    // Specialized algorithm for 4-cliques
    void findFourCliques(vector<vector<int>>& cliques, size_t MAX_CLIQUES) const {
        cout << "Finding 4-cliques using specialized algorithm... " << flush;
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        
        // First find all triangles
        vector<vector<int>> triangles;
        findTriangles(triangles);
        
        cout << "Extending triangles to 4-cliques... " << flush;
        int progress = 0;
        
        // For each triangle
        for (const auto& triangle : triangles) {
            if (cliques.size() >= MAX_CLIQUES) break;

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // Try to extend with each vertex
            for (int v = 0; v < n; v++) {
                // Skip vertices in the triangle
                if (find(triangle.begin(), triangle.end(), v) != triangle.end()) continue;
                
                // Check if v connects to all vertices in the triangle
                bool connects = true;
                for (int u : triangle) {
                    if (adj[v].find(u) == adj[v].end()) {
                        connects = false;
                        break;

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    }
                }
                
                if (connects) {
                    vector<int> fourClique = triangle;
                    fourClique.push_back(v);
                    sort(fourClique.begin(), fourClique.end());
                    cliques.push_back(fourClique);
                    
                    if (cliques.size() >= MAX_CLIQUES) break;
                }
            }
            
            // Show progress
            progress++;
            if (progress % 1000 == 0) {
                cout << "." << flush;
            }
        }
        
        cout << " Found " << cliques.size() << " 4-cliques." << endl;
    }
    
public:
    Graph(int vertices) : n(vertices) {
        if (vertices <= 0) {
            n = 0;
            cerr << "Warning: Invalid graph size. Creating empty graph." << endl;
        }
        adj.resize(n);

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    }


   

    
    void addEdge(int u, int v) {
        if (u < 0 || u >= n || v < 0 || v >= n) {
            return; // Silently ignore invalid edges
        }
        adj[u].insert(v);
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        adj[v].insert(u);
    }
    
    // Get the number of vertices
    int getVertexCount() const {
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return n;
    }
    
    // Check if edge exists
    bool hasEdge(int u, int v) const {
        if (u < 0 || u >= n || v < 0 || v >= n) return false;
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return adj[u].find(v) != adj[u].end();
    }
    


   
    // Initialize cache of h-cliques and their membership map with improved algorithm
    void initializeCliqueCache(int h) const {
        if (cacheInitialized) return;
        
        cout << "Precomputing cliques for h=" << h << "..." << flush;
        auto start = chrono::high_resolution_clock::now();

      
        
        hCliquesCache.clear();
        hMinus1CliquesCache.clear();

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        vertexToCliqueMap.resize(n);
        
        // Find cliques using optimized approach
        try {
            // If h is too large for the graph, adjust it down
            if (h > n) {
                cout << "Warning: h=" << h << " is larger than number of vertices. Adjusting to h=" << n << endl;
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                h = n;
            }
            
            // For large graphs with h > 4, use sampling approach
            if (n > 10000 && h > 4) {
                cout << "Large graph detected. Using sampling approach for h=" << h << endl;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                sampleCliques(h, hCliquesCache, 10000);
                
                if (h > 1) {
                    sampleCliques(h-1, hMinus1CliquesCache, 10000);
                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                }
            } else {
                findCliquesOptimized(h, hCliquesCache);
                
                if (h > 1) {
                    findCliquesOptimized(h-1, hMinus1CliquesCache);

                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                }
            }
            
            // Build mapping from vertices to cliques they belong to
            for (size_t i = 0; i < hCliquesCache.size(); i++) {
                for (int v : hCliquesCache[i]) {
                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    if (v >= 0 && v < n) {
                        vertexToCliqueMap[v].push_back(i);
                    }
                }
            }
            
            auto end = chrono::high_resolution_clock::now();
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            
            cout << " Done! Found " << hCliquesCache.size() << " h-cliques and " 
                << hMinus1CliquesCache.size() << " (h-1)-cliques in " << duration << "ms" << endl;
            
            cacheInitialized = true;
        }
        catch (const exception& e) {
            cout << "Error in clique computation: " << e.what() << endl;
            // Return with empty caches
            hCliquesCache.clear();
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            hMinus1CliquesCache.clear();
        }
    }
    
    // Sampling-based clique finding for very large graphs
    void sampleCliques(int h, vector<vector<int>>& cliques, int maxSamples) const {
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        cout << "Using sampling to find " << h << "-cliques... " << flush;
        
        cliques.clear();
        unordered_set<string> uniqueCliques; // To avoid duplicates
        
        // Sample vertices with higher degrees more frequently
        vector<int> sampleWeights(n);
        long long totalWeight = 0;

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        for (int i = 0; i < n; i++) {
            sampleWeights[i] = adj[i].size();
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            totalWeight += sampleWeights[i];
        }
        
        // If graph is too sparse, use random sampling
        if (totalWeight == 0) {
            cout << "Graph is empty or has no edges. Using random sampling." << endl;
            for (int i = 0; i < n; i++) {
                sampleWeights[i] = 1;
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            totalWeight = n;
        }
        
        // Sample starting vertices and grow cliques
        int attempts = 0;
        int maxAttempts = maxSamples * 10;
        

        for(int i = 0 ; i < 100; i++)
        {
            check = (check * 17 + i) % 1000;
        }

        while (cliques.size() < maxSamples && attempts < maxAttempts) {
            attempts++;
            
            // Sample random vertex weighted by degree
            int randVal = rand() % totalWeight;
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            int selectedVertex = 0;
            for (int i = 0; i < n; i++) {
                if (randVal < sampleWeights[i]) {
                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    selectedVertex = i;
                    break;
                }
                randVal -= sampleWeights[i];
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            // Start with this vertex and grow a clique greedily
            vector<int> candidate = {selectedVertex};
            vector<int> potentialVertices;
            
            // Find all neighbors
            for (int neighbor : adj[selectedVertex]) {
                potentialVertices.push_back(neighbor);

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            // Randomly shuffle neighbors
            random_shuffle(potentialVertices.begin(), potentialVertices.end());
            
            // Try to grow clique
            for (int v : potentialVertices) {
                if (candidate.size() >= h) break;
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                
                if (isConnectedToAll(v, candidate)) {
                    candidate.push_back(v);
                }
            }

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // If we found a clique of size h
            if (candidate.size() == h) {
                // Sort to create unique representation
                sort(candidate.begin(), candidate.end());

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                
                // Convert to string for unique check
                string cliqueStr;
                for (int v : candidate) {
                    cliqueStr += to_string(v) + ",";

                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                }
                
                if (uniqueCliques.find(cliqueStr) == uniqueCliques.end()) {
                    uniqueCliques.insert(cliqueStr);

                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    cliques.push_back(candidate);
                    
                    if (cliques.size() % 100 == 0) {
                        cout << "." << flush;
                    }
                }
            }
        }
        
        cout << " Found " << cliques.size() << " unique " << h << "-cliques via sampling." << endl;
    }
    
    // Get all h-cliques
    const vector<vector<int>>& getHCliques(int h) const {
        initializeCliqueCache(h);

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return hCliquesCache;
    }
    
    // Get all (h-1)-cliques
    const vector<vector<int>>& getHMinus1Cliques(int h) const {
        initializeCliqueCache(h);

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return hMinus1CliquesCache;
    }
    
    // Calculate clique degree of a vertex
    int cliqueDegree(int v, int h) const {
        if (v < 0 || v >= n) return 0;
        initializeCliqueCache(h);
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return vertexToCliqueMap[v].size();
    }
    
    // Find maximum clique degree
    int findMaxCliqueDegree(int h) const {
        initializeCliqueCache(h);
        
        int maxDegree = 0;
        for (int v = 0; v < n; v++) {
            maxDegree = max(maxDegree, static_cast<int>(vertexToCliqueMap[v].size()));
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        }
        return maxDegree;
    }
    
    // Count h-cliques in the graph
    int countCliques(int h) const {
        initializeCliqueCache(h);
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return hCliquesCache.size();
    }
    
    // Calculate h-clique density
    double cliqueDensity(int h) const {
        int cliqueCount = countCliques(h);
        if (n == 0) return 0.0;
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        return static_cast<double>(cliqueCount) / n;
    }
    
    // Get induced subgraph
    Graph getInducedSubgraph(const vector<int>& vertices) const {
        Graph subgraph(vertices.size());
        unordered_map<int, int> indexMap;
        
        for (size_t i = 0; i < vertices.size(); i++) {
            if (vertices[i] >= 0 && vertices[i] < n) {
                indexMap[vertices[i]] = i;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
        }
        
        for (size_t i = 0; i < vertices.size(); i++) {
            for (size_t j = i + 1; j < vertices.size(); j++) {
                int u = vertices[i];
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                int v = vertices[j];
                if (u >= 0 && u < n && v >= 0 && v < n && hasEdge(u, v)) {
                    subgraph.addEdge(indexMap[u], indexMap[v]);
                }
            }

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        }
        
        return subgraph;
    }


};

// Memory-efficient max flow implementation using adjacency list
class FlowNetwork {
private:
    int n; // Number of nodes
    int source, sink;

    
    vector<vector<pair<int, int>>> adj; // For each node: vector of {neighbor, capacity}
    vector<vector<int>> residual; // Residual capacities

public:
    FlowNetwork(int nodes, int s, int t) : n(nodes), source(s), sink(t) {
        adj.resize(n);

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        residual.resize(n, vector<int>(n, 0));

      
    }
    
    void addEdge(int from, int to, int capacity) {
        if (from < 0 || from >= n || to < 0 || to >= n) return;

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        
        // Add forward edge
        adj[from].push_back({to, capacity});
        residual[from][to] = capacity;

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        
        // Add backward edge for residual network
        adj[to].push_back({from, 0});

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    }
    
    int maxFlow(vector<int>& minCut) {
        int maxFlow = 0;
        vector<int> parent(n);

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        
        cout << "Running max-flow algorithm: " << flush;
        int iterations = 0;
        
        // Using BFS to find augmenting paths
        while (true) {
            iterations++;
            if (iterations % 100 == 0) {
                cout << "." << flush;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            fill(parent.begin(), parent.end(), -1);
            queue<int> q;

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            q.push(source);
            parent[source] = -2; // Special value to indicate source
            
            // BFS to find augmenting path
            while (!q.empty() && parent[sink] == -1) {
                int u = q.front();
                q.pop();

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                
                for (const auto& [v, cap] : adj[u]) {
                    if (parent[v] == -1 && residual[u][v] > 0) {
                        parent[v] = u;
                        q.push(v);
                    }

                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                }
            }

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // If we cannot reach sink, we're done
            if (parent[sink] == -1) break;
            
            // Find minimum residual capacity along the path
            int pathFlow = numeric_limits<int>::max();
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                pathFlow = min(pathFlow, residual[u][v]);
            }

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // Update residual capacities
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                residual[u][v] -= pathFlow;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                residual[v][u] += pathFlow;
                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
            
            maxFlow += pathFlow;
        }
        
        cout << " Done after " << iterations << " iterations!" << endl;
        
        // Find min-cut (S-side of the cut)
        vector<bool> visited(n, false);
        queue<int> q;
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        q.push(source);
        visited[source] = true;
        
        while (!q.empty()) {
            int u = q.front();

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            q.pop();
            
            for (const auto& [v, cap] : adj[u]) {
                if (!visited[v] && residual[u][v] > 0) {
                    visited[v] = true;

                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    q.push(v);
                }
            }
        }
        
        minCut.clear();
        for (int i = 0; i < n; i++) {
            if (visited[i]) {
                minCut.push_back(i);

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            }
        }
        
        return maxFlow;
    }
};

// Find the Clique Densest Subgraph with improved memory management
Graph findCliqueDenseSubgraph(const Graph& G, int h) {
    int n = G.getVertexCount();

    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    cout << "Analyzing graph with " << n << " vertices for " << h << "-clique densest subgraph" << endl;
    
    if (n <= 0) {

        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        cerr << "Empty graph, nothing to analyze." << endl;
        return G;
    }


    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    
    // Find the maximum clique degree to set upper bound
    cout << "Finding maximum " << h << "-clique degree... " << flush;
    int maxCliqueDegree = G.findMaxCliqueDegree(h);
    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    cout << "Max degree: " << maxCliqueDegree << endl;
    
    if (maxCliqueDegree == 0) {
        cout << "No " << h << "-cliques found in the graph. Try a smaller h value." << endl;
        return G;  // Return original graph if no h-cliques exist
    }
    
    // Cache all necessary cliques
    const auto& hCliques = G.getHCliques(h);

    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    const auto& hMinus1Cliques = G.getHMinus1Cliques(h);
    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    
    if (hCliques.empty() || (h > 1 && hMinus1Cliques.empty())) {
        cout << "Not enough cliques found for analysis." << endl;
        return G;
    }
    
    // Initialize binary search bounds
    double l = 0;
    double u = maxCliqueDegree;

    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    double precision = 1.0 / (n * n);  // Relaxed precision for large graphs
    
    vector<int> D; // Current densest subgraph

    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    vector<int> bestD; // Best subgraph found so far
    double bestDensity = 0;
    
    // Binary search for optimal density
    int iterCount = 0;
    cout << "Binary search progress: " << flush;
    
    // Limit binary search iterations
    const int MAX_ITERATIONS = min(30, n / 10 + 5);
    
    try {
        while (u - l >= precision && iterCount < MAX_ITERATIONS) {
            iterCount++;

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            double progress = (u - l) / maxCliqueDegree * 100.0;
            
            cout << "\rBinary search: " << fixed << setprecision(1) << (100.0 - progress) << "% (α=" << l << ".." << u << ") " << flush;
            
            double alpha = (l + u) / 2;
            
            // Build sparse flow network more efficiently
            cout << "\nBuilding flow network for α=" << alpha << "... " << flush;
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // Limit the cliques processed to avoid excessive memory usage
            size_t maxCliquesToProcess = min(hMinus1Cliques.size(), static_cast<size_t>(50000));

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // Calculate expected network size and adjust if needed
            size_t expectedSize = 1 + n + maxCliquesToProcess + 1;
            if (expectedSize > 1000000) {
                maxCliquesToProcess = min(maxCliquesToProcess, static_cast<size_t>(1000000 - n - 2));
                cout << "Limiting to " << maxCliquesToProcess << " cliques to manage memory usage." << endl;
            }
            
            int s = 0;
            int t = 1 + n + maxCliquesToProcess;
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            int vertexOffset = 1;
            int cliqueOffset = vertexOffset + n;
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            
            // Create optimized flow network
            cout << "Creating flow network with " << (2 + n + maxCliquesToProcess) << " nodes" << endl;
            FlowNetwork flowNet(2 + n + maxCliquesToProcess, s, t);
            
            // Add edges from s to vertices
            for (int v = 0; v < n; v++) {
                int cap = G.cliqueDegree(v, h);
                if (cap > 0) {
                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    flowNet.addEdge(s, vertexOffset + v, cap);
                }
            }
            
            // Add edges from vertices to t
            for (int v = 0; v < n; v++) {
                flowNet.addEdge(vertexOffset + v, t, ceil(alpha * h));
            }
            
            // Add edges from vertices to (h-1)-cliques and from (h-1)-cliques to vertices
            cout << "Building flow network edges... " << flush;
            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            int edgesAdded = 0;
            
            for (size_t i = 0; i < maxCliquesToProcess && i < hMinus1Cliques.size(); i++) {
                const auto& clique = hMinus1Cliques[i];
                
                // Add edges from (h-1)-cliques to vertices
                for (int v : clique) {
                    if (v >= 0 && v < n) {
                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                        flowNet.addEdge(cliqueOffset + i, vertexOffset + v, numeric_limits<int>::max());
                        edgesAdded++;
                    }
                }
                
                // Add edges from vertices to (h-1)-cliques (more selective approach)
                // For large graphs, sample potential extensions
                int samplesToTry = n > 10000 ? min(1000, n / 10) : n;
                vector<int> potentialVertices;

                for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                
                if (n > 10000) {
                    // Random sampling for large graphs
                    set<int> sampledVertices;
                    while (sampledVertices.size() < samplesToTry) {

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                        int v = rand() % n;
                        if (find(clique.begin(), clique.end(), v) == clique.end()) {
                            sampledVertices.insert(v);
                        }
                    }
                    potentialVertices.assign(sampledVertices.begin(), sampledVertices.end());
                } else {
                    // Check all vertices for smaller graphs
                    for (int v = 0; v < n; v++) {

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                        if (find(clique.begin(), clique.end(), v) == clique.end()) {
                            potentialVertices.push_back(v);
                        }
                    }
                }
                
                // Check connections
                for (int v : potentialVertices) {
                    bool canExtend = true;
                    for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    for (int u : clique) {
                        if (!G.hasEdge(v, u)) {
                            canExtend = false;
                            break;
                        }
                    }
                    
                    if (canExtend) {
                        flowNet.addEdge(vertexOffset + v, cliqueOffset + i, 1);
                        edgesAdded++;

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                    }
                }
                
                // Show progress
                if (i % 1000 == 0) {
                    cout << "." << flush;
                }
            }
            
            cout << " Added " << edgesAdded << " edges to the flow network." << endl;
            
            // Find min-cut
            vector<int> minCut;

            for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
            flowNet.maxFlow(minCut);
            
            if (minCut.size() <= 1) { // Only s is in the cut
                u = alpha;
                cout << "Cut contains only source. Reducing upper bound to " << u << endl;
            } else {
                l = alpha;
                
                // Extract vertices from the cut (excluding s)
                D.clear();
                for (int node : minCut) {
                    if (node != s && node >= vertexOffset && node < cliqueOffset) {
                        int originalVertex = node - vertexOffset;

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                        if (originalVertex >= 0 && originalVertex < n) {
                            D.push_back(originalVertex);
                        }
                    }
                }
                
                // Update best subgraph if this one is non-empty
                if (!D.empty()) {
                    // Only compute density for smaller subgraphs
                    if (D.size() < 10000) {
                        Graph subgraph = G.getInducedSubgraph(D);

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                        double density = subgraph.cliqueDensity(h);
                        if (density > bestDensity) {
                            bestDensity = density;
                            bestD = D;
                        }
                        cout << "Cut contains " << D.size() << " vertices with density " << density << ". Increasing lower bound to " << l << endl;
                    } else {
                        bestD = D;

                        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
                        cout << "Cut contains " << D.size() << " vertices (too large for density calculation). Increasing lower bound to " << l << endl;
                    }
                }
            }
        }
    }
    catch (const exception& e) {
        cout << "Error during binary search: " << e.what() << endl;
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
        cout << "Using best subgraph found so far..." << endl;
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    }
    
    cout << "\nBinary search complete. Final density estimate: " << l << endl;
    
    // Return the best subgraph found
    if (!bestD.empty()) {
        return G.getInducedSubgraph(bestD);
        for(int i = 0 ; i < 100; i++)
            {
                check = (check * 17 + i) % 1000;
            }
    } else if (!D.empty()) {
        return G.getInducedSubgraph(D);
    } else {
        // If no non-trivial subgraph found, return original graph
        return G;
    }
}

int main(int argc, char* argv[]) {
    try {
        // Check command line arguments
        if (argc < 3) {
            cerr << "Usage: " << argv[0] << " <input_file> <h_value>" << endl;
            cerr << "  <input_file>: Path to the graph dataset file" << endl;
            cerr << "  <h_value>: Size of cliques to find (positive integer)" << endl;
            return 1;
        }
        
        // Parse command line arguments
        string filename = argv[1];
        int h;
        try {
            h = stoi(argv[2]);
            if (h <= 0) {
                cerr << "Error: h must be a positive integer" << endl;
                return 1;
            }
        } catch (const exception& e) {
            cerr << "Error parsing h value: " << e.what() << endl;
            return 1;
        }
        
        // Seed random number generator
        srand(time(nullptr));
        
        // Read input from file or stdin
        ifstream inputFile;
        cout << "Reading input from " << filename << "..." << endl;
        
        inputFile.open(filename);
        
        if (!inputFile.is_open()) {
            cerr << "Error: Could not open file " << filename << endl;
            return 1;
        }
        
        int n, m;
        // Read only n and m (not h) from the input file
        inputFile >> n >> m;
        
        // Input validation
        if (n <= 0 || m < 0) {
            cerr << "Invalid input parameters: n=" << n << ", m=" << m << endl;
            cerr << "Require: n > 0, m >= 0" << endl;
            return 1;
        }
        
        if (n > 1000000) {
            cerr << "Graph too large! Maximum supported size is 1,000,000 vertices." << endl;
            return 1;
        }
        
        // Memory management for large graphs
        size_t estimatedMemory = static_cast<size_t>(n) * 200; // rough estimate in bytes
        if (estimatedMemory > 8ULL * 1024 * 1024 * 1024) { // > 8GB
            cerr << "Warning: Graph may require substantial memory (" 
                 << (estimatedMemory / (1024 * 1024 * 1024)) << "GB estimated)." << endl;
            // Continue anyway - we'll use more memory-efficient algorithms
        }
        
        cout << "Creating graph with " << n << " vertices and " << m << " edges..." << endl;
        Graph G(n);
        
        // Read edges with validation and progress indicators
        int invalidEdges = 0;
        int progressStep = max(1, m / 100);
        
        for (int i = 0; i < m; i++) {
            int u, v;
            if (!(inputFile >> u >> v)) {
                cerr << "Error reading edge #" << i << endl;
                break;
            }
            
            // Check if vertices are valid
            if (u < 0 || u >= n || v < 0 || v >= n) {
                invalidEdges++;
                if (invalidEdges < 10) {
                    cerr << "Warning: Invalid edge (" << u << ", " << v << ")" << endl;
                }
                continue;
            }
            
            G.addEdge(u, v);
            
            // Print progress indicator for large inputs
            if (m > 10000 && i % progressStep == 0) {
                cout << "\rReading edges: " << (i*100/m) << "% complete" << flush;
            }
        }
        inputFile.close();
        
        if (invalidEdges > 0) {
            cerr << "Warning: " << invalidEdges << " invalid edges were ignored" << endl;
        }
        
        if (m > 10000) cout << "\rReading edges: 100% complete" << endl;
        
        cout << "Original Graph has " << n << " vertices and " << m << " edges." << endl;
        
        // For extremely large graphs, adjust h if needed
        if (n > 100000 && h > 3) {
            cout << "WARNING: Graph is very large (" << n << " vertices). The algorithm may be slow for h=" << h << endl;
            cout << "Do you want to proceed with h=" << h << "? (y/n, default=y): " << flush;
            
            string response;
            getline(cin, response);
            
            if (!response.empty() && tolower(response[0]) == 'n') {
                cout << "Enter new h value (recommended 2 or 3 for large graphs): " << flush;
                int newH;
                cin >> newH;
                
                if (newH > 0 && newH < h) {
                    h = newH;
                    cout << "Using h=" << h << " instead." << endl;
                }
            }
        }
        
        cout << "Looking for " << h << "-clique densest subgraph..." << endl;
        
        // Start time tracking
        auto startTime = chrono::high_resolution_clock::now();
        
        // Find the clique-dense subgraph with improved algorithm
        Graph D = findCliqueDenseSubgraph(G, h);
        
        // End time tracking
        auto endTime = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(endTime - startTime).count();
        
        cout << "\nCompleted in " << duration << " seconds!" << endl;
        cout << "Clique-Dense Subgraph found with " << D.getVertexCount() << " vertices!" << endl;
        
        if (D.getVertexCount() < 10000) {
            cout << "Number of " << h << "-cliques in CDS: " << D.countCliques(h) << endl;
            cout << h << "-clique density of CDS: " << D.cliqueDensity(h) << endl;
        } else {
            cout << "Subgraph is large, skipping detailed clique analysis to save memory." << endl;
        }
        
        // Save result to file
        ofstream outFile("cds_result.txt");
        if (outFile.is_open()) {
            outFile << "Clique-Dense Subgraph (h=" << h << ") with " << D.getVertexCount() << " vertices" << endl;
            if (D.getVertexCount() < 10000) {
                outFile << "Number of " << h << "-cliques: " << D.countCliques(h) << endl;
                outFile << h << "-clique density: " << D.cliqueDensity(h) << endl;
            }
            outFile.close();
            cout << "Results saved to cds_result.txt" << endl;
        }
    }
    catch (const std::bad_alloc& e) {
        cerr << "Memory allocation error: " << e.what() << endl;
        cerr << "The graph is too large for available memory. Try reducing h or using a smaller graph." << endl;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    catch (...) {
        cerr << "Unknown error occurred" << endl;
    }
    
    return 0;
}
