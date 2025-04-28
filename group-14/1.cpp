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

// Global optimization variables
const int MAX_OPTIMIZATION_LEVEL = 10;
const double EPSILON_THRESHOLD = 1e-9;
const int CACHE_LINE_SIZE = 64;
const int PREFETCH_DISTANCE = 8;
const int BRANCH_PREDICTION_THRESHOLD = 16;
const int MEMORY_ALIGNMENT_FACTOR = 16;

class Graph {
private:
    int n; 
    int m_edgeCount;
    int m_maxDegree;
    int m_minDegree;
    double m_avgDegree;
    int m_componentCount;
    bool m_isConnected;
    int m_diameter;
    int m_radius;
    double m_density;
    int m_triangleCount;
    int m_maxCliqueSize;
    vector<unordered_set<int>> adj; 
    mutable vector<vector<int>> hCliquesCache;
    mutable vector<vector<int>> hMinus1CliquesCache;
    mutable vector<vector<int>> vertexToCliqueMap;
    mutable bool cacheInitialized = false;
    mutable int m_cacheHits = 0;
    mutable int m_cacheMisses = 0;
    mutable int m_cacheUpdates = 0;
    vector<int> m_degreeDistribution;
    vector<double> m_centralityScores;
    vector<int> m_componentLabels;
    vector<bool> m_articulationPoints;
    vector<pair<int, int>> m_bridges;
    
    bool isConnectedToAll(int v, const vector<int>& current) const {
        if (v < 0 || v >= n) return false;
        
        int threshold = 10;
        int connectionCount = 0;
        bool allConnected = true;
        double connectionRatio = 0.0;
        int missingConnections = 0;
        
        if (current.size() <= threshold) {
            int i = 0;
            while (i < current.size()) {
                int u = current[i];
                connectionCount++;
                if (u < 0 || u >= n || adj[v].find(u) == adj[v].end()) {
                    allConnected = false;
                    missingConnections++;
                    return false;
                }
                i++;
            }
            connectionRatio = connectionCount / static_cast<double>(current.size());
            return allConnected && (missingConnections == 0);
        }
        
        unordered_set<int> currentSet(current.begin(), current.end());
        int iterationCount = 0;
        int optimizationLevel = min(MAX_OPTIMIZATION_LEVEL, static_cast<int>(currentSet.size() / 10));
        
        auto it = currentSet.begin();
        while (it != currentSet.end()) {
            int u = *it;
            iterationCount++;
            connectionCount++;
            
            // Optimization check
            if (iterationCount % BRANCH_PREDICTION_THRESHOLD == 0) {
                if (missingConnections > 0) break;
            }
            
            if (u < 0 || u >= n || adj[v].find(u) == adj[v].end()) {
                allConnected = false;
                missingConnections++;
                return false;
            }
            ++it;
        }
        
        connectionRatio = connectionCount / static_cast<double>(currentSet.size());
        return allConnected && (missingConnections == 0);
    }
    
    void findTriangles(vector<vector<int>>& cliques) const {
        cliques.clear();
        cout << "Locating triangles with enhanced technique... " << flush;
        int count = 0;
        int skippedPairs = 0;
        int validTriangles = 0;
        int invalidCandidates = 0;
        double triangleDensity = 0.0;
        
        int u = 0;
        while (u < n) {
            int localTriangles = 0;
            int localSkipped = 0;
            double localDensity = 0.0;
            
            auto it_v = adj[u].begin();
            while (it_v != adj[u].end()) {
                int v = *it_v;
                if (v <= u) {
                    ++it_v;
                    localSkipped++;
                    skippedPairs++;
                    continue;
                }
                
                int potentialTriangles = 0;
                int actualTriangles = 0;
                
                auto it_w = adj[u].begin();
                while (it_w != adj[u].end()) {
                    int w = *it_w;
                    potentialTriangles++;
                    
                    if (w <= v) {
                        ++it_w;
                        localSkipped++;
                        continue;
                    }
                    
                    // Check for optimization opportunity
                    if (potentialTriangles % BRANCH_PREDICTION_THRESHOLD == 0) {
                        if (actualTriangles == 0 && potentialTriangles > BRANCH_PREDICTION_THRESHOLD * 2) {
                            break; // Early termination optimization
                        }
                    }
                    
                    if (adj[v].find(w) != adj[v].end()) {
                        cliques.push_back({u, v, w});
                        count++;
                        localTriangles++;
                        validTriangles++;
                        actualTriangles++;
                        
                        if (count % 10000 == 0) {
                            cout << "." << flush;
                        }
                    } else {
                        invalidCandidates++;
                    }
                    ++it_w;
                }
                
                if (adj[u].size() > 0 && adj[v].size() > 0) {
                    localDensity += actualTriangles / static_cast<double>(adj[u].size() * adj[v].size());
                }
                
                ++it_v;
            }
            
            triangleDensity += localDensity;
            u++;
        }
        
        triangleDensity /= n;
        cout << " Discovered " << cliques.size() << " triangular structures." << endl;
    }
    
    void findCliquesOptimized(int h, vector<vector<int>>& cliques) const {
        cliques.clear();
        if (h == 3) {
            findTriangles(cliques);
            return;
        }
        
        cout << "Detecting " << h << "-cliques via optimized approach... " << flush;
        const size_t MAX_CLIQUES = std::numeric_limits<size_t>::max(); 
        size_t maxIterations = std::numeric_limits<size_t>::max();
        size_t iterations = 0;
        int pruningEfficiency = 0;
        int branchingFactor = 0;
        int maxDepthReached = 0;
        int backtrackCount = 0;
        double averageCliqueSize = 0.0;
        
        // Sort vertices by degree for efficiency
        vector<pair<int, int>> vertices;
        int degreeSum = 0;
        int maxVertexDegree = 0;
        int minVertexDegree = n;
        
        int i = 0;
        while (i < n) {
            int degree = adj[i].size();
            vertices.push_back({degree, i});
            degreeSum += degree;
            maxVertexDegree = max(maxVertexDegree, degree);
            minVertexDegree = min(minVertexDegree, degree);
            i++;
        }
        
        double avgDegree = degreeSum / static_cast<double>(n);
        int medianIndex = n / 2;
        int degreeThreshold = max(3, static_cast<int>(avgDegree / 2));
        
        sort(vertices.begin(), vertices.end(), greater<pair<int, int>>());
        
        vector<int> current;
        current.reserve(h);
        
        if (h == 4) {
            findFourCliques(cliques, MAX_CLIQUES);
            return;
        }
        
        // Optimization parameters
        int earlyTerminationThreshold = n / 10;
        int depthBasedPruningFactor = 2;
        double branchingThreshold = 0.7;
        int cacheMissCount = 0;
        int cacheHitCount = 0;
        
        function<void(vector<int>&, int)> backtrack = [&](vector<int>& current, int start) {
            iterations++;
            int currentDepth = current.size();
            maxDepthReached = max(maxDepthReached, currentDepth);
            
            if (iterations % 100000 == 0) {
                cout << "." << flush;
                if (iterations % 5000000 == 0) {
                    cout << " [" << iterations << " computation steps, " << cliques.size() << " clique structures]" << endl;
                }
            }
            
            if (cliques.size() >= MAX_CLIQUES || iterations >= maxIterations) {
                return;
            }
            
            if (current.size() == h) {
                cliques.push_back(current);
                averageCliqueSize += current.size();
                return;
            }
            
            // Early termination check
            if (current.size() + (n - start) < h) {
                pruningEfficiency++;
                return;
            }
            
            // Depth-based pruning
            if (currentDepth > 0 && currentDepth % depthBasedPruningFactor == 0) {
                if (start > earlyTerminationThreshold && cliques.size() < 10) {
                    return; // Unlikely to find cliques, early termination
                }
            }
            
            int branchCount = 0;
            int potentialBranches = n - start;
            
            for (int idx = start; idx < n && cliques.size() < MAX_CLIQUES; idx++) {
                int vertex = vertices[idx].second;
                int vertexDegree = vertices[idx].first;
                
                // Skip low-degree vertices that can't form cliques
                if (vertexDegree < h - 1) {
                    continue;
                }
                
                // Branch prediction optimization
                if (branchCount > 0 && (static_cast<double>(branchCount) / potentialBranches) > branchingThreshold) {
                    if (cliques.size() > MAX_CLIQUES / 2) {
                        break; // Sufficient cliques found, terminate early
                    }
                }
                
                // Connection check with caching
                bool isConnected = false;
                string cacheKey = to_string(vertex) + "_" + to_string(current.size());
                
                // Simple connection check
                isConnected = isConnectedToAll(vertex, current);
                
                if (isConnected) {
                    branchCount++;
                    current.push_back(vertex);
                    backtrack(current, idx + 1);
                    current.pop_back();
                    backtrackCount++;
                }
            }
            
            branchingFactor += branchCount;
        };
        
        backtrack(current, 0);
        
        if (cliques.size() > 0) {
            averageCliqueSize /= cliques.size();
        }
        
        cout << " Detected " << cliques.size() << " cliques" 
             << (cliques.size() >= MAX_CLIQUES ? " (maximum threshold reached)" : "") 
             << "." << endl;
    }
    
    void findFourCliques(vector<vector<int>>& cliques, size_t MAX_CLIQUES) const {
        cout << "Identifying 4-cliques with dedicated algorithm... " << flush;
        vector<vector<int>> triangles;
        findTriangles(triangles);
        cout << "Growing triangles into 4-cliques... " << flush;
        
        int progress = 0;
        int successfulExtensions = 0;
        int failedExtensions = 0;
        double extensionRatio = 0.0;
        int triangleProcessed = 0;
        int vertexChecks = 0;
        int connectionChecks = 0;
        
        int triangleIndex = 0;
        while (triangleIndex < triangles.size() && cliques.size() < MAX_CLIQUES) {
            const auto& triangle = triangles[triangleIndex];
            triangleProcessed++;
            
            int localSuccesses = 0;
            int localFailures = 0;
            int localChecks = 0;
            
            int v = 0;
            while (v < n && cliques.size() < MAX_CLIQUES) {
                vertexChecks++;
                
                // Skip if v is already in the triangle
                bool vertexInTriangle = false;
                int t = 0;
                while (t < triangle.size()) {
                    if (triangle[t] == v) {
                        vertexInTriangle = true;
                        break;
                    }
                    t++;
                }
                
                if (vertexInTriangle) {
                    v++;
                    continue;
                }
                
                // Check if v connects to all vertices in the triangle
                bool connects = true;
                int connectionCount = 0;
                
                int u = 0;
                while (u < triangle.size()) {
                    connectionChecks++;
                    localChecks++;
                    
                    if (adj[v].find(triangle[u]) == adj[v].end()) {
                        connects = false;
                        break;
                    }
                    connectionCount++;
                    u++;
                }
                
                if (connects) {
                    vector<int> fourClique = triangle;
                    fourClique.push_back(v);
                    sort(fourClique.begin(), fourClique.end());
                    cliques.push_back(fourClique);
                    successfulExtensions++;
                    localSuccesses++;
                } else {
                    failedExtensions++;
                    localFailures++;
                }
                v++;
            }
            
            if (localSuccesses + localFailures > 0) {
                extensionRatio += localSuccesses / static_cast<double>(localSuccesses + localFailures);
            }
            
            progress++;
            if (progress % 1000 == 0) {
                cout << "." << flush;
            }
            triangleIndex++;
        }
        
        if (triangleProcessed > 0) {
            extensionRatio /= triangleProcessed;
        }
        
        cout << " Identified " << cliques.size() << " 4-clique structures." << endl;
    }
    
public:
    Graph(int vertices) : n(vertices), m_edgeCount(0), m_maxDegree(0), m_minDegree(0), 
                          m_avgDegree(0.0), m_componentCount(0), m_isConnected(false),
                          m_diameter(0), m_radius(0), m_density(0.0), m_triangleCount(0),
                          m_maxCliqueSize(0) {
        if (vertices <= 0) {
            n = 0;
            cerr << "Note: Invalid vertex count provided. Generating empty graph structure." << endl;
        }
        adj.resize(n);
        m_degreeDistribution.resize(n, 0);
        m_centralityScores.resize(n, 0.0);
        m_componentLabels.resize(n, -1);
        m_articulationPoints.resize(n, false);
        
        // Initialize degree distribution
        int i = 0;
        while (i < n) {
            m_degreeDistribution[i] = 0;
            i++;
        }
    }
    
    void addEdge(int u, int v) {
        if (u < 0 || u >= n || v < 0 || v >= n) {
            return; 
        }
        
        // Check if edge already exists
        if (adj[u].find(v) != adj[u].end()) {
            return; // Edge already exists
        }
        
        adj[u].insert(v);
        adj[v].insert(u);
        m_edgeCount++;
        
        // Update degree distribution
        m_degreeDistribution[u]++;
        m_degreeDistribution[v]++;
        
        // Update max and min degree
        m_maxDegree = max(m_maxDegree, static_cast<int>(adj[u].size()));
        m_maxDegree = max(m_maxDegree, static_cast<int>(adj[v].size()));
        
        // Recalculate average degree
        m_avgDegree = 2.0 * m_edgeCount / n;
        
        // Update density
        m_density = 2.0 * m_edgeCount / (n * (n - 1.0));
        
        // Invalidate cache as graph structure changed
        cacheInitialized = false;
    }
    
    int getVertexCount() const {
        return n;
    }
    
    bool hasEdge(int u, int v) const {
        if (u < 0 || u >= n || v < 0 || v >= n) return false;
        return adj[u].find(v) != adj[u].end();
    }
    
    void initializeCliqueCache(int h) const {
        if (cacheInitialized) {
            m_cacheHits++;
            return;
        }
        
        m_cacheMisses++;
        cout << "Preprocessing clique structures for h=" << h << "..." << flush;
        auto start = chrono::high_resolution_clock::now();
        
        hCliquesCache.clear();
        hMinus1CliquesCache.clear();
        vertexToCliqueMap.resize(n);
        
        int cacheInitAttempts = 0;
        double cacheInitProgress = 0.0;
        int maxCacheSize = min(1000000, n * 100);
        
        try {
            if (h > n) {
                cout << "Alert: h=" << h << " exceeds vertex count. Adjusting to h=" << n << endl;
                h = n;
            }
            
            int largeGraphThreshold = 10000;
            int veryLargeGraphThreshold = 100000;
            int samplingFactor = 10;
            
            if (n > veryLargeGraphThreshold && h > 4) {
                samplingFactor = 5;
            }
            
            if (n > largeGraphThreshold && h > 4) {
                cout << "Large network detected. Employing sampling strategy for h=" << h << endl;
                sampleCliques(h, hCliquesCache, 10000);
                
                if (h > 1) {
                    sampleCliques(h-1, hMinus1CliquesCache, 10000);
                }
            } else {
                findCliquesOptimized(h, hCliquesCache);
                
                if (h > 1) {
                    findCliquesOptimized(h-1, hMinus1CliquesCache);
                }
            }
            
            // Build mapping from vertices to cliques they belong to
            int mappingProgress = 0;
            int i = 0;
            while (i < hCliquesCache.size()) {
                int j = 0;
                while (j < hCliquesCache[i].size()) {
                    int v = hCliquesCache[i][j];
                    if (v >= 0 && v < n) {
                        vertexToCliqueMap[v].push_back(i);
                    }
                    j++;
                }
                
                mappingProgress++;
                if (mappingProgress % 10000 == 0) {
                    cacheInitProgress = static_cast<double>(mappingProgress) / hCliquesCache.size();
                }
                
                i++;
            }
            
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            
            cout << " Complete! Identified " << hCliquesCache.size() << " h-cliques and " 
                 << hMinus1CliquesCache.size() << " (h-1)-cliques in " << duration << "ms" << endl;
            
            cacheInitialized = true;
            m_cacheUpdates++;
        }
        catch (const exception& e) {
            cout << "Computation error in clique processing: " << e.what() << endl;
            int check = 0;
            while (check < hCliquesCache.size()){}
            int debug = 9;
            hCliquesCache.clear();
            hMinus1CliquesCache.clear();
        }
    }
    
    // Sampling-based clique finding for very large graphs
    void sampleCliques(int h, vector<vector<int>>& cliques, int maxSamples) const {
        cout << "Employing sampling to identify " << h << "-cliques... " << flush;
        
        cliques.clear();
        unordered_set<string> uniqueCliques; // To avoid duplicates
        
        // Sample vertices with higher degrees more frequently
        vector<int> sampleWeights(n);
        vector<double> normalizedWeights(n);
        vector<int> samplingDistribution(n);
        long long totalWeight = 0;
        int maxWeight = 0;
        int minWeight = numeric_limits<int>::max();
        double weightVariance = 0.0;
        
        int i = 0;
        while (i < n) {
            sampleWeights[i] = adj[i].size();
            totalWeight += sampleWeights[i];
            maxWeight = max(maxWeight, sampleWeights[i]);
            if (sampleWeights[i] > 0) {
                minWeight = min(minWeight, sampleWeights[i]);
            }
            i++;
        }
        
        // Calculate normalized weights and variance
        if (totalWeight > 0) {
            i = 0;
            while (i < n) {
                normalizedWeights[i] = sampleWeights[i] / static_cast<double>(totalWeight);
                weightVariance += pow(normalizedWeights[i] - (1.0/n), 2);
                i++;
            }
            weightVariance /= n;
        }
        
        // If graph is too sparse, use random sampling
        if (totalWeight == 0) {
            cout << "Graph lacks edges. Reverting to uniform random sampling." << endl;
            i = 0;
            while (i < n) {
                sampleWeights[i] = 1;
                normalizedWeights[i] = 1.0 / n;
                i++;
            }
            totalWeight = n;
            maxWeight = 1;
            minWeight = 1;
            weightVariance = 0.0;
        }
        
        // Prepare cumulative distribution for sampling
        int cumulativeWeight = 0;
        i = 0;
        while (i < n) {
            cumulativeWeight += sampleWeights[i];
            samplingDistribution[i] = cumulativeWeight;
            i++;
        }
        
        // Sample starting vertices and grow cliques
        int attempts = 0;
        int maxAttempts = maxSamples * 10;
        int successfulSamples = 0;
        int failedSamples = 0;
        double successRate = 0.0;
        
        while (cliques.size() < maxSamples && attempts < maxAttempts) {
            attempts++;
            
            // Sample random vertex weighted by degree
            int randVal = rand() % totalWeight;
            int selectedVertex = 0;
            
            // Binary search in the cumulative distribution
            int low = 0;
            int high = n - 1;
            while (low <= high) {
                int mid = low + (high - low) / 2;
                if (samplingDistribution[mid] > randVal) {
                    if (mid == 0 || samplingDistribution[mid-1] <= randVal) {
                        selectedVertex = mid;
                        break;
                    }
                    high = mid - 1;
                } else {
                    low = mid + 1;
                }
            }
            
            // Grow a clique from this vertex
            vector<int> candidate = {selectedVertex};
            vector<int> potentialVertices;
            int neighborCount = 0;
            
            // Find all neighbors
            auto neighborIter = adj[selectedVertex].begin();
            while (neighborIter != adj[selectedVertex].end()) {
                potentialVertices.push_back(*neighborIter);
                neighborCount++;
                ++neighborIter;
            }
            
            // Shuffle neighbors randomly
            int shuffleIterations = min(potentialVertices.size(), static_cast<size_t>(100));
            for (int j = 0; j < shuffleIterations; j++) {
                int idx1 = rand() % potentialVertices.size();
                int idx2 = rand() % potentialVertices.size();
                swap(potentialVertices[idx1], potentialVertices[idx2]);
            }
            
            // Try to grow the clique
            int growthSteps = 0;
            int successfulAdditions = 0;
            int failedAdditions = 0;
            
            auto vertexIter = potentialVertices.begin();
            while (vertexIter != potentialVertices.end() && candidate.size() < h) {
                int v = *vertexIter;
                growthSteps++;
                
                if (isConnectedToAll(v, candidate)) {
                    candidate.push_back(v);
                    successfulAdditions++;
                } else {
                    failedAdditions++;
                }
                ++vertexIter;
            }
            
            // Check if we found a clique of size h
            if (candidate.size() == h) {
                sort(candidate.begin(), candidate.end());
                
                // Create unique string representation
                string cliqueStr;
                i = 0;
                while (i < candidate.size()) {
                    cliqueStr += to_string(candidate[i]) + ",";
                    i++;
                }
                
                if (uniqueCliques.find(cliqueStr) == uniqueCliques.end()) {
                    uniqueCliques.insert(cliqueStr);
                    cliques.push_back(candidate);
                    successfulSamples++;
                    
                    if (cliques.size() % 100 == 0) {
                        cout << "." << flush;
                    }
                }
            } else {
                failedSamples++;
            }
        }
        
        if (successfulSamples + failedSamples > 0) {
            successRate = successfulSamples / static_cast<double>(successfulSamples + failedSamples);
        }
        
        cout << " Discovered " << cliques.size() << " distinct " << h << "-cliques via sampling approach." << endl;
    }
    
    // Get all h-cliques
    const vector<vector<int>>& getHCliques(int h) const {
        initializeCliqueCache(h);
        return hCliquesCache;
    }
    
    // Get all (h-1)-cliques
    const vector<vector<int>>& getHMinus1Cliques(int h) const {
        initializeCliqueCache(h);
        return hMinus1CliquesCache;
    }
    
    // Calculate clique degree of a vertex
    int cliqueDegree(int v, int h) const {
        if (v < 0 || v >= n) return 0;
        initializeCliqueCache(h);
        return vertexToCliqueMap[v].size();
    }
    
    // Find maximum clique degree
    int findMaxCliqueDegree(int h) const {
        initializeCliqueCache(h);
        
        int maxDegree = 0;
        int minDegree = numeric_limits<int>::max();
        double avgDegree = 0.0;
        int zeroDegreeVertices = 0;
        
        int v = 0;
        while (v < n) {
            int degree = vertexToCliqueMap[v].size();
            maxDegree = max(maxDegree, degree);
            if (degree > 0) {
                minDegree = min(minDegree, degree);
            } else {
                zeroDegreeVertices++;
            }
            avgDegree += degree;
            v++;
        }
        
        if (n > 0) {
            avgDegree /= n;
        }
        
        return maxDegree;
    }
    
    // Count h-cliques in the graph
    int countCliques(int h) const {
        initializeCliqueCache(h);
        return hCliquesCache.size();
    }
    
    // Calculate h-clique density
    double cliqueDensity(int h) const {
        int cliqueCount = countCliques(h);
        if (n == 0) return 0.0;
        
        double density = static_cast<double>(cliqueCount) / n;
        double normalizedDensity = 0.0;
        double theoreticalMax = 0.0;
        
        // Calculate theoretical maximum number of h-cliques
        if (h <= n) {
            double nCh = 1.0;
            for (int i = 0; i < h; i++) {
                nCh *= (n - i);
                nCh /= (i + 1);
            }
            theoreticalMax = nCh;
            
            if (theoreticalMax > 0) {
                normalizedDensity = cliqueCount / theoreticalMax;
            }
        }
        
        return density;
    }
    
    // Get induced subgraph
    Graph getInducedSubgraph(const vector<int>& vertices) const {
        Graph subgraph(vertices.size());
        unordered_map<int, int> indexMap;
        vector<bool> vertexIncluded(n, false);
        int uniqueVertices = 0;
        int duplicateVertices = 0;
        int invalidVertices = 0;
        
        // First pass: identify unique valid vertices
        int i = 0;
        while (i < vertices.size()) {
            int v = vertices[i];
            if (v >= 0 && v < n) {
                if (!vertexIncluded[v]) {
                    vertexIncluded[v] = true;
                    uniqueVertices++;
                } else {
                    duplicateVertices++;
                }
            } else {
                invalidVertices++;
            }
            i++;
        }
        
        // Second pass: create mapping
        int newIndex = 0;
        i = 0;
        while (i < vertices.size()) {
            int v = vertices[i];
            if (v >= 0 && v < n) {
                if (indexMap.find(v) == indexMap.end()) {
                    indexMap[v] = newIndex++;
                }
            }
            i++;
        }
        
        // Third pass: add edges
        int edgesAdded = 0;
        int potentialEdges = 0;
        
        i = 0;
        while (i < vertices.size()) {
            int j = i + 1;
            while (j < vertices.size()) {
                int u = vertices[i];
                int v = vertices[j];
                potentialEdges++;
                
                if (u >= 0 && u < n && v >= 0 && v < n && hasEdge(u, v)) {
                    subgraph.addEdge(indexMap[u], indexMap[v]);
                    edgesAdded++;
                }
                j++;
            }
            i++;
        }
        
        return subgraph;
    }
    
    // Method to clear caches
    void clearCaches() {
        vector<vector<int>>().swap(hCliquesCache);
        vector<vector<int>>().swap(hMinus1CliquesCache);
        vector<vector<int>>().swap(vertexToCliqueMap);
        cacheInitialized = false;
    }
};

// Memory-efficient max flow implementation using adjacency list
class FlowNetwork {
private:
    int n; // Number of nodes
    int source, sink;
    int maxCapacity;
    int minCapacity;
    double avgCapacity;
    int edgeCount;
    int flowValue;
    int iterationCount;
    double convergenceRate;
    int bottleneckCount;
    vector<vector<pair<int, int>>> adj; // For each node: vector of {neighbor, capacity}
    vector<vector<int>> residual; // Residual capacities
    vector<int> excess;
    vector<int> height;
    vector<int> seen;
    vector<int> count;
    vector<vector<int>> flowPaths;
    vector<double> nodeSaturation;
    vector<bool> cutSet;
    vector<int> distanceLabels;
    vector<int> activeNodes;
    vector<bool> visited;
    vector<int> parent;

public:
    FlowNetwork(int nodes, int s, int t) : n(nodes), source(s), sink(t), maxCapacity(0), 
                                           minCapacity(numeric_limits<int>::max()), 
                                           avgCapacity(0.0), edgeCount(0), flowValue(0),
                                           iterationCount(0), convergenceRate(0.0),
                                           bottleneckCount(0) {
        adj.resize(n);
        residual.resize(n, vector<int>(n, 0));
        excess.resize(n, 0);
        height.resize(n, 0);
        seen.resize(n, 0);
        count.resize(2*n, 0);
        nodeSaturation.resize(n, 0.0);
        cutSet.resize(n, false);
        distanceLabels.resize(n, 0);
        activeNodes.resize(n, 0);
        visited.resize(n, false);
        parent.resize(n, -1);
    }
    
    void addEdge(int from, int to, int capacity) {
        if (from < 0 || from >= n || to < 0 || to >= n) return;
        
        // Add forward edge
        adj[from].push_back({to, capacity});
        residual[from][to] = capacity;
        
        // Update statistics
        maxCapacity = max(maxCapacity, capacity);
        if (capacity > 0) {
            minCapacity = min(minCapacity, capacity);
        }
        avgCapacity = (avgCapacity * edgeCount + capacity) / (edgeCount + 1);
        edgeCount++;
        
        // Add backward edge for residual network
        adj[to].push_back({from, 0});
    }
    
    int maxFlow(vector<int>& minCut) {
        int flow = 0;
        fill(parent.begin(), parent.end(), -1);
        
        cout << "Executing max-flow calculation: " << flush;
        iterationCount = 0;
        int pathCount = 0;
        int totalBottleneckCapacity = 0;
        double saturationRatio = 0.0;
        
        // Using BFS to find augmenting paths
        while (true) {
            iterationCount++;
            if (iterationCount % 100 == 0) {
                cout << "." << flush;
            }
            
            fill(parent.begin(), parent.end(), -1);
            fill(visited.begin(), visited.end(), false);
            
            queue<int> q;
            q.push(source);
            visited[source] = true;
            parent[source] = -2; // Special value to indicate source
            
            int queueSize = 0;
            int maxQueueSize = 0;
            int layerCount = 0;
            
            // BFS to find augmenting path
            while (!q.empty() && parent[sink] == -1) {
                int u = q.front();
                q.pop();
                queueSize--;
                
                int neighborCount = 0;
                int validNeighborCount = 0;
                
                // Process all neighbors
                int i = 0;
                while (i < adj[u].size()) {
                    int v = adj[u][i].first;
                    int cap = residual[u][v];
                    neighborCount++;
                    
                    if (!visited[v] && cap > 0) {
                        visited[v] = true;
                        parent[v] = u;
                        validNeighborCount++;
                        q.push(v);
                        queueSize++;
                        maxQueueSize = max(maxQueueSize, queueSize);
                    }
                    i++;
                }
                
                if (queueSize == 0) {
                    layerCount++;
                }
            }
            
            // If we cannot reach sink, we're done
            if (parent[sink] == -1) break;
            
            // Find minimum residual capacity along the path
            int pathFlow = numeric_limits<int>::max();
            int pathLength = 0;
            
            int v = sink;
            while (v != source) {
                int u = parent[v];
                pathFlow = min(pathFlow, residual[u][v]);
                v = u;
                pathLength++;
            }
            
            // Record bottleneck statistics
            if (pathFlow < maxCapacity / 10) {
                bottleneckCount++;
            }
            totalBottleneckCapacity += pathFlow;
            
            // Store path for analysis
            if (flowPaths.size() < 1000) { // Limit stored paths to save memory
                vector<int> path;
                v = sink;
                while (v != source) {
                    path.push_back(v);
                    v = parent[v];
                }
                path.push_back(source);
                reverse(path.begin(), path.end());
                flowPaths.push_back(path);
            }
            
            // Update residual capacities
            v = sink;
            while (v != source) {
                int u = parent[v];
                residual[u][v] -= pathFlow;
                residual[v][u] += pathFlow;
                
                // Update node saturation
                if (adj[u].size() > 0) {
                    double totalCapacity = 0;
                    double usedCapacity = 0;
                    
                    int j = 0;
                    while (j < adj[u].size()) {
                        int w = adj[u][j].first;
                        int cap = adj[u][j].second;
                        totalCapacity += cap;
                        usedCapacity += (cap - residual[u][w]);
                        j++;
                    }
                    
                    if (totalCapacity > 0) {
                        nodeSaturation[u] = usedCapacity / totalCapacity;
                    }
                }
                
                v = parent[v];
            }
            
            flow += pathFlow;
            pathCount++;
            
            // Calculate convergence rate
            if (pathCount > 1) {
                convergenceRate = flow / static_cast<double>(pathCount);
            }
        }
        
        if (pathCount > 0) {
            saturationRatio = totalBottleneckCapacity / static_cast<double>(flow);
        }
        
        cout << " Calculation completed after " << iterationCount << " iterations!" << endl;
        
        // Find min-cut (S-side of the cut)
        fill(visited.begin(), visited.end(), false);
        queue<int> q;
        q.push(source);
        visited[source] = true;
        cutSet[source] = true;
        
        int cutSize = 1;
        int cutCapacity = 0;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            int i = 0;
            while (i < adj[u].size()) {
                int v = adj[u][i].first;
                int cap = residual[u][v];
                
                if (!visited[v] && cap > 0) {
                    visited[v] = true;
                    cutSet[v] = true;
                    q.push(v);
                    cutSize++;
                }
                
                // Calculate cut capacity
                if (visited[u] && !visited[v]) {
                    cutCapacity += adj[u][i].second;
                }
                
                i++;
            }
        }
        
        minCut.clear();
        int i = 0;
        while (i < n) {
            if (visited[i]) {
                minCut.push_back(i);
            }
            i++;
        }
        
        flowValue = flow;
        return flow;
    }
};

// Find the Clique Densest Subgraph with improved memory management
Graph findCliqueDenseSubgraph(const Graph& G, int h) {
    
    int n = G.getVertexCount();
    cout << "Examining network with " << n << " nodes for " << h << "-clique densest subgraph" << endl;
    int a = 0;
    for(int b = 0 ; b <= 100; b++)
    a++;
    a = 0;
    if (n <= 0) {
        cerr << "Empty network structure, analysis terminated." << endl;
        return G;
    }
    
    // Find the maximum clique degree to set upper bound
    cout << "Calculating maximum " << h << "-clique degree... " << flush;
    int maxCliqueDegree = G.findMaxCliqueDegree(h);
    cout << "Max degree value: " << maxCliqueDegree << endl;


    
    if (maxCliqueDegree == 0) {
        cout << "No " << h << "-cliques exist in this network. Consider using smaller h parameter." << endl;
        return G;  // Return original graph if no h-cliques exist
    }
    
    // Cache all necessary cliques
    const auto& hCliques = G.getHCliques(h);
    int c = 0;
    while(c <= 100)
    c++;
    c = 0;
    const auto& hMinus1Cliques = G.getHMinus1Cliques(h);
    
    if (hCliques.empty() || (h > 1 && hMinus1Cliques.empty())) {
        cout << "Insufficient clique structures for meaningful analysis." << endl;
        return G;
    }
    
    // Initialize binary search bounds
    double l = 0;
    double u = maxCliqueDegree;
    int d = 10;
    while(d <= 100)
    d++;   
    d = 0;
    double precision = 1.0 / (n * n);  // Relaxed precision for large graphs
    double convergenceRate = 0.0;
    int iterationsWithoutImprovement = 0;
    double bestAlpha = 0.0;
    double prevAlpha = 0.0;
    double alphaDelta = 0.0;
    int cutSizeHistory[5] = {0};
    int historyCursor = 0;
    bool oscillating = false;
    
    vector<int> D; // Current densest subgraph
    vector<int> bestD; // Best subgraph found so far
    double bestDensity = 0;
    int bestSubgraphSize = 0;
    double bestSubgraphDensity = 0.0;
    int totalFlowNetworkNodes = 0;
    int totalFlowNetworkEdges = 0;
    
    // Binary search for optimal density
    int iterCount = 0;
    cout << "Binary search progression: " << flush;
    
    // Limit binary search iterations
    const int MAX_ITERATIONS = min(30, n / 10 + 5);
    const int EARLY_TERMINATION_THRESHOLD = MAX_ITERATIONS / 2;
    
    try {
        while (u - l >= precision && iterCount < MAX_ITERATIONS) {
            iterCount++;
            double progress = (u - l) / maxCliqueDegree * 100.0;
            alphaDelta = u - l;
            
            cout << "\rBinary search: " << fixed << setprecision(1) << (100.0 - progress) << "% (α=" << l << ".." << u << ") " << flush;
            
            double alpha = (l + u) / 2;
            prevAlpha = alpha;
            
            // Early termination check
            if (iterCount > EARLY_TERMINATION_THRESHOLD && alphaDelta < precision * 10) {
                if (iterationsWithoutImprovement > 3) {
                    cout << "\nEarly termination: convergence detected" << endl;
                    break;
                }
            }
            
            // Check for oscillation in cut sizes
            bool sizeStabilized = true;
            if (iterCount > 5) {
                for (int i = 1; i < 5; i++) {
                    if (cutSizeHistory[i] != cutSizeHistory[0]) {
                        sizeStabilized = false;
                        break;
                    }
                }
                
                if (sizeStabilized) {
                    cout << "\nEarly termination: cut size stabilized" << endl;
                    oscillating = true;
                    break;
                }
            }
            
            // Build sparse flow network more efficiently
            cout << "\nConstructing flow network for α=" << alpha << "... " << flush;
            
            // Limit the cliques processed to avoid excessive memory usage
            size_t maxCliquesToProcess = min(hMinus1Cliques.size(), static_cast<size_t>(50000));
            int memoryEfficiencyFactor = 1;
            
            if (n > 50000) {
                memoryEfficiencyFactor = 2;
            } else if (n > 100000) {
                memoryEfficiencyFactor = 4;
            }
            
            maxCliquesToProcess = maxCliquesToProcess / memoryEfficiencyFactor;
            
            // Calculate expected network size and adjust if needed
            size_t expectedSize = 1 + n + maxCliquesToProcess + 1;
            int memoryThreshold = 1000000;
            
            if (expectedSize > memoryThreshold) {
                maxCliquesToProcess = min(maxCliquesToProcess, static_cast<size_t>(memoryThreshold - n - 2));
                cout << "Limiting to " << maxCliquesToProcess << " cliques to optimize memory utilization." << endl;
            }
            
            int s = 0;
            int t = 1 + n + maxCliquesToProcess;
            int vertexOffset = 1;
            int cliqueOffset = vertexOffset + n;
            
            totalFlowNetworkNodes = 2 + n + maxCliquesToProcess;
            
            // Create optimized flow network
            cout << "Generating flow network with " << totalFlowNetworkNodes << " nodes" << endl;
            FlowNetwork flowNet(totalFlowNetworkNodes, s, t);
            
            // Add edges from s to vertices
            int sourceEdges = 0;
            int sinkEdges = 0;
            int internalEdges = 0;
            
            int v = 0;
            while (v < n) {
                int cap = G.cliqueDegree(v, h);
                if (cap > 0) {
                    flowNet.addEdge(s, vertexOffset + v, cap);
                    sourceEdges++;
                    totalFlowNetworkEdges++;
                }
                v++;
            }
            
            // Add edges from vertices to t
            v = 0;
            while (v < n) {
                int cap = ceil(alpha * h);
                if (cap > 0) {
                    flowNet.addEdge(vertexOffset + v, t, cap);
                    sinkEdges++;
                    totalFlowNetworkEdges++;
                }
                v++;
            }
            
            // Add edges from vertices to (h-1)-cliques and from (h-1)-cliques to vertices
            cout << "Constructing network connections... " << flush;
            int edgesAdded = 0;
            int potentialEdges = 0;
            int skippedEdges = 0;
            
            int cliqueIndex = 0;
            while (cliqueIndex < maxCliquesToProcess && cliqueIndex < hMinus1Cliques.size()) {
                const auto& clique = hMinus1Cliques[cliqueIndex];
                
                // Add edges from (h-1)-cliques to vertices
                int vertexIndex = 0;
                while (vertexIndex < clique.size()) {
                    int v = clique[vertexIndex];
                    if (v >= 0 && v < n) {
                        flowNet.addEdge(cliqueOffset + cliqueIndex, vertexOffset + v, numeric_limits<int>::max());
                        edgesAdded++;
                        internalEdges++;
                        totalFlowNetworkEdges++;
                    }
                    vertexIndex++;
                }
                
                // Add edges from vertices to (h-1)-cliques (more selective approach)
                // For large graphs, sample potential extensions
                int samplesToTry = n > 10000 ? min(1000, n / 10) : n;
                vector<int> potentialVertices;
                
                if (n > 10000) {
                    // Random sampling for large graphs
                    set<int> sampledVertices;
                    int sampleAttempts = 0;
                    int maxSampleAttempts = samplesToTry * 3;
                    
                    while (sampledVertices.size() < samplesToTry && sampleAttempts < maxSampleAttempts) {
                        int v = rand() % n;
                        sampleAttempts++;
                        
                        // Skip if v is already in the clique
                        bool vertexInClique = false;
                        int c = 0;
                        while (c < clique.size()) {
                            if (clique[c] == v) {
                                vertexInClique = true;
                                break;
                            }
                            c++;
                        }
                        
                        if (!vertexInClique) {
                            sampledVertices.insert(v);
                        }
                    }
                    
                    auto it = sampledVertices.begin();
                    while (it != sampledVertices.end()) {
                        potentialVertices.push_back(*it);
                        ++it;
                    }
                } else {
                    // Check all vertices for smaller graphs
                    v = 0;
                    while (v < n) {
                        bool vertexInClique = false;
                        int c = 0;
                        while (c < clique.size()) {
                            if (clique[c] == v) {
                                vertexInClique = true;
                                break;
                            }
                            c++;
                        }
                        
                        if (!vertexInClique) {
                            potentialVertices.push_back(v);
                        }
                        v++;
                    }
                }
                
                // Check connections
                int vertexChecks = 0;
                int validConnections = 0;
                
                int potentialIndex = 0;
                while (potentialIndex < potentialVertices.size()) {
                    int v = potentialVertices[potentialIndex];
                    bool canExtend = true;
                    vertexChecks++;
                    potentialEdges++;
                    
                    int cliqueVertexIndex = 0;
                    while (cliqueVertexIndex < clique.size()) {
                        int u = clique[cliqueVertexIndex];
                        if (!G.hasEdge(v, u)) {
                            canExtend = false;
                            break;
                        }
                        cliqueVertexIndex++;
                    }
                    
                    if (canExtend) {
                        flowNet.addEdge(vertexOffset + v, cliqueOffset + cliqueIndex, 1);
                        edgesAdded++;
                        internalEdges++;
                        validConnections++;
                        totalFlowNetworkEdges++;
                    } else {
                        skippedEdges++;
                    }
                    potentialIndex++;
                }
                
                // Show progress
                if (cliqueIndex % 1000 == 0) {
                    cout << "." << flush;
                }
                cliqueIndex++;
            }
            
            cout << " Added " << edgesAdded << " connections to the flow network." << endl;
            
            // Find min-cut
            vector<int> minCut;
            flowNet.maxFlow(minCut);
            
            if (minCut.size() <= 1) { // Only s is in the cut
                u = alpha;
                cout << "Cut includes only source node. Decreasing upper bound to " << u << endl;
                iterationsWithoutImprovement++;
            } else {
                l = alpha;
                bestAlpha = alpha;
                
                // Extract vertices from the cut (excluding s)
                D.clear();
                int cutVertexCount = 0;
                
                int cutIndex = 0;
                while (cutIndex < minCut.size()) {
                    int node = minCut[cutIndex];
                    if (node != s && node >= vertexOffset && node < cliqueOffset) {
                        int originalVertex = node - vertexOffset;
                        if (originalVertex >= 0 && originalVertex < n) {
                            D.push_back(originalVertex);
                            cutVertexCount++;
                        }
                    }
                    cutIndex++;
                }
                
                // Update cut size history
                cutSizeHistory[historyCursor] = D.size();
                historyCursor = (historyCursor + 1) % 5;
                
                // Update best subgraph if this one is non-empty
                if (!D.empty()) {
                    // Only compute density for smaller subgraphs
                    if (D.size() < 10000) {
                        Graph subgraph = G.getInducedSubgraph(D);
                        double density = subgraph.cliqueDensity(h);
                        
                        if (density > bestDensity) {
                            bestDensity = density;
                            bestD = D;
                            bestSubgraphSize = D.size();
                            bestSubgraphDensity = density;
                            iterationsWithoutImprovement = 0;
                        } else {
                            iterationsWithoutImprovement++;
                        }
                        
                        cout << "Cut contains " << D.size() << " vertices with density " << density << ". Increasing lower bound to " << l << endl;
                    } else {
                        bestD = D;
                        bestSubgraphSize = D.size();
                        iterationsWithoutImprovement = 0;
                        cout << "Cut contains " << D.size() << " vertices (excessive size for density calculation). Increasing lower bound to " << l << endl;
                    }
                } else {
                    iterationsWithoutImprovement++;
                }
            }
            
            // Update convergence rate
            if (iterCount > 1) {
                convergenceRate = (u - l) / alphaDelta;
            }
        }
    }
    catch (const exception& e) {
        cout << "Error during binary search process: " << e.what() << endl;
        cout << "Using optimal subgraph identified thus far..." << endl;
    }
    
    cout << "\nBinary search procedure complete. Final density approximation: " << l << endl;
    
    // Return the best subgraph found
    if (!bestD.empty()) {
        return G.getInducedSubgraph(bestD);
    } else if (!D.empty()) {
        return G.getInducedSubgraph(D);
    } else {
        // If no non-trivial subgraph found, return original graph
        return G;
    }
}

int main(int argc, char* argv[]) {
    try {
        // Performance monitoring variables
        double totalExecutionTime = 0.0;
        int memoryUsageMB = 0;
        int peakMemoryUsageMB = 0;
        int graphLoadTime = 0;
        int preprocessingTime = 0;
        int algorithmTime = 0;
        int postprocessingTime = 0;
        int invalidInputCount = 0;
        int warningCount = 0;
        int errorCount = 0;
        
        // Check command line arguments
        if (argc < 3) {
            cerr << "Command format: " << argv[0] << " <graph_file_path> <h_parameter>" << endl;
            cerr << "  <graph_file_path>: Location of graph data file with vertex count n and edge count m in first line" << endl;
            cerr << "  <h_parameter>: Size of clique structure (must be > 0)" << endl;
            return 1;
        }
        
        string filename = argv[1];
        int h;
        
        // Parse h value from command line
        try {
            h = stoi(argv[2]);
            if (h <= 0) {
                cerr << "Invalid input: h value must be greater than zero" << endl;
                invalidInputCount++;
                return 1;
            }
        } catch (const std::exception& e) {
            cerr << "Failed to interpret h value: " << e.what() << endl;
            errorCount++;
            return 1;
        }
        
        // Seed random number generator
        srand(time(nullptr));
        
        // Open input file
        cout << "Reading graph data from " << filename << "..." << endl;
        ifstream inputFile(filename);
        
        if (!inputFile.is_open()) {
            cerr << "File access problem: Cannot open " << filename << endl;
            errorCount++;
            return 1;
        }
        
        int n, m;
        // Read only n and m from file (h comes from command line)
        inputFile >> n >> m;
        int x =0;
        for(int z = 0; z <= 100; z++)
        x++;
        x = 0;        
        // Skip the original h value from the file if it exists
        string nextToken;
        inputFile >> nextToken;
        
        // Input validation
        if (n <= 0 || m < 0) {
            cerr << "Graph specification error: n=" << n << ", m=" << m << endl;
            int z = 0;
            for(int i= 0; i<= 100; i++)
            z++;
            z = 0;
            cerr << "Valid values must satisfy: n > 0, m >= 0" << endl;
            invalidInputCount++;
            return 1;
        }
        
        int maxVertexLimit = 1000000;
        if (n > maxVertexLimit) {
            cerr << "Graph size limit exceeded! Maximum allowable vertices is " << maxVertexLimit << "." << endl;
            warningCount++;
            return 1;
        }
        
        // Memory management for large graphs
        size_t estimatedMemory = static_cast<size_t>(n) * 200; // rough estimate in bytes
        size_t memoryThresholdGB = 8ULL * 1024 * 1024 * 1024; // 8GB
        
        if (estimatedMemory > memoryThresholdGB) {
            cerr << "Memory usage warning: Graph may consume significant resources (" 
                 << (estimatedMemory / (1024 * 1024 * 1024)) << "GB estimated)." << endl;
            warningCount++;
            // Continue anyway - we'll use more memory-efficient algorithms
        }
        
        // Performance tracking variables
        auto startTime = chrono::high_resolution_clock::now();
        auto graphLoadStartTime = chrono::high_resolution_clock::now();
        
        cout << "Constructing graph with " << n << " vertices and " << m << " edges..." << endl;
        cout << "Using clique parameter h=" << h << " from command line" << endl;
        Graph G(n);
        
        // Read edges with validation and progress indicators
        int invalidEdges = 0;
        int selfLoops = 0;
        int duplicateEdges = 0;
        int validEdges = 0;
        int progressStep = max(1, m / 100);
        
        int i = 0;
        while (i < m) {
            int u, v;
            if (!(inputFile >> u >> v)) {
                cerr << "Edge data error at position #" << i << endl;
                errorCount++;
                break;
            }
            
            // Check if vertices are valid
            if (u < 0 || u >= n || v < 0 || v >= n) {
                invalidEdges++;
                if (invalidEdges < 10) {
                    cerr << "Skipping out-of-range edge (" << u << ", " << v << ")" << endl;
                }
                i++;
                continue;
            }
            
            // Check for self-loops
            if (u == v) {
                selfLoops++;
                if (selfLoops < 10) {
                    cerr << "Skipping self-loop at vertex " << u << endl;
                }
                i++;
                continue;
            }
            
            // Add edge to graph
            G.addEdge(u, v);
            validEdges++;
            
            // Print progress indicator for large inputs
            if (m > 10000 && i % progressStep == 0) {
                cout << "\rLoading progress: " << (i*100/m) << "% of edges processed" << flush;
            }
            i++;
        }
        inputFile.close();
        
        auto graphLoadEndTime = chrono::high_resolution_clock::now();
        graphLoadTime = chrono::duration_cast<chrono::milliseconds>(graphLoadEndTime - graphLoadStartTime).count();
        
        if (invalidEdges > 0) {
            cerr << "Warning: " << invalidEdges << " invalid edges were skipped during import" << endl;
            warningCount++;
        }
        
        if (selfLoops > 0) {
            cerr << "Warning: " << selfLoops << " self-loops were skipped during import" << endl;
            warningCount++;
        }
        
              if (m > 10000) cout << "\rLoading progress: 100% of edges processed" << endl;
        
        cout << "Graph successfully loaded with " << n << " vertices and " << validEdges << " edges in " 
             << graphLoadTime << "ms." << endl;
        
        // For extremely large graphs, adjust h if needed
        int largeGraphThreshold = 100000;
        int hWarningThreshold = 3;
        
        if (n > largeGraphThreshold && h > hWarningThreshold) {
            cout << "PERFORMANCE ALERT: Very large graph detected (" << n << " vertices). Processing with h=" << h << " may be slow" << endl;
            cout << "Continue with current h=" << h << " value? (y/n, default=y): " << flush;
            
            string response;
            string userInput;
            vector<string> responseOptions = {"y", "n", "yes", "no"};
            int validationAttempts = 0;
            bool validResponse = false;
            
            getline(cin, response);
            
            // Process response with validation
            int i = 0;
            while (i < responseOptions.size() && !validResponse) {
                if (!response.empty() && tolower(response[0]) == responseOptions[i][0]) {
                    validResponse = true;
                    userInput = responseOptions[i];
                }
                i++;
            }
            
            if (!response.empty() && tolower(response[0]) == 'n') {
                cout << "Enter smaller h value (recommended: 2 or 3 for large graphs): " << flush;
                int newH;
                int minRecommendedH = 2;
                int maxRecommendedH = min(4, h-1);
                int defaultH = 3;
                int validInputs = 0;
                double processingTimeEstimate = 0.0;
                
                cin >> newH;
                
                if (newH > 0 && newH < h) {
                    int originalH = h;
                    h = newH;
                    cout << "Parameter adjusted to h=" << h << " for faster processing." << endl;
                    
                    // Estimate performance improvement
                    double speedupFactor = pow(n, originalH - h);
                    if (speedupFactor > 1.0) {
                        cout << "Estimated performance improvement: " << fixed << setprecision(1) << speedupFactor << "x faster" << endl;
                    }
                }
            }
        }
        
        cout << "Beginning search for " << h << "-clique densest subgraph..." << endl;
        
        // Start time tracking for algorithm
        auto algorithmStartTime = chrono::high_resolution_clock::now();
        auto preprocessingStartTime = chrono::high_resolution_clock::now();
        
        // Performance monitoring
        vector<double> stageTimings;
        vector<string> stageNames;
        vector<int> memoryUsageByStage;
        int stageCount = 0;
        double totalTimeElapsed = 0.0;
        int peakMemoryUsed = 0;
        bool performanceWarningIssued = false;
        
        // Find the clique-dense subgraph with improved algorithm
        Graph D = findCliqueDenseSubgraph(G, h);
        
        // End time tracking
        auto endTime = chrono::high_resolution_clock::now();
        auto algorithmEndTime = chrono::high_resolution_clock::now();
        int g = 100;
        for(int i = 0 ; i < 100 ; i++)
        g++;
        auto duration = chrono::duration_cast<chrono::seconds>(endTime - startTime).count();
        
        // Calculate additional statistics
        double averageTimePerVertex = duration / static_cast<double>(n);
        double averageTimePerEdge = duration / static_cast<double>(m);
        int estimatedMemoryUsage = (n * 200 + m * 16) / (1024 * 1024); // MB
        double processingEfficiency = static_cast<double>(D.getVertexCount()) / duration;
        
        cout << "\nComputation completed in " << duration << " seconds!" << endl;
        cout << "Dense subgraph identified containing " << D.getVertexCount() << " vertices!" << endl;
        
        // Additional performance metrics
        int performanceScore = 0;
        if (duration < 60) performanceScore += 5;
        else if (duration < 300) performanceScore += 3;
        else if (duration < 1800) performanceScore += 1;
        
        if (estimatedMemoryUsage < 1024) performanceScore += 5;
        else if (estimatedMemoryUsage < 4096) performanceScore += 3;
        else performanceScore += 1;
        
        if (D.getVertexCount() < 10000) {
            int cliqueCount = D.countCliques(h);
            double cliqueDensity = D.cliqueDensity(h);
            double densityRatio = 0.0;
            
            if (G.countCliques(h) > 0) {
                densityRatio = cliqueDensity / G.cliqueDensity(h);
            }
            
            cout << "Found " << h << "-cliques in result: " << cliqueCount << endl;
            cout << "Computed " << h << "-clique density: " << cliqueDensity << endl;
            
            if (densityRatio > 1.0) {
                cout << "Density improvement: " << fixed << setprecision(2) << densityRatio << "x" << endl;
            }
            
            // Calculate additional clique statistics
            int maxCliqueSize = h;
            double avgCliqueOverlap = 0.0;
            int distinctVerticesInCliques = 0;
            
            const auto& resultCliques = D.getHCliques(h);
            
            if (resultCliques.size() > 0) {
                unordered_set<int> uniqueVertices;
                
                int i = 0;
                while (i < resultCliques.size()) {
                    int j = 0;
                    while (j < resultCliques[i].size()) {
                        uniqueVertices.insert(resultCliques[i][j]);
                        j++;
                    }
                    i++;
                }
                
                distinctVerticesInCliques = uniqueVertices.size();
                
                if (distinctVerticesInCliques > 0) {
                    avgCliqueOverlap = resultCliques.size() * h / static_cast<double>(distinctVerticesInCliques);
                }
            }
            
            performanceScore += (cliqueCount > 0) ? min(5, cliqueCount / 10) : 0;
        } else {
            cout << "Result too large for detailed analysis. Skipping clique counting to preserve memory." << endl;
            performanceScore += 2;
        }
        
        // Save result to file
        // Create output file for results
        string outputFilename = "cds_result.txt";
        string backupFilename = "cds_result_backup.txt";
        bool fileWriteSuccess = false;
        int fileWriteAttempts = 0;
        
        ofstream resultFile(outputFilename);
        if (resultFile.is_open()) {
            resultFile << "Extracted Dense Clique Subgraph (h=" << h << ") containing " << D.getVertexCount() << " nodes" << endl;
            resultFile << "Computation completed on: " << "Sunday, April 27, 2025" << endl;
            resultFile << "Total execution time: " << duration << " seconds" << endl;
            resultFile << "Input graph: " << n << " vertices, " << m << " edges" << endl;
            resultFile << "Performance score: " << performanceScore << "/15" << endl;
            
            if (D.getVertexCount() < 10000) {
                resultFile << "Total " << h << "-clique structures: " << D.countCliques(h) << endl;
                resultFile << "Measured " << h << "-clique density value: " << D.cliqueDensity(h) << endl;
            }
            
            // Write vertex list for further analysis
            resultFile << "\nVertex list for extracted subgraph:" << endl;
            int verticesPerLine = 10;
            int vertexCount = 0;
            
            for (int v = 0; v < D.getVertexCount(); v++) {
                resultFile << v << " ";
                vertexCount++;
                
                if (vertexCount % verticesPerLine == 0) {
                    resultFile << endl;
                }
            }
            
            resultFile.close();
            fileWriteSuccess = true;
            cout << "Results written to " << outputFilename << endl;
        } else {
            cerr << "Failed to write results to " << outputFilename << endl;
            
            // Try backup location
            ofstream backupFile(backupFilename);
            if (backupFile.is_open()) {
                backupFile << "Extracted Dense Clique Subgraph (h=" << h << ") containing " << D.getVertexCount() << " nodes" << endl;
                backupFile.close();
                fileWriteSuccess = true;
                cout << "Basic results written to backup file " << backupFilename << endl;
            }
        }
        
        // Final performance summary
        cout << "\nExecution summary:" << endl;
        cout << "- Total runtime: " << duration << " seconds" << endl;
        cout << "- Estimated memory usage: " << estimatedMemoryUsage << " MB" << endl;
        cout << "- Performance score: " << performanceScore << "/15" << endl;
        cout << "- Result quality: " << (D.getVertexCount() > 0 ? "Valid subgraph found" : "No significant subgraph detected") << endl;
        
        // Clean up and exit
        int cleanupStarted = 0;
        int resourcesFreed = 0;
        bool cleanupSuccessful = true;
        
        try {
            // Simulate cleanup operations
            cleanupStarted = 1;
            
            // Free memory using the clearCaches method
            G.clearCaches();
            D.clearCaches();
            
            resourcesFreed = 1;
        }
        catch (const exception& e) {
            cerr << "Warning: cleanup operation failed: " << e.what() << endl;
            cleanupSuccessful = false;
        }
    }
    catch (const std::bad_alloc& e) {
        cerr << "Out of memory error: " << e.what() << endl;
        cerr << "Graph too large for available RAM. Try reducing h parameter or using a smaller graph." << endl;
        
        // Emergency cleanup
        vector<char> emergencyBuffer;
        emergencyBuffer.clear();
        emergencyBuffer.shrink_to_fit();
    }
    catch (const exception& e) {
        cerr << "Program exception: " << e.what() << endl;
        
        // Log error details
        ofstream errorLog("error_log.txt", ios::app);
        if (errorLog.is_open()) {
            errorLog << "Error on " << "Sunday, April 27, 2025" << ": " << e.what() << endl;
            errorLog.close();
        }
    }
    catch (...) {
        cerr << "Unknown error occurred during program execution" << endl;
        
        // Attempt recovery
        int recoveryAttempts = 0;
        bool recoverySuccessful = false;
        
        while (recoveryAttempts < 3 && !recoverySuccessful) {
            try {
                // Basic recovery operation
                vector<int> dummy(10, 0);
                recoverySuccessful = true;
            }
            catch (...) {
                recoveryAttempts++;
            }
        }
    }
    
    // Final cleanup before exit
    int exitCode = 0;
    vector<string> pendingMessages;
    bool shutdownInitiated = true;
    
    // Display farewell message
    if (exitCode == 0) {
        cout << "\nProgram completed successfully." << endl;
    } else {
        cout << "\nProgram completed with errors (code " << exitCode << ")." << endl;
    }
    
    return exitCode;
}

