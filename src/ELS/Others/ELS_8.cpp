#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
using namespace std;
using namespace std::chrono;

// Global aggregators
int maxCliqueSize = 0;
long long totalMaximalCliques = 0;
vector<int> cliqueSizeDistribution; // index = clique size; value = count

// Merge thread-local aggregator into globals (sequential version)
void mergeAggregators(int localMaxClique, long long localCliqueCount, const vector<int>& localDist) {
    maxCliqueSize = max(maxCliqueSize, localMaxClique);
    totalMaximalCliques += localCliqueCount;
    if(localDist.size() > cliqueSizeDistribution.size())
        cliqueSizeDistribution.resize(localDist.size(), 0);
    for (size_t i = 0; i < localDist.size(); ++i)
        cliqueSizeDistribution[i] += localDist[i];
}

// Add an edge into our graph (0-indexed vertices)
void addEdge(int u, int v, vector<vector<int>>& adj) {
    if(u >= (int)adj.size() || v >= (int)adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    // Since we eventually sort and deduplicate, simply push_back
    adj[u].push_back(v);
    adj[v].push_back(u);
}

// Two-pointer intersection of two sorted vectors
vector<int> intersectVectors(const vector<int>& A, const vector<int>& B) {
    vector<int> res;
    res.reserve(min(A.size(), B.size()));
    auto itA = A.begin(), itB = B.begin();
    while(itA != A.end() && itB != B.end()){
        if(*itA < *itB)
            ++itA;
        else if(*itB < *itA)
            ++itB;
        else {
            res.push_back(*itA);
            ++itA; ++itB;
        }
    }
    return res;
}

// Two-pointer set difference: A \ B; both sorted.
vector<int> diffVectors(const vector<int>& A, const vector<int>& B) {
    vector<int> res;
    res.reserve(A.size());
    auto itA = A.begin(), itB = B.begin();
    while(itA != A.end()){
        if(itB == B.end() || *itA < *itB){
            res.push_back(*itA);
            ++itA;
        } else if(*itA == *itB) {
            ++itA; ++itB;
        } else {
            ++itB;
        }
    }
    return res;
}

/**
 * Recursive Bron–KerboschPivot on sorted vectors.
 * P, R, X are sorted vectors of vertices.
 */
void BronKerboschPivot(vector<int> P, vector<int> R, vector<int> X,
                         const vector<vector<int>>& adj,
                         int &localMax, long long &localCount, vector<int>& localDist) {
    if(P.empty() && X.empty()){
        int cliqueSize = R.size();
        if(cliqueSize > 1){
            localMax = max(localMax, cliqueSize);
            localCount++;
            if(localDist.size() <= (size_t)cliqueSize)
                localDist.resize(cliqueSize+1, 0);
            localDist[cliqueSize]++;
        }
        return;
    }
    // Choose a pivot from P∪X arbitrarily (here, the first element)
    vector<int> unionPX = P;
    unionPX.insert(unionPX.end(), X.begin(), X.end());
    sort(unionPX.begin(), unionPX.end());
    int bestPivot = unionPX.front();
    size_t bestDegree = 0;
    // Use intersection count with P as pivot degree.
    for (int candidate : unionPX) {
        vector<int> common = intersectVectors(P, adj[candidate]);
        if(common.size() > bestDegree){
            bestDegree = common.size();
            bestPivot = candidate;
        }
    }
    // Compute P \ N(bestPivot)
    vector<int> diffP = diffVectors(P, adj[bestPivot]);
    // For each vertex in diffP, recurse with new sets.
    for (int v : diffP) {
        vector<int> newR = R;
        newR.push_back(v);
        vector<int> newP = intersectVectors(P, adj[v]);
        vector<int> newX = intersectVectors(X, adj[v]);
        BronKerboschPivot(newP, newR, newX, adj, localMax, localCount, localDist);
        // Remove v from P and add it to X (maintain sorted order)
        P.erase(remove(P.begin(), P.end(), v), P.end());
        X.push_back(v);
        sort(X.begin(), X.end());
    }
}

/* Compute degeneracy ordering by repeatedly removing a vertex with smallest degree */
vector<int> degeneracyOrder(const vector<vector<int>>& adj) {
    int n = adj.size();
    vector<int> d(n);
    for (int i = 0; i < n; i++) {
        d[i] = adj[i].size();
    }
    int maxDeg = *max_element(d.begin(), d.end());
    vector<vector<int>> bucket(maxDeg + 1);
    for (int i = 0; i < n; i++) {
        bucket[d[i]].push_back(i);
    }
    vector<bool> removed(n, false);
    vector<int> order;
    order.reserve(n);
    for (int k = 0; k < n; k++) {
        int currDeg = 0;
        while (currDeg <= maxDeg && bucket[currDeg].empty())
            currDeg++;
        if (currDeg > maxDeg)
            break;
        int u = bucket[currDeg].back();
        bucket[currDeg].pop_back();
        removed[u] = true;
        order.push_back(u);
        for (int v : adj[u]) {
            if (!removed[v]) {
                int oldDeg = d[v];
                auto &bkt = bucket[oldDeg];
                auto it = find(bkt.begin(), bkt.end(), v);
                if (it != bkt.end())
                    bkt.erase(it);
                d[v]--;
                bucket[d[v]].push_back(v);
            }
        }
    }
    reverse(order.begin(), order.end());
    return order;
}

// Outer BronKerbosch using degeneracy ordering (sequential version)
void BronKerboschDegeneracy(const vector<vector<int>>& adj) {
    int n = adj.size();
    vector<int> ordering = degeneracyOrder(adj);
    vector<int> pos(n);
    for (int i = 0; i < n; i++) {
        pos[ordering[i]] = i;
    }
    // Reset global aggregators.
    maxCliqueSize = 0;
    totalMaximalCliques = 0;
    cliqueSizeDistribution.clear();

    // Process each vertex sequentially.
    for (int i = 0; i < n; i++){
        cout << "[DEBUG] Processing node " << ordering[i]
             << " (loop index " << i << ")" << endl;
        int vi = ordering[i];
        vector<int> P, X, R;
        // Build P: neighbors with pos > pos[vi]
        for (int w : adj[vi]) {
            if (pos[w] > pos[vi])
                P.push_back(w);
        }
        // Build X: neighbors with pos < pos[vi]
        for (int w : adj[vi]) {
            if (pos[w] < pos[vi])
                X.push_back(w);
        }
        sort(P.begin(), P.end());
        sort(X.begin(), X.end());
        R.push_back(vi);

        int localMaxClique = 0;
        long long localCliqueCount = 0;
        vector<int> localDist;
        BronKerboschPivot(P, R, X, adj, localMaxClique, localCliqueCount, localDist);
        mergeAggregators(localMaxClique, localCliqueCount, localDist);
    }
}

int main(int argc, char* argv[]) {
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " [input file]" << endl;
        return 1;
    }
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    if(!infile){
        cerr << "Failed to open input file." << endl;
        return 1;
    }

    string line;
    vector<pair<int,int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            edgeList.push_back({u,v});
            maxVertex = max(maxVertex, max(u,v));
        }
    }
    infile.close();

    // Build and optimize adjacency list.
    vector<vector<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    // Sort and deduplicate neighbor lists.
    for(auto &neighbors : adj){
        sort(neighbors.begin(), neighbors.end());
        neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    auto start = high_resolution_clock::now();
    BronKerboschDegeneracy(adj);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start).count();

    outfile << "Largest size of the clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << "\n";
    outfile << "Execution time (ms): " << duration << "\n";
    outfile << "Distribution of different size cliques:\n";
    for(size_t i = 0; i < cliqueSizeDistribution.size(); i++){
        if(cliqueSizeDistribution[i] > 0)
            outfile << "Size " << i << ": " << cliqueSizeDistribution[i] << "\n";
    }
    outfile.close();
    return 0;
}