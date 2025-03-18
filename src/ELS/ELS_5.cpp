#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <thread>
#include <mutex>
#include <future>
#include <cstdint>
using namespace std;
using namespace std::chrono;

// ---------------- Bitset Class ----------------
class Bitset {
public:
    int size;
    vector<uint64_t> bits;
    Bitset(int n, bool fill = false): size(n), bits((n + 63) / 64, fill ? ~0ULL : 0ULL) {
        if(fill && n % 64 != 0){
            bits.back() &= ((uint64_t)1 << (n % 64)) - 1;
        }
    }
    Bitset(): size(0) {}
    
    bool empty() const {
        for(auto b : bits)
            if(b != 0) return false;
        return true;
    }
    
    void set(int i) {
        int idx = i / 64, pos = i % 64;
        bits[idx] |= (1ULL << pos);
    }
    
    void reset(int i) {
        int idx = i / 64, pos = i % 64;
        bits[idx] &= ~(1ULL << pos);
    }
    
    bool test(int i) const {
        int idx = i / 64, pos = i % 64;
        return bits[idx] & (1ULL << pos);
    }
    
    // Bitwise AND (intersection)
    Bitset intersect(const Bitset &other) const {
        Bitset result(size);
        for(size_t i = 0; i < bits.size(); i++){
            result.bits[i] = bits[i] & other.bits[i];
        }
        return result;
    }
    // Bitwise difference: A \ other.
    Bitset difference(const Bitset &other) const {
        Bitset result(size);
        for(size_t i = 0; i < bits.size(); i++){
            result.bits[i] = bits[i] & ~(other.bits[i]);
        }
        return result;
    }
    // OR with other (in-place union)
    void orWith(const Bitset &other) {
        for(size_t i = 0; i < bits.size(); i++){
            bits[i] |= other.bits[i];
        }
    }
    // Count bits set.
    int count() const {
        int cnt = 0;
        for(auto b : bits)
            cnt += __builtin_popcountll(b);
        return cnt;
    }
    // Return all indices with bit set.
    vector<int> getElements() const {
        vector<int> elems;
        for (int i = 0; i < size; i++){
            if(test(i))
                elems.push_back(i);
        }
        return elems;
    }
};

// ---------------- Global Clique Statistics ----------------
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
unordered_map<int, int> cliqueSizeDistribution;
mutex clique_mutex;  // Protects updates

// ---------------- Helper Functions ----------------
void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

// ---------------- Bron–Kerbosch with Bitset (In–Place backtracking) ----------------
int choosePivot(const Bitset &unionSet, const Bitset &P, const vector<Bitset> &bitAdj) {
    int pivot = -1, maxCount = -1;
    vector<int> candidates = unionSet.getElements();
    for (int candidate : candidates) {
        int cnt = (bitAdj[candidate].intersect(P)).count();
        if(cnt > maxCount){
            maxCount = cnt;
            pivot = candidate;
        }
    }
    return pivot;
}

void BronKerboschPivotBit(Bitset P, Bitset R, Bitset X, const vector<Bitset>& bitAdj) {
    if(P.empty() && X.empty()){
        int cliqueSize = R.count();
        if(cliqueSize > 1) {
            lock_guard<mutex> lock(clique_mutex);
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    }
    Bitset unionPX = P;
    unionPX.orWith(X);
    int pivot = choosePivot(unionPX, P, bitAdj);
    if(pivot == -1){
        vector<int> elems = unionPX.getElements();
        if(!elems.empty())
            pivot = elems[0];
        else
            return;
    }
    Bitset diff = P.difference(bitAdj[pivot]);
    vector<int> diffElements = diff.getElements();
    for (int v : diffElements) {
        Bitset newR = R;
        newR.set(v);
        Bitset newP = P.intersect(bitAdj[v]);
        Bitset newX = X.intersect(bitAdj[v]);
        BronKerboschPivotBit(newP, newR, newX, bitAdj);
        P.reset(v);
        X.set(v);
    }
}

// ---------------- Main Function ----------------
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    auto start = high_resolution_clock::now();
    
    // Read input file and build graph.
    ifstream ip("wiki-Vote.txt");
    if(!ip){
        cerr << "Failed to open file." << endl;
        return 1;
    }
    
    unordered_set<int> nodes;
    int a, b;
    int maxVertex = 0;
    string line;
    while(getline(ip, line)){
        if(line.empty() || line[0]=='#')
            continue;
        istringstream ss(line);
        ss >> a >> b;
        if(!ss.fail()){
            maxVertex = max(maxVertex, max(a, b));
            nodes.insert(a);
            nodes.insert(b);
        }
    }
    ip.clear();
    ip.seekg(0, ios::beg);
    
    vector<unordered_set<int>> adj(maxVertex + 1);
    while(getline(ip, line)){
        if(line.empty() || line[0]=='#')
            continue;
        istringstream ss(line);
        ss >> a >> b;
        if(!ss.fail()){
            addEdge(a, b, adj);
            nodes.insert(a);
            nodes.insert(b);
        }
    }
    ip.close();
    cout << "Adjacency list created." << endl;
    
    // --- Compute degeneracy ordering (unchanged) ---
    // (We still use the unordered_set based method here.)
    unordered_map<int, unordered_set<int>> tempadj;
    for (int i = 0; i < adj.size(); i++){
        tempadj[i] = adj[i];
    }
    unordered_map<int, int> deg;
    unordered_map<int, unordered_set<int>> deglist;
    for (int node : nodes) {
        int d = (node < adj.size()) ? adj[node].size() : 0;
        deg[node] = d;
        deglist[d].insert(node);
    }
    cout << "Degree list created." << endl;
    int global_maxdeg = 0;
    for (auto &p : deg)
        global_maxdeg = max(global_maxdeg, p.second);
    vector<int> degenlist;
    degenlist.reserve(nodes.size());
    int minDegree = 0;
    while(degenlist.size() < nodes.size()){
        while(minDegree <= global_maxdeg && deglist[minDegree].empty())
            minDegree++;
        if(minDegree > global_maxdeg)
            break;
        int u = *deglist[minDegree].begin();
        deglist[minDegree].erase(u);
        degenlist.push_back(u);
        for (auto v : tempadj[u]){
            tempadj[v].erase(u);
            int curr = deg[v];
            deglist[curr].erase(v);
            deg[v]--;
            int newdeg = deg[v];
            deglist[newdeg].insert(v);
            if(newdeg < minDegree)
                minDegree = newdeg;
        }
    }
    cout << "Degeneracy ordering created." << endl;
    
    // Build lookup: pos[v] is the position of v in degenlist.
    unordered_map<int, int> pos;
    for (int i = 0; i < degenlist.size(); i++){
        pos[degenlist[i]] = i;
    }
    
    // Build Bitset representation for each vertex's neighbors.
    int nVertices = maxVertex + 1;
    vector<Bitset> bitAdj;
    bitAdj.reserve(nVertices);
    for (int i = 0; i < nVertices; i++){
        Bitset bs(nVertices, false);
        for (int neigh : adj[i])
            bs.set(neigh);
        bitAdj.push_back(bs);
    }
    
    // ---------------- Process nodes in batches using a thread pool (via std::async) ----------------
    // ---------------- Process nodes in batches using a thread pool (via std::async) ----------------
int totalNodes = degenlist.size();
int batchSize = 10000;  // Tune as needed.
vector<future<void>> batchFutures;
for (int batchStart = 0; batchStart < totalNodes; batchStart += batchSize) {
    int batchEnd = min(batchStart + batchSize, totalNodes);
    // Launch one async task per batch.
    batchFutures.push_back(async(launch::async, [batchStart, batchEnd, &degenlist, &pos, &bitAdj, nVertices, &adj](){
        for (int i = batchStart; i < batchEnd; i++){
            int vi = degenlist[i];
            // Create Bitset representations for P and X.
            Bitset P(nVertices, false), X(nVertices, false), R(nVertices, false);
            // For each neighbor of vi, use degeneracy ordering (using pos)
            for (int neighbor : adj[vi]) { // Using original 'adj'
                if(pos.find(neighbor) != pos.end()){
                    if (pos[neighbor] > i)
                        P.set(neighbor);
                    else if (pos[neighbor] < i)
                        X.set(neighbor);
                }
            }
            R.set(vi);
            BronKerboschPivotBit(P, R, X, bitAdj);
        }
    }));
}
for(auto &f : batchFutures)
    f.get();
    
    // ---------------- Output results ----------------
    ofstream op("output.txt");
    if(!op){
        cerr << "Failed to open output file." << endl;
        return 1;
    }
    op << "\nNo. of maximal cliques: " << totalMaximalCliques << "\n";
    op << "\nClique size distribution:\n";
    for (auto &p : cliqueSizeDistribution){
        op << "Clique size: " << p.first << ", Frequency: " << p.second << "\n";
    }
    op << "Maximum clique size: " << maxCliqueSize << "\n";
    op.close();
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start);
    cout << "Number of nodes: " << nodes.size() << "\n";
    cout << "No. of elements in degeneracy ordering: " << degenlist.size() << "\n";
    cout << "Time taken: " << duration.count() << " milliseconds" << "\n";
    
    return 0;
}