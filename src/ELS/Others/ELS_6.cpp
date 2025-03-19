#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
using namespace std;
using namespace std::chrono;

void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

unordered_set<int> setintersect(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for(int a : A) {
        if(B.find(a) != B.end())
            res.insert(a);
    }
    return res;
}

unordered_set<int> setdiff(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for(int a : A) {
        if(B.find(a) == B.end())
            res.insert(a);
    }
    return res;
}

int maxCliqueSize = 0;
int totalMaximalCliques = 0;
unordered_map<int, int> cliqueSizeDistribution;

void BronKerboschPivot(unordered_set<int> P, unordered_set<int> R, unordered_set<int> X, 
                         const vector<unordered_set<int>>& adj) {
    unordered_set<int> unionPX = P;
    unionPX.insert(X.begin(), X.end());
    if(unionPX.empty()){
        int cliqueSize = R.size();
        if(cliqueSize > 1) {
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    }
    int u = *unionPX.begin();
    unordered_set<int> diff = setdiff(P, adj[u]);
    for(int v : diff) {
        unordered_set<int> newP = setintersect(P, adj[v]);
        unordered_set<int> newX = setintersect(X, adj[v]);
        unordered_set<int> newR = R;
        newR.insert(v);
        BronKerboschPivot(newP, newR, newX, adj);
        P.erase(v);
        X.insert(v);
    }
}

int main(int argc, char* argv[]){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    auto startTime = high_resolution_clock::now();
    
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " [input file]" << endl;
        return 1;
    }
    
    ifstream ip(argv[1]);
    if(!ip){
        cerr << "Failed to open file: " << argv[1] << endl;
        return 1;
    }
    
    // Read file to determine max vertex and collect nodes.
    unordered_set<int> nodes;
    int a, b, maxVertex = 0;
    string line;
    while(getline(ip, line)){
        if(line.empty() || line[0] == '#')
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
    
    // Build the graph.
    vector<unordered_set<int>> adj(maxVertex + 1);
    while(getline(ip, line)){
        if(line.empty() || line[0] == '#')
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
    
    // Compute degeneracy ordering using a bucket-based method.
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
    int global_maxdeg = 0;
    for(auto &p : deg){
        global_maxdeg = max(global_maxdeg, p.second);
    }
    vector<int> degenlist;
    degenlist.reserve(nodes.size());
    int minDegree = 0;
    while(degenlist.size() < nodes.size()){
        while(minDegree <= global_maxdeg && deglist[minDegree].empty()){
            minDegree++;
        }
        if(minDegree > global_maxdeg)
            break; 
        int u = *deglist[minDegree].begin();
        deglist[minDegree].erase(u);
        degenlist.push_back(u);
        for(auto v : tempadj[u]){
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
    
    // Build position lookup.
    unordered_map<int, int> pos;
    for(int i = 0; i < degenlist.size(); i++){
        pos[degenlist[i]] = i;
    }
    
    // Sequential outer loop for Bron–Kerbosch pivot.
    for (int i = 0; i < (int)degenlist.size(); i++){
        int vi = degenlist[i];
        unordered_set<int> P, X, R;
        for (int w : adj[vi]) {
            if(pos.find(w) != pos.end()){
                if(pos[w] > i)
                    P.insert(w);
                else if(pos[w] < i)
                    X.insert(w);
            }
        }
        R.insert(vi);
        BronKerboschPivot(P, R, X, adj);
    }
    
    // Output the required final information.
    ofstream op("output.txt");
    if(!op){
        cerr << "Failed to open output file." << endl;
        return 1;
    }
    op << "No. of maximal cliques: " << totalMaximalCliques << "\n";
    op << "Maximum clique size: " << maxCliqueSize << "\n";
    op << "Clique size distribution:\n";
    for(auto &p : cliqueSizeDistribution){
        op << "Clique size: " << p.first << ", Frequency: " << p.second << "\n";
    }
    op.close();
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - startTime);
    cout << "Time taken: " << duration.count() << " milliseconds" << "\n";
    return 0;
}