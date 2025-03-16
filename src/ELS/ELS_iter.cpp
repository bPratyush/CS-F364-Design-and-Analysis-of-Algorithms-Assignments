#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <stack>
using namespace std;
using namespace std::chrono;

// Eppstein, Löffler & Strash (2010) algorithm for finding all maximal cliques in an undirected graph

void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

unordered_set<int> setintersect(const unordered_set<int>& A, const unordered_set<int>& B){
    unordered_set<int> res;
    for (int a : A) {
        if(B.find(a) != B.end())
            res.insert(a);
    }
    return res;
}

unordered_set<int> setdiff(const unordered_set<int>& A, const unordered_set<int>& B) {
    unordered_set<int> res;
    for (int a : A) {
        if(B.find(a) == B.end())
            res.insert(a);
    }
    return res;
}

int maxCliqueSize = 0;
int totalMaximalCliques = 0;
unordered_map<int, int> cliqueSizeDistribution;

// --- Iterative Bron–Kerbosch Pivot using an explicit stack ---
struct Frame {
    unordered_set<int> P;
    unordered_set<int> R;
    unordered_set<int> X;
    // Candidates computed from P \ N(u) (the "diff" list)
    vector<int> diff;
    // Current index in diff
    int idx;
    // candidate that caused this frame to be created (for parent's update)
    // For the initial frame, candidate = -1.
    int candidate;
};

void IterativeBronKerboschPivot(const vector<unordered_set<int>>& adj,
                                const unordered_set<int>& initP,
                                const unordered_set<int>& initR,
                                const unordered_set<int>& initX) {
    stack<Frame> st;
    Frame init;
    init.P = initP;
    init.R = initR;
    init.X = initX;
    init.idx = 0;
    init.candidate = -1;
    st.push(init);

    while(!st.empty()){
        Frame &current = st.top();
        // Compute union(P, X)
        unordered_set<int> unionPX = current.P;
        unionPX.insert(current.X.begin(), current.X.end());
        // Base case: if union is empty, report current.R as a maximal clique.
        if(unionPX.empty()){
            int cliqueSize = current.R.size();
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
            // Pop this frame and update parent's candidate.
            int finishedCandidate = current.candidate;
            st.pop();
            if(!st.empty() && finishedCandidate != -1) {
                st.top().P.erase(finishedCandidate);
                st.top().X.insert(finishedCandidate);
            }
            continue;
        }
        // If we have not computed the candidate list (diff) yet, do so.
        if(current.diff.empty() && current.idx == 0) {
            // Choose a pivot u arbitrarily from union(P,X)
            int u = *unionPX.begin();
            unordered_set<int> diffSet = setdiff(current.P, adj[u]);
            // Convert to vector (order does not matter)
            for (int v : diffSet)
                current.diff.push_back(v);
        }
        // If there is an untried candidate in the diff list, process it.
        if(current.idx < current.diff.size()){
            int v = current.diff[current.idx];
            current.idx++;
            // Create new frame for the recursive call with v added.
            Frame next;
            next.R = current.R;
            next.R.insert(v);
            next.P = setintersect(current.P, adj[v]);
            next.X = setintersect(current.X, adj[v]);
            next.idx = 0;
            next.diff.clear();
            next.candidate = v;
            st.push(next);
        } else {
            // All candidates in current frame are exhausted; pop and update parent's sets.
            int finishedCandidate = current.candidate;
            st.pop();
            if(!st.empty() && finishedCandidate != -1) {
                st.top().P.erase(finishedCandidate);
                st.top().X.insert(finishedCandidate);
            }
        }
    } 
}

// Compute degeneracy ordering by repeatedly removing the vertex with smallest degree.
vector<int> degeneracyorder(const vector<unordered_set<int>>& adj){
    int n = adj.size();
    vector<bool> used(n, false);
    vector<int> degree(n, 0);
    for(int i = 0; i < n; i++)
        degree[i] = adj[i].size();
    vector<int> ordering;
    for(int k = 0; k < n; k++){
        int u = -1;
        int minDeg = INT_MAX;
        for(int i = 0; i < n; i++){
            if(!used[i] && degree[i] < minDeg){
                minDeg = degree[i];
                u = i;
            }
        }
        if(u == -1) break;
        used[u] = true;
        ordering.push_back(u);
        for(int w: adj[u]){
            if(!used[w])
                degree[w]--;
        }
    }
    reverse(ordering.begin(), ordering.end());
    return ordering;
}
void BronKerboschDegeneracy(const vector<unordered_set<int>>& adj) {
    int n = adj.size();
    vector<int> ordering = degeneracyorder(adj);
    vector<int> pos(n, 0);
    for (int i = 0; i < int(ordering.size()); i++){
        pos[ordering[i]] = i;
    }
    // For each vertex v in degeneracy order, use it as the starting point.
    for (int i = 0; i < int(ordering.size()); i++){
        int v = ordering[i];
        unordered_set<int> P;
        for (int w : adj[v]) {
            if(pos[w] > pos[v])
                P.insert(w);
        }
        unordered_set<int> X;
        for (int w : adj[v]) {
            if(pos[w] < pos[v])
                X.insert(w);
        }
        unordered_set<int> R;
        R.insert(v);
        IterativeBronKerboschPivot(adj, P, R, X);
    }
}

int main(int argc, char* argv[]) {
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " [input file]" << "\n";
        return 1;
    }
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0]=='#') continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            edgeList.push_back({u, v});
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();

    // Build adjacency list based on maximum vertex found.
    vector<unordered_set<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);

    //debug printing for large graphs.
    for (int i = 0; i < int(adj.size()); ++i) {
        cout << "Vertex " << i << ": ";
        if(adj[i].empty()){
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for (int neighbor : neighbors)
                cout << neighbor << " ";
        }
        cout << "\n";
    }

    auto start = high_resolution_clock::now();
    BronKerboschDegeneracy(adj);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    outfile << "Largest size of the clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << "\n";
    outfile << "Execution time (ms): " << duration.count() << "\n";
    outfile << "Distribution of different size cliques:\n";
    for(const auto& pair : cliqueSizeDistribution)
        outfile << "Size " << pair.first << ": " << pair.second << "\n";
    outfile.close();

    return 0;
}