#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <map>
#include <thread>
using namespace std;
using namespace std::chrono;

int maxCliqueSize = 0;
int totalMaximalCliques = 0;
map<int, int> cliqueSizeDistribution;
mutex clique_mutex;  // protects updates below

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

//
// EXPAND: same as TTT algorithm but now Q (the current clique) is passed as a parameter.
// A recursion depth parameter is added; if depth is below THRESHOLD, recursive calls are spawned in parallel.
// (You may adjust THRESHOLD as needed.)
//
const int THRESHOLD = 2;  // control thread spawn depth

void EXPAND(unordered_set<int> SUBG, unordered_set<int> CAND, vector<unordered_set<int>>& adj, vector<int> Q, int depth = 0) {
    if(SUBG.empty()){
        int cliqueSize = Q.size();
        if(cliqueSize > 1) {
            lock_guard<mutex> lock(clique_mutex);
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    } else {
        int u = -1;
        int maxCount = 0;
        for (int x : SUBG) {
            int cnt = 0;
            for (int y : CAND) {
                if(adj[x].find(y) != adj[x].end())
                    ++cnt;
            }
            if(cnt > maxCount) {
                maxCount = cnt;
                u = x;
            }
        }
        unordered_set<int> Extu;
        for (int v : CAND) {
            if(adj[u].find(v) == adj[u].end())
                Extu.insert(v);
        }
        unordered_set<int> FINI;        
        vector<thread> threads;
        while(!Extu.empty()){
            int q = *Extu.begin();
            Q.push_back(q);
            unordered_set<int> SUBGq;
            for (int v : SUBG) {
                if(adj[q].find(v) != adj[q].end())
                    SUBGq.insert(v);
            }
            unordered_set<int> CANDq;
            for (int v : CAND) {
                if(adj[q].find(v) != adj[q].end())
                    CANDq.insert(v);
            }
            if(depth < THRESHOLD) {
                threads.emplace_back(EXPAND, SUBGq, CANDq, ref(adj), Q, depth+1);
            } else {
                EXPAND(SUBGq, CANDq, adj, Q, depth+1);
            }
            CAND.erase(q);
            FINI.insert(q);
            Q.pop_back();
            Extu.clear();
            for (int v : CAND) {
                if(adj[u].find(v) == adj[u].end())
                    Extu.insert(v);
            }
        }
        for(auto &thr : threads)
            thr.join();
    }
}

void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    for(int i = 0; i < V; i++)
        Vset.insert(i);
    vector<int> Q;  // initially empty
    EXPAND(Vset, Vset, adj, Q, 0);
}

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    if(argc < 2) {
        cerr << "Usage: " << argv[0] << " input_file" << endl;
        return 1;
    }
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int,int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0]=='#')
            continue;
        istringstream iss(line);
        int u, v;
        if(iss >> u >> v){
            edgeList.push_back({u, v});
            maxVertex = max(maxVertex, max(u, v));
        }
    }
    infile.close();
    vector<unordered_set<int>> adj(maxVertex + 1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    // Print adjacency list (unchanged)
    for(int i = 0; i < adj.size(); ++i){
        cout << "Vertex " << i << ": ";
        if(adj[i].empty()){
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for(int neighbor : neighbors){
                cout << neighbor << " ";
            }
        }
        cout << endl;
    }
    
    auto start = high_resolution_clock::now();
    CLIQUES(adj, maxVertex + 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;
    for(const auto& pair : cliqueSizeDistribution)
        outfile << "Size " << pair.first << ": " << pair.second << endl;
    outfile.close();
    return 0;
}