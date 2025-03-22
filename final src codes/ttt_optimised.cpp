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
using namespace std;
using namespace std::chrono;

// Tomita, Tanata & Takahashi (2006) algorithm for finding all maximal cliques in an undirected graph - Optimised

vector<int> Q;   
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
map<int, int> cliqueSizeDistribution;

void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

void EXPAND(unordered_set<int> SUBG, unordered_set<int> CAND, vector<unordered_set<int>>& adj) {
    if(SUBG.empty()) {
        int cliqueSize = Q.size();
        if(cliqueSize > 1) {
            maxCliqueSize = max(maxCliqueSize, cliqueSize);
            totalMaximalCliques++;
            cliqueSizeDistribution[cliqueSize]++;
        }
        return;
    } 
    int u = -1, maxCount = 0;
    for(auto it = SUBG.begin(); it != SUBG.end(); ++it) {
        int cnt = 0;
        for(auto jt = CAND.begin(); jt != CAND.end(); ++jt) {
            if(adj[*it].find(*jt) != adj[*it].end())
                ++cnt;
        }
        if(cnt > maxCount) {
            maxCount = cnt;
            u = *it;
        }
    }
    unordered_set<int> Extu;
    Extu.reserve(CAND.size());
    for(auto it = CAND.begin(); it != CAND.end(); ++it) {
        if(adj[u].find(*it) == adj[u].end())
            Extu.insert(*it);
    }
    while(!Extu.empty()){
        int q = *Extu.begin();
        Q.push_back(q);
        unordered_set<int> SUBGq;
        SUBGq.reserve(SUBG.size());
        for(auto it = SUBG.begin(); it != SUBG.end(); ++it) {
            if(adj[q].find(*it) != adj[q].end())
                SUBGq.insert(*it);
        }
        unordered_set<int> CANDq;
        CANDq.reserve(CAND.size());
        for(auto it = CAND.begin(); it != CAND.end(); ++it) {
            if(adj[q].find(*it) != adj[q].end())
                CANDq.insert(*it);
        }
        EXPAND(std::move(SUBGq), std::move(CANDq), adj);
        CAND.erase(q);
        Extu.erase(q);
        Q.pop_back();
    }
}
    
void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    vector<int> ordering;
    for (int i = 0; i < V; i++){
        Vset.insert(i);
        ordering.push_back(i);
    }
    for(auto it = ordering.begin(); it != ordering.end(); ++it) {
        cout << "[DEBUG] Processing node " << *it 
             << " (loop index " << distance(ordering.begin(), it) << ")" << endl;
    }
    EXPAND(Vset, Vset, adj);
}
    
int main(int argc, char* argv[]) {
    struct rlimit rl;
    rlim_t stack_size = 512 * 1024 * 1024;
    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result != 0) {
        cerr << "Error getting stack limit: " << strerror(errno) << endl;
        return 1;
    }
    if (rl.rlim_cur < stack_size) {
        rl.rlim_cur = stack_size;
        if (rl.rlim_max < rl.rlim_cur) {
            rl.rlim_max = rl.rlim_cur; 
        }

        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0) {
            cerr << "Error setting stack limit: " << strerror(errno) << endl;
        } else {
            cerr << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << endl;
        }
    }

    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    
    if(argc < 2){
        cerr << "Usage: " << argv[0] << " [input file]" << endl;
        return 1;
    }
    
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int, int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0] == '#')
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
    
    auto start = high_resolution_clock::now();
    CLIQUES(adj, maxVertex + 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    outfile << "Largest size of the clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << "\n";
    outfile << "Execution time (ms): " << duration.count() << "\n";
    outfile << "Distribution of different size cliques:\n";
    for(auto it = cliqueSizeDistribution.begin(); it != cliqueSizeDistribution.end(); ++it) {
        outfile << "Size " << it->first << ": " << it->second << "\n";
    }
    outfile.close();
    return 0;
}