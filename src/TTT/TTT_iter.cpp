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
#include <stack>
using namespace std;
using namespace std::chrono;

//Tomita, Tanata & Takahashi (2006) algorithm for finding all maximal cliques in an undirected graph

vector<int> Q;  // (only used in the recursive version; the iterative version stores clique R per frame)
int maxCliqueSize = 0;
int totalMaximalCliques = 0;
map<int,int> cliqueSizeDistribution;

void addEdge(int u, int v, vector<unordered_set<int>>& adj) {
    if(u >= adj.size() || v >= adj.size()){
        int newSize = max(u, v) + 1;
        adj.resize(newSize);
    }
    // Treat edge as undirected
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
        if(B.find(a)==B.end())
            res.insert(a);
    }
    return res;
}

struct Frame {
    unordered_set<int> SUBG;   
    unordered_set<int> CAND;    
    unordered_set<int> R;       
    vector<int> diff;           
    int idx;                    
    int pivot;                 
    int candidate;              
};

void IterativeEXPAND(const vector<unordered_set<int>>& adj,
                     const unordered_set<int>& initSUBG,
                     const unordered_set<int>& initCAND,
                     const unordered_set<int>& initR) {
    stack<Frame> st;
    Frame init;
    init.SUBG = initSUBG;
    init.CAND = initCAND;
    init.R = initR;
    init.idx = 0;
    init.diff.clear();
    init.pivot = -1;
    init.candidate = -1;
    st.push(init);

    while(!st.empty()){
        Frame &cur = st.top();
        // If SUBG is empty then current R is a maximal clique.
        if(cur.SUBG.empty()){
            int sizeR = cur.R.size();
            maxCliqueSize = max(maxCliqueSize, sizeR);
            totalMaximalCliques++;
            cliqueSizeDistribution[sizeR]++;
            int finishedCandidate = cur.candidate;  // candidate that led to this frame
            st.pop();
            if(!st.empty() && finishedCandidate != -1) {
                // In the parent's frame, remove the candidate that was completed.
                st.top().CAND.erase(finishedCandidate);
                // Also, force recomputation of diff in parent's frame by clearing it.
                st.top().diff.clear();
                st.top().idx = 0;
            }
            continue;
        }
        // If we have not computed the candidate list for this frame, compute it.
        if(cur.diff.empty() && cur.idx == 0) {
            // Choose a pivot arbitrarily from SUBG (or SUBG ∪ CAND; here we choose from SUBG)
            cur.pivot = *cur.SUBG.begin();
            // The candidate list is: all v in CAND that are NOT adjacent to pivot.
            unordered_set<int> diffSet = setdiff(cur.CAND, adj[cur.pivot]);
            for (int v : diffSet)
                cur.diff.push_back(v);
        }
        // If there remain candidates in the diff list, process the next candidate.
        if(cur.idx < cur.diff.size()){
            int v = cur.diff[cur.idx];
            cur.idx++;
            // Create a new frame corresponding to adding v to R.
            Frame next;
            next.R = cur.R;
            next.R.insert(v);
            next.SUBG = unordered_set<int>();
            for (int w : cur.SUBG)
                if(adj[v].find(w) != adj[v].end())
                    next.SUBG.insert(w);
            next.CAND = unordered_set<int>();
            for (int w : cur.CAND)
                if(adj[v].find(w) != adj[v].end())
                    next.CAND.insert(w);
            next.diff.clear();
            next.idx = 0;
            next.pivot = -1;
            next.candidate = v;
            st.push(next);
        } else {
            // All candidates in the current frame have been processed.
            int finishedCandidate = cur.candidate;
            st.pop();
            if(!st.empty() && finishedCandidate != -1) {
                st.top().CAND.erase(finishedCandidate);
                st.top().diff.clear();
                st.top().idx = 0;
            }
        }
    }
}

void CLIQUES(vector<unordered_set<int>>& adj, int V) {
    unordered_set<int> Vset;
    for (int i = 0; i < V; i++)
        Vset.insert(i);
    IterativeEXPAND(adj, Vset, Vset, unordered_set<int>());
}

int main(int argc, char* argv[]) {
    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    string line;
    vector<pair<int,int>> edgeList;
    int maxVertex = 0;
    while(getline(infile, line)){
        if(line.empty() || line[0]=='#')
            continue;
        istringstream iss(line);
        int u,v;
        if(iss >> u >> v){
            edgeList.push_back({u,v});
            maxVertex = max(maxVertex, max(u,v));
        }
    }
    infile.close();
    
    // Build the adjacency list based on maximum vertex.
    vector<unordered_set<int>> adj(maxVertex+1);
    for(auto &edge : edgeList)
        addEdge(edge.first, edge.second, adj);
    
    //Debug printing for large graphs
    for(int i = 0; i < (int)adj.size(); ++i){
        cout << "Vertex " << i << ": ";
        if(adj[i].empty()){
            cout << "No neighbors";
        } else {
            vector<int> neighbors(adj[i].begin(), adj[i].end());
            sort(neighbors.begin(), neighbors.end());
            for (int neighbor : neighbors)
                cout << neighbor << " ";
        }
        cout << endl;
        cout.flush();
    }
    
    auto start = high_resolution_clock::now();
    // Use maxVertex+1 as the vertex count
    CLIQUES(adj, maxVertex+1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    outfile << "Largest size of the clique: " << maxCliqueSize << endl;
    outfile << "Total number of maximal cliques: " << totalMaximalCliques << endl;
    outfile << "Execution time (ms): " << duration.count() << endl;
    outfile << "Distribution of different size cliques:" << endl;
    for (const auto& pair : cliqueSizeDistribution)
        outfile << "Size " << pair.first << ": " << pair.second << endl;
    outfile.close();
    
    return 0;
}