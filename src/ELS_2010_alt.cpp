#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
using namespace std;

//Eppstein, Löffler & Strash (2010) algorithm for finding all maximal cliques in an undirected graph
void addEdge(int u,int v,vector<unordered_set<int>>& adj) {
    if (u>=adj.size()||v>=adj.size()){
        int newSize=max(u,v)+1;
        adj.resize(newSize);
    }
    adj[u].insert(v);
    adj[v].insert(u);
}

unordered_set<int> setintersect(const unordered_set<int>& A,const unordered_set<int>& B){
    unordered_set<int> res;
    for(int a:A){
        if(B.find(a)!=B.end()) res.insert(a);
    }
    return res;
}

unordered_set<int> setdiff(const unordered_set<int>& A,const unordered_set<int>& B) {
    unordered_set<int> res;
    for(int a:A){
        if(B.find(a)==B.end()) res.insert(a);
    }
    return res;
}

void BronKerboschPivot(unordered_set<int> P,unordered_set<int> R,unordered_set<int> X,const vector<unordered_set<int> >& adj) {
    unordered_set<int> unionPX=P;
    //if P union X = empty then report R as a maximal clique
    unionPX.insert(X.begin(),X.end());
    if (unionPX.empty()){
        cout << "Maximal Clique: ";
        for (int v:R) cout <<v<< " ";
        cout << endl;
        return;
    }
    //Choose a pivot u from P union X (arbitrarily).
    int u = *unionPX.begin();
    //For each vertex v belonging to P \ gamma(u)
    unordered_set<int> diff=setdiff(P,adj[u]);
    for(int v:diff){
        unordered_set<int> newP=setintersect(P,adj[v]);
        unordered_set<int> newX=setintersect(X,adj[v]);
        unordered_set<int> newR=R;
        newR.insert(v);
        BronKerboschPivot(newP, newR, newX, adj);
        P.erase(v);
        X.insert(v);
    }
}

/*Degeneracy ordering, can be computed by a simple greedy strategy of repeatedly removing
a vertex with smallest degree (and its incident edges) from the graph until it is empty*/
vector<int> degeneracyorder(const vector<unordered_set<int> >& adj){
    int n=adj.size();
    vector<bool> used(n,false);
    vector<int> degree(n,0);
    for(int i=0;i<n;i++) degree[i]=adj[i].size();
    vector<int> ordering;
    for(int k=0;k<n;k++) {
        int u=-1;
        int minDeg=INT_MAX;
        for(int i=0;i<n;i++) {
            if(!used[i]&&degree[i]<minDeg) {
                minDeg=degree[i];
                u=i;
            }
        }
        if(u==-1) break;
        used[u]=true;
        ordering.push_back(u);
        for(int w:adj[u]){
            if (!used[w]) degree[w]--;
        }
    }
    reverse(ordering.begin(), ordering.end());
    return ordering;
}

void BronKerboschDegeneracy(const vector<unordered_set<int> >&adj) {
    int n=adj.size();
    vector<int> ordering=degeneracyorder(adj);
    vector<int> pos(n,0);
    for(int i=0;i<n;i++) pos[ordering[i]] = i;
    // Process each vertex in the degeneracy ordering.
    for (int i=0;i<n;i++) {
        int vi=ordering[i];
        // P <- gamma(vi) intersection {vi+1, ..., vn-1}
        unordered_set<int> P;
        for(int w:adj[vi]){
            if(pos[w]>pos[vi]) P.insert(w);
        }
        // X <- gamma(vi) intersection {v0, ..., vi-1}
        unordered_set<int> X;
        for(int w:adj[vi]){
            if(pos[w]<pos[vi]) X.insert(w);
        }
        // R is initialized to {vi}
        unordered_set<int> R;
        R.insert(vi);
        BronKerboschPivot(P, R, X, adj);
    }
}

int main(int argc, char* argv[]) {
    ifstream infile(argv[1]);
    string line;
    // Skip comment lines
    while (getline(infile, line)) {
        if (line[0] != '#') break;
    }
    // Read number of vertices and edges
    istringstream iss(line);
    int V, E;
    iss >> V >> E;
    vector<unordered_set<int>> adj(V);
    int u, v;
    while (infile >> u >> v) {
        addEdge(u, v, adj);
    }
    infile.close();
    BronKerboschDegeneracy(adj);
    return 0;
}