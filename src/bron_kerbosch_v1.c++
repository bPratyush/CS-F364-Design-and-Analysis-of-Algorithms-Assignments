#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
using namespace std;

void addEdge(int u,int v,vector<unordered_set<int> > &adj) {
    adj[u].insert(v);
    adj[v].insert(u);
}

vector<int> degeneracyOrdering(int V, const vector<unordered_set<int> > &adj) {
    vector<int>order;
    vector<int>deg(V);
    vector<unordered_set<int> >bucket(V);
    for(int v=0;v<V;v++){
        deg[v]=adj[v].size();
        bucket[deg[v]].insert(v);
    }
    for(int d=0;d<V;d++){
        while(!bucket[d].empty()){
            int v=*bucket[d].begin();
            bucket[d].erase(v);
            order.push_back(v);
            for (int u:adj[v]) {
                if (deg[u] > d) {
                    bucket[deg[u]].erase(u);
                    deg[u]--;
                    bucket[deg[u]].insert(u);
                }
            }
        }
    }
    return order;
}

void bronKerbosch2(unordered_set<int> R,unordered_set<int> P,unordered_set<int> X,const vector<unordered_set<int> > &adj) {
    if(P.empty()&&X.empty()){
        cout << "Maximal Clique: ";
        for(int v:R) cout << v << " ";
        cout << endl;
        return;
    }
    unordered_set<int> unionPX=P;
    unionPX.insert(X.begin(),X.end());
    int pivot=*unionPX.begin();
    unordered_set<int>candi;
    for (int v:P){
        if(adj[pivot].find(v)==adj[pivot].end()) candi.insert(v);
    }
    for(int v:candi){
        unordered_set<int> newR=R;
        newR.insert(v);
        unordered_set<int>newP,newX;
        for(int n:adj[v]){
            if(P.find(n)!=P.end()) newP.insert(n);
            if(X.find(n)!=X.end()) newX.insert(n);
        }
        bronKerbosch2(newR,newP,newX,adj);
        P.erase(v);
        X.insert(v);
    }
}

void bronKerbosch3(int V,const vector<unordered_set<int> >&adj) {
    vector<int> order=degeneracyOrdering(V, adj);
    unordered_set<int> P,X;
    for(int v=0;v<V;v++) P.insert(v);
    for(int v:order){
        unordered_set<int> newP,newX;
        for(int n:adj[v]){
            if(P.find(n)!=P.end()) newP.insert(n);
            if(X.find(n)!=X.end()) newX.insert(n);
        }
        unordered_set<int> newR;
        newR.insert(v);
        bronKerbosch2(newR,newP,newX,adj);
        P.erase(v);
        X.insert(v);
    }
}

int main() {
    //Example Graph
    int V=5;
    vector<unordered_set<int> >adj(V);
    addEdge(0, 1, adj);
    addEdge(0, 2, adj);
    addEdge(1, 2, adj);
    addEdge(1, 3, adj);
    addEdge(3, 4, adj);
    cout << "Maximal cliques in the graph:\n";
    bronKerbosch3(V, adj);
    return 0;
}