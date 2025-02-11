#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <cstdlib>
using namespace std;
//Eppstein, Löffler & Strash (2010) algorithm for finding all maximal cliques in an undirected graph
void addEdge(int u,int v,vector<unordered_set<int> > &adj){
    adj[u].insert(v);
    adj[v].insert(u);
}

vector<int> degeneracyOrdering(int V,vector<unordered_set<int> > &adj){
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

void bronKerbosch2(unordered_set<int> R,unordered_set<int> P,unordered_set<int> X,vector<unordered_set<int> > &adj){
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

void bronKerbosch3(int V,vector<unordered_set<int> >&adj){
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

int main(int argc, char* argv[]){
    //Input Format: First line contains number of vertices V and next lines contain edges in (u,v) format
    ifstream infile(argv[1]);
    int V;
    infile >>V;
    vector<unordered_set<int> > adj(V);
    int u,v;
    while (infile>>u>>v) addEdge(u,v,adj);
    infile.close();
    cout << "Maximal cliques in the graph:" << endl;
    bronKerbosch3(V, adj);
    return 0;
}