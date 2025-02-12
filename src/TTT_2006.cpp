#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <cstdlib>
using namespace std;
//Tomita, Tanata & Takahashi (2006) algorithm for finding all maximal cliques in an undirected graph
void addEdge(int u,int v,vector<unordered_set<int> > &adj){
    adj[u].insert(v);
    adj[v].insert(u);
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
    int pivot=-1,maxcnt=-1;
    for(int u:unionPX){
        int cnt=0;
        for(int v:P){
            if(adj[u].count(v) == 0) cnt++;
        }
        if(cnt>maxcnt){
            maxcnt=cnt;
            pivot=u;
        }
    }
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

int main(int argc, char* argv[]){
    //Input Format: First line contains number of vertices V and next lines contain edges in (u,v) format
    ifstream infile(argv[1]);
    int V;
    infile >>V;
    vector<unordered_set<int> > adj(V);
    int u,v;
    while (infile>>u>>v) addEdge(u,v,adj);
    infile.close();
    unordered_set<int> R, P, X;
    for (int i = 0; i < V; i++)
    P.insert(i);
    cout << "Maximal cliques in the graph:" << endl;
    bronKerbosch2(R,P,X, adj);
    return 0;
}