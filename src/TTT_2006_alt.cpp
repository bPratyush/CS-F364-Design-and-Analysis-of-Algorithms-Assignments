#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

vector<int> Q;

void EXPAND(unordered_set<int> SUBG,unordered_set<int> CAND,const vector<unordered_set<int> >& adj){
    if(SUBG.empty()){
        cout << "Clique: ";
        for(int v:Q) cout << v << " ";
        cout << endl;
        return;
    }
    else{
        // u := a vertex in SUBG that maximises |CAND ∩ gamma(u)|
        int u=-1;
        int maxCount=0;
        for(int x:SUBG){
            int cnt=0;
            for(int y:CAND)
                if(adj[x].find(y)!=adj[x].end()) ++cnt;
            if(cnt>maxCount){
                maxCount=cnt;
                u=x;
            }
        }
        // let Extu = CAND - gamma(u);
        unordered_set<int> Extu;
        for (int v : CAND){
            if (adj[u].find(v)==adj[u].end()) Extu.insert(v);
        }
        // FINI := empty
        unordered_set<int> FINI;
        while (!Extu.empty()){
            int q=*Extu.begin();
            cout << q << ", ";
            // Q := Q union {q}
            Q.push_back(q);
            // SUBGq = SUBG ∩ gamma(q)
            unordered_set<int> SUBGq;
            for(int v:SUBG){
                if(adj[q].find(v)!=adj[q].end()) SUBGq.insert(v);
            }
            // CANDq = CAND ∩ gamma(q)
            unordered_set<int> CANDq;
            for (int v:CAND){
                if (adj[q].find(v)!=adj[q].end()) CANDq.insert(v);
            }
            // EXPAND(SUBGq, CANDq)
            EXPAND(SUBGq,CANDq,adj);
            // CAND := CAND - {q} and FINI := FINI union {q};
            CAND.erase(q);
            FINI.insert(q);
            cout << "back, ";
            // Q := Q - {q}
            Q.pop_back();
            // Recompute Extu = CAND - gamma(u)
            Extu.clear();
            for (int v:CAND){
                if (adj[u].find(v)==adj[u].end()) Extu.insert(v);
            }
        }
    }
}

void CLIQUES(const vector<unordered_set<int> >& adj, int V) {
    unordered_set<int> Vset;
    for (int i = 0; i < V; i++) Vset.insert(i);
    EXPAND(Vset, Vset, adj);
}

void addEdge(int u,int v,vector<unordered_set<int> > &adj){
    adj[u].insert(v);
    adj[v].insert(u);
}

int main(int argc, char* argv[]){
    ifstream infile(argv[1]);
    int V;
    infile >>V;
    vector<unordered_set<int> > adj(V);
    int u,v;
    while (infile>>u>>v) addEdge(u,v,adj);
    infile.close();
    cout << "Maximal cliques in the graph:" << endl;
    CLIQUES(adj, V);
    return 0;
}