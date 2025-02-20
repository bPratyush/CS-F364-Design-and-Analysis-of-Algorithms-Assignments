#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <cstdlib>
using namespace std;
//Chiba & Nishizeki (1985) algorithm for finding all maximal cliques in a graph
vector<int> newToOld;

void printclique(const unordered_set<int>& C){
    cout << "Clique: ";
    for(int v:C) cout<<newToOld[v]<<" ";
    cout << endl;
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

void UPDATE(int i,unordered_set<int>& C,int n,const vector<unordered_set<int> >& adj,vector<int>& S,vector<int>& T){
    if(i==n+1){
        printclique(C);
        return;
    }
    unordered_set<int>diff=setdiff(C,adj[i]);
    if(!diff.empty()) UPDATE(i+1,C,n,adj,S,T);
    for(int y=1;y<=n;y++) {
        if(C.find(y)==C.end()&&y!=i){
            unordered_set<int> inter=setintersect(adj[y],diff);
            T[y]=inter.size();
        }
    }
    unordered_set<int> cap=setintersect(C,adj[i]);
    for(int x:cap){
        unordered_set<int> Nx=setdiff(adj[x],C);
        if(Nx.find(i)!=Nx.end()) Nx.erase(i);
        for(int y:Nx){
            unordered_set<int> inter=setintersect(adj[y],diff);
            T[y]+=inter.size();
        }
    }
    for(int x:diff){
        unordered_set<int> Nx=setdiff(adj[x],C);
        for(int y:Nx){
            S[y]++;
        }
    }
    bool FLAG=true;
    int capSize=cap.size(); 
    for(int y:adj[i]){
        if(C.find(y)==C.end()&&y<i&&T[y]==capSize){
            FLAG=false;
            break;
        }
    }
    vector<int> diffVec(diff.begin(),diff.end());
    sort(diffVec.begin(),diffVec.end()); 
    for(size_t k=0;k<diffVec.size();k++) {
        int j=diffVec[k];
        unordered_set<int> Nj=setdiff(adj[j],C);
        for(int y:Nj){
            if(y<i&&T[y]==capSize){
                if(adj[y].find(j)!=adj[y].end()) S[y]--;
                else if(k==0&&S[y]>0) FLAG=false;
            }
        }
    }
    if(!cap.empty()&&!diffVec.empty()){
        int jp=diffVec.back();
        for(int y:cap) {
            if(y<i&&T[y]==(int)adj[i].size()&&S[y]==0){
                if(jp<y||jp<i-1){
                    FLAG=false;
                    break;
                }
            }
        }
    }
    for(int x:cap){
        for(int y:adj[x]){
            if(C.find(y)==C.end()&&y!=i) T[y]=0;
        }
    }
    for(int x:diff) {
        for(int y:adj[x]){
            if(C.find(y)==C.end()) S[y]=0;
        }
    }
    if(FLAG){
        unordered_set<int> SAVE=diff;
        unordered_set<int> newC=cap;
        newC.insert(i);
        unordered_set<int> oldC=C;
        C=newC;
        UPDATE(i+1,C,n,adj,S,T);
        unordered_set<int> recovered;
        for(int x:C){
            if(x!=i) recovered.insert(x);
        }
        recovered.insert(SAVE.begin(),SAVE.end());
        C=recovered;
    }
}

int main(int argc,char* argv[]){
    ifstream infile(argv[1]);
    int n;
    infile >> n;
    vector<pair<int,int> > edges;
    int u, v;
    while(infile >> u >> v) edges.push_back(pair<int,int>(u+1,v+1));
    infile.close();
    int maxVertex=n;
    vector<unordered_set<int> > origAdj(maxVertex+1);
    for (size_t i=0;i<edges.size();i++) {
        pair<int,int> e=edges[i];
        origAdj[e.first].insert(e.second);
        origAdj[e.second].insert(e.first);
    }
    vector<pair<int,int> > degList; 
    for (int i=1;i<=n;i++) {
        degList.push_back(pair<int,int>((int)origAdj[i].size(),i));
    }
    sort(degList.begin(), degList.end());
    vector<int> oldToNew(n+1,0);
    newToOld.resize(n+1,0);
    for (int i=0;i<n;i++) {
        int oldId=degList[i].second;
        int newId=i+1; 
        oldToNew[oldId]=newId;
        newToOld[newId]=oldId;
    }
    vector<unordered_set<int> > adj(n+1);
    for (int i=1;i<=n;i++){
        int newU=oldToNew[i];
        for (auto it=origAdj[i].begin();it!=origAdj[i].end();++it){
            int neighbor=*it;
            if (neighbor<=n){
                int newV=oldToNew[neighbor];
                adj[newU].insert(newV);
            }
        }
    }
    vector<int> S(n+1,0),T(n+1,0);
    unordered_set<int> C;
    C.insert(1);
    UPDATE(2, C, n, adj, S, T);
    return 0;
}