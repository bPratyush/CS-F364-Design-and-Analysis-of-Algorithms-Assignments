/*
 Finding Densest Subgraph according to a given notion of the chosen density

 This version fixes the issue for disconnected graphs by computing connected components
 and then running the densest subgraph procedure on each component separately.
*/

#include <bits/stdc++.h>
using namespace std;
 
typedef long long ll;
const ll INF = 1e15;
 
struct Edge {
    int u, v;
    ll cap, flow;
};
 
// Dinic max flow implementation
struct Dinic {
    int n, s, t;
    vector<Edge> edges;
    vector<vector<int>> adj;
    vector<int> level, ptr;
    
    Dinic(int n, int s, int t) : n(n), s(s), t(t) {
        adj.resize(n);
        level.resize(n);
        ptr.resize(n);
    }
    
    void addEdge(int u, int v, ll cap) {
        edges.push_back({u, v, cap, 0});
        edges.push_back({v, u, 0, 0});
        int m = edges.size();
        adj[u].push_back(m-2);
        adj[v].push_back(m-1);
    }
    
    bool bfs() {
        fill(level.begin(), level.end(), -1);
        level[s] = 0;
        queue<int> q;
        q.push(s);
        while(!q.empty()){
            int u = q.front();
            q.pop();
            for (int idx : adj[u]) {
                Edge &e = edges[idx];
                if(level[e.v] < 0 && e.flow < e.cap) {
                    level[e.v] = level[u] + 1;
                    q.push(e.v);
                }
            }
        }
        return level[t] >= 0;
    }
    
    ll dfs(int u, ll pushed) {
        if(pushed == 0)
            return 0;
        if(u == t)
            return pushed;
        for (int &cid = ptr[u]; cid < (int) adj[u].size(); cid++) {
            int idx = adj[u][cid];
            Edge &e = edges[idx];
            if(level[u] + 1 != level[e.v] || e.flow >= e.cap)
                continue;
            ll tr = dfs(e.v, min(pushed, e.cap - e.flow));
            if(tr == 0)
                continue;
            e.flow += tr;
            edges[idx ^ 1].flow -= tr;
            return tr;
        }
        return 0;
    }
 
    ll maxFlow() {
        ll flow = 0;
        while(bfs()){
            fill(ptr.begin(), ptr.end(), 0);
            while (ll pushed = dfs(s, INF))
                flow += pushed;
        }
        return flow;
    }
    
    void minCut(vector<bool> &cut) {
        cut.assign(n, false);
        queue<int> q;
        q.push(s);
        cut[s] = true;
        while (!q.empty()){
            int u = q.front();
            q.pop();
            for (int idx : adj[u]) {
                Edge &e = edges[idx];
                if(e.flow < e.cap && !cut[e.v]){
                    cut[e.v] = true;
                    q.push(e.v);
                }
            }
        }
    }
};
 
// Densest Subgraph Solver for edge-based density (h = 2)
struct DensestSubgraphSolver {
    int n, m; // n: number of vertices; m: number of edges in this (sub)graph
    vector<vector<int>> graph; 
    vector<int> deg;
    
    DensestSubgraphSolver(int n, int m) : n(n), m(m) {
        graph.resize(n);
        deg.assign(n, 0);
    }
    
    void addEdge(int u, int v) {
        graph[u].push_back(v);
        graph[v].push_back(u);
        deg[u]++;
        deg[v]++;
    }
 
    pair<Dinic, int> buildFlowNetwork(double alpha) {
        int N = n + 2; 
        int s = 0, t = N - 1;
        Dinic dinic(N, s, t);
        
        for (int i = 0; i < n; i++) {
            dinic.addEdge(s, i+1, m);
            ll cap = (ll)(m + 2*alpha - deg[i] + 1e-8);
            if(cap < 0) cap = 0;
            dinic.addEdge(i+1, t, cap);
        }
 
        for (int u = 0; u < n; u++) {
            for (int v : graph[u]) {
                if(u < v) { 
                    dinic.addEdge(u+1, v+1, 1);
                    dinic.addEdge(v+1, u+1, 1);
                }
            }
        }
        return {dinic, t};
    }
 
    vector<int> computeDensestSubgraph() {
        double l = 0.0, u = 0.0;
        for (int d : deg)
            u = max(u, (double)d);
        
        double eps = 1e-7;
        vector<int> bestSubgraph;
 
        while(u - l >= eps) {
            double alpha = (l + u) / 2.0;
 
            auto [dinic, t] = buildFlowNetwork(alpha);
            ll flow = dinic.maxFlow();
            vector<bool> cut;
            dinic.minCut(cut);
 
            vector<int> S;
            for (int i = 1; i <= n; i++) 
                if(cut[i])
                    S.push_back(i-1);
 
            if(S.empty())
                u = alpha;
            else { 
                l = alpha; 
                bestSubgraph = S; 
            }
        }
        return bestSubgraph;
    }
};
 
 
// solve for densest subgraph on each component and output the best one.
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
 
    ifstream fin("input7.txt");
    if (!fin.is_open()) {
        cerr << "Error opening input file\n";
        return 1;
    }
 
    int n, m;
    fin >> n >> m; 
    vector<vector<int>> fullGraph(n);
    vector<pair<int,int>> inputEdges;
    for (int i = 0; i < m; i++){
        int u, v;
        fin >> u >> v; 
        fullGraph[u].push_back(v);
        fullGraph[v].push_back(u);
        inputEdges.push_back({u, v});
    }
 
    // Compute connected components of the full graph.
    vector<int> comp(n, -1);
    int compCount = 0;
    function<void(int, int)> dfs = [&](int v, int c) {
        comp[v] = c;
        for (int w : fullGraph[v])
            if(comp[w] == -1)
                dfs(w, c);
    };
    for (int i = 0; i < n; i++){
        if(comp[i] == -1){
            dfs(i, compCount++);
        }
    }
 
    // Group vertices by component.
    vector<vector<int>> comps(compCount);
    for (int i = 0; i < n; i++){
        comps[comp[i]].push_back(i);
    }
 
    vector<int> bestGlobalSub;
    double bestGlobalDensity = -1.0;
 
    for (int c = 0; c < compCount; c++){
        if(comps[c].empty())
            continue;
 
        unordered_map<int,int> mapping;
        for (int i = 0; i < (int)comps[c].size(); i++){
            mapping[comps[c][i]] = i;
        }
 
        int comp_n = comps[c].size();
        int comp_m = 0;
        DensestSubgraphSolver solver(comp_n, 0);
        solver.graph.clear();
        solver.graph.resize(comp_n);
        solver.deg.assign(comp_n, 0);
 
        for (int u : comps[c]){
            for (int v : fullGraph[u]){
                if(mapping.count(v)) {
                    int u_new = mapping[u], v_new = mapping[v];
                    if(u_new < v_new) {
                        solver.addEdge(u_new, v_new);
                        comp_m++;
                    }
                }
            }
        }
        solver.m = comp_m;
 
        // Compute densest subgraph for this component.
        vector<int> compSub = solver.computeDensestSubgraph();
        int edgeCount = 0;
        set<int> compSubSet(compSub.begin(), compSub.end());
        for (int u : compSub) {
            for (int v : solver.graph[u]) {
                if (compSubSet.count(v))
                    edgeCount++;
            }
        }
        edgeCount /= 2;
        double density = compSub.empty() ? 0.0 : (double)edgeCount / compSub.size();
 
        if(density > bestGlobalDensity){
            bestGlobalDensity = density;
            bestGlobalSub = compSub;
            for (int &x : bestGlobalSub)
                x = comps[c][x];
        }
    }
 
    if(bestGlobalSub.empty()){
        cout << "No dense subgraph found.\n";
        return 0;
    }
 
    cout << "Densest Subgraph Vertices (" << bestGlobalSub.size() << " vertices):\n";
    for (int v : bestGlobalSub)
        cout << v << " ";
    cout << "\nEdge-density of subgraph: " << bestGlobalDensity << "\n";
    
    return 0;
}
