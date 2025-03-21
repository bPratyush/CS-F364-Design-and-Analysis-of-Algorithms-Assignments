//// filepath: /Users/bpratyush/Documents/GitHub/CS-F364-Design-and-Analysis-of-Algorithms-Assignments/src/CT/degeneracycalc.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>
using namespace std;
using ll = long long;

unordered_map<ll, unordered_set<ll>> adj;

// Compute the degeneracy of a graph following the standard bucket-based approach
ll compute_degeneracy(unordered_set<ll> &nodes) {
    // deg[node] holds the current degree of node
    unordered_map<ll, ll> deg;
    // deglist[d] holds all nodes of degree d
    unordered_map<ll, unordered_set<ll>> deglist;

    // Initialize degrees
    for (auto node : nodes) {
        ll d = adj[node].size();
        deg[node] = d;
        deglist[d].insert(node);
    }

    ll max_degeneracy = 0;
    ll minDegree = 0;

    // Copy of adjacency for safe neighbor updates
    unordered_map<ll, unordered_set<ll>> tempadj = adj;

    // We remove nodes one by one, always taking from the smallest degree bucket
    while (!nodes.empty()) {
        // Advance minDegree until we find a non-empty bucket
        while (minDegree < (ll)nodes.size() && deglist[minDegree].empty()) {
            minDegree++;
        }
        // If we've exceeded the range, we are done
        if (minDegree >= (ll)nodes.size()) {
            break;
        }

        // Pick a node from the current smallest degree bucket
        ll u = *deglist[minDegree].begin();
        deglist[minDegree].erase(u);
        nodes.erase(u);

        // Record the degeneracy so far
        max_degeneracy = max(max_degeneracy, minDegree);

        // Decrement the degree of each neighbor
        for (auto v : tempadj[u]) {
            tempadj[v].erase(u);
            deglist[deg[v]].erase(v);
            deg[v]--;
            deglist[deg[v]].insert(v);
        }
    }

    return max_degeneracy;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream ip("email-Enron.txt");
    if (!ip.is_open()) {
        cerr << "Error: unable to open file.\n";
        return 1;
    }

    unordered_set<ll> nodes;
    string line;

    // Read edges
    while (getline(ip, line)) {
        if (line.empty() || line[0] == '#') continue;
        ll a, b;
        istringstream ss(line);
        if (ss >> a >> b) {
            adj[a].insert(b);
            adj[b].insert(a);
            nodes.insert(a);
            nodes.insert(b);
        }
    }
    ip.close();

    cout << "Adjacency list created.\n";

    auto start = chrono::high_resolution_clock::now();
    ll degeneracy = compute_degeneracy(nodes);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    cout << "Degeneracy of the graph: " << degeneracy << "\n";
    cout << "Time taken (ms): " << duration << "\n";

    return 0;
}