#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

int n, m;
vector<vector<int>> adj;  // adjacency list
vector<int> deg;          // degrees for vertex ordering
vector<int> order;        // order of vertices by degree
vector<int> S, T;         // helper arrays for tests
int cliqueCount = 0;
int maxCliqueSize = 0;

void printClique(const vector<int>& C) {
    cliqueCount++;
    int sz = C.size();
    maxCliqueSize = max(maxCliqueSize, sz);
    cout << "Clique " << cliqueCount << ": ";
    for (auto v : C)
        cout << v << " ";
    cout << endl;
}

// UPDATE procedure from the pseudocode
void UPDATE(int i, vector<int>& C) {
    if (i == n + 1) {
        printClique(C);
        return;
    }

    // Step 1: recursive call if C is already a clique of G_i
    bool subset = true;
    for (auto v : C) {
        if (find(adj[i].begin(), adj[i].end(), v) == adj[i].end()) {
            subset = false;
            break;
        }
    }
    if (subset)
        UPDATE(i + 1, C);

    // Prepare for tests
    // Compute T[y] = |N(y) ∩ C ∩ N(i)|
    fill(T.begin(), T.end(), 0);
    for (auto x : C) {
        if (find(adj[i].begin(), adj[i].end(), x) != adj[i].end()) {
            for (auto y : adj[x]) {
                if (y != i && find(C.begin(), C.end(), y) == C.end()) {
                    T[y]++;
                }
            }
        }
    }

    // Compute S[y] = |N(y) ∩ (C - N(i))|
    fill(S.begin(), S.end(), 0);
    for (auto x : C) {
        if (find(adj[i].begin(), adj[i].end(), x) == adj[i].end()) {
            for (auto y : adj[x]) {
                if (find(C.begin(), C.end(), y) == C.end())
                    S[y]++;
            }
        }
    }

    // Maximality test
    bool FLAG = true;
    for (auto y : adj[i]) {
        if (find(C.begin(), C.end(), y) == C.end() && y < i && T[y] == count_if(C.begin(), C.end(), [&](int v) {
            return find(adj[i].begin(), adj[i].end(), v) != adj[i].end();
        })) {
            FLAG = false;
            break;
        }
    }

    // Lexicographical test
    vector<int> CNi;
    for (auto v : C) {
        if (find(adj[i].begin(), adj[i].end(), v) != adj[i].end()) {
            CNi.push_back(v);
        }
    }
    sort(C.begin(), C.end());
    int p = C.size() - CNi.size();

    // Simplified for illustration, full test per paper can be added if needed
    if (!FLAG) {
        return;
    }

    // Reinitialize S and T (already done by fill())

    if (FLAG) {
        vector<int> SAVE;
        for (auto v : C) {
            if (find(adj[i].begin(), adj[i].end(), v) == adj[i].end()) {
                SAVE.push_back(v);
            }
        }

        // C = (C ∩ N(i)) ∪ {i}
        vector<int> newC;
        for (auto v : C) {
            if (find(adj[i].begin(), adj[i].end(), v) != adj[i].end()) {
                newC.push_back(v);
            }
        }
        newC.push_back(i);
        UPDATE(i + 1, newC);

        // C = old C after recursive call (implicit by copy in newC)
    }
}

void CLIQUE() {
    // Order vertices by non-decreasing degree
    order.resize(n + 1);
    for (int i = 1; i <= n; i++)
        order[i] = i;

    sort(order.begin() + 1, order.end(), [&](int u, int v) {
        return deg[u] < deg[v];
    });

    S.resize(n + 1, 0);
    T.resize(n + 1, 0);

    vector<int> C;  // start with empty clique
    UPDATE(2, C);
}

int main() {
    cout << "Enter number of vertices and edges: ";
    cin >> n >> m;

    adj.resize(n + 1);
    deg.resize(n + 1, 0);

    cout << "Enter edges (u v):" << endl;
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
        deg[u]++;
        deg[v]++;
    }

    CLIQUE();

    cout << "\nTotal maximal cliques found: " << cliqueCount << endl;
    cout << "Size of the largest clique: " << maxCliqueSize << endl;

    return 0;
}