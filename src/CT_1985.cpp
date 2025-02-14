#include <bits/stdc++.h>
using namespace std;

vector<int> T_, S_;
set<vector<int>> printedCliques; 

void printClique(const unordered_set<int> &C) {
    vector<int> v(C.begin(), C.end());
    sort(v.begin(), v.end());

    if (printedCliques.count(v) == 0) {
        cout << "Clique: ";
        for (int x : v) cout << x << " ";
        cout << "\n";
        
       
        printedCliques.insert(v);
    }
}

unordered_set<int> setdiff(const unordered_set<int> &A, const unordered_set<int> &B) {
    unordered_set<int> r;
    for (int x : A) if (!B.count(x)) r.insert(x);
    return r;
}

unordered_set<int> setintersect(const unordered_set<int> &A, const unordered_set<int> &B) {
    unordered_set<int> r;
    if (A.size() < B.size()) {
        for (int x : A) if (B.count(x)) r.insert(x);
    } else {
        for (int x : B) if (A.count(x)) r.insert(x);
    }
    return r;
}

void UPDATE(int i, unordered_set<int> &C, int n, const vector<unordered_set<int>> &adj) {
    if (i == n) {
        printClique(C);
        return;
    }

    unordered_set<int> diff = setdiff(C, adj[i]);
    if (!diff.empty()) UPDATE(i + 1, C, n, adj);

    for (int y = 0; y < n; y++) {
        if (!C.count(y)) {
            T_[y] = 0;
            S_[y] = 0;
        }
    }

    unordered_set<int> cap = setintersect(C, adj[i]);

    for (int x : cap) {
        for (int y : adj[x]) {
            if (y != i && !C.count(y)) T_[y]++;
        }
    }
    for (int x : diff) {
        for (int y : adj[x]) {
            if (!C.count(y)) S_[y]++;
        }
    }

    bool FLAG = true;
    int cs = (int)cap.size();
    for (int y : adj[i]) {
        if (!C.count(y) && y < i && T_[y] == cs) {
            FLAG = false;
            break;
        }
    }

    vector<int> dv(diff.begin(), diff.end());
    sort(dv.begin(), dv.end());

    for (int j : dv) {
        if (!FLAG) break;
        for (int y : adj[j]) {
            if (y < i && !C.count(y) && T_[y] == cs) {
                if (adj[y].count(j)) S_[y]--;
                else {
                    FLAG = false;
                    break;
                }
            }
        }
    }

    if (!cap.empty() && !dv.empty()) {
        int jb = dv.back();
        for (int y : cap) {
            if (y < i && (int)adj[i].size() == T_[y] && S_[y] == 0) {
                if (jb < y || jb < i - 1) {
                    FLAG = false;
                    break;
                }
            }
        }
    }

    for (int x : cap) {
        for (int y : adj[x]) {
            if (!C.count(y) && y != i) T_[y] = 0;
        }
    }
    for (int x : diff) {
        for (int y : adj[x]) {
            if (!C.count(y)) S_[y] = 0;
        }
    }

    if (FLAG) {
        unordered_set<int> sv = diff;
        unordered_set<int> nc = cap;
        nc.insert(i);
        unordered_set<int> oc = C;
        C = nc;
        UPDATE(i + 1, C, n, adj);
        C = oc;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<unordered_set<int>> adj(n);
    int u, v;
    while (cin >> u >> v) {
        if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
            adj[u].insert(v);
            adj[v].insert(u);
        }
    }

    T_.resize(n, 0);
    S_.resize(n, 0);
    unordered_set<int> C;
    C.insert(0);
    UPDATE(1, C, n, adj);

    return 0;
}
