#include <bits/stdc++.h>
using namespace std;

// Global arrays to track S[y], T[y] in the style of Chiba–Nishizeki ( as named by them)
vector<int> T_, S_;

// tracking cliques have been printed to avoid duplicates.
set<vector<int>> printedCliques;

/* An extra utility function for printing cliques */
void printClique(const unordered_set<int> &C) {
    vector<int> v(C.begin(), C.end());
    sort(v.begin(), v.end());
    if (printedCliques.count(v) == 0) {
        cout << "Clique: ";
        for (int x : v) {
            cout << x << " ";
        }
        cout << "\n";
        printedCliques.insert(v);
    }
}

/* An extra utility function for finding set difference */
unordered_set<int> setdiff(const unordered_set<int> &A, const unordered_set<int> &B) {
    unordered_set<int> r;
    for (int x : A) {
        if (!B.count(x)) {
            r.insert(x);
        }
    }
    return r;
}

/* An extra utility function for finding set intersection */
unordered_set<int> setintersect(const unordered_set<int> &A, const unordered_set<int> &B) {
    unordered_set<int> r;
    if (A.size() < B.size()) {
        for (int x : A) {
            if (B.count(x)) {
                r.insert(x);
            }
        }
    } else {
        for (int x : B) {
            if (A.count(x)) {
                r.insert(x);
            }
        }
    }
    return r;
}

/*
  Variables parallel to the pseudo code
  pseudocode | mycode
    C            C
    N           adj
    T            T
    S           S
*/

void UPDATE(int i, unordered_set<int> &C, int n, const vector<unordered_set<int>> &adj) {
    if (i == n) { // if i = n + 1
        printClique(C); // then print out a new clique
        return;
    }
    vector<int> oldT = T_;
    vector<int> oldS = S_;
    /* You'll understand below check line 168 so as to why I am creating two dummy arrays*/
    
    unordered_set<int> diff = setdiff(C, adj[i]);
    if (!diff.empty()) { 
      /*
        1: if C N(i) then UPDATE (i + 1, C);
      */
        UPDATE(i + 1, C, n, adj);
        T_ = oldT;
        S_ = oldS;
    }


   /* setting T_[y] = 0, S_[y] = 0 to start fresh
    for the "add i" checks.*/
    for (int y = 0; y < n; y++) {
        if (!C.count(y)) {
            T_[y] = 0;
            S_[y] = 0;
        }
    }

    
    unordered_set<int> cap = setintersect(C, adj[i]);

  
    for (int x : cap) { // 2: for each vertex x in (C ∩ N(i))
        for (int y : adj[x]) {
            if (y != i && !C.count(y)) { // do for each vertex y in (N(x)-C-{i})
                T_[y]++; // do T[y] := T[y] + 1;
            }
        }
    }
    for (int x : diff) { // for each vertex x C-N(i)
        for (int y : adj[x]) {
            if (!C.count(y)) { // do For each vertex y N(x)- C
                S_[y]++;  // do S[y] := S[y] + 1;
            }
        }
    }

    bool FLAG = true; // FLAG := true;
    int cs = (int)cap.size(); // size of (C ∩ N(i))
    for (int y : adj[i]) {
        if (!C.count(y) && y < i && T_[y] == cs) { // if there exists a vertex y in (N(i)-C) such that y<i and T[y] = |C ∩ N(i)|
            FLAG = false; // then FLAG := false {(Cf3 N(i))U{i} is not a clique of Gi } 
            break;
        }
    }

   // sort all the vertices in C- N(i) in ascending order j1 <j2 <...<jp, where p= |C-N(i)|;
  
    vector<int> dv(diff.begin(), diff.end());
    sort(dv.begin(), dv.end());

    // If we haven't already failed, check if we can reduce S_[y].
    // If that fails, set FLAG = false.
    for (int j : dv) {
        if (!FLAG) break;
        for (int y : adj[j]) {
            if (y < i && !C.count(y) && T_[y] == cs) {
                 // do for each vertex y in N(jk)-C such that y<i and T[y]= |C ∩ N(i)|
                if (adj[y].count(j)) { // do if y>=jk
                    S_[y]--; // then S[y]:= S[y]- {alter S[y] to S(y)}
                } else {            
                    FLAG = false;
                    break;
                    /*
                    else if (j is the first vertex which satisfies y<jk)
                    then {S[y] = S(y)}
                    if (S[y] + k -1 = p) and (y >= jk-1) {jo = 0}
		    then FLAG :=false; { C is not lexico, largest}
                    */
                }
            }
        }
    }

    if (!cap.empty() && !dv.empty()) { // 7: if C ∩ N(i) is not equal to Phi  
        int jb = dv.back(); 
        for (int y : cap) {
            if (y < i && (int)adj[i].size() == T_[y] && S_[y] == 0) {
                if (jb < y || jb < i - 1) {
                    FLAG = false;
                    break;
                    /* 
                    I am doing the following part of the pseudocode in this logic
                    
                    then for each vertex yC:Ct_J{i} such that y<i, T[y]=ICON(i) and S[y]=0
                 {access y from the adjacency list of a vertex in C fq N(i)}
                 do if jp <y then FLAG := false {C is not lexico, largest.}
                 else if jp < i-1 then FLAG :=false; { C is not lexico, largest.}     
                 */
                }
            }
        }
    }

    // 8 and 9 have been managed with the help of using dummy arrays. 
    // chiba's code tried to manually "reset" T_[y] and S_[y]
    // here for neighbors of cap & diff. But that conflicts with recursion,
    // because if I skip i or add i, I want to backtrack them differently.
    // Instead, I rely on the "back up / restore" approach at the start/end
    // of UPDATE(...). So I removed those loops.

    if (FLAG) {
        /* the code below does the 10: part of the pseudocode 
        */
        unordered_set<int> oc = C;
        unordered_set<int> nc = cap;
        nc.insert(i);      
        C = nc;
        UPDATE(i + 1, C, n, adj);
        C = oc;
    }

    // Finally, I restore T_ and S_ so that when I return to the caller,
    // I haven't permanently changed them. This is essential for not
    // corrupting other branches in the recursion.
    T_ = oldT;
    S_ = oldS;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<unordered_set<int>> adj(n);

    // Read edges until EOF
    int u, v;
    while (cin >> u >> v) {
        if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
            adj[u].insert(v);
            adj[v].insert(u);
        }
    }

    // Initialize T_ and S_ = 0 for all vertices
    T_.resize(n, 0);
    S_.resize(n, 0);

    unordered_set<int> C;
    C.insert(0);

    // Start recursion from i = 1
    UPDATE(1, C, n, adj);

    return 0;
}
