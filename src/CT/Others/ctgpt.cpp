//// filepath: /Users/bpratyush/Documents/GitHub/CS-F364-Design-and-Analysis-of-Algorithms-Assignments/src/CT/CLIQUE.cpp
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <numeric>
using namespace std;
using namespace std::chrono;

static vector<int> S, T;          // S[y], T[y] arrays from pseudocode
static bool FLAG;                 // to track if the clique is valid & lexicographically largest
static int n;                     // number of vertices
static vector<unordered_set<int>> adj;  // adjacency list
static vector<int> numbering;     // order of vertices by non-decreasing degree
static vector<int> degree; 
// Helper: Print the clique
void printClique(const vector<int>& C, ofstream &outfile) {
    outfile << "Found clique of size " << C.size() << ": ";
    for (auto &v : C) outfile << v << " ";
    outfile << endl;
}

// UPDATE(i,C) from your pseudocode
void UPDATE(int i, vector<int>& C, ofstream &outfile) {
    // If we've passed the last vertex, print the clique
    if (i == n) {
        printClique(C, outfile);
        return;
    }

    // 1: If (C - N(i)) is not empty => call UPDATE(i+1, C)
    // numbering[i] is the actual vertex 'cur'
    int cur = numbering[i];
    bool hasOutside = false;
    for (auto &vx : C) {
        if (adj[cur].find(vx) == adj[cur].end()) {
            hasOutside = true;
            break;
        }
    }
    if (!hasOutside) {
        UPDATE(i + 1, C, outfile);
    }

    // Prepare T[y] = |N(y) ∩ C ∩ N(cur)| and S[y] = |N(y) ∩ (C - N(cur))|
    // (2 & 3 in pseudocode)
    // We assume V - C - {cur} = all other vertices not in C or cur
    // Pseudocode references each y in that set:
    vector<bool> inC(n, false);
    for (auto &vx : C) {
        inC[vx] = true;
    }
    for (int y = 0; y < n; y++) {
        if (!inC[y] && y != cur) {
            T[y] = 0;
            S[y] = 0;
        }
    }
    // 2: Increase T[y] for each y in N(x)-C-{cur}, x in C ∩ N(cur)
    for (auto &x : C) {
        if (adj[cur].find(x) != adj[cur].end()) {
            // x ∈ C ∩ N(cur)
            for (auto &y : adj[x]) {
                if (!inC[y] && y != cur) {
                    T[y]++;
                }
            }
        }
    }
    // 3: Increase S[y] for each y in N(x)-C, x in C - N(cur)
    for (auto &x : C) {
        if (adj[cur].find(x) == adj[cur].end()) {
            // x ∈ C - N(cur)
            for (auto &y : adj[x]) {
                if (!inC[y]) {
                    S[y]++;
                }
            }
        }
    }

    // 4: Maximality test
    // if ∃ y ∈ N(cur)-C s.t. y<i and T[y] = |C ∩ N(cur)| => FLAG=false
    FLAG = true;
    int sizeCcapNcur = 0;
    for (auto &vx : C) {
        if (adj[cur].find(vx) != adj[cur].end()) {
            sizeCcapNcur++;
        }
    }
    // N(cur)-C = all v in adj[cur] not in C
    for (auto &y : adj[cur]) {
        if (!inC[y]) {
            // y < i? (pseudocode uses y<i, but we have 0-based vs 1-based indexing.)
            // We'll interpret as: numbering[pos[y]] < numbering[i]? Or just y < numbering[i]?
            // Below is a simple check consistent with pseudocode's "vertex label" notion:
            if (y < cur && T[y] == sizeCcapNcur) {
                FLAG = false;
                break;
            }
        }
    }

    // 5: sort all vertices in C - N(cur) in ascending order j1<j2<...<jp
    // (We have them in the iteration above, but let's do a fresh vector to follow logic.)
    vector<int> cMinusNcur;
    for (auto &vx : C) {
        if (adj[cur].find(vx) == adj[cur].end()) {
            cMinusNcur.push_back(vx);
        }
    }
    sort(cMinusNcur.begin(), cMinusNcur.end());
    int p = (int)cMinusNcur.size();

    // 6,7: Lexicographic test (cases S[y]>=1 and S[y]=0)
    // This is quite detailed in the pseudocode; we do a simplified check:
    for (int k = 0; k < p && FLAG; k++) {
        int jk = cMinusNcur[k];
        // for each y in N(jk)-C where y < cur and T[y]=|C ∩ N(cur)|
        // attempt to see if we should decrement S[y] or break FLAG
        for (auto &y : adj[jk]) {
            if (!inC[y] && y != cur && y < cur) {
                if (T[y] == sizeCcapNcur) {
                    // if y >= jk => S[y]-- else do extended checks
                    if (y >= jk) {
                        S[y] = max(0, S[y] - 1);
                    } else {
                        // Check if (S[y]+k == p) etc. => FLAG=false
                        // (Simplified version below; adapt as needed to fully match the pseudocode)
                        if ((S[y] + (k)) == p) {
                            FLAG = false;
                            break;
                        }
                    }
                }
            }
        }
    }
    // 7 continued: if (C ∩ N(cur)) not empty => more checks for y<i, T[y]=sizeCcapNcur, S[y]=0
    if (sizeCcapNcur > 0 && FLAG) {
        for (int y = 0; y < n; y++) {
            if (!inC[y] && y != cur && y < cur) {
                if (T[y] == sizeCcapNcur && S[y] == 0) {
                    // compare y with the last in cMinusNcur => if out of lex range => FLAG=false
                    if (!cMinusNcur.empty() && cMinusNcur.back() < y) {
                        FLAG = false;
                        break;
                    }
                    // another condition with i => if cMinusNcur.back() < i-1 => FLAG=false
                    if (cMinusNcur.empty() || cMinusNcur.back() < (i - 1)) {
                        FLAG = false;
                        break;
                    }
                }
            }
            if (!FLAG) break;
        }
    }

    // 8 & 9: Reinitialize S[y], T[y] after usage
    for (auto &x : C) {
        if (adj[cur].find(x) != adj[cur].end()) {
            for (auto &y : adj[x]) {
                if (!inC[y] && y != cur) {
                    T[y] = 0;
                }
            }
        }
    }
    for (auto &x : C) {
        if (adj[cur].find(x) == adj[cur].end()) {
            for (auto &y : adj[x]) {
                if (!inC[y]) {
                    S[y] = 0;
                }
            }
        }
    }

    // 10: If FLAG => (C ∩ N(cur)) union {cur} => new clique
    if (FLAG) {
        // SAVE = C - N(cur)
        vector<int> SAVE;
        for (auto &vx : C) {
            if (adj[cur].find(vx) == adj[cur].end()) {
                SAVE.push_back(vx);
            }
        }
        // C = (C ∩ N(cur)) ∪ {cur}
        vector<int> newC;
        for (auto &vx : C) {
            if (adj[cur].find(vx) != adj[cur].end()) {
                newC.push_back(vx);
            }
        }
        newC.push_back(cur);

        UPDATE(i + 1, newC, outfile);

        // restore C = (C - {cur}) ∪ SAVE
        // effectively revert the changes
        C.clear();
        for (auto &vx : newC) {
            if (vx != cur) C.push_back(vx);
        }
        for (auto &sv : SAVE) {
            C.push_back(sv);
        }
    }
}

// The main CLIQUE driver
void CLIQUE(ofstream &outfile) {
    // Sort vertices by non-decreasing degree
    // (like d(1) <= d(2) <= ... per pseudocode)
    iota(numbering.begin(), numbering.end(), 0);
    sort(numbering.begin(), numbering.end(), [&](int a, int b){
        // compare degree first; if tie, compare actual IDs
        if (degree[a] == degree[b]) return a < b;
        return degree[a] < degree[b];
    });

    // Initialize S, T to zero
    S.assign(n, 0);
    T.assign(n, 0);

    // Start with C = {1} i.e. numbering[0] in 0-based
    vector<int> C;
    C.push_back(numbering[0]);
    UPDATE(1, C, outfile);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    ifstream infile(argv[1]);
    ofstream outfile("output.txt");
    if (!infile.is_open()) {
        outfile << "Could not open input file." << endl;
        return 0;
    }

    vector<pair<int,int>> edgeList;
    int maxVertex = -1;
    {
        string line;
        while (getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            istringstream iss(line);
            int u, v;
            if (iss >> u >> v) {
                edgeList.push_back({u, v});
                maxVertex = max(maxVertex, max(u, v));
            }
        }
    }
    infile.close();

    // Build adjacency
    n = maxVertex + 1;
    adj.assign(n, {});
    degree.assign(n, 0);
    numbering.assign(n, 0);

    for (auto &e : edgeList) {
        int u = e.first, v = e.second;
        if(u<n && v<n) {
            adj[u].insert(v);
            adj[v].insert(u);
            degree[u]++;
            degree[v]++;
        }
    }

    auto start = high_resolution_clock::now();
    CLIQUE(outfile);
    auto end = high_resolution_clock::now();
    auto dur = duration_cast<milliseconds>(end - start).count();

    outfile << "Execution time (ms): " << dur << endl;
    outfile.close();
    return 0;
}