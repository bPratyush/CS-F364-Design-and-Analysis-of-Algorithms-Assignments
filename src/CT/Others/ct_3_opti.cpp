#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <set>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <sys/resource.h>
#include <utility> // for std::move

using namespace std;
using namespace chrono;

// ------------------------------------------------------------
// Global variables
// ------------------------------------------------------------
static vector<int> T_, S_;
// Instead of an unordered_map<int, int>, use a direct vector for clique-size counts.
static vector<int> cliqueSizeCount;
static int numCliques    = 0;
static int maxCliqueSize = 0;
static vector<unordered_set<int>> adj;
static int n             = 0;
static bool FLAG         = true;
// Store printed cliques in a set of sorted vectors.
static set<vector<int>> printedCliques;

// ------------------------------------------------------------
// Updates maxCliqueSize and ensures cliqueSizeCount has enough size
// ------------------------------------------------------------
static inline void ensureCliqueSize(int sz) {
    if (sz > maxCliqueSize) {
        maxCliqueSize = sz;
    }
    if (static_cast<int>(cliqueSizeCount.size()) <= sz) {
        cliqueSizeCount.resize(sz + 1, 0);
    }
}

// ------------------------------------------------------------
// Convert an unordered_set clique into a sorted vector and record it
// ------------------------------------------------------------
static inline void recordClique(const unordered_set<int>& C) {
    // Convert to sorted vector
    vector<int> cliqueVec;
    cliqueVec.reserve(C.size());
    for (int x : C) {
        cliqueVec.push_back(x);
    }
    sort(cliqueVec.begin(), cliqueVec.end());

    int sz = static_cast<int>(cliqueVec.size());
    ensureCliqueSize(sz);

    // Insert if not already present and print the updated clique count.
    if (!printedCliques.count(cliqueVec)) {
        ++cliqueSizeCount[sz];
        printedCliques.insert(std::move(cliqueVec));
        ++numCliques;
        cout << "[DEBUG] Total cliques so far: " << numCliques << endl;
    }
}

// ------------------------------------------------------------
// Recursive function exploring node i with the current clique C
// ------------------------------------------------------------
static void UPDATE(int i, unordered_set<int>& C) {
    // cout << "[DEBUG] Processing node " << i << ", current clique: ";
    // for (int x : C)
    //     cout << x << " ";
    // cout << endl;
    
    if (i == n) {
        // cout << "[DEBUG] Final clique at node " << i << ": ";
        // for (int x : C)
        //     cout << x << " ";
        // cout << endl;
        recordClique(C);
        return;
    }

    // diff = elements in C that are NOT neighbors of i
    const auto& neighborsOfI = adj[i];
    unordered_set<int> diff;
    diff.reserve(C.size());
    for (int x : C) {
        if (!neighborsOfI.count(x)) {
            diff.insert(x);
        }
    }
    if (!diff.empty()) {
       
        UPDATE(i + 1, C);
    }

    // cap = elements in C that ARE neighbors of i
    unordered_set<int> cap;
    cap.reserve(C.size());
    for (int x : C) {
        if (neighborsOfI.count(x)) {
            cap.insert(x);
        }
    }


    // Update T_ and S_
    for (int x : cap) {
        for (int y : adj[x]) {
            if ((y != i) && !C.count(y)) {
                ++T_[y];
            }
        }
    }
    for (int x : diff) {
        for (int y : adj[x]) {
            if (!C.count(y)) {
                ++S_[y];
            }
        }
    }

    FLAG = true;
    const int capSize = static_cast<int>(cap.size());

    // Preliminary maximality check
    for (int y : neighborsOfI) {
        if (!C.count(y) && (y < i) && (T_[y] == capSize)) {
            
            FLAG = false;
            break;
        }
    }

    // Sort diff for lexicographic tests
    vector<int> dv;
    dv.reserve(diff.size());
    for (int x : diff) {
        dv.emplace_back(x);
    }
    sort(dv.begin(), dv.end());

    // Lexicographic checks - Step 6.
    const int p_int = static_cast<int>(dv.size());
    for (int k = 0; k < p_int && FLAG; k++) {
        int j_k = dv[k]; // current vertex from sorted diff (C - N(i))
        
        for (int y : adj[j_k]) {
            if ((y < i) && !C.count(y) && (T_[y] == capSize)) {
                
                if (y >= j_k) {
                    
                    S_[y] = std::max(0, S_[y] - 1);
                   
                } else {
                    if ((k == 0 || y < dv[k - 1]) && ((S_[y] + k) == p_int) && (y >= (j_k - 1))) {
                        
                        FLAG = false;
                        break;
                    }
                }
            }
        }
    }

    // Lexico test - Step 7.
    int jp = dv.empty() ? 0 : dv.back();

    if (!cap.empty()) {
        for (int y = 0; y < i && FLAG; y++) {
            if (!C.count(y) && (T_[y] == capSize) && (S_[y] == 0)) {
                
                if (jp < y) {
                    
                    FLAG = false;
                    break;
                }
            }
        }
    } else {
        if (jp < (i - 1)) {

            FLAG = false;
        }
    }

    // Reset T_, S_ for nodes in adjacency of cap
    for (int x : cap) {
        for (int y : adj[x]) {
            if ((y != i) && !C.count(y)) {
                T_[y] = 0;
            }
        }
    }
    for (int x : cap) {
        for (int y : adj[x]) {
            if (!C.count(y)) {
                S_[y] = 0;
            }
        }
    }

    // Recurse if feasible: add i into cap -> new clique => update C, then restore C
    if (FLAG) {
        auto diffSave = std::move(diff);
        cap.insert(i);

        auto oldC = std::move(C);
        C = std::move(cap);

        UPDATE(i + 1, C);

        // Restore C
        C = std::move(oldC);
        for (int val : diffSave) {
            C.insert(val);
        }
    } 
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main(int argc, char* argv[]) {
    
    struct rlimit rl;
    rlim_t stack_size = 512 * 1024 * 1024;
    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result != 0) {
        cerr << "Error getting stack limit: " << strerror(errno) << endl;
        return 1;
    }
    if (rl.rlim_cur < stack_size) {
        rl.rlim_cur = stack_size;
        if (rl.rlim_max < rl.rlim_cur) {
            rl.rlim_max = rl.rlim_cur; 
        }

        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0) {
            cerr << "Error setting stack limit: " << strerror(errno) << endl;
        } else {
            cerr << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << endl;
        }
    }
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [input file]" << "\n";
        return 1;
    }

    // Read input edges
    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error: Unable to open input file " << argv[1] << "\n";
        return 1;
    }

    int u, v;
    int maxVertex = -1;
    vector<pair<int,int>> edges;
    edges.reserve(200000);
    string line;
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        istringstream iss(line);
        if (iss >> u >> v) {
            maxVertex = std::max(maxVertex, std::max(u, v));
            edges.emplace_back(u, v);
        }
    }
    infile.close();

    n = maxVertex + 1;
    adj.clear();
    adj.resize(n);

    // Build adjacency list
    for (auto &p : edges) {
        int a = p.first, b = p.second;
        if (a != b && a >= 0 && a < n && b >= 0 && b < n) {
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }

    // Initialize T_ and S_
    T_.assign(n, 0);
    S_.assign(n, 0);

    // Start with the clique {0}
    unordered_set<int> C;
    C.insert(0);

    auto startTime = high_resolution_clock::now();
    UPDATE(1, C);
    auto stopTime = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stopTime - startTime).count();

    // Output results
    ofstream outfile("output_3.txt");
    if (!outfile) {
        cerr << "Error: Unable to open output_3.txt for writing\n";
        return 1;
    }

    outfile << "Largest size of a clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << numCliques << "\n";
    outfile << "Execution time (ms): " << duration << "\n";
    outfile << "Distribution of clique sizes:\n";
    for (int i = 1; i <= maxCliqueSize; ++i) {
        int countVal = (i < (int)cliqueSizeCount.size()) ? cliqueSizeCount[i] : 0;
        outfile << "Size " << i << ": " << countVal << "\n";
    }
    outfile.close();
    return 0;
}