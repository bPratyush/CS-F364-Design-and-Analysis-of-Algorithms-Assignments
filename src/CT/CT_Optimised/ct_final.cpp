#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <sys/resource.h>
#include <utility> 
#include <set>
using namespace std;
using namespace chrono;

// Chiba and Nishizeki (1985) algorithm for finding all maximal cliques in an undirected graph

static vector<int> T_, S_;
static vector<int> cliqueSizeCount;
static int numCliques = 0;
static int maxCliqueSize = 0;
static vector<unordered_set<int>> adj;
static int n = 0;
static bool FLAG = true; 

static inline void ensureCliqueSize(int sz) {
    if (sz > maxCliqueSize)
        maxCliqueSize = sz;
    if ((int)cliqueSizeCount.size() <= sz)
        cliqueSizeCount.resize(sz + 1, 0);
}

static inline int getMax(const unordered_set<int>& C) {
    int mx = -1;
    for (int x : C)
        if (x > mx)
            mx = x;
    return mx;
}

static inline void recordClique(const unordered_set<int>& C) {
    vector<int> cliqueVec;
    cliqueVec.reserve(C.size());
    for (int x : C)
        cliqueVec.push_back(x);
    sort(cliqueVec.begin(), cliqueVec.end());
    int sz = static_cast<int>(cliqueVec.size());
    if (sz < 2)
        return;
    ensureCliqueSize(sz);

    ++cliqueSizeCount[sz];
    ++numCliques;
    cout << "[DEBUG] Total cliques so far: " << numCliques << endl;
}
 
static void UPDATE(int i, unordered_set<int>& C) {
    if (!C.empty()) {
        int m = getMax(C);
        if (i <= m) {
            UPDATE(m + 1, C);
            return;
        }
    }
    
    if (i == n) {
        recordClique(C);
        return;
    }

    const auto& neighborsOfI = adj[i];
    unordered_set<int> diff;
    diff.reserve(C.size());
    for (int x : C) {
        if (!neighborsOfI.count(x))
            diff.insert(x);
    }

    if (!diff.empty())
        UPDATE(i + 1, C);
    
    unordered_set<int> cap;
    cap.reserve(C.size());
    for (int x : C) {
        if (neighborsOfI.count(x))
            cap.insert(x);
    }

    for (int x : cap) {
        for (int y : adj[x]) {
            if ((y != i) && !C.count(y))
                ++T_[y];
        }
    }
    for (int x : diff) { 
        for (int y : adj[x]) {
            if (!C.count(y))
                ++S_[y];
        }
    }
    
    FLAG = true;
    const int capSize = static_cast<int>(cap.size());

    for (int y : neighborsOfI) {
        if (!C.count(y) && (y < i) && (T_[y] == capSize)) {
            FLAG = false;
            break;
        }
    }

    vector<int> dv;
    dv.reserve(diff.size());
    for (int x : diff)
        dv.push_back(x);
    sort(dv.begin(), dv.end());

    const int p_int = static_cast<int>(dv.size());
    for (int k = 0; k < p_int && FLAG; k++) {
        int j_k = dv[k]; 
        for (int y : adj[j_k]) {
            if ((y < i) && !C.count(y) && (T_[y] == capSize)) {
                if (y >= j_k) {
                    S_[y]--;
                } else { 
                    if (k == 0 || y >= dv[k - 1]) {
                        if ((S_[y] + k) == p_int) {
                            FLAG = false;
                            break;  
                        }
                    }
                }
            }
        }
    }

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
        if (jp < (i - 1))
            FLAG = false;
    }
    
    for (int x : cap) {
        for (int y : adj[x]) {
            if ((y != i) && !C.count(y))
                T_[y] = 0;
        }
    }
    for (int x : diff) {  
        for (int y : adj[x]) {
            if (!C.count(y))
                S_[y] = 0;
        }
    }
    
    if (FLAG) {
        auto diffSave = std::move(diff);
        cap.insert(i);
        auto oldC = std::move(C);
        C = std::move(cap);
        UPDATE(i + 1, C);
        C = std::move(oldC);
        for (int val : diffSave)
            C.insert(val);
    }
}

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
        if (rl.rlim_max < rl.rlim_cur)
            rl.rlim_max = rl.rlim_cur;
        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0)
            cerr << "Error setting stack limit: " << strerror(errno) << endl;
        else
            cerr << "Stack size increased to " << (stack_size / (1024 * 1024)) << " MB" << endl;
    }
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " [input file]" << "\n";
        return 1;
    }

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
        if (line.empty() || line[0] == '#')
            continue;
        istringstream iss(line);
        if (iss >> u >> v) {
            maxVertex = max(maxVertex, max(u, v));
            edges.emplace_back(u, v);
        }
    }
    infile.close();
    n = maxVertex + 1;
    adj.clear();
    adj.resize(n);

    for (auto &p : edges) {
        int a = p.first, b = p.second;
        if (a != b && a >= 0 && a < n && b >= 0 && b < n) {
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }

    T_.assign(n, 0);
    S_.assign(n, 0);
    
    vector<pair<int,int>> degreeIndex;
degreeIndex.reserve(n);
for (int i = 0; i < n; i++) {
    degreeIndex.push_back({static_cast<int>(adj[i].size()), i});
}

sort(degreeIndex.begin(), degreeIndex.end());
vector<int> oldToNew(n), newToOld(n);
for (int newIndex = 0; newIndex < n; newIndex++) {
    int oldIndex = degreeIndex[newIndex].second;
    oldToNew[oldIndex] = newIndex;
    newToOld[newIndex] = oldIndex;
}
vector<unordered_set<int>> newAdj(n);
for (int i = 0; i < n; i++) {
    int new_i = oldToNew[i];
    for (int neighbor : adj[i]) {
        int newNeighbor = oldToNew[neighbor];
        if (new_i != newNeighbor)
            newAdj[new_i].insert(newNeighbor);
    }
}
adj = std::move(newAdj);
int stidx = 0;
while(stidx < degreeIndex.size() && degreeIndex[stidx].first == 0) {
    stidx++;
}
    unordered_set<int> C;
    C.insert(0);
    
    auto startTime = high_resolution_clock::now();
    UPDATE(1, C);
    auto stopTime = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stopTime - startTime).count();
    ofstream outfile("output.txt");
    if (!outfile) {
        cerr << "Error: Unable to open output.txt for writing\n";
        return 1;
    }
    outfile << "Largest size of a clique: " << maxCliqueSize << "\n";
    outfile << "Total number of maximal cliques: " << numCliques << "\n";
    outfile << "Execution time (ms): " << duration << "\n";
    outfile << "Distribution of clique sizes:\n";
    for (int i = 2; i <= maxCliqueSize; ++i) {
        int countVal = (i < (int)cliqueSizeCount.size()) ? cliqueSizeCount[i] : 0;
        outfile << "Size " << i << ": " << countVal << "\n";
    }
    outfile.close();
    return 0;
}