#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include<map>
#include <set>
#include <chrono>

using namespace std;

unordered_map<int, unordered_set<int>> graph; // Adjacency list
int v,e;
set<int> P; // Candidate vertices
int totalCliques =0;
map<int,int> sizeCounts;

vector<int> computeDegeneracyOrder(unordered_map<int, unordered_set<int>>& graph, int v, int e) {
    vector<int> order;
    order.reserve(v);

    vector<int> degree(v, 0);
    vector<bool> removed(v, false);
    vector<set<int>> buckets(v);

    // Initialize degrees and place vertices in buckets
    for (int i = 0; i < v; i++) {
        degree[i] = graph[i].size();  // Degree is 0 if the vertex has no neighbors
        buckets[degree[i]].insert(i);  // Place vertex in the bucket of its degree
    }

    int remainingVertices = v;  // Track the number of vertices left to process

    while (remainingVertices > 0) {  // Keep processing until all vertices are processed
        for (int d = 0; d < v; d++) {
            while (!buckets[d].empty()) {
                int u = *buckets[d].begin();
                buckets[d].erase(u);

                if (removed[u]) continue;

                removed[u] = true;
                order.push_back(u);
                remainingVertices--;  // Decrease count of unprocessed vertices

                // Update neighbors' degrees and re-bucket them
                for (int neighbor : graph[u]) {
                    if (!removed[neighbor]) {
                        int oldDegree = degree[neighbor];
                        buckets[oldDegree].erase(neighbor);  // Remove from old bucket

                        degree[neighbor]--;
                        buckets[degree[neighbor]].insert(neighbor);  // Move to new bucket
                    }
                }
            }
        }
    }

    return order;
}


// Bron-Kerbosch algorithm with pivoting
void BronKerbosch(set<int> R, set<int> P, set<int> X) {
    if (P.empty() && X.empty()) {
        // Found a maximal clique, write to file
        // outfile << "Clique:";
        // for (int v : R) outfile << " " << v;
        // outfile << endl;
        totalCliques++;
        sizeCounts[R.size()]++;
        return;
    }

    // Choose pivot (first element in P ∪ X)
    int u = *P.begin();
    for (int v : X) {
        u = v;
        break;
    }

    // Iterate over non-neighbors of pivot
    set<int> P_copy = P;
    for (int v : P_copy) {
        if (graph[u].count(v)) continue;

        set<int> newR = R, newP, newX;
        newR.insert(v);

        // Intersect P and X with neighbors of v
        for (int w : P) if (graph[v].count(w)) newP.insert(w);
        for (int w : X) if (graph[v].count(w)) newX.insert(w);

        BronKerbosch(newR, newP, newX);

        // Remove v from P and add to X
        P.erase(v);
        X.insert(v);
    }
}

// Load graph from file
void loadGraph(const string &filename) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }

    string line;
    bool firstLine = true;  // To read the first line with v and e

    while (getline(infile, line)) {
        if (line[0] == '#') continue;  // Skip comments

        if (firstLine) {
            // Read the first line for vertices and edges
            sscanf(line.c_str(), "%d %d", &v, &e);
            cout << "Vertices: " << v << ", Edges: " << e << endl;
            firstLine = false;
            continue;
        }

        int u, w;
        if (sscanf(line.c_str(), "%d %d", &u, &w) == 2) {
            graph[u].insert(w);
            graph[w].insert(u);  // Undirected graph: add both directions
        }
       
    }
    cout << "Graph loaded" << endl;

    infile.close();
}

void BronKerboschDegen(set<int> R, set<int> P, set<int> X)
{
    cout << "degen called" << endl;
    vector<int> order = computeDegeneracyOrder(graph, v, e);
    // cout << "Degeneracy order is:-" << endl;
    // for(int i=0; i<v;i++)
    // {
    //     cout<<order[i]<<" ";
    // }
    cout << "degen order done" << endl;
    cout<<endl;
    for(int i=0; i<v;i++)
    {
        int currentVertex = order[i];
        set<int> RForPivot = {currentVertex};
        set<int> PForPivot, XForPivot;
        //set P is intersection of neighbors of currentVeretx and the set of elements that come after currentVertex in the degeneracy ordering
        int j = i+1;
        while(j<v )
        {
            if( (graph[currentVertex].find(order[j]) != graph[currentVertex].end()))
            {PForPivot.insert(order[j]);}
            j++;
        }
        //set X is intersection of neighbors of currentVertex and the set of elements that come before currentVertex in the degener
        j = i-1;
        while(j>=0 )
        {
            if((graph[currentVertex].find(order[j]) != graph[currentVertex].end()))
            {XForPivot.insert(order[j]);}
            j--;
        }
        BronKerbosch(RForPivot, PForPivot, XForPivot);
    }
}

// Main function
int main() {
    // Read the graph
    auto start = std::chrono::high_resolution_clock::now();  // Start timer
    loadGraph("datenron.txt");
    cout<<"Version 1: Bron-Kerbosch with Pivot and Degenaracy"<<endl;

    // Open output file
    // ofstream outfile("output.txt");
    // if (!outfile) {
    //     cerr << "Error: Unable to open output file!" << endl;
    //     return 1;
    // }

    // Initialize P with all vertices
    for (const auto &pair : graph) {
        P.insert(pair.first);
    }

    // Run Bron-Kerbosch
    BronKerboschDegen(set<int>(), P, set<int>() );
    cout << "Total Cliques: " << totalCliques << endl;
    cout << "Size of Cliques: " << endl;
    for (auto it = sizeCounts.begin(); it != sizeCounts.end(); it++) {
        cout << it->first << " : " << it->second << endl;
    }
    auto stop = std::chrono::high_resolution_clock::now();   // Stop timer
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    
    std::cout << "Execution Time: " << duration.count() << " ms" << std::endl;
    return 0;
}