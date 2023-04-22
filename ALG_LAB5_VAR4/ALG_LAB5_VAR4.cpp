#include <iostream> // ввод-вывод
#include <vector> // шаблон класса для контейнеров последовательностей
#include <queue> // шаблон очереди
#include <stack> // шаблон стека
#include <fstream>
#include <chrono>

using namespace std;


//A class representing a graph.

class Graph {
public:
    /**
    * Returns the degree of vertex u in the graph.
    * The degree of a vertex is the number of edges incident to it.
    *
    * @param u the vertex to get the degree of
    * @return the degree of vertex u
    */
    int get_degree(int u) const;

    /**
    * Runs Kruskal's algorithm on the graph to find a minimum spanning tree.
    * Kruskal's algorithm works by repeatedly adding the shortest edge that does
    * not create a cycle in the current tree.
    *
    * @return a vector of pairs representing the edges in the minimum spanning tree
    *         Each pair contains the two vertices that the edge connects, and the weight of the edge.
    */
    vector<pair<int, int>> kruskal();
    
      /*  *Constructor for the Graph class.
        *
        * @param vertices The number of vertices in the graph.
        * @param directed Whether the graph is directed(true) or undirected(false).Default
        is false.
        * */
        Graph(int vertices, bool directed = false);
   /* /
        *Adds an edge to the graph.
        *
        * @param u The starting vertex of the edge.
        * @param v The ending vertex of the edge.
        * @param weight The weight of the edge.Default is 1.
        * /*/
        void add_edge(int u, int v, int weight);
   /* /
        *Returns the adjacency matrix of the graph.
        *
        *@return The adjacency matrix of the graph.
        * /*/
        vector<vector<int>> get_adj_matrix();
   /* /
        *Returns the adjacency list of the graph.
        *
        *@return The adjacency list of the graph.
        * /*/
        vector<vector<pair<int, int>>> get_adj_list();
   /* /
        *Returns the edges of the graph.
        *
        *@return The edges of the graph.
        * /*/
        vector<pair<int, pair<int, int>>> get_edges();
   /* /
        *Returns the incidence matrix of the graph.
        *
        *@return The incidence matrix of the graph.
        * /*/
        vector<vector<int>> get_inc_matrix();
    /*/
        *Prints the adjacency matrix of the graph.
        * /*/
        void print_adj_matrix() const;
   /* /
        *Prints the incidence matrix of the graph.
        * /*/
        void print_inc_matrix() const;
    /*/
        *Prints the adjacency list of the graph.
        * /*/
        void print_adj_list() const;
    void print_edge_list() const;
private:
    int vertices_; // The number of vertices in the graph.
    bool directed_; // Whether the graph is directed or undirected.
    vector<vector<int>> adj_matrix_; // The adjacency matrix of the graph.
    vector<vector<pair<int, int>>> adj_list_; // The adjacency list of the graph.
    vector<pair<int, pair<int, int>>> edges_; // The edges of the graph.
    vector<vector<int>> inc_matrix_; // The incidence matrix of the graph.
}; 
///
//*Constructor for the Graph class.
//* @param vertices The number of vertices in the graph.
//* @param directed A boolean indicating whether the graph is directed or undirected.
//* /
Graph::Graph(int vertices, bool directed) {
    vertices_ = vertices;
    directed_ = directed;
    adj_matrix_ = vector<vector<int>>(vertices, vector<int>(vertices, 0));
    adj_list_ = vector<vector<pair<int, int>>>(vertices);
}
///
//*Adds an edge between two vertices in the graph.
//* @param u The index of the first vertex.
//* @param v The index of the second vertex.
//* @param weight The weight of the edge.
//* /
void Graph::add_edge(int u, int v, int weight) {
    adj_matrix_[u][v] = weight;
    adj_list_[u].push_back(make_pair(v, weight));
    edges_.push_back(make_pair(weight, make_pair(u, v)));
    if (directed_) {
        adj_matrix_[u][v] = weight;
    }
    else {
        adj_matrix_[u][v] = adj_matrix_[v][u] = weight;
    }
}
///
//*Returns the adjacency matrix of the graph.
//* @return A 2D vector representing the adjacency matrix.
//* /
vector<vector<int>> Graph::get_adj_matrix() {
    return adj_matrix_;
}
///
//*Returns the adjacency list of the graph.
//* @return A 2D vector representing the adjacency list.
//* /
vector<vector<pair<int, int>>> Graph::get_adj_list() {
    return adj_list_;
}
///
//*Returns the edges of the graph.
//* @return A vector of pairs representing the edges.
//* /
vector<pair<int, pair<int, int>>> Graph::get_edges() {
    return edges_;
}

///
//*Prints the adjacency matrix of the graph.
//* /
void Graph::print_adj_matrix() const {
    std::cout << "    ";
    for (int i = 0; i < vertices_; i++) { 
        std::cout << "V" << i << " "; 
        if (i < 10) cout << " ";
    }
    std::cout << endl;
    for (int i = 0; i < vertices_; i++) {
        cout << "V" << i << " ";
        if (i < 10) cout << " ";
        for (int j = 0; j < vertices_; j++) {
            cout << adj_matrix_[i][j] << "  ";
            if (adj_matrix_[i][j] < 10) cout << " ";
        }
        cout << endl;
    }
    cout << endl;
}


///
//*Prints the adjacency list of the graph.
//* /
void Graph::print_adj_list() const {
    for (int i = 0; i < vertices_; i++) {
        cout << i << ": ";
        for (int j = 0; j < vertices_; j++) {
            if (adj_matrix_[i][j] != 0) {
                cout << j << "(" << adj_matrix_[i][j] << ") ";
            }
        }
        cout << endl;
    }
    cout << endl;
}

int Graph::get_degree(int u) const {
    int degree = 0;
    for (int v = 0; v < adj_matrix_.size(); v++) {
        if (adj_matrix_[u][v] != 0) {
            degree++;
        }
    }
    return degree;
}

void Graph::print_edge_list() const {
    for (auto edge : edges_) {
        cout << edge.second.first << " -> " << edge.second.second;
        if (edge.first != 0) {
            cout << " (" << edge.first << ")";
        }
        std::cout << endl;
    }
    std::cout << endl;
}

///
//
//Finds the minimum spanning tree of the graph using Kruskal's algorithm.
//
//@return A vector of pairs representing the edges in the minimum spanning tree.
//* /
vector<pair<int, int>> Graph::kruskal() {
    // Create a vector to store the edges in the minimum spanning tree
    vector<pair<int, int>> mst;// Create a vector to store the parent of each vertex in the union-find data structure
vector<int> parent(vertices_);

// Initialize each vertex's parent to itself
for (int i = 0; i < vertices_; i++) {
    parent[i] = i;
}

// Sort the edges by weight
sort(edges_.begin(), edges_.end(), [](const pair<int ,pair<int, int>>& a, const pair<int, pair<int, int>>& b) {
    return abs(a.first) < abs(b.first);
    });

// Iterate over the edges in order of increasing weight
for (auto edge : edges_) {
    int u = edge.second.first;
    int v = edge.second.second;

    // Find the parents of u and v in the union-find data structure
    while (parent[u] != u) {
        u = parent[u];
    }
    while (parent[v] != v) {
        v = parent[v];
    }

    // If the parents are different, add the edge to the minimum spanning tree
    if (u != v) {
        mst.push_back(edge.second);

        // Merge the two trees by setting one tree's parent to the other tree's parent
        parent[u] = v;
    }
}

return mst;
}

/**
 * Generate a random graph with the given number of vertices and minimum edges per vertex.
 *
 * @param num_vertices: The number of vertices in the graph.
 * @param min_edges_per_vertex: The minimum number of edges each vertex should have.
 * @return A Graph object representing the generated graph.
 */
Graph generate_graph(int num_vertices, int min_edges_per_vertex) {
    Graph g(num_vertices);

    // Create a connected graph with num_vertices-1 edges using BFS
    queue<int> q;
    vector<bool> visited(num_vertices, false);
    visited[0] = true;
    q.push(0);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        int num_edges_added = 0;
        for (int v = 0; v < num_vertices; v++) {
            if (u == v || g.get_adj_matrix()[u][v] != 0) {
                continue;
            }
            int weight = rand() % 20 + 1;
            g.add_edge(u, v, weight);
            num_edges_added++;
            if (num_edges_added == num_vertices - 1) {
                break;
            }
        }
        for (int v = 0; v < num_vertices; v++) {
            if (g.get_adj_matrix()[u][v] != 0 && !visited[v]) {
                visited[v] = true;
                q.push(v);
            }
        }
    }

    // Add additional edges until each vertex has at least min_edges_per_vertex
    for (int u = 0; u < num_vertices; u++) {
        while (g.get_degree(u) < min_edges_per_vertex + 1) {
            int v = rand() % num_vertices;
            if (u == v || g.get_adj_matrix()[u][v] != 0) {
                continue;
            }
            int weight = rand() % 20 + 1;
            g.add_edge(u, v, weight);
        }
    }

    return g;
}



int main() {
    srand(time(0));
    int num_vertices[] = { 10, 20, 50, 100 };
    int min_edges_per_vertex[] = { 3, 4, 10, 20 };
    int num_graphs = 5;
    int j = 0;
    
    Graph g(4);
    g.add_edge(0, 1, 2);
    g.add_edge(0, 2, 3);
    g.add_edge(0, 3, 4);
    g.add_edge(1, 2, 1);
    g.add_edge(2, 3, 5);

    // Compute the minimum spanning tree using Kruskal's algorithm
    vector<pair<int, int>> mst = g.kruskal();

    // Print the edges of the minimum spanning tree
    cout << "Edges in the minimum spanning tree:" << endl;
    for (auto edge : mst) {
        cout << edge.first << " -- " << edge.second << endl;
    }
    ofstream fout("output.txt");
    for (j = 0; j < 4; j++) {

        for (int i = 0; i < num_graphs; i++) {

            Graph g = generate_graph(num_vertices[j], min_edges_per_vertex[j]);

            std::cout << "Adjacency matrix:" << endl;
            g.print_adj_matrix();

            /*std::cout << "Incidence matrix:" << endl;
            g.print_inc_matrix();*/

            std::cout << "Adjacency list:" << endl;
            g.print_adj_list();

           /* std::cout << "Edge list:" << endl;
            g.print_edge_list();*/

            std::cout << "Graph " << i + 1 << " with " << num_vertices[j] << " vertices " << endl;
            fout << "Graph " << i + 1 << " with " << num_vertices[j] << " vertices " << endl;
            chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
            vector<pair<int, int>> krusk = g.kruskal();
            chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
            chrono::duration<double, milli> milli_diff = end - start;

            std::cout << "Time kruskal: " << milli_diff.count() << " milli sec" << endl;
            fout << "Time kruskal: " << milli_diff.count() << " milli sec" << endl;
            std::cout << endl << endl << endl;
            fout << endl << endl << endl;
        }
    }
    
    fout.close();
    return EXIT_SUCCESS;
}