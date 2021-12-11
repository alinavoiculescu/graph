#include <bits/stdc++.h>
using namespace std;

ifstream fin("biconex.in");
ofstream fout("biconex.out");

const int INF = 0x3f3f3f3f;
const int MAX = 5000; //500000

void CountingSort(vector<int>&);

class Graph
{
///PRIVATE
    int numberOfNodes;
    int numberOfEdges;
    bool directed;
    vector<vector<int>> adjacencyList;
    vector<pair<int, int>> adjacencyListWeightedGraph[MAX];
    vector<vector<int>> adjacencyMatrixWeightedGraph;
    vector<pair<int,int>> adjacencyListWithEdgesNumber[MAX];

    void DFS(int, bool[]);

    void _BiconnectedComponents(int, int, int[], int[], stack<int>&, bool[], vector<vector<int>>&);

    void _StronglyConnectedComponents(int, int&, int[], int[], stack<int>&, bool[], bool[], vector<vector<int>>&);

    void _TopologicalSort(int, bool[], vector<int>&);

    void _CriticalConnections(int, bool[], int[], int[], int[], vector<vector<int>>&);

    int FindRoot(int, int[]);

    void UnionRoots(int, int, int[], int[]);

    bool MaxFlowBFS(int, int, vector<vector<int>>, int[]);


///PUBLIC
    public:
        Graph();
        Graph(int, int, bool);
        Graph(int, bool);
        void Build(istream&);
        void BuildFromVector(vector<vector<int>>);
        void BuildWeightedGraph(istream&);
        void BuildAdjacencyMatrixWeightedGraph(istream&);
        void BuildGraphWithEdgesNumber(istream&);
        int GetNumberOfNodes();
        int GetNumberOfEdges();

        void BFS(int, int[]);

        int NumberOfConnectedComponents();

        vector<vector<int>> BiconnectedComponents(int);

        vector<vector<int>> StronglyConnectedComponents();

        vector<int> TopologicalSort();

        //Checks if a graph exists
        bool HavelHakimi(vector<int>&);

        vector<vector<int>> CriticalConnections();

        //Prim's Minimum Spanning Tree
        void MinimumSpanningTree(int, ostream&);

        void Dijkstra(int, ostream&);

        void BellmanFord(int, ostream&);

        void DisjointSetForests(istream&, ostream&);

        //Ford-Flurkenson && Edmonds-Karp
        int MaxFlow(int, int);

        int TreeDiameter();

        vector<vector<int>> RoyFloyd();

        vector<int> EulerianPath(int);
};

//LeetCode Critical Connections Problem
class Solution {
public:
    vector<vector<int>> CriticalConnections(int n, vector<vector<int>>& connections)
    {
        Graph G(n, connections.size(), 0);
        G.BuildFromVector(connections);

        return G.CriticalConnections();
    }
};

int main()
{




/*  ----------------------------EulerianCycle----------------------------
        int numberOfNodes, numberOfEdges, startNode = 0;
        vector<int> path;

        fin >> numberOfNodes >> numberOfEdges;

        Graph G(numberOfNodes, numberOfEdges, 0);
        G.BuildGraphWithEdgesNumber(fin);

        path = G.EulerianPath(startNode);

        for (unsigned int i = 0; i < path.size(); i++)
            fout << path[i] + 1 << " ";
    ----------------------------EulerianCycle----------------------------*/



/*  ----------------------------RoyFloyd----------------------------
    int numberOfNodes;
    fin >> numberOfNodes;

    Graph G(numberOfNodes, 1);
    G.BuildAdjacencyMatrixWeightedGraph(fin);

    vector<vector<int>> distance = G.RoyFloyd();

    for (int i = 0; i < numberOfNodes; i++)
    {
        for (int j = 0; j < numberOfNodes; j++)
        {
            fout << distance[i][j] << " ";
        }
        fout << "\n";
    }
    ----------------------------RoyFloyd----------------------------*/


/*  ----------------------------TreeDiameter----------------------------
    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes;
    numberOfEdges = numberOfNodes - 1;

    Graph G(numberOfNodes, numberOfEdges, 0);
    G.Build(fin);

    fout << G.TreeDiameter();
    ----------------------------TreeDiameter----------------------------*/



/*  ----------------------------MaxFlow----------------------------
    int numberOfNodes, numberOfEdges, source, destination;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1);
    G.BuildWeightedGraph(fin);

    source = 0;
    destination = G.GetNumberOfNodes() - 1;

    fout << G.MaxFlow(source, destination);
    ----------------------------MaxFlow----------------------------*/



/*  ----------------------------Disjoint----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 0);

    G.DisjointSetForests(fin, fout);

    ----------------------------Disjoint----------------------------*/



/*  ----------------------------BellmanFord----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1);
    G.BuildWeightedGraph(fin);

    startNode = 0;
    G.BellmanFord(startNode, fout);

    ----------------------------BellmanFord----------------------------*/



/*  ----------------------------Dijkstra----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1);
    G.BuildWeightedGraph(fin);

    startNode = 0;
    G.Dijkstra(startNode, fout);

    ----------------------------Dijkstra----------------------------*/



/*  ----------------------------APM----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 0);
    G.BuildWeightedGraph(fin);

    startNode = 0;
    G.MinimumSpanningTree(startNode, fout);

    ----------------------------APM----------------------------*/



/*  ----------------------------CCs----------------------------

    Solution S;
    int n = 4;
    vector<vector<int>> connections = {{0,1},{1,2},{2,0},{1,3}}, criticalConnections;

    criticalConnections = S.CriticalConnections(n, connections);

    for (unsigned int i = 0; i < criticalConnections.size(); i++)
        cout << "Critical connection " << i + 1 << ": " << criticalConnections[i][0] << " - " << criticalConnections[i][1] << "\n";

    ----------------------------CCs----------------------------*/



/*  ----------------------------HH----------------------------

    Graph G;
    vector<int> degrees = {5,5,5,3,2,2,2};

    if (G.HavelHakimi(degrees))
        cout << "Yes";
    else
       cout << "No";

    ----------------------------HH----------------------------*/



/*  ----------------------------TS----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1);
    G.Build(fin);

    vector<int> topologicalSort = G.TopologicalSort();

    reverse(topologicalSort.begin(), topologicalSort.end());

    for (auto node : topologicalSort)
        fout << node + 1 << " ";

    ----------------------------TS----------------------------*/


/*  ----------------------------SCCs----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1);
    G.Build(fin);

    vector<vector<int>> stronglyConnectedComponents = G.StronglyConnectedComponents();
    fout << stronglyConnectedComponents.size() << "\n";

    for (unsigned int i = 0; i < stronglyConnectedComponents.size(); i++)
    {
        vector<int> stronglyConnectedComponent = stronglyConnectedComponents[i];

        for (unsigned int j = 0; j < stronglyConnectedComponent.size(); j++)
        {
            int node = stronglyConnectedComponent[j] + 1;
            fout << node << " ";
        }
        fout << "\n";
    }

    ----------------------------SCCs----------------------------*/


/*  ----------------------------BCCs----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 0);
    G.Build(fin);

    vector<vector<int>> biconnectedComponents = G.BiconnectedComponents(0);
    fout << biconnectedComponents.size() << "\n";

    for (unsigned int i = 0; i < biconnectedComponents.size(); i++)
    {
        vector<int> biconnectedComponent = biconnectedComponents[i];

        for (unsigned int j = 0; j < biconnectedComponent.size(); j++)
        {
            int node = biconnectedComponent[j] + 1;
            fout << node << " ";
        }
        fout << "\n";
    }

    ----------------------------BCCs----------------------------*/


/*  ----------------------------DFS----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 0);
    G.Build(fin);

    fout << G.numberOfConnectedComponents();

    ----------------------------DFS----------------------------*/

/*  ----------------------------BFS----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    fin >> numberOfNodes >> numberOfEdges >> startNode;

    Graph G(numberOfNodes, numberOfEdges, 1);
    G.Build(fin);

    int distance[numberOfNodes];
    startNode--;
    G.BFS(startNode, distance);

    for (int i = 0; i < numberOfNodes; i++)
        fout << distance[i] << " ";

    ----------------------------BFS----------------------------*/

    fin.close();
    fout.close();

    return 0;
}

Graph::Graph() : numberOfNodes(0), numberOfEdges(0), directed(0) { }

Graph::Graph(int _numberOfNodes, int _numberOfEdges, bool _directed) : numberOfNodes(_numberOfNodes), numberOfEdges(_numberOfEdges), directed(_directed)
{
    for (int i = 0; i < _numberOfNodes; i++)
        adjacencyList.push_back( {} );
        //adjacencyList.push_back(vector<int>());
}

Graph::Graph(int _numberOfNodes, bool _directed) : numberOfNodes(_numberOfNodes), numberOfEdges(INF), directed(_directed)
{
    for (int i = 0; i < _numberOfNodes; i++)
        adjacencyMatrixWeightedGraph.push_back( {} );
        //adjacencyMatrixWeightedGraph.push_back(vector<int>());
}

void Graph::Build(istream &fin)
{
    for (int i = 0; i < numberOfEdges; i++)
    {
        int firstNode, secondNode;
        fin >> firstNode >> secondNode;
        firstNode--;
        secondNode--;

        adjacencyList[firstNode].push_back(secondNode);
        if (!directed)
            adjacencyList[secondNode].push_back(firstNode);
    }
}

void Graph::BuildAdjacencyMatrixWeightedGraph(istream &fin)
{
    for (int i = 0; i < numberOfNodes; i++)
    {
        adjacencyMatrixWeightedGraph[i].resize(numberOfNodes);

        for (int j = 0; j < numberOfNodes; j++)
            fin >> adjacencyMatrixWeightedGraph[i][j];
    }
}

void Graph::BuildFromVector(vector<vector<int>> edges)
{
    for(unsigned int i = 0; i < edges.size(); i++)
    {
        int firstNode, secondNode;
        firstNode = edges[i][0];
        secondNode = edges[i][1];

        adjacencyList[firstNode].push_back(secondNode);
        if (!directed)
            adjacencyList[secondNode].push_back(firstNode);

    }
}

void Graph::BuildWeightedGraph(istream &fin)
{
    for (int i = 0; i < numberOfEdges; i++)
    {
        int firstNode, secondNode, weight;
        fin >> firstNode >> secondNode >> weight;
        firstNode--;
        secondNode--;

        adjacencyListWeightedGraph[firstNode].push_back(make_pair(secondNode, weight));

        if (!directed)
            adjacencyListWeightedGraph[secondNode].push_back(make_pair(firstNode, weight));
    }
}

void Graph::BuildGraphWithEdgesNumber(istream &fin)
{
    for (int i = 0; i < numberOfEdges; i++)
    {
        int firstNode, secondNode;
        fin >> firstNode >> secondNode;
        firstNode--;
        secondNode--;

        adjacencyListWithEdgesNumber[firstNode].push_back(make_pair(secondNode, i));

        if (!directed)
            adjacencyListWithEdgesNumber[secondNode].push_back(make_pair(firstNode, i));
    }
}

int Graph::GetNumberOfNodes()
{
    return numberOfNodes;
}

int Graph::GetNumberOfEdges()
{
    return numberOfEdges;
}

void Graph::BFS(int startNode, int distance[])
{
    bool visited[numberOfNodes];
    queue<int> queueBFS;

    for(int i = 0; i < numberOfNodes; i++)
    {
        visited[i] = 0;
        distance[i] = -1;
    }

    visited[startNode] = 1;
    distance[startNode] = 0;
    queueBFS.push(startNode);

    while (!queueBFS.empty())
    {
        int currentNode = queueBFS.front();

        for (unsigned int i = 0; i < adjacencyList[currentNode].size(); i++)
        {
            int adjacentNode = adjacencyList[currentNode][i];

            if (!visited[adjacentNode])
            {
                visited[adjacentNode] = 1;
                distance[adjacentNode] = distance[currentNode] + 1;
                queueBFS.push(adjacentNode);
            }
        }
        queueBFS.pop();
    }
}

void Graph::DFS(int startNode, bool visited[])
{
    visited[startNode] = 1;

    for (unsigned int i = 0; i < adjacencyList[startNode].size(); i++) //parcurg vecinii nodului
    {
        int nextNode = adjacencyList[startNode][i];
        if (!visited[nextNode])
            DFS(nextNode, visited);
    }
}

int Graph::NumberOfConnectedComponents()
{
    bool visited[numberOfNodes];
    int _numberOfConnectedComponents = 0;

    for (int i = 0; i < numberOfNodes; i++)
        visited[i] = 0;

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (!visited[i])
        {
            DFS(i, visited);
            _numberOfConnectedComponents++;
        }
    }

    return _numberOfConnectedComponents;
}

void Graph::_BiconnectedComponents(int node, int id, int ids[], int low[], stack<int> &stackBiconnectedComponents, bool visited[], vector<vector<int>> &biconnectedComponents)
{
    stackBiconnectedComponents.push(node);
    visited[node] = 1;
    ids[node] = low[node] = id++;

    for (unsigned int i = 0; i < adjacencyList[node].size(); i++)
    {
        int adjacentNode = adjacencyList[node][i];

        if (visited[adjacentNode])
            low[node] = min(low[node], ids[adjacentNode]);
        else
        {
            _BiconnectedComponents(adjacentNode, id, ids, low, stackBiconnectedComponents, visited, biconnectedComponents);

            low[node] = min(low[node], low[adjacentNode]);

            if (low[adjacentNode] >= ids[node])
            {
                vector<int> biconnectedComponent;
                int nodeBiconnectedComponent;

                do {
                    nodeBiconnectedComponent = stackBiconnectedComponents.top();
                    biconnectedComponent.push_back(nodeBiconnectedComponent);
                    stackBiconnectedComponents.pop();
                    } while (nodeBiconnectedComponent != adjacentNode);
                biconnectedComponent.push_back(node);

                biconnectedComponents.push_back(biconnectedComponent);
            }
        }
    }
}

void Graph::_StronglyConnectedComponents(int node, int &id, int ids[], int low[], stack<int> &stackStronglyConnectedComponents, bool onStack[], bool visited[], vector<vector<int>> &stronglyConnectedComponents)
{
    stackStronglyConnectedComponents.push(node);
    onStack[node] = 1;
    visited[node] = 1;
    ids[node] = low[node] = id++;

    for (unsigned int i = 0; i < adjacencyList[node].size(); i++)
    {
        int adjacentNode = adjacencyList[node][i];

        if (!visited[adjacentNode])
        {
            _StronglyConnectedComponents(adjacentNode, id, ids, low, stackStronglyConnectedComponents, onStack, visited, stronglyConnectedComponents);

            low[node] = min(low[node], low[adjacentNode]);
        }
        else if (onStack[adjacentNode])
            low[node] = min(low[node], low[adjacentNode]);
    }

    if (ids[node] == low[node])
    {
        vector<int> stronglyConnectedComponent;
        int nodeStronglyConnectedComponent;

        do {
            nodeStronglyConnectedComponent = stackStronglyConnectedComponents.top();
            stronglyConnectedComponent.push_back(nodeStronglyConnectedComponent);
            stackStronglyConnectedComponents.pop();
            onStack[nodeStronglyConnectedComponent] = 0;
            } while (nodeStronglyConnectedComponent != node);

        stronglyConnectedComponents.push_back(stronglyConnectedComponent);
    }
}

vector<vector<int>> Graph::BiconnectedComponents(int node)
{
    int id = 0, ids[numberOfNodes], low[numberOfNodes];
    bool visited[numberOfNodes];
    vector<vector<int>> biconnectedComponents;
    stack<int> stackBiconnectedComponents;

    for (int i = 0; i < numberOfNodes; i++)
        visited[i] = 0;

    _BiconnectedComponents(node, id, ids, low, stackBiconnectedComponents, visited, biconnectedComponents);

    return biconnectedComponents;
}

vector<vector<int>> Graph::StronglyConnectedComponents()
{
    int id = 0, ids[numberOfNodes], low[numberOfNodes];
    bool onStack[numberOfNodes], visited[numberOfNodes];
    vector<vector<int>> stronglyConnectedComponents;
    stack<int> stackStronglyConnectedComponents;

    for (int i = 0; i < numberOfNodes; i++)
    {
        onStack[i] = 0;
        visited[i] = 0;
    }

    for (int i = 0; i < numberOfNodes; i++)
        if (!visited[i])
            _StronglyConnectedComponents(i, id, ids, low, stackStronglyConnectedComponents, onStack, visited, stronglyConnectedComponents);

    return stronglyConnectedComponents;
}

void Graph::_TopologicalSort(int startNode, bool visited[], vector<int> &topologicalSort)
{
    visited[startNode] = 1;

    for (unsigned int i = 0; i < adjacencyList[startNode].size(); i++) //parcurg vecinii nodului
    {
        int nextNode = adjacencyList[startNode][i];
        if (!visited[nextNode])
            _TopologicalSort(nextNode, visited, topologicalSort);
    }
    topologicalSort.push_back(startNode);
}


vector<int> Graph::TopologicalSort()
{
    vector<int> topologicalSort;
    bool visited[numberOfNodes];

    for(int i = 0; i < numberOfNodes; i++)
        visited[i] = 0;

    for(int i = 0; i < numberOfNodes; i++)
        if(!visited[i])
            _TopologicalSort(i, visited, topologicalSort);

    return topologicalSort;
}

void CountingSort(vector<int>& toSort)
{
    vector<int> countingSort(MAX, 0);
    int maxValue = -1;

    for (unsigned int i = 0; i < toSort.size(); i++)
    {
        countingSort[toSort[i]]++;

        if (toSort[i] > maxValue)
            maxValue = toSort[i];
    }

    int i = 0;

    for (int j = maxValue; j >= 0; j--)
        while (countingSort[j] != 0)
        {
            toSort[i++] = j;
            countingSort[j]--;
        }
}

bool Graph::HavelHakimi(vector<int> &degrees)
{
    numberOfNodes = degrees.size();

    if (numberOfNodes < 1 || degrees[0] == 0)
        return 1;

    if (accumulate(degrees.begin(), degrees.end(), 0) % 2)
        return 0;

    CountingSort(degrees);

    if (degrees[0] > numberOfNodes - 1)
        return 0;

    int element = degrees[0];
    degrees.erase(degrees.begin() + 0);

    for (int i = 0; i < element; i++)
    {
        if(degrees[i] > 0)
            degrees[i]--;
        else
            return 0;
    }

    return HavelHakimi(degrees);
}

void Graph::_CriticalConnections(int node, bool visited[], int disc[], int low[], int parent[], vector<vector<int>> &criticalConnections)
{
    static int time = 0;

    visited[node] = 1;
    disc[node] = low[node] = ++time;

    for (unsigned int i = 0; i < adjacencyList[node].size(); i++)
    {
        int adjacentNode = adjacencyList[node][i];

        if (!visited[adjacentNode])
        {
            parent[adjacentNode] = node;
            _CriticalConnections(adjacentNode, visited, disc, low, parent, criticalConnections);

            low[node] = min(low[node], low[adjacentNode]);

            if(low[adjacentNode] > disc[node])
            {
                vector<int> criticalConnection;

                criticalConnection.push_back(node);
                criticalConnection.push_back(adjacentNode);

                criticalConnections.push_back(criticalConnection);
            }
        }
        else if (adjacentNode != parent[node])
            low[node] = min(low[node], disc[adjacentNode]);
    }
}

vector<vector<int>> Graph::CriticalConnections()
{
    vector<vector<int>> criticalConnections;
    bool visited[numberOfNodes];
    int disc[numberOfNodes], low[numberOfNodes], parent[numberOfNodes];

    for (int i = 0; i < numberOfNodes; i++)
    {
        parent[i] = -1;
        visited[i] = 0;
    }

    for (int i = 0; i < numberOfNodes; i++)
        if (!visited[i])
            _CriticalConnections(i, visited, disc, low, parent, criticalConnections);

    return criticalConnections;
}

void Graph::MinimumSpanningTree(int startNode, ostream &fout)
{
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> minHeap;
    int totalCost = 0, minimumSpanningTreeEdges = 0;
    vector<int> key(numberOfNodes, INF);
    vector<int> parent(numberOfNodes, -1);
    vector<bool> inMST(numberOfNodes, 0);

    minHeap.push(make_pair(0, startNode));
    key[startNode] = 0;

    while(!minHeap.empty())
    {
        int node = minHeap.top().second;
        minHeap.pop();

        if(!inMST[node])
        {
            inMST[node] = 1;

            for (auto i : adjacencyListWeightedGraph[node])
            {
                int adjacentNode = i.first;
                int weight = i.second;

                if (!inMST[adjacentNode] && weight < key[adjacentNode])
                {
                    if (key[adjacentNode] != INF)
                        totalCost -= key[adjacentNode];
                    else
                        minimumSpanningTreeEdges++;

                    key[adjacentNode] = weight;
                    parent[adjacentNode] = node;
                    minHeap.push(make_pair(key[adjacentNode], adjacentNode));

                    totalCost += key[adjacentNode];
                }
            }
        }
    }

    fout << totalCost << "\n" << minimumSpanningTreeEdges << "\n";

    for (int i = 1; i < numberOfNodes; i++)
        fout << parent[i] + 1 << " " << i + 1 << "\n";
}

void Graph::Dijkstra(int startNode, ostream &fout)
{
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> minHeap;
    vector<int> key(numberOfNodes, INF);
    vector<bool> inDijkstra(numberOfNodes, 0);

    minHeap.push(make_pair(0, startNode));
    key[startNode] = 0;

    while(!minHeap.empty())
    {
        int node = minHeap.top().second;
        minHeap.pop();

        if(!inDijkstra[node])
        {
            inDijkstra[node] = 1;

            for (auto i : adjacencyListWeightedGraph[node])
            {
                int adjacentNode = i.first;
                int adjacentNodeWeight = i.second;

                if (key[node] + adjacentNodeWeight < key[adjacentNode])
                {
                    key[adjacentNode] = key[node] + adjacentNodeWeight;
                    minHeap.push(make_pair(key[adjacentNode], adjacentNode));
                }
            }
        }
    }

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (i != startNode)
        {
            if (key[i] == INF)
            {
                key[i] = 0;
            }
            fout << key[i] << " ";
        }
    }
}

void Graph::BellmanFord(int startNode, ostream &fout)
{
    vector<int> key(numberOfNodes, INF);
    queue<int> queueBF;
    long long int iterations = 0;
    bool negativeCycle = 0;

    key[startNode] = 0;
    queueBF.push(startNode);

    while (!queueBF.empty())
    {
        if (iterations > (long long int)(numberOfNodes - 1) * numberOfEdges)
        {
            negativeCycle = 1;
            break;
        }

        int node = queueBF.front();
        queueBF.pop();

        for (auto i : adjacencyListWeightedGraph[node])
        {
            int adjacentNode = i.first;
            int adjacentNodeWeight = i.second;

            if (key[node] + adjacentNodeWeight < key[adjacentNode])
            {

                key[adjacentNode] = key[node] + adjacentNodeWeight;
                queueBF.push(adjacentNode);
            }
        }
        iterations++;
    }

    if (negativeCycle)
        fout << "Ciclu negativ!";
    else
        for (int i = 0; i < numberOfNodes; i++)
            if (i != startNode)
            {
                if (key[i] == INF)
                    key[i] = 0;
                fout << key[i] << " ";
            }
}

int Graph::FindRoot(int node, int parent[])
{
    if (parent[node] == node)
        return node;

    return FindRoot(parent[node], parent);
}

void Graph::UnionRoots(int node1, int node2, int parent[], int height[])
{
    if (node1 != node2)
    {
        if (height[node1] > height[node2])
        {
            height[node1] += height[node2];
            parent[node2] = node1;
        }
        else
        {
            height[node2] += height[node1];
            parent[node1] = node2;
        }
    }
}

void Graph::DisjointSetForests(istream &fin, ostream &fout)
{
    int parent[numberOfNodes], height[numberOfNodes];

    for (int i = 0; i < numberOfNodes; i++)
    {
        parent[i] = i;
        height[i] = 1;
    }

    for (int i = 0; i < numberOfEdges; i++)
    {
        int code, node1, node2;
        fin >> code >> node1 >> node2;

        node1--;
        node2--;

        int node1Root, node2Root;
        node1Root = FindRoot(node1, parent);
        node2Root = FindRoot(node2, parent);

        if (code == 1)
            UnionRoots(node1Root, node2Root, parent, height);
        else if (code == 2)
        {
            if (node1Root == node2Root)
                fout << "DA\n";
            else
                fout << "NU\n";
        }
    }
}

bool Graph::MaxFlowBFS(int source, int destination, vector<vector<int>> residualGraph, int parent[])
{
    bool visited[numberOfNodes] = {0};
    queue<int> queueMaxFlowBFS;

    queueMaxFlowBFS.push(source);
    visited[source] = 1;
    parent[source] = -1;

    while (!queueMaxFlowBFS.empty())
    {
        int node = queueMaxFlowBFS.front();
        queueMaxFlowBFS.pop();

        for (int i = 0; i < numberOfNodes; i++)
        {
            if (!visited[i] && residualGraph[node][i] > 0)
            {
                if (i == destination)
                {
                    parent[i] = node;
                    return 1;
                }

                queueMaxFlowBFS.push(i);
                parent[i] = node;
                visited[i] = 1;
            }
        }
    }

    return 0;
}

int Graph::MaxFlow(int source, int destination)
{
    vector<vector<int>> residualGraph(numberOfNodes, vector<int> (numberOfNodes, 0));
    int parent[numberOfNodes];
    int maxFlow = 0;

    for (int i = 0; i < numberOfNodes; i++)
        for (auto j : adjacencyListWeightedGraph[i])
        {
            int adjacentNode = j.first;
            residualGraph[i][adjacentNode] = j.second;
        }

    while(MaxFlowBFS(source, destination, residualGraph, parent))
    {
        int pathFlow = INF;

        for(int i = destination; i != source; i = parent[i])
        {
            int node = parent[i];
            pathFlow = min(pathFlow, residualGraph[node][i]);
        }

        for(int i = destination; i != source; i = parent[i])
        {
            int node = parent[i];
            residualGraph[node][i] -= pathFlow;
            residualGraph[i][node] += pathFlow;
        }

        maxFlow += pathFlow;
    }

    return maxFlow;
}

int Graph::TreeDiameter()
{
    int distance[numberOfNodes], node1 = 0, node2 = 0, diameter = 0;
    BFS(node1, distance);

    for (int i = 1; i < numberOfNodes; i++)
        if (distance[i] > distance[node2])
            node2 = i;

    BFS(node2, distance);

    for (int i = 0; i < numberOfNodes; i++)
        diameter = max(diameter, distance[i]);

    return diameter + 1;
}

vector<vector<int>> Graph::RoyFloyd()
{
    vector<vector<int>> distance(numberOfNodes, vector<int> (numberOfNodes, 0));

    for (int i = 0; i < numberOfNodes; i++)
        for (int j = 0; j < numberOfNodes; j++)
            distance[i][j] = adjacencyMatrixWeightedGraph[i][j];

    for (int k = 0; k < numberOfNodes; k++)
    {
        for (int i = 0; i < numberOfNodes; i++)
        {
            for (int j = 0; j < numberOfNodes; j++)
            {
                if (i != j &&
                    distance[i][k] && distance[k][j] &&
                   (distance[i][j] > distance[i][k] + distance[k][j] || !distance[i][j]))
                {
                    distance[i][j] = distance[i][k] + distance[k][j];
                }
            }
        }
    }

    return distance;
}

vector<int> Graph::EulerianPath(int startNode)
{
    vector<bool> visited(numberOfEdges, 0);
    stack<int> stackEuler;
    vector<int> path;

    for(int i = 0; i < numberOfNodes; i++)
        if (adjacencyListWithEdgesNumber[i].size() % 2 == 1)
            return {-2};

    stackEuler.push(startNode);

    while (!stackEuler.empty())
    {
        int node = stackEuler.top();

        if (!adjacencyListWithEdgesNumber[node].empty())
        {
            int edgeNumber = adjacencyListWithEdgesNumber[node].back().second;
            int adjacentNode = adjacencyListWithEdgesNumber[node].back().first;

            adjacencyListWithEdgesNumber[node].pop_back();

            if (!visited[edgeNumber])
            {
                visited[edgeNumber] = 1;
                stackEuler.push(adjacentNode);
            }
        }
        else
        {
            stackEuler.pop();
            path.push_back(node);
        }
    }

    path.pop_back();

    return path;
}
