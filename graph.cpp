#include <bits/stdc++.h>
using namespace std;

ifstream fin("biconex.in");
ofstream fout("biconex.out");

const int INF = 0x3f3f3f3f;
const int MAX = 5000; //500000 eulerian cycle
                      //100000 apm, tree diameter

void CountingSort(vector<int>&);

class Graph
{
///PRIVATE
    int numberOfNodes;
    int numberOfEdges;
    bool directed;
    bool fromOne;  //if first node is 1 => true

    vector<vector<int>> adjacencyList;
    vector<pair<int, int>> adjacencyListWeightedGraph[MAX];
    vector<vector<int>> adjacencyMatrixWeightedGraph;
    vector<pair<int,int>> adjacencyListWithEdgesNumber[MAX];

    void functionBiconnectedComponents(int, int, int[], int[], stack<int>&, bool[], vector<vector<int>>&);

    void functionStronglyConnectedComponents(int, int&, int[], int[], stack<int>&, bool[], bool[], vector<vector<int>>&);

    void functionTopologicalSort(int, bool[], vector<int>&);

    void functionCriticalConnections(int, bool[], int[], int[], int[], vector<vector<int>>&);

    int FindRoot(int, int[]);

    void UnionRoots(int, int, int[], int[]);

    bool MaxFlowBFS(int, int, vector<vector<int>>, int[]);

    bool MaximumBipartiteMatchingBFS(vector<vector<int>>&, vector<int>&, vector<int>&, vector<int>&);

    bool MaximumBipartiteMatchingDFS(int, vector<vector<int>>&, vector<int>&, vector<int>&, vector<int>&);


///PUBLIC
    public:
        Graph();
        Graph(int, int, bool, bool);
        Graph(int, bool, bool);

        ~Graph() { }

        void Build(istream&);
        void BuildFromVector(vector<vector<int>>);
        void BuildWeightedGraph(istream&);
        void BuildAdjacencyMatrixWeightedGraph(istream&);
        void BuildGraphWithEdgesNumber(istream&);
        void BuildTransposeWeightedGraph(istream&);

        int GetNumberOfNodes();
        int GetNumberOfEdges();
        bool IsDirected();
        bool IsFromOne();

        void BFS(int, int[]);

        void DFS(int, bool[]);

        int NumberOfConnectedComponents();

        vector<vector<int>> BiconnectedComponents(int);

        vector<vector<int>> StronglyConnectedComponents();

        vector<int> TopologicalSort();

        //Checks if a graph exists
        bool HavelHakimi(vector<int>&);

        vector<vector<int>> CriticalConnections();

        //Prim's Minimum Spanning Tree
        vector<int> MinimumSpanningTree(int, int&, int&);

        vector<int> Dijkstra(int);

        vector<int> BellmanFord(int, bool&);

        void DisjointSetForests(istream&, ostream&);

        //Ford-Fulkerson && Edmonds-Karp
        int MaxFlow(int, int);

        int TreeDiameter();

        vector<vector<int>> RoyFloyd();

        vector<int> EulerianPath(int);

        int HamiltonianCycleCost();

        //Hopcroft-Karp
        vector<int> MaximumBipartiteMatching(int, int, int, vector<vector<int>>);
};

//LeetCode Critical Connections Problem
class Solution {
public:
    vector<vector<int>> CriticalConnections(int n, vector<vector<int>>& connections)
    {
        Graph G(n, connections.size(), 0, 0);
        G.BuildFromVector(connections);

        return G.CriticalConnections();
    }
};

int main()
{
/*  ----------------------------MaximumMatchingBipartiteGraph----------------------------
    int numberOfNodesLeft, numberOfNodesRight, numberOfEdges;
    vector<vector<int>> adjacencyListBipartiteGraph;
    vector<int> left;

    fin >> numberOfNodesLeft >> numberOfNodesRight >> numberOfEdges;

    for (int i = 0; i < numberOfNodesLeft + 1; i++)
        adjacencyListBipartiteGraph.push_back( {} );

    for (int i = 0; i < numberOfEdges; i++)
    {
        int node1, node2;
        fin >> node1 >> node2;

        adjacencyListBipartiteGraph[node1].push_back(node2);
    }

    Graph G(numberOfNodesLeft + numberOfNodesRight, numberOfEdges, 0, 1);

    left = G.MaximumBipartiteMatching(numberOfNodesLeft, numberOfNodesRight, numberOfEdges, adjacencyListBipartiteGraph);

    fout << left[0] << "\n";
    for (int i = 1; i < numberOfNodesLeft + 1; i++)
    {
        if (left[i])
            fout << i << " " << left[i] << "\n";
    }
    ----------------------------MaximumMatchingBipartiteGraph----------------------------*/



/*  ----------------------------HamiltonianCycle----------------------------
    int numberOfNodes, numberOfEdges, minimumCost;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1, 0);
    G.BuildTransposeWeightedGraph(fin);

    minimumCost = G.HamiltonianCycleCost();

    if (minimumCost != -1)
        fout << minimumCost;
    else
        fout << "Nu exista solutie";
    ----------------------------HamiltonianCycle----------------------------*/



/*  ----------------------------EulerianCycle----------------------------
        int numberOfNodes, numberOfEdges, startNode = 0;
        vector<int> path;

        fin >> numberOfNodes >> numberOfEdges;

        Graph G(numberOfNodes, numberOfEdges, 0, 1);
        G.BuildGraphWithEdgesNumber(fin);

        path = G.EulerianPath(startNode);

        for (unsigned int i = 0; i < path.size(); i++)
            fout << path[i] + 1 << " ";
    ----------------------------EulerianCycle----------------------------*/



/*  ----------------------------RoyFloyd----------------------------
    int numberOfNodes;
    fin >> numberOfNodes;

    Graph G(numberOfNodes, 1, 1);
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

    Graph G(numberOfNodes, numberOfEdges, 0, 1);
    G.Build(fin);

    fout << G.TreeDiameter();

    ----------------------------TreeDiameter----------------------------*/



/*  ----------------------------MaxFlow----------------------------

    int numberOfNodes, numberOfEdges, source, destination;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1, 1);
    G.BuildWeightedGraph(fin);

    source = 0;
    destination = G.GetNumberOfNodes() - 1;

    fout << G.MaxFlow(source, destination);

    ----------------------------MaxFlow----------------------------*/



/*  ----------------------------Disjoint----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 0, 1);

    G.DisjointSetForests(fin, fout);

    ----------------------------Disjoint----------------------------*/



/*  ----------------------------BellmanFord----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    bool negativeCycle = 0;
    vector<int> key(numberOfNodes);

    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1, 1);
    G.BuildWeightedGraph(fin);

    startNode = 0;
    key = G.BellmanFord(startNode, negativeCycle);

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

    ----------------------------BellmanFord----------------------------*/



/*  ----------------------------Dijkstra----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    vector<int> key(numberOfNodes);
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1, 1);
    G.BuildWeightedGraph(fin);

    startNode = 0;
    key = G.Dijkstra(startNode);

    for (int i = 0; i < numberOfNodes; i++)
    {
        if (i != startNode)
        {
            if (key[i] == INF)
                key[i] = 0;

            fout << key[i] << " ";
        }
    }

    ----------------------------Dijkstra----------------------------*/



/*  ----------------------------APM----------------------------

    int numberOfNodes, numberOfEdges, startNode, totalCost = 0, minimumSpanningTreeEdges = 0, fromOne = 1;
    vector<int> parent(numberOfNodes);

    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 0, fromOne);
    G.BuildWeightedGraph(fin);

    startNode = 0;
    parent = G.MinimumSpanningTree(startNode, totalCost, minimumSpanningTreeEdges);

    fout << totalCost << "\n" << minimumSpanningTreeEdges << "\n";

    if (fromOne)
        for (int i = 1; i < numberOfNodes; i++)
            fout << parent[i] + 1 << " " << i + 1 << "\n";
    else
        for (int i = 1; i < numberOfNodes; i++)
            fout << parent[i] << " " << i << "\n";

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

    Graph G(numberOfNodes, numberOfEdges, 1, 1);
    G.Build(fin);

    vector<int> topologicalSort = G.TopologicalSort();

    reverse(topologicalSort.begin(), topologicalSort.end());

    for (auto node : topologicalSort)
        fout << node + 1 << " ";

    ----------------------------TS----------------------------*/


/*  ----------------------------SCCs----------------------------

    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph G(numberOfNodes, numberOfEdges, 1, 1);
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

    Graph G(numberOfNodes, numberOfEdges, 0, 1);
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

    Graph G(numberOfNodes, numberOfEdges, 0, 1);
    G.Build(fin);

    fout << G.NumberOfConnectedComponents();

    ----------------------------DFS----------------------------*/

/*  ----------------------------BFS----------------------------

    int numberOfNodes, numberOfEdges, startNode;
    fin >> numberOfNodes >> numberOfEdges >> startNode;

    Graph G(numberOfNodes, numberOfEdges, 1, 1);
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

Graph::Graph() : numberOfNodes(0), numberOfEdges(0), directed(0), fromOne(0) { }

Graph::Graph(int _numberOfNodes, int _numberOfEdges, bool _directed, bool _fromOne) : numberOfNodes(_numberOfNodes), numberOfEdges(_numberOfEdges), directed(_directed), fromOne(_fromOne)
{
    for (int i = 0; i < _numberOfNodes; i++)
        adjacencyList.push_back( {} );
        //adjacencyList.push_back(vector<int>());
}

Graph::Graph(int _numberOfNodes, bool _directed, bool _fromOne) : numberOfNodes(_numberOfNodes), numberOfEdges(INF), directed(_directed), fromOne(_fromOne)
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

        if (fromOne)
        {
            firstNode--;
            secondNode--;
        }

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

        if (fromOne)
        {
            firstNode--;
            secondNode--;
        }

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

        if (fromOne)
        {
            firstNode--;
            secondNode--;
        }

        adjacencyListWithEdgesNumber[firstNode].push_back(make_pair(secondNode, i));

        if (!directed)
            adjacencyListWithEdgesNumber[secondNode].push_back(make_pair(firstNode, i));
    }
}

void Graph::BuildTransposeWeightedGraph(istream &fin)
{
    for (int i = 0; i < numberOfEdges; i++)
    {
        int firstNode, secondNode, weight;
        fin >> firstNode >> secondNode >> weight;

        if (fromOne)
        {
            firstNode--;
            secondNode--;
        }

        adjacencyListWeightedGraph[secondNode].push_back(make_pair(firstNode, weight));

        if (!directed)
            adjacencyListWeightedGraph[firstNode].push_back(make_pair(secondNode, weight));
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

bool Graph::IsDirected()
{
    return directed;
}

bool Graph::IsFromOne()
{
    return fromOne;
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

///     SHOW BFS ORDER FROM START NODE:
//        if (fromOne)
//            cout << currentNode + 1 << " ";
//        else
//            cout << currentNode << " ";

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

///     SHOW DFS ORDER FROM START NODE:
//        if (fromOne)
//            cout << startNode + 1 << " ";
//        else
//            cout << startNode << " ";

    for (unsigned int i = 0; i < adjacencyList[startNode].size(); i++)
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

void Graph::functionBiconnectedComponents(int node, int id, int ids[], int low[], stack<int> &stackBiconnectedComponents, bool visited[], vector<vector<int>> &biconnectedComponents)
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
            functionBiconnectedComponents(adjacentNode, id, ids, low, stackBiconnectedComponents, visited, biconnectedComponents);

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

void Graph::functionStronglyConnectedComponents(int node, int &id, int ids[], int low[], stack<int> &stackStronglyConnectedComponents, bool onStack[], bool visited[], vector<vector<int>> &stronglyConnectedComponents)
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
            functionStronglyConnectedComponents(adjacentNode, id, ids, low, stackStronglyConnectedComponents, onStack, visited, stronglyConnectedComponents);

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

    functionBiconnectedComponents(node, id, ids, low, stackBiconnectedComponents, visited, biconnectedComponents);

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
            functionStronglyConnectedComponents(i, id, ids, low, stackStronglyConnectedComponents, onStack, visited, stronglyConnectedComponents);

    return stronglyConnectedComponents;
}

void Graph::functionTopologicalSort(int startNode, bool visited[], vector<int> &topologicalSort)
{
    visited[startNode] = 1;

    for (unsigned int i = 0; i < adjacencyList[startNode].size(); i++)
    {
        int nextNode = adjacencyList[startNode][i];
        if (!visited[nextNode])
            functionTopologicalSort(nextNode, visited, topologicalSort);
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
            functionTopologicalSort(i, visited, topologicalSort);

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

void Graph::functionCriticalConnections(int node, bool visited[], int disc[], int low[], int parent[], vector<vector<int>> &criticalConnections)
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
            functionCriticalConnections(adjacentNode, visited, disc, low, parent, criticalConnections);

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
            functionCriticalConnections(i, visited, disc, low, parent, criticalConnections);

    return criticalConnections;
}

vector<int> Graph::MinimumSpanningTree(int startNode, int &totalCost, int &minimumSpanningTreeEdges)
{
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> minHeap;
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

    return parent;
}

vector<int> Graph::Dijkstra(int startNode)
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

    return key;
}

vector<int> Graph::BellmanFord(int startNode, bool &negativeCycle)
{
    vector<int> key(numberOfNodes, INF);
    queue<int> queueBF;
    long long int iterations = 0;

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

    return key;
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

int Graph::HamiltonianCycleCost()
{
    ///NOTE: Adjacency list (adjacencyListWeightedGraph) is for the Transpose Graph
    //       In this way, we know every node throw which we can reach a node x
    int minCost[1 << numberOfNodes][numberOfNodes], minimumCost = INF;

    for (int i = 0; i < (1 << numberOfNodes); i++)
        for (int j = 0; j < numberOfNodes ; j++)
            minCost[i][j] = INF;

    minCost[1][0] = 0;

    //for each sequence (written in bits)
    //                   e.g. sequence {2,3} is 01100
    //                        node 2 is 00100 (1<<2)
    for (int i = 0; i < (1 << numberOfNodes); i++)
    {
        //for each node we want to add in the hamiltonian path
        for (int j = 0; j < numberOfNodes; j++)
        {
            //if the node j is in sequence i
            if (i & (1 << j))
            {
                for (auto k : adjacencyListWeightedGraph[j])
                {
                    int adjacentNode = k.first;
                    int weight = k.second;

                    if (i & (1 << adjacentNode))
                        minCost[i][j] = min(minCost[i][j], minCost[i ^ (1 << j)][adjacentNode] + weight);
                }
            }
        }
    }

    for (auto i : adjacencyListWeightedGraph[0])
    {
        int adjacentNode = i.first;
        int weight = i.second;
        minimumCost = min(minimumCost, minCost[(1 << numberOfNodes) - 1][adjacentNode] + weight);
    }

    if (minimumCost == INF)
        return -1;

    return minimumCost;
}

bool Graph::MaximumBipartiteMatchingBFS(vector<vector<int>> &adjacencyListBipartiteGraph, vector<int> &left, vector<int> &right, vector<int> &distance)
{
    queue<int> queueMaximumMatching;

    for (unsigned int i = 1; i < left.size(); i++)
    {
        if (!left[i])
        {
            distance[i] = 0;
            queueMaximumMatching.push(i);
        }
        else
        {
            distance[i] = INF;
        }
    }

    distance[0] = INF;

    while(!queueMaximumMatching.empty())
    {
        int nodeLeft = queueMaximumMatching.front();
        queueMaximumMatching.pop();

        if (distance[nodeLeft] < distance[0])
        {
            for (auto nodeRight : adjacencyListBipartiteGraph[nodeLeft])
            {
                if (distance[right[nodeRight]] == INF)
                {
                    distance[right[nodeRight]] = distance[nodeLeft] + 1;
                    queueMaximumMatching.push(right[nodeRight]);
                }
            }
        }
    }

    if (distance[0] == INF)
        return 0;

    return 1;
}

bool Graph::MaximumBipartiteMatchingDFS(int startNodeLeft, vector<vector<int>> &adjacencyListBipartiteGraph, vector<int> &left, vector<int> &right, vector<int> &distance)
{
    if (!startNodeLeft)
        return 1;

    for (auto nodeRight : adjacencyListBipartiteGraph[startNodeLeft])
    {
        if (distance[right[nodeRight]] == distance[startNodeLeft] + 1)
        {
            if (MaximumBipartiteMatchingDFS(right[nodeRight], adjacencyListBipartiteGraph, left, right, distance))
            {
                right[nodeRight] = startNodeLeft;
                left[startNodeLeft] = nodeRight;
                return 1;
            }
        }
    }

    distance[startNodeLeft] = INF;

    return 0;
}

vector<int> Graph::MaximumBipartiteMatching(int numberOfNodesLeft, int numberOfNodesRight, int numberOfEdges, vector<vector<int>> adjacencyListBipartiteGraph)
{
    vector<int> left(numberOfNodesLeft + 1, 0), right(numberOfNodesRight + 1, 0), distance(numberOfNodesLeft + 1, 0);
    int maxNumberOfNodesMatched = 0;

    while (MaximumBipartiteMatchingBFS(adjacencyListBipartiteGraph, left, right, distance))
    {
        for (int i = 1; i < numberOfNodesLeft + 1; i++)
        {
            if (!left[i] && MaximumBipartiteMatchingDFS(i, adjacencyListBipartiteGraph, left, right, distance))
                maxNumberOfNodesMatched++;
        }
    }

    left[0] = maxNumberOfNodesMatched;

    return left;
}
