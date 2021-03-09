#include <bits/stdc++.h>
//#include <gurobi_c++.h>

using namespace std;

#define mp(a, b) make_pair(a, b)
#define vb vector<bool>
#define vi vector<int>
#define ii pair<int, int>
#define EPS 1e-4

int mode;

struct Edge 
{
    int u, v; 
    double len;
    int id;
    Edge(int u = 0, int v = 0, double len = 0.0, int id = 0) : u(u), v(v), len(len), id(id) {}
    bool operator< (Edge& e)
    {
        return len < e.len;
    }
};

// Union find used for cycle detection efficiently
struct UnionFind
{
    vi pset;   // who is the pointer of the set
    UnionFind(int n) 
    {
        pset.assign(n, 0);
        for(int i = 0; i < n; ++i) pset[i] = i;
    }
    int findSet(int i)
    {
        return (pset[i] == i ? i : (pset[i] = findSet(pset[i])));
    }
    bool isSameSet(int i, int j)
    {
        return findSet(i) == findSet(j);
    }
    void unionSet(int i, int j)     // make set i point to j
    {
        if(isSameSet(i, j)) return;
        pset[findSet(i)] = findSet(j);
    }
};

struct AdjInfo
{
    int v;
    double len;
    int id;
    AdjInfo(int v, double len, int id) : v(v), len(len), id(id) {}
};

/*
The file is read in the following format
n m
list of m with two nodes and length. (u v c)
combination(n, 2) lines of requirement values (r01, r02, r03, ..., r0n, r12, r13, ..., rn-1n)
*/

int n; // number of vertices
int m; // number of edges
int K; // max size of each group
vector<Edge> edges; // edges given
vector<vector<double>> req;  // requirement values
vector<vector<AdjInfo>> adjList; // adjacency list

// Return the neighbor of node u for a given edge
inline int getNeighbor(int u, Edge& e)
{
    return (e.u == u ? e.v : e.u);
}

bool inline eq(double a, double b)
{
    return abs(a-b) < EPS;
}
bool inline leq(double a, double b)
{
    return a < b || abs(a-b) < EPS;
}
bool inline lt(double a, double b)
{
    return a < b && abs(a-b) > EPS;
}

// Stores candidate solution (tree)
vector<vector<double>> dist;
vector<bool> seen;
struct Solution
{
    vector<vector<AdjInfo>> adj;
    set<int> usedEdges;
    double objective;
    
    Solution()
    {
        adj.resize(n, vector<AdjInfo>());
        objective = 0;
    }
    
    void clear()
    {
        for(int i = 0; i < n; ++i)
        {
            adj[i].clear();
        }
        usedEdges.clear();
        objective = 0;
    }

    void computeObjectiveFun()
    {
        int cur;
        for(int i = 0; i < n; ++i)
        {
            fill(dist[i].begin(), dist[i].end(), DBL_MAX);
        }
        this->objective = 0;
        // BFS for each node to compute the distances
        for(int node = 0; node < n; ++node)
        {
            dist[node][node] = 0;
            queue<int> q;
            q.push(node);
            while(q.size())
            {
                cur = q.front();
                q.pop();
                for(AdjInfo& ainfo: this->adj[cur])
                {
                    if(dist[node][ainfo.v] == DBL_MAX)
                    {
                        dist[node][ainfo.v] = ainfo.len + dist[node][cur];
                        q.push(ainfo.v);
                    }
                }
            }
            for(int v = node+1; v < n; v++)
            {
                objective += dist[node][v]*req[node][v];
            }
        }
    }

    inline bool hasEdge(int idx)
    {
        return usedEdges.find(idx) != usedEdges.end();
    }

    // Removes edge from the solution (doesn't recompute anything)
    void removeEdge(Edge& edge)
    {
        //this->usedEdge[edge.id] = false;
        usedEdges.erase(edge.id);
        for(auto it = this->adj[edge.u].begin(); it !=  this->adj[edge.u].end(); ++it)
        {
            if(it->id == edge.id)
            {
                this->adj[edge.u].erase(it);
                break;
            }
        }
        for(auto it = this->adj[edge.v].begin(); it !=  this->adj[edge.v].end(); ++it)
        {
            if(it->id == edge.id)
            {
                this->adj[edge.v].erase(it);
                break;
            }
        }
    }

};


// printing functions for debugging only purpose
inline void print(Edge& e)
{
    printf("(%d, %d, %.2f, %d)\n", e.u, e.v, e.len, e.id);
}
void print(vector<Edge>& edges)
{
    int cnt = 0;
    for(Edge& e: edges)
    {
        printf("%d: ", cnt++);
        print(e);
    }
    putchar('\n');
}
void print(Solution& s)
{
    printf("Edges used:\n");
    for(auto it = s.usedEdges.begin(); it != s.usedEdges.end(); ++it)
    {
        print(edges[*it]);
    }
    putchar('\n');
    printf("Objective value = %.2f\n", s.objective);
    putchar('\n');
}

struct Info
{
    int id;
    bool legit;
    Info(int id = -1, int legit = false) : id(id), legit(legit) {}
};

struct Group 
{
    vector<Info> node;
};

vector<Group> group;

/* Find best solution with minimal path trees using Dijkstra */
void findInitialSolution(Solution& best)
{
    best.objective = DBL_MAX;
    int cur;
    // Generate Min Path Tree solution
    for(int i = 0; i < n; ++i)
    {
        // perform Dijkstra in the node i
        vector<double> dist(n, DBL_MAX);
        vector<int> uEdge(n, -1);
        cur = i;
        dist[cur] = 0.0;
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        pq.push(mp(dist[cur], cur));
        double w;
        while(pq.size())
        {
            cur = pq.top().second;
            w = pq.top().first;
            pq.pop();
            if(lt(dist[cur], w))
                continue;
            for(AdjInfo& ainfo : adjList[cur])
            {
                if(dist[ainfo.v] > dist[cur] + ainfo.len)
                {
                    dist[ainfo.v] = dist[cur] + ainfo.len;
                    uEdge[ainfo.v] = ainfo.id;
                    pq.push(mp(dist[ainfo.v], ainfo.v));
                }
            }
        }
        // construct Solution for minimum path tree from node
        Solution sol;
        Edge e;
        int cnt = 0;
        for(int& edgeID : uEdge)
        {
            if(edgeID > -1)
            {
                //sol.usedEdge[edgeID] = true;
                sol.usedEdges.insert(edgeID);
                e = edges[edgeID];
                sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                cnt++;
            }
        }
        assert(cnt == n-1);
        sol.computeObjectiveFun();
        if(lt(sol.objective, best.objective))
        {
            best = sol;
        }
    }
    printf("Best Obj = %.10f\n", best.objective);
}

int computeSubtree(int cur, vi& subtreeSize, vi& dad, Solution& best)
{
    subtreeSize[cur] = 1;
    for(AdjInfo& ainfo : best.adj[cur])
    {
        if(dad[ainfo.v] > -1 || ainfo.v == 0)
            continue;
        dad[ainfo.v] = cur;
        computeSubtree(ainfo.v, subtreeSize, dad, best);
        subtreeSize[cur] += subtreeSize[ainfo.v];
    }
    return subtreeSize[cur];
}

void DFS(int cur, int invalid, Group& gp, Solution& best, vi& subtreeSize)
{
    seen[cur] = true;
    gp.node.push_back(Info(cur, true));
    for(AdjInfo& ainfo : best.adj[cur])
    {
        if(ainfo.v == invalid || seen[ainfo.v]) 
            continue;
        subtreeSize[ainfo.v] = 0;
        DFS(ainfo.v, invalid, gp, best, subtreeSize);
    }
}

void divideGroups(Solution& best)
{
    vector<int> dad(n, -1);
    vector<int> subtreeSize(n, 0);
    computeSubtree(0, subtreeSize, dad, best);
    fill(seen.begin(), seen.end(), false);
    queue<int> q;
    for(int i = 0; i < n; ++i)
    {
        if(subtreeSize[i] == K)
            q.push(i);
    }
    int cur, tmpPar;
    while(q.size())
    {
        while(q.size())
        {
            cur = q.front();
            q.pop();
            Group gp;
            DFS(cur, dad[cur], gp, best, subtreeSize);
            group.push_back(gp);
            tmpPar = dad[cur];
            while(tmpPar != -1)
            {
                subtreeSize[tmpPar] -= subtreeSize[cur];
                tmpPar = dad[tmpPar];
            }
            subtreeSize[cur] = 0;
        }
        for(int i = 0; i < n; ++i)
        {
            if(!seen[i] && subtreeSize[i] == K)
                q.push(i);
        }
    }
    set<ii> order;
    for(int i = 0; i < n; ++i)
    {
        if(!seen[i])
        {
            if(subtreeSize[i] >= K)
            {
                order.insert({subtreeSize[i], i});
            }
        }
    }
    while(order.size())
    {
        cur = order.begin()->second;
        order.erase(order.begin());
        if(seen[cur])
            continue;
        vector<ii> subtrees;
        for(AdjInfo& ainfo : best.adj[cur])
        {
            if(subtreeSize[ainfo.v] == 0 || subtreeSize[ainfo.v] > subtreeSize[cur])
                continue;
            subtrees.push_back({subtreeSize[ainfo.v], ainfo.v});
            assert(subtreeSize[ainfo.v] < K);
        }
        sort(subtrees.begin(), subtrees.end());
        if((int) subtrees.size() == 1 || subtreeSize[cur] == K)
        {
            assert(subtreeSize[cur] == K);
            Group gp;
            DFS(cur, dad[cur], gp, best, subtreeSize);
            group.push_back(gp);
            tmpPar = cur;
            while(tmpPar != -1)
            {
                auto it = order.find({subtreeSize[tmpPar], tmpPar});
                subtreeSize[tmpPar] -= K;
                if(it != order.end())
                {
                    order.erase(it);
                    if(subtreeSize[tmpPar] >= K)
                    {
                        order.insert({subtreeSize[tmpPar], tmpPar});
                    }
                }
                tmpPar = dad[tmpPar];
            }
            continue;
        }
        assert((int) subtrees.size() > 1 && subtreeSize[cur] > K);
        int pos, ipos, numNodes, rmNodes;
        pos = 0;
        rmNodes = 0;
        while(subtreeSize[cur] - rmNodes > K)
        {
            ipos = pos;
            numNodes = subtrees[ipos].first + 1;
            pos++;
            while(pos < (int) subtrees.size() && numNodes + subtrees[pos].first <= K)
            {
                numNodes += subtrees[pos].first;
                pos++;
            }
            if(pos-ipos == 1)
            {
                break;
            }
            rmNodes += numNodes-1;
            Group gp;
            for(int i = ipos; i < pos; ++i)
            {
                DFS(subtrees[i].second, dad[subtrees[i].second], gp, best, subtreeSize);
            }
            gp.node.push_back(Info(dad[subtrees[ipos].second], false));
            group.push_back(gp);
        }
        int revPos = (int) subtrees.size()-1;
        while(revPos >= pos && subtreeSize[cur] - rmNodes > K)
        {
            rmNodes += subtrees[revPos].first;
            Group gp;
            DFS(subtrees[revPos].second, dad[subtrees[revPos].second], gp, best, subtreeSize);
            group.push_back(gp);
            revPos--;
        }
        if(subtreeSize[cur] - rmNodes == K)
        {
            rmNodes = subtreeSize[cur];
            Group gp;
            DFS(cur, dad[cur], gp, best, subtreeSize);
            group.push_back(gp);
        }
        tmpPar = cur;
        while(tmpPar != -1)
        {
            auto it = order.find({subtreeSize[tmpPar], tmpPar});
            subtreeSize[tmpPar] -= rmNodes;
            if(it != order.end())
            {
                order.erase(it);
                if(subtreeSize[tmpPar] >= K)
                {
                    order.insert({subtreeSize[tmpPar], tmpPar});
                }
            }
            tmpPar = dad[tmpPar];
        }
        for(int i = 0; i < n; ++i)
        {
            if(!seen[i])
            {
                if(subtreeSize[i] >= K)
                {
                    order.insert({subtreeSize[i], i});
                }
            }
        }
    }
    if(!seen[0])
    {
        Group gp;
        DFS(0, -1, gp, best, subtreeSize);
        group.push_back(gp);
        subtreeSize[0] = 0;
        putchar('\n');
    }
    int totNodes = 0;
    for(Group& gp : group)
    {
        printf("group:");
        for(Info& node : gp.node)
        {
            printf(" %d", node.id);
            totNodes++;
        }
        putchar('\n');
    }
    printf("total Nodes = %d\n", totNodes);
    printf("total Groups = %d\n", (int) group.size());
}

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        printf("usage: ./mathheuristic K < inputFile\n");
        return -1;
    }
    K = atoi(argv[1]);
    printf("K = %d\n", K);
    cin >> n >> m;
    edges.resize(m);
    adjList.resize(n, vector<AdjInfo>());
    for(int i = 0; i < m; ++i)
    {
        cin >> edges[i].u >> edges[i].v >> edges[i].len;
        edges[i].id = i;
        adjList[edges[i].u].push_back(AdjInfo(edges[i].v, edges[i].len, edges[i].id));
        adjList[edges[i].v].push_back(AdjInfo(edges[i].u, edges[i].len, edges[i].id));
    }
    //getIdxFlow.resize(n, vector<int>(2*m));
    seen.resize(n);
    dist.resize(n, vector<double>(n));
    req.resize(n, vector<double>(n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            cin >> req[i][j];
            req[j][i] = req[i][j];
        }
    }
    Solution best;
    findInitialSolution(best);
    divideGroups(best);
    return 0;
}