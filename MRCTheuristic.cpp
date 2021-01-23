/*
This code is my implementation from:
Masone, A., Nenni, M.E., Sforza, A. et al. The Minimum Routing Cost Tree Problem. 
Soft Comput 23, 2947–2957 (2019). https://doi.org/10.1007/s00500-018-3557-3
*/

#include <bits/stdc++.h>

using namespace std;

#define mp(a, b) make_pair(a, b)
#define EPS 1e-5

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

struct AdjInfo
{
    int v;
    double len;
    int id;
    AdjInfo(int v, double len, int id) : v(v), len(len), id(id) {}
};

struct Level
{
    int node, dIn;
    Level(int node = 0, int dIn = 0) : node(node), dIn(dIn) {}
};

int n; // number of vertices
int m; // number of edges
vector<Edge> edges; // edges given
vector<vector<AdjInfo>> adjList;
vector<vector<int>> n_;
vector<vector<int>> getID;
vector<vector<Level>> L;
vector<int> nodeLevel;
vector<int> subtreeSize;
vector<bool> seen;
int lmax;

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
        this->objective = 0;
        fill(seen.begin(), seen.end(), false);
        this->countSubTreeSz(0);
        Edge* e;
        int subSz;
        for(auto it = usedEdges.begin(); it != usedEdges.end(); ++it)
        {
            e = &edges[*it];
            subSz = min(subtreeSize[e->u], subtreeSize[e->v]);
            this->objective += e->len*(n-subSz)*subSz;
        }
    }

    int countSubTreeSz(int cur)
    {
        seen[cur] = true;
        int cnt = 1;
        for(AdjInfo& ainfo : adj[cur])
        {
            if(seen[ainfo.v]) 
                continue;
            cnt += countSubTreeSz(ainfo.v);
        }
        subtreeSize[cur] = cnt;
        return cnt;
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

bool comp(const Level& l1, const Level& l2)
{
    return l1.dIn > l2.dIn;
}

int DFS(int cur, vector<int> subTree, vector<bool>& seen, vector<int>* adj)
{
    seen[cur] = true;
    int sz, tmp;
    sz = 1;
    for(int& v : adj[cur])
    {
        if(seen[v])
            continue;
        tmp = DFS(v, subTree, seen, adj);
        n_[cur][v] += tmp;
        sz += tmp;
    }
    subTree[cur] = sz;
    return sz;
}

int mode;

void Phase1()
{
    int cur;
    n_.assign(n, vector<int>(n, 0));
    getID.assign(n, vector<int>(n, -1));
    set<int> V; // V'
    for(int i = 0; i < n; ++i)
    {
        V.insert(i);
        // perform Dijkstra in the node (i)
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

        // Mark edges from minimum path spanning tree and compute nij with DFS
        vector<int> adj[n];
        Edge e;
        for(int i = 0; i < n; ++i)
        {
            if(uEdge[i] > -1)
            {
                e = edges[uEdge[i]];
                adj[getNeighbor(i, e)].push_back(i);
                getID[getNeighbor(i, e)][i] = e.id;
            }
        }
        vector<int> subTree(n, 0);
        vector<bool> seen(n, false);
        DFS(i, subTree, seen, adj);
    }

    
    int l = 0;
    vector<int> dOut(n, 0);
    L.resize(n+1, vector<Level>());
    nodeLevel.assign(n, 0);
    vector<int> dIn(n, 0);
    // Compute dIn of all nodes
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            dIn[i] += n_[j][i];
        }
    }

    // separate the nodes in levels
    while(V.size())
    {
        l++;
        int dmin = INT_MAX;
        for(auto it = V.begin(); it != V.end(); ++it)
        {
            dOut[*it] = 0;
            for(auto it2 = V.begin(); it2 != V.end(); ++it2)
            {
                dOut[*it] += n_[*it][*it2];
            }
            dmin = min(dmin, dOut[*it]);
        }
        for(auto it = V.begin(); it != V.end(); ++it)
        {
            if(dOut[*it] == dmin)
            {
                L[l].push_back(Level(*it, dIn[*it]));
                nodeLevel[*it] = l;
            }
        }
        for(Level& lvl : L[l])
        {
            V.erase(lvl.node);
        }       
    }
    lmax = l;
}

void Phase2(Solution& sol)
{
    // sort nodes at each level by ascending inDegree
    for(int i = 1; i <= lmax; ++i)
    {
        sort(L[i].begin(), L[i].end(), comp);
        for(int j = 0; j < (int) L[i].size()-1; ++j)
        {
            assert(L[i][j].dIn >= L[i][j+1].dIn);
        }
    }
    int k, l1, l2, currentNode, edgeToAdd;
    double cmpVal;
    l1 = lmax;
    set<int> T;
    while((int) T.size() < n)
    {
        if(mode == 0)
            cmpVal = DBL_MAX;
        else if (mode == 1) 
            cmpVal = 0;
        else
            cmpVal = DBL_MAX;
        l2 = l1;
        do
        {
            edgeToAdd = -1;
            for(Level& lvl : L[l2])
            {
                k = lvl.node;
                if(T.find(k) != T.end())
                    continue;
                currentNode = -1;
                if(T.empty())
                {
                    T.insert(k);
                }
                else
                {
                    for(auto it = T.begin(); it != T.end(); ++it)
                    {
                        if(n_[*it][k] <= 0)
                            continue;
                        // Separated for each mode
                        if(mode == 0)
                        {
                            if(lt(edges[getID[*it][k]].len, cmpVal))
                            {
                                currentNode = k;
                                edgeToAdd = getID[*it][k];
                                cmpVal = edges[getID[*it][k]].len;
                            }
                        }
                        else if(mode == 1)
                        {
                            if(lt(cmpVal, n_[*it][k]))
                            {
                                currentNode = k;
                                edgeToAdd = getID[*it][k];
                                cmpVal = n_[*it][k];
                            }
                        }
                        else
                        {
                            if(lt((double) edges[getID[*it][k]].len/n_[*it][k], cmpVal))
                            {
                                currentNode = k;
                                edgeToAdd = getID[*it][k];
                                cmpVal = (double) edges[getID[*it][k]].len/n_[*it][k];
                            }
                        }
                    }
                    if(currentNode != -1)
                        break;
                }
            }
            l2--;
        } while (currentNode == -1);
        if(edgeToAdd > -1)
        {
            T.insert(currentNode);
            Edge* e = &edges[edgeToAdd];
            sol.usedEdges.insert(e->id);
            sol.adj[e->u].push_back(AdjInfo(e->v, e->len, e->id));
            sol.adj[e->v].push_back(AdjInfo(e->u, e->len, e->id));
        }
        bool done = true;
        for(Level& lvl : L[l1])
        {
            if(T.find(lvl.node) == T.end())
            {
                done = false;
            }
        }
        if(done)
            l1--;
    }
    sol.computeObjectiveFun();
    printf("%.10f\n", sol.objective);
}

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        printf("usage: ./heuristic mode < inputFile\n");
        return -1;
    }
    mode = atoi(argv[1]);
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
    // ignore requirement values (is always 1)
    seen.resize(n);
    subtreeSize.resize(n);
    Phase1();
    Solution sol;
    Phase2(sol);
}
