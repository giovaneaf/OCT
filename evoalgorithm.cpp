#include <bits/stdc++.h>

using namespace std;

#define vb vector<bool>
#define vi vector<int>
#define vvi vector<vi>

struct Edge 
{
    int u, v, len, id;
    Edge(int u = 0, int v = 0, int len = 0, int id = 0) : u(u), v(v), len(len), id(id) {}
};

/*
The file is read in the following format
n m
list of m with two nodes and length. (u v c)
combination(n, 2) lines of requirement values (r01, r02, r03, ..., r0n, r12, r13, ..., rn-1n)
*/

int n; // number of vertices
int m; // number of edges
vector<Edge> edges; // edges given
vvi req; // requirement values

// seed used to generate random numbers and shuffle array
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//unsigned seed = 4116840420;

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
    int v, len, id;
    AdjInfo(int v, int len, int id) : v(v), len(len), id(id) {}
};

// Stores candidate solution (tree)
struct Solution
{
    vector<vector<AdjInfo>> adj;
    vvi dist;
    vb usedEdge;
    int objective;
    float fitness;
    Solution()
    {
        adj.resize(n, vector<AdjInfo>());
        dist.resize(n, vi(n));
        usedEdge.resize(m, false);
        objective = 0;
        fitness = 0.0;
    }
};

// printing functions for debugging only purpose
inline void print(Edge& e)
{
    printf("(%d, %d, %d, %d)\n", e.u, e.v, e.len, e.id);
}
void print(vector<Edge>& edges)
{
    for(Edge& e: edges)
    {
        print(e);
    }
    putchar('\n');
}
void print(Solution& s)
{
    printf("Edges used:\n");
    for(int i = 0; i < m; ++i)
    {
        if(s.usedEdge[i])
        {
            print(edges[i]);
        }
    }
    putchar('\n');
    printf("Objective value = %d\n", s.objective);
    printf("Fitness value = %f\n", s.fitness);
    putchar('\n');
}
void print(vb& usedEdges)
{
    for(int i = 0; i < m; ++i)
    {
        if(usedEdges[i])
        {
            print(edges[i]);
        }
    }
}

// evolutionary/genetic algorithm
struct Evolutionary
{
    vector<Solution> solutions;
    Evolutionary(int popSize)
    {
        solutions.resize(popSize);
        genInitialPop(popSize);
        mutateInserting(solutions[0]);
    }

    /* Generate popSize initial solutions (trees) by shuffling the edges
	   and inserting the edges like Kruskal Algorithm */
    void genInitialPop(int popSize)
    {
        vector<Edge> cpy = edges;
        int minObj, maxObj, numForests, idx;
        minObj = INT_MAX;
        maxObj = INT_MIN;
        for(int i = 0; i < popSize; ++i)
        {
            shuffle(begin(cpy), end(cpy), default_random_engine(seed));
            UnionFind uf(n);
            numForests = n;
            Solution sol;
            for(Edge& e: cpy)
            {
                if(!uf.isSameSet(e.u, e.v))
                {
                    uf.unionSet(e.u, e.v);
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    sol.usedEdge[e.id] = true;
                    numForests--;
                }
                if(numForests == 1) // If the tree is done
                {
                    computeObjectiveFun(sol);
                    minObj = min(minObj, sol.objective);
                    maxObj = max(maxObj, sol.objective);
                    break;
                }
            }
            solutions[i] = sol;
        }
        for(int i = 0; i < popSize; ++i)
        {
            if(minObj == maxObj)
                solutions[i].fitness = 1.0;
            else
                solutions[i].fitness = 1.0 - ((float) solutions[i].objective - minObj)/(maxObj - minObj);
        }
    }
	/* Input: Adjacency list of the tree
	   Output: Objective function value */
    void computeObjectiveFun(Solution& s)
    {
        int cur;
        for(int i = 0; i < n; ++i)
        {
            fill(s.dist[i].begin(), s.dist[i].end(), INT_MAX);
        }
        s.objective = 0;
        // BFS for each node to compute the distances
        for(int node = 0; node < n; ++node)
        {
            s.dist[node][node] = 0;
            queue<int> q;
            q.push(node);
            while(q.size())
            {
                cur = q.front();
                q.pop();
                for(AdjInfo& ainfo: s.adj[cur])
                {
                    if(s.dist[node][ainfo.v] == INT_MAX)
                    {
                        s.dist[node][ainfo.v] = ainfo.len + s.dist[node][cur];
                        q.push(ainfo.v);
                    }
                }
            }
            for(int v = node+1; v < n; v++)
            {
                s.objective += s.dist[node][v]*req[node][v];
            }
        }
    }

    // Mutate when inserting a new edge in the solution - O(n^2)
    void mutateInserting(Solution& s)
    {
        vi possibleEdges(m-(n-1));
        int idx = 0;
        for(int i = 0; i < m; ++i)
        {
            if(!s.usedEdge[i])
            {
                possibleEdges[idx++] = i;
            }
        }
        int rdint = possibleEdges[rand()%(m-(n-1))];
        Edge edge = edges[rdint];
        // Find cycle in graph with BFS
        vi dad(n, -1);
        vi eIdx(n, -1);
        queue<int> q;
        q.push(edge.u);
        int cur;
        while(q.size())
        {
            cur = q.front();
            q.pop();
            for(AdjInfo& ainfo: s.adj[cur])
            {
                if(dad[ainfo.v] == -1)
                {
                    dad[ainfo.v] = cur;
                    eIdx[ainfo.v] = ainfo.id;
                    q.push(ainfo.v);
                }
                if(ainfo.v == edge.v)
                    break;
            }
        }
        while(q.size()) // Assert queue is empty
            q.pop();
        // insert chosen edge in tree
        s.usedEdge[rdint] = true;
        s.adj[edge.u].push_back(AdjInfo(edge.v, edge.len, edge.id));
        s.adj[edge.v].push_back(AdjInfo(edge.u, edge.len, edge.id));
        cur = edge.v;
        int e1, e2, rmEdge, minObj, curObj, curNode;
        e1 = e2 = rmEdge = edge.id;
        minObj = curObj = s.objective;
        // traverse all edges in cycle
        while(cur != edge.u)
        {
            e1 = e2;
            e2 = eIdx[cur];
            vb seen(n, false);
            q.push(cur);
            // find all nodes that distances are outdated
            list<int> nodesToUpdate;
            while(q.size())
            {
                curNode = q.front();
                seen[curNode] = true;
                nodesToUpdate.push_back(curNode);
                q.pop();
                for(AdjInfo& ainfo: s.adj[curNode])
                {
                    if(seen[ainfo.v] || ainfo.id == e1 || ainfo.id == e2)
                        continue;
                    q.push(ainfo.v);
                    seen[ainfo.v] = true;
                }
            }
            // update the distances of the values doing BFS and updating the objective function
            for(int& node: nodesToUpdate)
            {
                for(int i = 0; i < n; ++i)
                    curObj -= s.dist[node][i]*req[node][i];
                fill(s.dist[node].begin(), s.dist[node].end(), INT_MAX);
                s.dist[node][node] = 0;
                q.push(node);
                
                while(q.size())
                {
                    curNode = q.front();
                    q.pop();
                    for(AdjInfo& ainfo: s.adj[curNode])
                    {
                        if(ainfo.id == e2)
                            continue;
                        
                        if(s.dist[node][ainfo.v] == INT_MAX)
                        {
                            s.dist[node][ainfo.v] = s.dist[ainfo.v][node] = ainfo.len + s.dist[node][curNode];
                            q.push(ainfo.v);
                        }
                    }
                }
                for(int i = 0; i < n; ++i)
                    curObj += s.dist[node][i]*req[node][i];
            }
            // after updating check if the objective function is lower than the last seen
            if(curObj < minObj)
            {
                minObj = curObj;
                rmEdge = e2;
            }
            cur = dad[cur];
        }
        // remove the edge s.t. objective function is minimum
        edge = edges[rmEdge];
        s.usedEdge[rmEdge] = false;
        for(auto it = s.adj[edge.u].begin(); it !=  s.adj[edge.u].end(); ++it)
        {
            if(it->id == edge.id)
            {
                s.adj[edge.u].erase(it);
                break;
            }
        }
        for(auto it = s.adj[edge.v].begin(); it !=  s.adj[edge.v].end(); ++it)
        {
            if(it->id == edge.id)
            {
                s.adj[edge.v].erase(it);
                break;
            }
        }
        // call this to update the new distances correctly
        computeObjectiveFun(s);
        assert(minObj == s.objective);
    }

};

int main()
{
    srand(seed);
    cin >> n >> m;
    edges.resize(m);
    for(int i = 0; i < m; ++i)
    {
        cin >> edges[i].u >> edges[i].v >> edges[i].len;
        edges[i].id = i;
    }
    req.resize(n, vi(n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            cin >> req[i][j];
            req[j][i] = req[i][j];
        }
    }
    Evolutionary ev(5);
    return 0;
}