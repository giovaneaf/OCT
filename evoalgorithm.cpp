#include <bits/stdc++.h>
#include <gurobi_c++.h>

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

// seed used to generate random numbers
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//unsigned seed = 369140336;

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

// Return the neighbor of node u for a given edge
inline int getNeighbor(int u, Edge& e)
{
    return (e.u == u ? e.v : e.u);
}

// Removes edge from the solution (doesn't recompute anything)
void removeEdge(Edge& edge, Solution& s)
{
    s.usedEdge[edge.id] = false;
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
}

// IMPORTANT: Assertions should be removed when testing

// evolutionary/genetic algorithm
struct Evolutionary
{
    vector<Solution> solutions;
    Evolutionary(int popSize)
    {
        solutions.resize(popSize);
        genInitialPop(popSize);
        mutateInserting(solutions[0]);
        mutateRemoving(solutions[0]);
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
        // Selecting edge to insert
        vi possibleEdges(m-(n-1));
        int idx = 0;
        for(int i = 0; i < m; ++i)
        {
            if(!s.usedEdge[i])
            {
                possibleEdges[idx++] = i;
            }
        }
        int rdInt = possibleEdges[rand()%((int) possibleEdges.size())];
        Edge edge = edges[rdInt];
        // Find cycle in graph with BFS (path from one endpoint to the other)
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
                if(ainfo.v == edge.v)   // path found
                    break;
            }
        }
        while(q.size()) // Empty queue for later usage
            q.pop();
        // insert chosen edge in tree
        s.usedEdge[rdInt] = true;
        s.adj[edge.u].push_back(AdjInfo(edge.v, edge.len, edge.id));
        s.adj[edge.v].push_back(AdjInfo(edge.u, edge.len, edge.id));
        cur = edge.v;
        int e1, e2, rmEdge, minObj, curObj, curNode;
        e1 = e2 = rmEdge = edge.id;         // e2 represents the removed edge
        minObj = curObj = s.objective;
        // traverse all edges in cycle
        while(cur != edge.u)
        {
            e1 = e2;
            e2 = eIdx[cur];
            vb updated(n, false);
            q.push(cur);
            // find all nodes that distances are outdated
            while(q.size())
            {
                curNode = q.front();
                updated[curNode] = true;
                q.pop();
                for(AdjInfo& ainfo: s.adj[curNode])
                {
                    if(updated[ainfo.v] || ainfo.id == e1 || ainfo.id == e2)
                        continue;
                    q.push(ainfo.v);
                    updated[ainfo.v] = true;
                }
            }
            // update the distances of the values doing BFS and updating the objective function
            Edge newEdge = edges[e1];
            int neighbor = getNeighbor(cur, newEdge);
            curNode = cur;
            list<int> nodesToUpdate;
            for(int i = 0; i < n; ++i)
            {
                if(!updated[i])
                {
                    curObj -= s.dist[curNode][i]*req[curNode][i];
                    s.dist[curNode][i] = s.dist[i][curNode] = newEdge.len + s.dist[neighbor][i];
                    curObj += s.dist[curNode][i]*req[curNode][i];
                    nodesToUpdate.push_back(i);
                }
            }
            vb seen(n, false);
            q.push(curNode);
            // update remaining nodes
            while(q.size())
            {
                curNode = q.front();
                seen[curNode] = true;
                q.pop();
                for(AdjInfo& ainfo: s.adj[curNode])
                {
                    if(seen[ainfo.v] || ainfo.id == e1 || ainfo.id == e2)
                        continue;
                    q.push(ainfo.v);
                    seen[ainfo.v] = true;
                    for(int& i : nodesToUpdate)
                    {
                        curObj -= s.dist[ainfo.v][i]*req[ainfo.v][i];
                        s.dist[ainfo.v][i] = s.dist[i][ainfo.v] = ainfo.len + s.dist[curNode][i];
                        curObj += s.dist[ainfo.v][i]*req[ainfo.v][i];
                    }
                }
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
        assert(edge.id == rmEdge);
        removeEdge(edge, s);
        // call this to update the new distances correctly
        computeObjectiveFun(s);
        assert(minObj == s.objective);
    }

    // Mutate when considering to remove a random edge - O(m*n^2)
    void mutateRemoving(Solution& s)
    {
        // selecting edge to remove
        vi possibleEdges(n-1);
        int idx = 0;
        for(int i = 0; i < m; ++i)
        {
            if(s.usedEdge[i])
            {
                possibleEdges[idx++] = i;
            }
        }
        int rdInt = possibleEdges[rand()%((int) possibleEdges.size())];
        Edge edge = edges[rdInt];
        // remove chosen edge
        removeEdge(edge, s);
        // find nodes in one set when removing the chosen edge
        vi inA(n, 0);
        queue<int> q;
        q.push(edge.u);
        int curNode;
        int szA = 0;
        while(q.size())
        {
            curNode = q.front();
            inA[curNode] = 1;
            szA++;
            q.pop();
            for(AdjInfo& ainfo: s.adj[curNode])
            {
                assert(ainfo.id != edge.id); // assert edge was removed
                if(inA[ainfo.v])
                    continue;
                inA[ainfo.v] = 1;
                q.push(ainfo.v);
            }
        }
        // Find the best edge to add when removing the chosen edge
        int minObj, curObj, addEdge;
        minObj = curObj = s.objective;
        addEdge = rdInt;
        for(int i = 0; i < m; ++i)
        {
            // XOR is used to ensure the edge is connecting the two disconnected sets A and B
            if(!s.usedEdge[i] && (inA[edges[i].u]^inA[edges[i].v]))
            {
                curNode = edges[i].u;
                list<int> nodesToUpdate;
                for(int j = 0; j < n; ++j)
                {
                    if(inA[curNode]^inA[j]) // need to be updated
                    {
                        curObj -= s.dist[curNode][j]*req[curNode][j];
                        s.dist[curNode][j] = s.dist[j][curNode] = edges[i].len + s.dist[edges[i].v][j];
                        curObj += s.dist[curNode][j]*req[curNode][j];
                        nodesToUpdate.push_back(j);
                    }
                }
                vb seen(n, false);
                q.push(curNode);
                // update the distance values from one set to the other
                while(q.size())
                {
                    curNode = q.front();
                    seen[curNode] = true;
                    q.pop();
                    for(AdjInfo& ainfo: s.adj[curNode])
                    {
                        if(seen[ainfo.v])
                            continue;
                        assert((inA[curNode]^inA[ainfo.v]) == 0);
                        q.push(ainfo.v);
                        seen[ainfo.v] = true;
                        for(int& j : nodesToUpdate)
                        {
                            curObj -= s.dist[ainfo.v][j]*req[ainfo.v][j];
                            s.dist[ainfo.v][j] = s.dist[j][ainfo.v] = ainfo.len + s.dist[curNode][j];
                            curObj += s.dist[ainfo.v][j]*req[ainfo.v][j];
                        }
                    }
                }
                // assertion to check if curObj value is correct
                Solution ss = s;
                ss.usedEdge[i] = true;
                ss.adj[edges[i].u].push_back(AdjInfo(edges[i].v, edges[i].len, edges[i].id));
                ss.adj[edges[i].v].push_back(AdjInfo(edges[i].u, edges[i].len, edges[i].id));
                computeObjectiveFun(ss);
                assert(ss.objective == curObj);
                if(curObj < minObj)
                {
                    minObj = curObj;
                    addEdge = i;
                }
            }
        }
        // Insert the best edge in solution
        Edge bestEdge = edges[addEdge];
        s.usedEdge[addEdge] = true;
        s.adj[bestEdge.u].push_back(AdjInfo(bestEdge.v, bestEdge.len, bestEdge.id));
        s.adj[bestEdge.v].push_back(AdjInfo(bestEdge.u, bestEdge.len, bestEdge.id));
        // update solution objective function and distances
        s.objective = minObj;
        assert(inA[bestEdge.u]^inA[bestEdge.v]); // assert that the edge form a tree
        curNode = bestEdge.u;
        int neighbor = getNeighbor(curNode, bestEdge);
        for(int i = 0; i < n; ++i)
        {
            if(inA[curNode]^inA[i]) // if the values are updated by the edge
                s.dist[curNode][i] = s.dist[i][curNode] = bestEdge.len + s.dist[neighbor][i];
        }
        vb seen(n, false);
        q.push(curNode);
        while(q.size())
        {
            curNode = q.front();
            seen[curNode] = true;
            q.pop();
            for(AdjInfo& ainfo: s.adj[curNode])
            {
                if(seen[ainfo.v] || (inA[curNode]^inA[ainfo.v]))
                    continue;
                assert((inA[curNode]^inA[ainfo.v]) == 0);
                q.push(ainfo.v);
                seen[ainfo.v] = true;
                for(int i = 0; i < n; ++i)
                {
                    if(inA[curNode]^inA[i])
                    {
                        s.dist[ainfo.v][i] = s.dist[i][ainfo.v] = ainfo.len + s.dist[curNode][i];
                    }
                }
            }
        }
        // assertion to check if distances are correctly updated
        int tmp = 0;
        for(int i = 0; i < n; ++i)
            for(int j = i+1; j < n; ++j)
                tmp += s.dist[i][j]*req[i][j];
        assert(tmp == s.objective);
    }

};

int main()
{
printf("seed = %u\n", seed);
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