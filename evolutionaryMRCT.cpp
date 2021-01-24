#include <bits/stdc++.h>

using namespace std;

#define mp(a, b) make_pair(a, b)
#define vb vector<bool>
#define vi vector<int>
#define ii pair<int, int>
#define EPS 1e-3

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
vector<Edge> edges; // edges given
vector<vector<double>> req; // requirement values
vector<vector<AdjInfo>> adjList; // used for PTAS crossover

// seed used to generate random numbers
unsigned seed;
// seeds used for testing
unsigned seedVector[] = {280192806, 871237442, 2540188929, 107472404, 3957311442, 316851227, 619606212, 1078082709, 916212990, 698598169};
//Mersenne Twister: Good quality random number generator
std::mt19937 rng;
map<ii, Edge*> edgeMap;
map<int, list<vector<int>>> prufferCodes;

// Return the neighbor of node u for a given edge
inline int getNeighbor(int u, Edge& e)
{
    return (e.u == u ? e.v : e.u);
}

void createPruffer(vector<int>& pruffer, list<vector<int>>& trees, int cur, int K)
{
    if((int) pruffer.size() == cur)
    {
        trees.push_back(pruffer);
        return ;
    }
    for(int i = 0; i < K; ++i)
    {
        pruffer[cur] = i;
        createPruffer(pruffer, trees, cur+1, K);
    }
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
vector<int> subtreeSize;
vector<bool> seen;

struct Solution;

void neighbor(int i, Solution& sol);

struct Solution
{
    vector<vector<AdjInfo>> adj;
    //vector<vector<double>> dist;
    //vb usedEdge;
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

    void fillDist()
    {
        int cur;
        for(int i = 0; i < n; ++i)
        {
            fill(dist[i].begin(), dist[i].end(), DBL_MAX);
        }
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
        }
    }

    inline bool hasEdge(int idx)
    {
        return usedEdges.find(idx) != usedEdges.end();
    }

    void mutateInserting()
    {
        this->fillDist();
        // Selecting edge to insert
        vi possibleEdges(m-(n-1));
        int idx = 0;
        for(int i = 0; i < m; ++i)
        {
            if(!hasEdge(i))
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
            for(AdjInfo& ainfo: this->adj[cur])
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
        usedEdges.insert(rdInt);
        this->adj[edge.u].push_back(AdjInfo(edge.v, edge.len, edge.id));
        this->adj[edge.v].push_back(AdjInfo(edge.u, edge.len, edge.id));
        cur = edge.v;
        int e1, e2, rmEdge, curNode;
        double minObj, curObj;
        e1 = e2 = rmEdge = edge.id;         // e2 represents the removed edge
        minObj = curObj = this->objective;
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
                for(AdjInfo& ainfo: this->adj[curNode])
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
                    curObj -= dist[curNode][i]*req[curNode][i];
                    dist[curNode][i] = dist[i][curNode] = newEdge.len + dist[neighbor][i];
                    curObj += dist[curNode][i]*req[curNode][i];
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
                for(AdjInfo& ainfo: this->adj[curNode])
                {
                    if(seen[ainfo.v] || ainfo.id == e1 || ainfo.id == e2)
                        continue;
                    q.push(ainfo.v);
                    seen[ainfo.v] = true;
                    for(int& i : nodesToUpdate)
                    {
                        curObj -= dist[ainfo.v][i]*req[ainfo.v][i];
                        dist[ainfo.v][i] = dist[i][ainfo.v] = ainfo.len + dist[curNode][i];
                        curObj += dist[ainfo.v][i]*req[ainfo.v][i];
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
        this->removeEdge(edge);
        // call this to update the new distances correctly
        this->objective = minObj;
        assert(eq(minObj, this->objective));
    }

    // Mutate when considering to remove a random edge - O(m*n^2)
    void mutateRemoving()
    {
        // selecting edge to remove
        vi possibleEdges(n-1);
        int idx = 0;
        for(auto it = usedEdges.begin(); it != usedEdges.end(); ++it)
        {
            possibleEdges[idx++] = *it;
        }
        assert(idx == n-1);
        int rdInt = possibleEdges[rand()%((int) possibleEdges.size())];
        Edge edge = edges[rdInt];
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
            for(AdjInfo& ainfo: this->adj[curNode])
            {
                if(ainfo.id == edge.id)
                    continue;
                if(inA[ainfo.v])
                    continue;
                inA[ainfo.v] = 1;
                q.push(ainfo.v);
            }
        }
        // remove chosen edge
        this->removeEdge(edge);
        // pre compute the sums
        vector<double> sum(n, 0.0);
        vector<double> dst(n, DBL_MAX);
        for(int i = 0; i < n; ++i)
        {
            fill(dst.begin(), dst.end(), DBL_MAX);
            q.push(i);
            dst[i] = 0.0;
            while(q.size())
            {
                curNode = q.front();
                q.pop();
                for(AdjInfo& ainfo : adj[curNode])
                {
                    if(eq(dst[ainfo.v], DBL_MAX))
                    {
                        dst[ainfo.v] = dst[curNode] + ainfo.len;
                        sum[i] += dst[ainfo.v];
                        q.push(ainfo.v);
                    }
                }
            }
        }
        // Find the best edge to add when removing the chosen edge
        int addEdge;
        addEdge = rdInt;
        int a, b;
        a = edges[rdInt].u;
        b = edges[rdInt].v;
        if(inA[b]) 
            swap(a, b);
        double best, tmp;
        best = sum[a]*(n-szA) + sum[b]*szA + edges[rdInt].len*szA*(n-szA);
        for(int i = 0; i < m; ++i)
        {
            // XOR is used to ensure the edge is connecting the two disconnected sets A and B
            if(!hasEdge(i) && (inA[edges[i].u]^inA[edges[i].v]))
            {
                a = edges[i].u;
                b = edges[i].v;
                if(inA[b]) 
                    swap(a, b);
                tmp = sum[a]*(n-szA) + sum[b]*szA + edges[rdInt].len*szA*(n-szA);
                if(lt(tmp, best))
                {
                    best = tmp;
                    addEdge = i;
                }
            }            
        }
        // Insert the best edge in solution
        Edge bestEdge = edges[addEdge];
        usedEdges.insert(addEdge);
        this->adj[bestEdge.u].push_back(AdjInfo(bestEdge.v, bestEdge.len, bestEdge.id));
        this->adj[bestEdge.v].push_back(AdjInfo(bestEdge.u, bestEdge.len, bestEdge.id));
        // update solution objective function and distances
        this->computeObjectiveFun();
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

void neighbor(int i, Solution& sol)
{
    vi possibleEdges(n-1);
    int idx = 0;
    for(auto it = sol.usedEdges.begin(); it != sol.usedEdges.end(); ++it)
    {
        possibleEdges[idx++] = *it;
    }
    assert(idx == n-1);
    int rdInt = possibleEdges[i];
    Edge edge = edges[rdInt];
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
        for(AdjInfo& ainfo: sol.adj[curNode])
        {
            if(ainfo.id == edge.id)
                continue;
            if(inA[ainfo.v])
                continue;
            inA[ainfo.v] = 1;
            q.push(ainfo.v);
        }
    }
    // remove chosen edge
    sol.removeEdge(edge);
    // pre compute the sums
    vector<double> sum(n, 0.0);
    vector<double> dst(n, DBL_MAX);
    for(int i = 0; i < n; ++i)
    {
        fill(dst.begin(), dst.end(), DBL_MAX);
        q.push(i);
        dst[i] = 0.0;
        while(q.size())
        {
            curNode = q.front();
            q.pop();
            for(AdjInfo& ainfo : sol.adj[curNode])
            {
                if(eq(dst[ainfo.v], DBL_MAX))
                {
                    dst[ainfo.v] = dst[curNode] + ainfo.len;
                    sum[i] += dst[ainfo.v];
                    q.push(ainfo.v);
                }
            }
        }
    }
    // Find the best edge to add when removing the chosen edge
    int addEdge;
    addEdge = rdInt;
    int a, b;
    a = edges[rdInt].u;
    b = edges[rdInt].v;
    if(inA[b]) 
        swap(a, b);
    double best, tmp;
    best = sum[a]*(n-szA) + sum[b]*szA + edges[rdInt].len*szA*(n-szA);
    for(int i = 0; i < m; ++i)
    {
        // XOR is used to ensure the edge is connecting the two disconnected sets A and B
        if(!sol.hasEdge(i) && (inA[edges[i].u]^inA[edges[i].v]))
        {
            a = edges[i].u;
            b = edges[i].v;
            if(inA[b]) 
                swap(a, b);
            tmp = sum[a]*(n-szA) + sum[b]*szA + edges[rdInt].len*szA*(n-szA);
            if(lt(tmp, best))
            {
                best = tmp;
                addEdge = i;
            }
        }            
    }
    // Insert the best edge in solution
    Edge bestEdge = edges[addEdge];
    sol.usedEdges.insert(addEdge);
    sol.adj[bestEdge.u].push_back(AdjInfo(bestEdge.v, bestEdge.len, bestEdge.id));
    sol.adj[bestEdge.v].push_back(AdjInfo(bestEdge.u, bestEdge.len, bestEdge.id));
    // update solution objective function and distances
    sol.computeObjectiveFun();
}

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

void buildRandomSolution(vector<Edge>& edge, Solution& sol)
{
    int m = (int) edge.size();
    UnionFind uf(n);
    vector<int> adj[n];
    vector<int> idxEdge(n, 0);
    for(int i = 0; i < m; ++i)
    {
        adj[edge[i].u].push_back(i);
        adj[edge[i].v].push_back(i);
    }
    // randomly sort the edges
    set<int> nodes;
    set<int>::iterator it;
    for(int i = 0; i < n; ++i)
    {
        shuffle(begin(adj[i]), end(adj[i]), default_random_engine(seed));
        nodes.insert(i);
    }
    int cur = 0;
    int nEdges = 0;
    Edge e;
    
    while(nEdges < n-1)
    {
        it = nodes.lower_bound(cur);
        if(it == nodes.end())
            it = nodes.begin();
        cur = *it;
        if(idxEdge[cur] < (int) adj[cur].size())    // has edge to test
        {
            e = edge[adj[cur][idxEdge[cur]]];
            if(!uf.isSameSet(e.u, e.v))
            {
                uf.unionSet(e.u, e.v);
                nEdges++;
                sol.usedEdges.insert(e.id);
                sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            }
            idxEdge[cur]++;
        }
        else
        {
            nodes.erase(cur);
        }
        cur++;
    }
}

void buildProbGreedySolution(vector<Edge>& edge, vb& fixedEdge, Solution& sol)
{
    int m = (int) edge.size();
    UnionFind uf(n);
    // put fixed edges in solution
    vector<int> idx;
    double fitSum, minVal, maxVal, rngDbl, accVal;
    minVal = DBL_MAX;
    maxVal = 0;
    int solSize = 0;
    for(int i = 0; i < m; ++i)
    {
        if(fixedEdge[i])
        {
            uf.unionSet(edge[i].u, edge[i].v);
            sol.usedEdges.insert(edge[i].id);
            sol.adj[edge[i].u].push_back(AdjInfo(edge[i].v, edge[i].len, edge[i].id));
            sol.adj[edge[i].v].push_back(AdjInfo(edge[i].u, edge[i].len, edge[i].id));
            solSize++;
        }
        else
        {
            idx.push_back(i);
            minVal = min(minVal, edge[i].len);
            maxVal = max(maxVal, edge[i].len);
        }
    }
    // create fitness for each edge
    vector<double> fitness((int) idx.size());
    vb unavEdge((int) idx.size(), false);
    fitSum = 0.0;
    for(int i = 0; i < (int) idx.size(); ++i)
    {
        if(eq(minVal, maxVal))
            fitness[i] = 1.0;
        else
            fitness[i] = 1.0 - ((edge[idx[i]].len - minVal)/(maxVal - minVal)) + 0.4;
        fitSum += fitness[i];
        assert(leq(fitness[i], 1.4));
    }
    int chosen;
    Edge e;
    while(solSize < n-1)    // while solution isn't a tree
    {
        std::uniform_real_distribution<double> distrib(0.0, fitSum);
        rngDbl = distrib(rng);
        assert(leq(rngDbl, fitSum));
        accVal = 0.0;
        chosen = -1;
        for(int i = 0; i < (int) idx.size(); ++i)
        {
            if(unavEdge[i])  // edge unavailable
                continue;
            if(leq(rngDbl, accVal + fitness[i]))    // solution chosen
            {
                chosen = i;
                break;
            }
            accVal += fitness[i];
        }
        if(chosen == -1)    // handling possible no edge selection
        {
            for(int i = (int) idx.size(); i >= 0; ++i)
            {
                if(!unavEdge[i])
                {
                    chosen = i;
                    break;
                }
            }
        }
        fitSum -= fitness[chosen];
        e = edge[idx[chosen]];
        unavEdge[chosen] = true;
        // try inserting edge in solution
        if(!uf.isSameSet(e.u, e.v))
        {
            uf.unionSet(e.u, e.v);
            sol.usedEdges.insert(e.id);
            sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            solSize++;
        }
    }
}

void buildMinPathSolution(vector<Edge>& edge, Solution& sol)
{
    // generate adjacency list to perform Dijkstra
    vector<AdjInfo> adj[n];
    for(Edge& e : edge)
    {
        adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
        adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
    }
    // perform Dijkstra in the random node (cur)
    int cur = rand()%n;
    vector<double> dist(n, DBL_MAX);
    vector<int> uEdge(n, -1);
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
        for(AdjInfo& ainfo : adj[cur])
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
    Edge e;
    int cnt = 0;
    for(int& edgeID : uEdge)
    {
        if(edgeID > -1)
        {
            sol.usedEdges.insert(edgeID);
            cnt++;
            e = edges[edgeID];
            sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
        }
    }
    assert(cnt == n-1);
}

void buildPTASSolution(vector<Edge>& edge, Solution& sol)
{
    // generate adjacency list
    vector<AdjInfo> adj[n];
    for(Edge& e : edge)
    {
        adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
        adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
    }
    vector<int> reservoir(5);
    vector<bool> selected(n);
    vector<int> tmp(n);
    int K, i, cur, cnt;
    while(true) // break means found a valid solution
    {
        fill(selected.begin(), selected.end(), false);
        K = min(n-1, 2 + rand()%4);
        reservoir[0] = rand()%n;
        map<int, int> nToK;
        selected[reservoir[0]] = true;
        for(i = 1; i < K; ++i)
        {
            cnt = 0;
            for(AdjInfo& ainfo : adjList[reservoir[i-1]])
            {
                if(selected[ainfo.v])   // node already selected
                    continue;
                tmp[cnt++] = ainfo.v;
            }
            if(cnt == 0)
                break;
            reservoir[i] = tmp[rand()%cnt];
            selected[reservoir[i]] = true;   
        }
        if(i < K)     // invalid set
            continue;
        for(i = 0; i < K; ++i)
        {
            nToK[reservoir[i]] = i;
        }
        vector<double> dist(n, DBL_MAX);
        vector<int> uEdge(n, -1);
        vector<int> superNode(n, -1);
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        // multi-source Dijsktra
        for(i = 0; i < K; ++i)
        {
            dist[reservoir[i]] = 0.0;
            pq.push(mp(0.0, reservoir[i]));
            superNode[reservoir[i]] = i;
        }
        double w;
        while(pq.size())
        {
            cur = pq.top().second;
            w = pq.top().first;
            pq.pop();
            if(lt(dist[cur], w))
                continue;
            for(AdjInfo& ainfo : adj[cur])
            {
                if(selected[ainfo.v])
                    continue;
                if(dist[ainfo.v] > dist[cur] + ainfo.len)
                {
                    dist[ainfo.v] = dist[cur] + ainfo.len;
                    uEdge[ainfo.v] = ainfo.id;
                    superNode[ainfo.v] = superNode[cur];
                    pq.push(mp(dist[ainfo.v], ainfo.v));
                }
            }
        }
        // construct Solution for minimum path tree from node
        Edge e;
        for(int& edgeID : uEdge)
        {
            if(edgeID > -1)
            {
                sol.usedEdges.insert(edgeID);
                e = edges[edgeID];
                sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            }
        }
        if(K == 2)
        {
            auto it = edgeMap.find(mp(min(reservoir[0], reservoir[1]), max(reservoir[0], reservoir[1])));
            if(it == edgeMap.end())
            {
                // invalid tree
                sol.clear();
                continue;
            }
            e = *(it->second);
            sol.usedEdges.insert(e.id);
            sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            break;
        }
        vector<vector<double>> sumReq(K, vector<double>(K, 0.0));
        for(i = 0; i < n; ++i)
        {
            for(int j = i+1; j < n; ++j)
            {
                if(i == j || superNode[i] == superNode[j]) 
                    continue;
                sumReq[superNode[i]][superNode[j]] += req[i][j];
                sumReq[superNode[j]][superNode[i]] += req[i][j];
            }
        }
        
        list<vector<int>> trees = prufferCodes[K];
        vector<int> deg(K);
        double bestValue = DBL_MAX;
        double tmpValue;
        vector<vector<AdjInfo>> bestAdj;
        for(vector<int>& pruffer : trees)
        {
            fill(deg.begin(), deg.end(), 1);
            vector<vector<AdjInfo>> adjTmp(K, vector<AdjInfo>());
            for(int& node : pruffer)
            {
                deg[node]++;
            }
            set<int> s;
            for(i = 0; i < K; ++i)
            {
                if(deg[i] == 1)
                    s.insert(i);
            }
            int lowest, node;
            bool valid = true;
            for(i = 0; i < K-2; ++i)
            {
                assert((int) s.size());
                lowest = reservoir[*s.begin()];
                s.erase(s.begin());
                node = reservoir[pruffer[i]];
                auto it = edgeMap.find(mp(min(lowest, node), max(lowest, node)));
                if(it == edgeMap.end())
                {
                    // invalid tree, skip
                    valid = false;
                    break;
                }
                e = *(it->second);
                adjTmp[nToK[e.u]].push_back(AdjInfo(nToK[e.v], e.len, e.id));
                adjTmp[nToK[e.v]].push_back(AdjInfo(nToK[e.u], e.len, e.id));
                deg[pruffer[i]]--;
                if(deg[pruffer[i]] == 1)
                    s.insert(pruffer[i]);
            }
            if(valid)
            {
                lowest = reservoir[*s.begin()];
                s.erase(s.begin());
                node = reservoir[*s.begin()];
                auto it = edgeMap.find(mp(min(lowest, node), max(lowest, node)));
                if(it == edgeMap.end())
                {
                    // invalid tree, skip
                    valid = false;
                }
                if(valid)
                {
                    // test this tree
                    e = *(it->second);
                    adjTmp[nToK[e.u]].push_back(AdjInfo(nToK[e.v], e.len, e.id));
                    adjTmp[nToK[e.v]].push_back(AdjInfo(nToK[e.u], e.len, e.id));
                    tmpValue = 0.0;
                    vector<vector<double>> _dist(K, vector<double>(K, DBL_MAX));
                    for(i = 0; i < K; ++i)
                    {
                        _dist[i][i] = 0.0;
                        pq.push(mp(_dist[i][i], i));
                        double w;
                        while(pq.size())
                        {
                            cur = pq.top().second;
                            w = pq.top().first;
                            pq.pop();
                            if(lt(_dist[i][cur], w))
                                continue;
                            for(AdjInfo& ainfo : adjTmp[cur])
                            {
                                assert(selected[reservoir[ainfo.v]]);
                                if(_dist[i][ainfo.v] > _dist[i][cur] + ainfo.len)
                                {
                                    _dist[i][ainfo.v] = _dist[i][cur] + ainfo.len;
                                    pq.push(mp(_dist[i][ainfo.v], ainfo.v));
                                }
                            }
                        }
                        for(int j = 1; j < K; ++j)
                        {
                            tmpValue += _dist[i][j]*sumReq[i][j];
                        }
                    }
                    if(lt(tmpValue, bestValue))
                    {
                        bestValue = tmpValue;
                        bestAdj = adjTmp;
                    }
                }
            }
        }

        if((int) bestAdj.size() == K-1)
        {
            for(int i = 0; i < K; ++i)
            {
                for(AdjInfo& ainfo : bestAdj[i])
                {
                    e = edges[ainfo.id];
                    if(!sol.hasEdge(e.id))
                    {
                        sol.usedEdges.insert(e.id);
                        sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                        sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    }
                }
            }
            break;
        }
        else
        {
            sol.clear();
        }
    }
}

// evolutionary/genetic algorithm
struct Evolutionary
{
    vector<Solution> solutions;
    vector<Solution> offspring;
    vector<ii> wins;
    vector<double> fitness;
    vector<double> newFitness;
    int popSize;                // initial population size
    int numTour;                // number of tournaments per generation
    int numGen;                 // number of generations
    int numCrossover;           // number of crossovers
    int offspringSize;  
    int numMutations;
    Evolutionary(int popSize, int numGen, int numCrossover, int numMutations)
    {
        this->popSize = popSize;
        this->numGen = numGen;
        this->numCrossover = numCrossover;
        this->offspringSize = numCrossover+popSize;
        this->numTour = 10*offspringSize;
        this->numMutations = numMutations;
        solutions.resize(popSize);
        fitness.resize(popSize);
        newFitness.resize(popSize);
        offspring.resize(offspringSize);
        wins.resize(offspringSize);
    }

    Solution run()
    {
        if(mode == 0)
            genRandomPop();
        else if(mode == 1)
            genGreedyProbPop();
        else if(mode == 2)
            genMinPathPop();
        else if(mode == 3)
        {
            genPTASPop();
        }
        for(Solution& sol : solutions)
        {
            assert(lt(0, sol.objective) && lt(sol.objective, DBL_MAX));
        }
        int gen = 1;
        double maxObj, minObj;
        Solution best;
        best.objective = DBL_MAX;
        double fitSum;
        double rngDbl;
        double accVal;
        int rngInt;
        int notImproving = 0;
        double curBestVal = DBL_MAX;
        Solution* tmpBest;
        while(gen <= numGen && notImproving < 25)
        {
            printf("Generation = %d\n", gen);
            minObj = DBL_MAX;
            maxObj = 0;
            // find best solution
            tmpBest = &best;
            for(Solution& sol : solutions)
            {
                minObj = min(minObj, sol.objective);
                maxObj = max(maxObj, sol.objective);
                if(sol.objective < tmpBest->objective)      // update if solution is better
                {
                    tmpBest = &sol;
                }
            }
            best = *tmpBest;
            // Evaluate fitness ([0, 1] interval, greater is better)
            fitSum = 0;
            for(int i = 0; i < popSize; ++i)
            {
                if(abs(minObj - maxObj) < EPS)
                    fitness[i] = 1.0;
                else
                    fitness[i] = 1.0 - (solutions[i].objective - minObj)/(maxObj - minObj) + 0.2;
                fitSum += fitness[i];
                assert(leq(fitness[i], 1.2));
            }
            std::uniform_real_distribution<double> distrib(0.0, fitSum);
            // Crossover between parents
            int id1, id2;
            for(int i = 0; i < numCrossover; ++i)
            {
                id1 = id2 = popSize-1;
                rngDbl = distrib(rng);
                accVal = 0.0;
                for(int j = 0; j < popSize; ++j)
                {
                    if(leq(rngDbl, accVal + fitness[j]))    // parent chosen
                    {
                        id1 = j;
                        break;
                    }
                    accVal += fitness[j];
                }
                rngDbl = distrib(rng);
                accVal = 0.0;
                for(int j = 0; j < popSize; ++j)
                {
                    if(leq(rngDbl, accVal + fitness[j]))    // parent chosen
                    {
                        id2 = j;
                        break;
                    }
                    accVal += fitness[j];
                }
                offspring[i] = crossover(solutions[id1], solutions[id2]);
            }
            int idx = numCrossover;
            for(int i = 0; i < popSize; ++i)
            {
                offspring[idx++] = solutions[i];
            }
            Solution* solPtr;
            tmpBest = &best;
            for(int i = 0; i < numMutations; ++i)
            {
                solPtr = &offspring[rand()%offspringSize];
                rngInt = rand()%2;
                if(rngInt)
                {
                    solPtr->mutateInserting();
                }
                else
                {
                    solPtr->mutateRemoving();
                }
                if(solPtr->objective < tmpBest->objective)      // update if solution is better
                {
                    tmpBest = solPtr;
                }
            }
            best = *tmpBest;
            for(int i = 0; i < offspringSize; ++i)
            {
                wins[i] = mp(0, i);
            }
            tournamentSelection(offspring, wins);
            sort(wins.begin(), wins.end(), greater<ii>());
            for(int i = 0; i < min(offspringSize, popSize); ++i)
            {
                solutions[i] = offspring[wins[i].second];
            }
            if(lt(best.objective, curBestVal))
            {
                curBestVal = best.objective;
                notImproving = 0;
            }
            else
            {
                notImproving++;
            }
            gen++;
            printf("Best so far = %.10f\n", best.objective);
        }
        return best;
    }

    Solution crossover(Solution& s1, Solution& s2)
    {
        vector<Edge> avEdges;
        vb fixedEdge;
        map<int, int> mm;
        for(auto it = s1.usedEdges.begin(); it != s1.usedEdges.end(); ++it)
        {
            mm[*it]++;
        }
        for(auto it = s2.usedEdges.begin(); it != s2.usedEdges.end(); ++it)
        {
            mm[*it]++;
        }
        for(auto it = mm.begin(); it != mm.end(); ++it)
        {
            avEdges.push_back(edges[it->first]);
            if(it->second == 2)
            {
                fixedEdge.push_back(true);
            }
            else
            {
                fixedEdge.push_back(false);
            }
        }
        Solution sol;
        if((int) avEdges.size() == n-1)
        {
            sol = s1;
        }
        else
        {
            // Calling solver
            if(mode == 0)
            {
                // Calling random crossover
                buildRandomSolution(avEdges, sol);
                sol.computeObjectiveFun();
            }
            else if(mode == 1)
            {
                // Calling a greedy crossover selecting lowest cost edges with higher probability
                buildProbGreedySolution(avEdges, fixedEdge, sol);
                sol.computeObjectiveFun();
            }
            else if(mode == 2)
            {
                // Calling a greedy shortest path tree from a random node
                buildMinPathSolution(avEdges, sol);
                sol.computeObjectiveFun();
            }
            else
            {
                // Calling PTAS crossover
                buildPTASSolution(avEdges, sol);
                sol.computeObjectiveFun();
            }
        }
        return sol;
    }

    void tournamentSelection(vector<Solution>& offspring, vector<ii>& wins)
    {
        int i, rdInt;
        int best;
        int K = 2 + rand()%(popSize/10);;
        double p, curP, rdDbl;
        p = 0.9;
        vector<int> reservoir(K);
        std::uniform_real_distribution<double> rd01(0.0, 1.0);
        for(int tour = 0; tour < this->numTour; tour++)
        {
            // Reservoir Algorithm to sample K random solutions
            for(i = 0; i < K; ++i)
                reservoir[i] = i;
            for(; i < (int) wins.size(); ++i)
            {
                rdInt = rand()%(i+1);
                if(rdInt < K)
                {
                    reservoir[rdInt] = i;
                }
            }
            priority_queue<pair<double, int>> pq;
            for(i = 0; i < K; ++i)
            {
                pq.push(mp(-offspring[reservoir[i]].objective, i));
            }
            best = reservoir[K-1];
            curP = 0.8;
            while(pq.size())
            {
                rdDbl = rd01(rng);
                assert(leq(rdDbl, 1.0));
                if(leq(rdDbl, curP))
                {
                    best = reservoir[pq.top().second];
                    break;
                }
                pq.pop();
                curP *= p;
            }
            wins[best].first++;
        }
    }

    /* Generate popSize initial solutions (trees) by shuffling the edges
	   and inserting the edges like Kruskal Algorithm */
    void genRandomPop()
    {
        printf("RandomPop\n");
        vector<Edge> cpy = edges;
        int numForests;
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
                    //sol.usedEdge[e.id] = true;
                    sol.usedEdges.insert(e.id);
                    numForests--;
                }
                if(numForests == 1) // If the tree is done
                {
                    break;
                }
            }
            assert(numForests == 1);
            sol.computeObjectiveFun();
            solutions[i] = sol;
        }
    }

    /* Generate popSize initial solutions with minimal path trees of
       random vertices using Dijkstra */
    void genMinPathPop()
    {
        printf("MinPathPop\n");
        // generate adjacency list to perform Dijkstra
        vector<AdjInfo> adj[n];
        for(Edge& e : edges)
        {
            adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
        }
        int i, rdInt, cur;
        int reservoirSz = min(popSize, n);
        vector<int> reservoir(popSize);
        // Reservoir Algorithm to sample reservoirSz random solutions
        for(i = 0; i < reservoirSz; ++i)
            reservoir[i] = i;
        for(; i < n; ++i)
        {
            rdInt = rand()%(i+1);
            if(rdInt < reservoirSz)
            {
                reservoir[rdInt] = i;
            }
        }
        // Generate Min Path Tree solution
        for(i = 0; i < reservoirSz; ++i)
        {
            // perform Dijkstra in the node (reservoir[i])
            vector<double> dist(n, DBL_MAX);
            vector<int> uEdge(n, -1);
            cur = reservoir[i];
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
                for(AdjInfo& ainfo : adj[cur])
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
            solutions[i] = sol;
        }

        // Generate Min Path Tree solution with some edges removed with some probability
        for(i = reservoirSz; i < popSize; ++i)
        {
            // invalid edges
            vector<bool> tabu(m, false);
            cur = rand()%n;
            for(int j = 0; j < m; ++j)
            {
                if(solutions[cur].hasEdge(j))
                {
                    if((rand()%10) == 0)    // 10% of chance to remove an used edge from the sol
                        tabu[j] = true;
                }
            }
            // perform Dijkstra in the node (cur)
            vector<double> dist(n, DBL_MAX);
            vector<int> uEdge(n, -1);
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
                for(AdjInfo& ainfo : adj[cur])
                {
                    if(tabu[ainfo.id])
                        continue;
                    if(dist[ainfo.v] > dist[cur] + ainfo.len)
                    {
                        dist[ainfo.v] = dist[cur] + ainfo.len;
                        uEdge[ainfo.v] = ainfo.id;
                        pq.push(mp(dist[ainfo.v], ainfo.v));
                    }
                }
            }
            bool connected = true;
            for(int j = 0; j < n; ++j)
            {
                if(lt(dist[j], DBL_MAX))       // node reacheable
                    continue;
                connected = false;
                break;
            }
            if(!connected)                     // try again!
            {
                i--;
                continue;
            }
            // construct Solution for minimum path tree from node
            Solution sol;
            Edge e;
            int cnt = 0;
            for(int& edgeID : uEdge)
            {
                if(edgeID > -1)
                {
                    sol.usedEdges.insert(edgeID);
                    e = edges[edgeID];
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    cnt++;
                }
            }
            assert(cnt == n-1);
            sol.computeObjectiveFun();
            solutions[i] = sol;
        }  
    }

    /* Generate popSize initial solutions with highest probability of
       of chosing the lowest cost edges (MST-like) */
    void genGreedyProbPop()
    {
        printf("MST-like pop\n");

        for(int i = 0; i < popSize; ++i)
        {
            UnionFind uf(n);
            double fitSum, minVal, maxVal, rngDbl, accVal;
            minVal = DBL_MAX;
            maxVal = 0;
            int solSize = 0;
            for(int j = 0; j < m; ++j)
            {
                minVal = min(minVal, edges[j].len);
                maxVal = max(maxVal, edges[j].len);
            }
            // create fitness for each edge
            vector<double> fitness(m);
            vb unavEdge(m, false);
            fitSum = 0.0;
            for(int j = 0; j < m; ++j)
            {
                if(eq(minVal, maxVal))
                    fitness[j] = 1.0;
                else
                    fitness[j] = 1.0 - ((edges[j].len - minVal)/(maxVal - minVal)) + 0.4;
                fitSum += fitness[j];
                assert(leq(fitness[j], 1.4));
            }
            int chosen;
            Edge e;
            Solution sol;
            while(solSize < n-1)    // while solution isn't a tree
            {
                std::uniform_real_distribution<double> distrib(0.0, fitSum);
                rngDbl = distrib(rng);
                assert(leq(rngDbl, fitSum));
                accVal = 0.0;
                chosen = m-1;
                for(int j = 0; j < m; ++j)
                {
                    if(unavEdge[j])  // edge unavailable
                        continue;
                    if(leq(rngDbl, accVal + fitness[j]))    // solution chosen
                    {
                        chosen = j;
                        break;
                    }
                    accVal += fitness[j];
                }
                fitSum -= fitness[chosen];
                e = edges[chosen];
                unavEdge[chosen] = true;
                // try inserting edge in solution
                if(!uf.isSameSet(e.u, e.v))
                {
                    uf.unionSet(e.u, e.v);
                    sol.usedEdges.insert(e.id);
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    solSize++;
                }
            }
            sol.computeObjectiveFun();
            solutions[i] = sol;;
        }
    }

    void genPTASPop()
    {
        printf("PTAS pop\n");

        adjList.assign(n, vector<AdjInfo>());
        for(Edge& e : edges)
        {
            adjList[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            adjList[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            edgeMap[mp(min(e.u, e.v), max(e.u, e.v))] = &e;
        }
        int popIdx = 0;
        vector<int> reservoir(5);
        vector<bool> selected(n);
        vector<int> tmp(n);
        int K, i, cur, cnt;
        while(popIdx < popSize)
        {
            fill(selected.begin(), selected.end(), false);
            K = min(n-1, 2 + rand()%4);
            reservoir[0] = rand()%n;
            map<int, int> nToK;
            selected[reservoir[0]] = true;
            for(i = 1; i < K; ++i)
            {
                cnt = 0;
                for(AdjInfo& ainfo : adjList[reservoir[i-1]])
                {
                    if(selected[ainfo.v])   // node already selected
                        continue;
                    tmp[cnt++] = ainfo.v;
                }
                if(cnt == 0)
                    break;
                reservoir[i] = tmp[rand()%cnt];
                selected[reservoir[i]] = true;   
            }
            if(i < K)     // invalid set
                continue;
            for(i = 0; i < K; ++i)
            {
                nToK[reservoir[i]] = i;
            }
            vector<double> dist(n, DBL_MAX);
            vector<int> uEdge(n, -1);
            vector<int> superNode(n, -1);
            priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
            // multi-source Dijsktra
            for(i = 0; i < K; ++i)
            {
                dist[reservoir[i]] = 0.0;
                pq.push(mp(0.0, reservoir[i]));
                superNode[reservoir[i]] = i;
            }
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
                    if(selected[ainfo.v])
                        continue;
                    if(dist[ainfo.v] > dist[cur] + ainfo.len)
                    {
                        dist[ainfo.v] = dist[cur] + ainfo.len;
                        uEdge[ainfo.v] = ainfo.id;
                        superNode[ainfo.v] = superNode[cur];
                        pq.push(mp(dist[ainfo.v], ainfo.v));
                    }
                }
            }
            // construct Solution for minimum path tree from node
            Solution sol;
            Edge e;
            for(int& edgeID : uEdge)
            {
                if(edgeID > -1)
                {
                    sol.usedEdges.insert(edgeID);
                    e = edges[edgeID];
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                }
            }
            if(K == 2)
            {
                auto it = edgeMap.find(mp(min(reservoir[0], reservoir[1]), max(reservoir[0], reservoir[1])));
                if(it == edgeMap.end())
                {
                    // invalid tree
                    continue;
                }
                e = *(it->second);
                sol.usedEdges.insert(e.id);
                sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                sol.computeObjectiveFun();
                solutions[popIdx++] = sol;
                continue;
            }
            vector<vector<double>> sumReq(K, vector<double>(K, 0.0));
            for(i = 0; i < n; ++i)
            {
                for(int j = i+1; j < n; ++j)
                {
                    if(i == j || superNode[i] == superNode[j]) 
                        continue;
                    sumReq[superNode[i]][superNode[j]] += req[i][j];
                    sumReq[superNode[j]][superNode[i]] += req[i][j];
                }
            }
            
            list<vector<int>> trees = prufferCodes[K];
            vector<int> deg(K);
            double bestValue = DBL_MAX;
            double tmpValue;
            vector<vector<AdjInfo>> bestAdj;
            for(vector<int>& pruffer : trees)
            {
                fill(deg.begin(), deg.end(), 1);
                vector<vector<AdjInfo>> adjTmp(K, vector<AdjInfo>());
                for(int& node : pruffer)
                {
                    deg[node]++;
                }
                set<int> s;
                for(i = 0; i < K; ++i)
                {
                    if(deg[i] == 1)
                        s.insert(i);
                }
                int lowest, node;
                bool valid = true;
                for(i = 0; i < K-2; ++i)
                {
                    assert((int) s.size());
                    lowest = reservoir[*s.begin()];
                    s.erase(s.begin());
                    node = reservoir[pruffer[i]];
                    auto it = edgeMap.find(mp(min(lowest, node), max(lowest, node)));
                    if(it == edgeMap.end())
                    {
                        // invalid tree, skip
                        valid = false;
                        break;
                    }
                    e = *(it->second);
                    adjTmp[nToK[e.u]].push_back(AdjInfo(nToK[e.v], e.len, e.id));
                    adjTmp[nToK[e.v]].push_back(AdjInfo(nToK[e.u], e.len, e.id));
                    deg[pruffer[i]]--;
                    if(deg[pruffer[i]] == 1)
                        s.insert(pruffer[i]);
                }
                if(valid)
                {
                    lowest = reservoir[*s.begin()];
                    s.erase(s.begin());
                    node = reservoir[*s.begin()];
                    auto it = edgeMap.find(mp(min(lowest, node), max(lowest, node)));
                    if(it == edgeMap.end())
                    {
                        // invalid tree, skip
                        valid = false;
                    }
                    if(valid)
                    {
                        // test this tree
                        e = *(it->second);
                        adjTmp[nToK[e.u]].push_back(AdjInfo(nToK[e.v], e.len, e.id));
                        adjTmp[nToK[e.v]].push_back(AdjInfo(nToK[e.u], e.len, e.id));
                        tmpValue = 0.0;
                        vector<vector<double>> _dist(K, vector<double>(K, DBL_MAX));
                        for(i = 0; i < K; ++i)
                        {
                            _dist[i][i] = 0.0;
                            pq.push(mp(_dist[i][i], i));
                            double w;
                            while(pq.size())
                            {
                                cur = pq.top().second;
                                w = pq.top().first;
                                pq.pop();
                                if(lt(_dist[i][cur], w))
                                    continue;
                                for(AdjInfo& ainfo : adjTmp[cur])
                                {
                                    assert(selected[reservoir[ainfo.v]]);
                                    if(_dist[i][ainfo.v] > _dist[i][cur] + ainfo.len)
                                    {
                                        _dist[i][ainfo.v] = _dist[i][cur] + ainfo.len;
                                        pq.push(mp(_dist[i][ainfo.v], ainfo.v));
                                    }
                                }
                            }
                            for(int j = 1; j < K; ++j)
                            {
                                tmpValue += _dist[i][j]*sumReq[i][j];
                            }
                        }
                        if(lt(tmpValue, bestValue))
                        {
                            bestValue = tmpValue;
                            bestAdj = adjTmp;
                        }
                    }
                }
            }

            if((int) bestAdj.size() == K-1)
            {
                for(int i = 0; i < K; ++i)
                {
                    for(AdjInfo& ainfo : bestAdj[i])
                    {
                        e = edges[ainfo.id];
                        if(!sol.hasEdge(e.id))
                        {
                            sol.usedEdges.insert(e.id);
                            sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                            sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                        }
                    }
                }
                sol.computeObjectiveFun();
                assert(lt(0, sol.objective) && lt(sol.objective, DBL_MAX));   // ensure value is correct
                solutions[popIdx++] = sol;
            }
        }
    }

};

int main(int argc, char* argv[])
{
    if(argc != 5)
    {
        printf("usage: ./evolutionary popSize numGen numCrossovers numMutations < inputFile\n");
        return -1;
    }
    cin >> n >> m;
    edges.resize(m);
    for(int i = 0; i < m; ++i)
    {
        cin >> edges[i].u >> edges[i].v >> edges[i].len;
        edges[i].id = i;
    }
    subtreeSize.resize(n);
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
    ofstream log("log.txt", ios::app);
    log << fixed << setprecision(10);
    for(mode = 2; mode >= 2; mode--)
    {
        if(mode == 0)
        {
            printf("Random Mode Selected\n");
            log << "RANDOM\n";
        }
        else if(mode == 1)
        {
            printf("MST Mode Selected\n");
            log << "MST\n";
        }
        else if(mode == 2)
        {
            printf("Minimum Path Mode Selected\n");
            log << "Minimum Path\n";
        }
        else
        {
            printf("PTAS Mode Selected\n");
            log << "PTAS\n";
            for(int k = 3; k <= 5; k++)
            {
                list<vector<int>> lst;
                vector<int> tmp(k-2);
                createPruffer(tmp, lst, 0, k);
                prufferCodes[k] = lst;
            }
        }
        for(int seedid = 0; seedid < 10; ++seedid)
        {
            seed = seedVector[seedid];
            printf("seed = %u\n", seed);
            //Initialize seeds
            srand(seed);
            rng.seed(seed);
            Evolutionary ev(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
            chrono::steady_clock::time_point begin, end;
            begin = chrono::steady_clock::now();
            Solution best = ev.run();
            printf("Best Value Found = %.10f\n", best.objective);
            end = chrono::steady_clock::now();
            cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << endl;
            log << best.objective << endl;
        }
    }
    log.close();
    return 0;
}