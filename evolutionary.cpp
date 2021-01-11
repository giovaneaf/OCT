#include <bits/stdc++.h>
//#include <gurobi_c++.h>

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

// seed used to generate random numbers
unsigned seed;
// seeds used for testing
unsigned seedVector[] = {280192806, 871237442, 2540188929, 107472404, 3957311442, 316851227, 619606212, 1078082709, 916212990, 698598169};
//Mersenne Twister: Good quality random number generator
std::mt19937 rng;
map<ii, Edge*> edgeMap;

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

// IMPORTANT: Assertions should be removed when testing

// Stores candidate solution (tree)
struct Solution
{
    vector<vector<AdjInfo>> adj;
    vector<vector<double>> dist;
    vb usedEdge;
    double objective;
    
    Solution()
    {
        adj.resize(n, vector<AdjInfo>());
        dist.resize(n, vector<double>(n));
        usedEdge.resize(m, false);
        objective = 0;
    }
    
    /* Input: Adjacency list of the tree
	   Output: Objective function value */
    void computeObjectiveFun()
    {
        int cur;
        for(int i = 0; i < n; ++i)
        {
            fill(this->dist[i].begin(), this->dist[i].end(), DBL_MAX);
        }
        this->objective = 0;
        // BFS for each node to compute the distances
        for(int node = 0; node < n; ++node)
        {
            this->dist[node][node] = 0;
            queue<int> q;
            q.push(node);
            while(q.size())
            {
                cur = q.front();
                q.pop();
                for(AdjInfo& ainfo: this->adj[cur])
                {
                    if(this->dist[node][ainfo.v] == DBL_MAX)
                    {
                        this->dist[node][ainfo.v] = ainfo.len + this->dist[node][cur];
                        q.push(ainfo.v);
                    }
                }
            }
            for(int v = node+1; v < n; v++)
            {
                this->objective += this->dist[node][v]*req[node][v];
            }
        }
    }

    // Mutate when inserting a new edge in the solution - O(n^2)
    void mutateInserting()
    {
        // Selecting edge to insert
        vi possibleEdges(m-(n-1));
        int idx = 0;
        for(int i = 0; i < m; ++i)
        {
            if(!this->usedEdge[i])
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
        this->usedEdge[rdInt] = true;
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
                    curObj -= this->dist[curNode][i]*req[curNode][i];
                    this->dist[curNode][i] = this->dist[i][curNode] = newEdge.len + this->dist[neighbor][i];
                    curObj += this->dist[curNode][i]*req[curNode][i];
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
                        curObj -= this->dist[ainfo.v][i]*req[ainfo.v][i];
                        this->dist[ainfo.v][i] = this->dist[i][ainfo.v] = ainfo.len + this->dist[curNode][i];
                        curObj += this->dist[ainfo.v][i]*req[ainfo.v][i];
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
        this->computeObjectiveFun();
        assert(eq(minObj, this->objective));
    }

    // Mutate when considering to remove a random edge - O(m*n^2)
    void mutateRemoving()
    {
        // selecting edge to remove
        vi possibleEdges(n-1);
        int idx = 0;
        for(int i = 0; i < m; ++i)
        {
            if(this->usedEdge[i])
            {
                possibleEdges[idx++] = i;
            }
        }
        int rdInt = possibleEdges[rand()%((int) possibleEdges.size())];
        Edge edge = edges[rdInt];
        // remove chosen edge
        this->removeEdge(edge);
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
                assert(ainfo.id != edge.id); // assert edge was removed
                if(inA[ainfo.v])
                    continue;
                inA[ainfo.v] = 1;
                q.push(ainfo.v);
            }
        }
        // Find the best edge to add when removing the chosen edge
        double minObj, curObj;
        int addEdge;
        minObj = curObj = this->objective;
        addEdge = rdInt;
        for(int i = 0; i < m; ++i)
        {
            // XOR is used to ensure the edge is connecting the two disconnected sets A and B
            if(!this->usedEdge[i] && (inA[edges[i].u]^inA[edges[i].v]))
            {
                curNode = edges[i].u;
                list<int> nodesToUpdate;
                for(int j = 0; j < n; ++j)
                {
                    if(inA[curNode]^inA[j]) // need to be updated
                    {
                        curObj -= this->dist[curNode][j]*req[curNode][j];
                        this->dist[curNode][j] = this->dist[j][curNode] = edges[i].len + this->dist[edges[i].v][j];
                        curObj += this->dist[curNode][j]*req[curNode][j];
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
                    for(AdjInfo& ainfo: this->adj[curNode])
                    {
                        if(seen[ainfo.v])
                            continue;
                        assert((inA[curNode]^inA[ainfo.v]) == 0);
                        q.push(ainfo.v);
                        seen[ainfo.v] = true;
                        for(int& j : nodesToUpdate)
                        {
                            curObj -= this->dist[ainfo.v][j]*req[ainfo.v][j];
                            this->dist[ainfo.v][j] = this->dist[j][ainfo.v] = ainfo.len + this->dist[curNode][j];
                            curObj += this->dist[ainfo.v][j]*req[ainfo.v][j];
                        }
                    }
                }
                if(curObj < minObj)
                {
                    minObj = curObj;
                    addEdge = i;
                }
            }
        }
        // Insert the best edge in solution
        Edge bestEdge = edges[addEdge];
        this->usedEdge[addEdge] = true;
        this->adj[bestEdge.u].push_back(AdjInfo(bestEdge.v, bestEdge.len, bestEdge.id));
        this->adj[bestEdge.v].push_back(AdjInfo(bestEdge.u, bestEdge.len, bestEdge.id));
        // update solution objective function and distances
        this->objective = minObj;
        assert(inA[bestEdge.u]^inA[bestEdge.v]); // assert that the edge form a tree
        curNode = bestEdge.u;
        int neighbor = getNeighbor(curNode, bestEdge);
        for(int i = 0; i < n; ++i)
        {
            if(inA[curNode]^inA[i]) // if the values are updated by the edge
                this->dist[curNode][i] = this->dist[i][curNode] = bestEdge.len + this->dist[neighbor][i];
        }
        vb seen(n, false);
        q.push(curNode);
        while(q.size())
        {
            curNode = q.front();
            seen[curNode] = true;
            q.pop();
            for(AdjInfo& ainfo: this->adj[curNode])
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
                        this->dist[ainfo.v][i] = this->dist[i][ainfo.v] = ainfo.len + this->dist[curNode][i];
                    }
                }
            }
        }
        // assertion to check if distances are correctly updated
        double tmp = 0;
        for(int i = 0; i < n; ++i)
            for(int j = i+1; j < n; ++j)
                tmp += this->dist[i][j]*req[i][j];
        assert(eq(tmp, this->objective));
    }

    // Removes edge from the solution (doesn't recompute anything)
    void removeEdge(Edge& edge)
    {
        this->usedEdge[edge.id] = false;
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
    for(int i = 0; i < m; ++i)
    {
        if(s.usedEdge[i])
        {
            print(edges[i]);
        }
    }
    putchar('\n');
    printf("Objective value = %.2f\n", s.objective);
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


// Variables used for solver
/*static bool setupGurobi = true;
static int constrCnt = 0;
GRBEnv env = GRBEnv(true);
vector<vector<int>> getIdxFlow;
vector<double> O;

string getNewConstr()
{
    return "C" + to_string(constrCnt++);
}*/

/* 
Flow formulation retrieved from:
PhD Thesis - The Optimum Communication Spanning Tree Problem (2015)
Author: Carlos Luna-Mota
Advisor: Elena Fernández
*/
/*
Solution gurobiSolverFlow(vector<Edge>& avEdges, vb& fixedEdge, double timeLimit)
{
    assert((int) avEdges.size() == (int) fixedEdge.size());
    Solution sol;
    vector<vector<int>> Nmin, Nplus;
    vector<Edge> mEdges;
    double d = 0;
    int m = 2*(int) avEdges.size();
    mEdges.resize(m);
    Nmin.resize(n, vector<int>());
    Nplus.resize(n, vector<int>());
    int cnt = 0;
    for(Edge& e : avEdges)
    {
        mEdges[cnt] = {e.u, e.v, e.len, cnt};
        mEdges[cnt+1] = {e.v, e.u, e.len, cnt+1};
        Nmin[e.v].push_back(cnt);
        Nplus[e.u].push_back(cnt);
        Nplus[e.v].push_back(cnt+1);
        Nmin[e.u].push_back(cnt+1);
        d += e.len;
        cnt += 2;
    }
    try 
    {
        if(setupGurobi)
        {
            env.set("OutputFlag", "0");
            env.start();
            O.resize(n, 0.0);
            for(int u = 0; u < n; ++u)
            {
                for(int v = u+1; v < n; ++v)
                {
                    O[u] += req[u][v];
                }
            }
            setupGurobi = false;
        }
        int cnt = 0;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < m; ++v)
            {
                getIdxFlow[u][v] = cnt++;
            }
        }
        env.set("TimeLimit", to_string(timeLimit));
        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create binary variables
        GRBVar* x = model.addVars(m/2, GRB_BINARY);

        // Create continuous variables
        GRBVar* f = model.addVars(n*m, GRB_CONTINUOUS);

        GRBLinExpr obj, linexpr, linexpr2;;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < m; ++v)
            {
                obj.addTerms(&mEdges[v].len, &f[getIdxFlow[u][v]], 1);
                model.addConstr(f[getIdxFlow[u][v]] >= 0, getNewConstr());  //(h)
            }
        }
        model.setObjective(obj, GRB_MINIMIZE);              // (a)
        const double one = 1.0;
        const double negOne = -1.0;
        for(int i = 0; i < m/2; ++i)
        {
            linexpr.addTerms(&one, &x[i], 1);
        }
        model.addConstr(linexpr == n-1, getNewConstr());      // (b)
        linexpr.clear();
        double val;
        for(int o = 0; o < n; ++o)
        {
            for(int& inEdge : Nmin[o])
            {
                linexpr.addTerms(&one, &f[getIdxFlow[o][inEdge]], 1);
            }
            model.addConstr(linexpr == 0, getNewConstr());      // (c)
            linexpr.clear();
            for(int& outEdge : Nplus[o])
            {
                linexpr.addTerms(&one, &f[getIdxFlow[o][outEdge]], 1);
            }
            model.addConstr(linexpr == O[o], getNewConstr());   // (e)
            linexpr.clear();
            for(int j = 0; j < n; ++j)
            {
                if(o == j) continue;
                for(int& inEdge : Nmin[j])
                {
                    linexpr.addTerms(&one, &f[getIdxFlow[o][inEdge]], 1);
                }
                for(int& outEdge : Nplus[j])
                {
                    linexpr.addTerms(&negOne, &f[getIdxFlow[o][outEdge]], 1);
                }
                val = (j <= o ? 0.0 : req[o][j]);
                model.addConstr(linexpr == val, getNewConstr());    // (d)
                linexpr.clear();
            }
            cnt = 0;
            for(int j = 0; j < m/2; ++j)
            {
                model.addConstr(f[getIdxFlow[o][cnt]] + f[getIdxFlow[o][cnt+1]] <= O[o]*x[j]); // (f)
                cnt += 2;
            }
        }
        for(int i = 0; i < m/2; ++i)
        {
            if(fixedEdge[i])
            {
                model.addConstr(x[i] == 1, getNewConstr());
            }
        }
        // Optimize model
        model.optimize();
        // Set new solution to return
        for(int i = 0; i < m/2; i++)
        {
            if(x[i].get(GRB_DoubleAttr_X) > 0.99)
            {
                sol.usedEdge[avEdges[i].id] = true;
                sol.adj[avEdges[i].u].push_back(AdjInfo(avEdges[i].v, avEdges[i].len, avEdges[i].id));
                sol.adj[avEdges[i].v].push_back(AdjInfo(avEdges[i].u, avEdges[i].len, avEdges[i].id));
            }
        }
    }
    catch(GRBException e) 
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } 
    catch(...) 
    {
        cout << "Exception during optimization" << endl;
    }
    return sol;
}
*/
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
    for(int i = 0; i < n; ++i)
    {
        shuffle(begin(adj[i]), end(adj[i]), default_random_engine(seed));
    }
    int cur = 0;
    int nEdges = 0;
    Edge e;
    while(nEdges < n-1)
    {
        if(idxEdge[cur] < (int) adj[cur].size())    // has edge to test
        {
            e = edge[adj[cur][idxEdge[cur]]];
            if(!uf.isSameSet(e.u, e.v))
            {
                uf.unionSet(e.u, e.v);
                nEdges++;
                sol.usedEdge[e.id] = true;
                sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            }
            idxEdge[cur]++;
        }
        cur = ((cur == n-1) ? 0 : cur+1);
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
            sol.usedEdge[edge[i].id] = true;
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
            sol.usedEdge[e.id] = true;
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
    while(pq.size())
    {
        cur = pq.top().second;
        pq.pop();
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
            sol.usedEdge[edgeID] = true;
            cnt++;
            e = edges[edgeID];
            sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
        }
    }
    assert(cnt == n-1);
}

// evolutionary/genetic algorithm
struct Evolutionary
{
    vector<Solution> solutions;
    vector<Solution> offspring;
    vector<ii> wins;
    vector<double> fitness;
    int popSize;                // initial population size
    int numTour;                // number of tournaments per generation
    int numGen;                 // number of generations
    int numCrossover;           // number of crossovers
    int offspringSize;  
    Evolutionary(int popSize, int numGen, int numCrossover)
    {
        this->popSize = popSize;
        this->numGen = numGen;
        this->numCrossover = numCrossover;
        this->offspringSize = numCrossover+popSize;
        this->numTour = 10*offspringSize;
        solutions.resize(popSize);
        fitness.resize(popSize);
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
            sol.computeObjectiveFun();
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
        while(gen <= numGen && notImproving < 25)
        {
            printf("Generation = %d\n", gen);
            minObj = DBL_MAX;
            maxObj = 0;
            // find best solution
            for(Solution& sol : solutions)
            {
                minObj = min(minObj, sol.objective);
                maxObj = max(maxObj, sol.objective);
                if(sol.objective < best.objective)      // update if solution is better
                {
                    best = sol;
                }
            }
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
            int numMutations = numCrossover/2;
            Solution* solPtr;
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
                if(solPtr->objective < best.objective)      // update if solution is better
                {
                    best = *solPtr;
                }
            }
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
        bool equal = true;
        vector<Edge> avEdges;
        vb fixedEdge;
        for(int i = 0; i < m; ++i)
        {
            if((!s1.usedEdge[i]) && (!s2.usedEdge[i]))
                continue;
            // Edge used in at least one tree
            avEdges.push_back(edges[i]);
            if(s1.usedEdge[i] && s2.usedEdge[i])
            {
                fixedEdge.push_back(true);
            }
            else
            {
                equal = false;
                fixedEdge.push_back(false);
            }         
        }
        Solution sol;
        if(equal)
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
                    sol.usedEdge[e.id] = true;
                    numForests--;
                }
                if(numForests == 1) // If the tree is done
                {
                    break;
                }
            }
            assert(numForests == 1);
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
            while(pq.size())
            {
                cur = pq.top().second;
                pq.pop();
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
                    sol.usedEdge[edgeID] = true;
                    e = edges[edgeID];
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    cnt++;
                }
            }
            assert(cnt == n-1);
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
                if(solutions[cur].usedEdge[j])
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
            while(pq.size())
            {
                cur = pq.top().second;
                pq.pop();
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
                    sol.usedEdge[edgeID] = true;
                    e = edges[edgeID];
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    cnt++;
                }
            }
            assert(cnt == n-1);
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
                    sol.usedEdge[e.id] = true;
                    sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                    solSize++;
                }
            }
            solutions[i] = sol;;
        }
    }

    void genPTASPop()
    {
        vector<AdjInfo> adj[n];
        for(Edge& e : edges)
        {
            adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
            edgeMap[mp(min(e.u, e.v), max(e.u, e.v))] = &e;
        }
        int popIdx = 0;
        vector<int> reservoir(5);
        vector<bool> selected(n);
        int K, i, rdInt;
        while(popIdx < popSize)
        {
            fill(selected.begin(), selected.end(), false);
            K = min(n-1, 2 + rand()%4);
            for(i = 0; i < K; ++i)
                reservoir[i] = i;
            for(; i < n; ++i)
            {
                rdInt = rand()%(i+1);
                if(rdInt < K)
                {
                    reservoir[rdInt] = i;
                }
            }
            for(i = 0; i < K; ++i)
            {
                selected[reservoir[i]] = true;
            }
            shuffle(reservoir.begin(), reservoir.begin()+K, default_random_engine(seed));
            int root = -1;
            for(i = 0; i < K; ++i)
            {
                for(AdjInfo& ainfo : adj[reservoir[i]])
                {
                    if(!selected[ainfo.v])
                    {
                        root = i;
                        break;
                    }
                }
            }
            if(root == -1)
                continue;
            vector<double> dist(n, DBL_MAX);
            vector<int> uEdge(n, -1);
            int cur = reservoir[root];
            dist[cur] = 0.0;
            priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
            pq.push(mp(dist[cur], cur));
            while(pq.size())
            {
                cur = pq.top().second;
                pq.pop();
                for(AdjInfo& ainfo : adj[cur])
                {
                    if(selected[ainfo.v])
                        continue;
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
            for(int& edgeID : uEdge)
            {
                if(edgeID > -1)
                {
                    sol.usedEdge[edgeID] = true;
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
                sol.usedEdge[e.id] = true;
                sol.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                sol.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                solutions[popIdx++] = sol;
                continue;
            }
            
            // for K > 2
            /*for(int i = 0; i < K; ++i)
            {
                printf("%d ", reservoir[i]);
            }
            putchar('\n');
            printf("%d\n", root);*/
            vector<int> tmp(K-2);
            list<vector<int>> trees;
            createPruffer(tmp, trees, 0, K);
            vector<int> deg(K);
            Solution best, aux;
            best.objective = DBL_MAX;
            for(vector<int>& pruffer : trees)
            {
                fill(deg.begin(), deg.end(), 1);
                aux = sol;
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
                    aux.usedEdge[e.id] = true;
                    aux.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                    aux.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
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
                        aux.usedEdge[e.id] = true;
                        aux.adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
                        aux.adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
                        aux.computeObjectiveFun();
                        if(aux.objective < best.objective)
                        {
                            best = aux;
                        }
                    }
                }
            }

            if(lt(best.objective, DBL_MAX))
            {
                solutions[popIdx++] = best;
            }

        }
    }

};

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        printf("usage: ./evolutionary popSize numGen numCrossovers < inputFile\n");
        return -1;
    }
    cin >> n >> m;
    edges.resize(m);
    for(int i = 0; i < m; ++i)
    {
        cin >> edges[i].u >> edges[i].v >> edges[i].len;
        edges[i].id = i;
    }
    //getIdxFlow.resize(n, vector<int>(2*m));
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
        for(int seedid = 0; seedid < 10; ++seedid)
        {
            seed = seedVector[seedid];
            printf("seed = %u\n", seed);
            //Initialize seeds
            srand(seed);
            rng.seed(seed);
            Evolutionary ev(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
            chrono::steady_clock::time_point begin, end;
            begin = chrono::steady_clock::now();
            Solution best = ev.run();
            /*vector<Edge> vEdges;
            for(int i = 0; i < m; ++i)
            {
                if(best.usedEdge[i])
                    vEdges.push_back(edges[i]);
            }
            vector<bool> usedEdge(n-1, true);
            Solution validation = gurobiSolverFlow(vEdges, usedEdge, INT_MAX);
            validation.computeObjectiveFun();
            assert(eq(validation.objective, best.objective));*/
            printf("Best Value Found = %.10f\n", best.objective);
            end = chrono::steady_clock::now();
            cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << endl;
            log << best.objective << endl;
        }
    }
    log.close();
    return 0;
}

/*
// Variables used for Gurobi solver
vector<vector<double>> reqMin2;
vector<vector<vector<int>>> getIdx;
static bool setupGurobi = true;
static int constrCnt = 0;
GRBEnv env = GRBEnv(true);

Solution gurobiSolverAnc(vector<Edge>& avEdges, vb& fixedEdge)
{
    assert((int) avEdges.size() == (int) fixedEdge.size());
    Solution sol;
    vector<vector<int>> Nmin;
    vector<Edge> mEdges;
    double d = 0;
    int m = 2*(int) avEdges.size();
    mEdges.resize(m);
    Nmin.resize(n, vector<int>());
    int cnt = 0;
    for(Edge& e : avEdges)
    {
        mEdges[cnt] = {e.u, e.v, e.len, cnt};
        mEdges[cnt+1] = {e.v, e.u, e.len, cnt+1};
        Nmin[e.v].push_back(cnt);
        Nmin[e.u].push_back(cnt+1);
        d += e.len;
        cnt += 2;
    }
    try 
    {       

        if(setupGurobi)
        {
            env.set("OutputFlag", "0");
            env.start();
            // Computed needed values for formulation
            reqMin2.resize(n, vector<double>(n));
            for(int u = 0; u < n; ++u)
            {
                for(int v = u+1; v < n; ++v)
                {
                    reqMin2[u][v] = -2*req[u][v];
                }
            }
            getIdx.resize(n, vector<vector<int>>(n, vector<int>(n)));
            cnt = 0;
            for(int u = 0; u < n; ++u)
                for(int v = 0; v < n; ++v)
                    for(int w = 0; w < n; ++w)
                        getIdx[u][v][w] = cnt++;
            // Don't need to be computed anymore
            setupGurobi = false;
        }

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create binary variables
        GRBVar* x = model.addVars(m, GRB_BINARY);
        GRBVar* y = model.addVars(n*n, GRB_BINARY);
        GRBVar* z = model.addVars(n*n*n, GRB_BINARY);

        // Create integer variables
        GRBVar* eta = model.addVars(n, GRB_INTEGER);

        // Create continuous variables
        GRBVar* delta = model.addVars(n, GRB_CONTINUOUS);
        GRBVar* rho = model.addVars(n*n*n, GRB_CONTINUOUS);

        GRBLinExpr obj;
        for(int u = 0; u < n; ++u)
        {
            for(int v = u+1; v < n; ++v)
            {
                obj.addTerms(&req[u][v], &delta[u], 1);
                obj.addTerms(&req[u][v], &delta[v], 1);
                for(int w = 0; w < n; ++w)
                {
                    obj.addTerms(&reqMin2[u][v], &rho[getIdx[w][u][v]], 1);
                }
            }
        }
        model.setObjective(obj, GRB_MINIMIZE);              // (01)
        GRBLinExpr linexpr, linexpr2;

        int root = 0;
        const double one = 1.0;
        for(int& inEdge : Nmin[root])
        {
            linexpr.addTerms(&one, &x[inEdge], 1);
        }
        model.addConstr(linexpr == 0, getNewConstr());      // (02)
        model.addConstr(delta[root] == 0, getNewConstr());  // (04)
        model.addConstr(eta[root] == 0, getNewConstr());    // (08)
        linexpr.clear();

        for(int u = 1; u < n; ++u)
        {
            for(int& inEdge : Nmin[u])
            {
                linexpr.addTerms(&one, &x[inEdge], 1);
            }
            model.addConstr(linexpr == 1, getNewConstr());  // (03)
            model.addConstr(delta[u] >= 0, getNewConstr()); // (05)
            model.addConstr(eta[u] >= 0, getNewConstr());   // (09)
            linexpr.clear();
        }
        for(Edge& e : mEdges)
        {
            model.addConstr(delta[e.v] >= delta[e.u] + x[e.id]*e.len - (1 - x[e.id])*d, getNewConstr());            // (06)
            model.addConstr(delta[e.v] <= delta[e.u] + x[e.id]*e.len + (1 - x[e.id])*d, getNewConstr());            // (07)
            model.addConstr(eta[e.v] >= eta[e.u] + x[e.id] - (1 - x[e.id])*n, getNewConstr());                      // (10)
            model.addConstr(eta[e.v] <= eta[e.u] + x[e.id] + (1 - x[e.id])*n, getNewConstr());                      // (11)
            model.addConstr(y[getIdx[0][e.u][e.v]] >= x[e.id], getNewConstr());                                     // (14)
        }
        for(int i = 0; i < m/2; ++i)
        {
            // Constraint to ensure the fixed edges are in solution
            if(fixedEdge[i])
            {
                model.addConstr(x[2*i]+x[2*i+1] == 1, getNewConstr());
            }
        }
        for(int v = 0; v < n; ++v)
        {
            model.addConstr(y[getIdx[0][v][v]] == 1, getNewConstr());   //(13)
            for(int u = 0; u < n; ++u)
            {
                linexpr.addTerms(&one, &y[getIdx[0][u][v]], 1);
            }
            model.addConstr(linexpr == eta[v]+1);               //(12)
            linexpr.clear();
        }
        int uv, uw, vw, wuv, wu, wv;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < n; ++v)
            {
                if(u == v) continue;
                uv = getIdx[0][u][v];
                for(int w = 0; w < n; ++w)
                {
                    uw = getIdx[0][u][w];
                    vw = getIdx[0][v][w];
                    wuv = getIdx[w][u][v];
                    wu = getIdx[0][w][u];
                    wv = getIdx[0][w][v];
                    model.addConstr(y[uv] + y[vw] <= 1 + y[uw], getNewConstr());    //(15)
                    model.addConstr(2*z[wuv] <= y[wu] + y[wv], getNewConstr());     //(16)
                    model.addConstr(z[wuv]+1 >= y[wu] + y[wv], getNewConstr());     //(17)
                    for(int& eid : Nmin[w])
                    {
                        linexpr.addTerms(&mEdges[eid].len, &x[eid], 1);
                        linexpr2.addTerms(&mEdges[eid].len, &z[wuv], 1);
                    }
                    model.addConstr(rho[wuv] <= linexpr, getNewConstr());           //(18)
                    model.addConstr(rho[wuv] <= linexpr2, getNewConstr());          //(18)
                    linexpr.clear();
                    linexpr2.clear();
                }
            }
        }
        // Optimize model
        model.optimize();
        // Set new solution to return
        int idx = 0;
        for(int i = 0; i < m; i += 2)
        {
            if(x[i].get(GRB_DoubleAttr_X) > 0.99 || x[i+1].get(GRB_DoubleAttr_X) > 0.99)
            {
                sol.usedEdge[avEdges[idx].id] = true;
                sol.adj[avEdges[idx].u].push_back(AdjInfo(avEdges[idx].v, avEdges[idx].len, avEdges[idx].id));
                sol.adj[avEdges[idx].v].push_back(AdjInfo(avEdges[idx].u, avEdges[idx].len, avEdges[idx].id));
            }
            idx++;
        }
    }
    catch(GRBException e) 
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } 
    catch(...) 
    {
        cout << "Exception during optimization" << endl;
    }
    constrCnt = 0;
    return sol;
}
*/

/*
// hash to store solution for crossover (when solver is called)
struct Hash
{
    double estSearchTime = 0.0;     // time to query in map using EMA
    double estSolverTime = 0.0;     // time to call solver
    double alpha = 0.99;            // parameter for Exponential Moving Average (EMA)
    long long int nCalls = 0;       // number of calls
    double remPercentage = 0.05;    // remove 5% of the table when it's slow
    struct Entry
    {
        Solution sol;               // solver solutions
        long long int solverTime;   // solver time for this call
    };
    struct Order
    {
        double time;
        map<vector<int>, Entry>::iterator it;
        bool operator< (const Order& o) const
        {
            return time < o.time;
        }
    };
    map<vector<int>, Entry> table;
    multiset<Order> ms;              // used to remove fastest solver calls
    chrono::steady_clock::time_point begin, end;
    long long int elapsedTime;
    bool growing = true;
    double timeLimit = 30;           // time limit in seconds of gurobi

    void lookUp(vector<Edge>& avEdges, vb& fixedEdge, Solution& sol)
    {
        nCalls++;
        // key (used edges) to find in map
        vector<int> key((int) avEdges.size());
        int nFixedEdges = 0;
        for(int i = 0; i < (int) avEdges.size(); ++i)
        {
            key[i] = avEdges[i].id;
            if(fixedEdge[i]) 
                nFixedEdges++;
        }
        // search in map and update search time
        begin = chrono::steady_clock::now();
        map<vector<int>, Entry>::iterator it = table.find(key);
        end = chrono::steady_clock::now();
        elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        estSearchTime = alpha*estSearchTime + (1-alpha)*elapsedTime;
        if(leq(estSearchTime*3.0, estSolverTime))
            growing = true;
        else
            growing = false;        
        if(it == table.end())    // key not found insert in map
        {
            if(!growing)
            {
                table.erase(ms.begin()->it);
                ms.erase(ms.begin());
            }
            // call solver and update solver time
            Entry entry;
            begin = chrono::steady_clock::now();
            cout << "called solver with " << (int) avEdges.size() << " edges which " << nFixedEdges << " are already set\n";
            entry.sol = gurobiSolverFlow(avEdges, fixedEdge, timeLimit);
            end = chrono::steady_clock::now();
            elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
            std::cout << "Solver time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << endl;
            estSolverTime += (elapsedTime - estSolverTime)/nCalls;
            entry.sol.computeObjectiveFun();
            entry.solverTime = elapsedTime;
            table[key] = entry;
            sol = entry.sol;
            Order order;
            order.it = table.find(key);
            order.time = elapsedTime;
            ms.insert(order);
            return ;
        }
        // if key is in map
        sol = it->second.sol;
        estSolverTime += (it->second.solverTime - estSolverTime)/nCalls;
    }

};

Hash myHash;
*/

