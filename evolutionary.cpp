#include <bits/stdc++.h>
#include <gurobi_c++.h>

using namespace std;

#define vb vector<bool>
#define vi vector<int>
#define EPS 1e-9

struct Edge 
{
    int u, v; 
    double len;
    int id;
    Edge(int u = 0, int v = 0, double len = 0.0, int id = 0) : u(u), v(v), len(len), id(id) {}
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
//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
unsigned seed = 24438823;

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
        assert(minObj == this->objective);
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
                // assertion to check if curObj value is correct
                Solution ss = *this;
                ss.usedEdge[i] = true;
                ss.adj[edges[i].u].push_back(AdjInfo(edges[i].v, edges[i].len, edges[i].id));
                ss.adj[edges[i].v].push_back(AdjInfo(edges[i].u, edges[i].len, edges[i].id));
                ss.computeObjectiveFun();
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
        int tmp = 0;
        for(int i = 0; i < n; ++i)
            for(int j = i+1; j < n; ++j)
                tmp += this->dist[i][j]*req[i][j];
        assert(tmp == this->objective);
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

// Variable used for Gurobi solver
vector<vector<double>> reqMin2;
vector<vector<vector<int>>> getIdx;
static bool computeValues = true;
static int constrCnt = 0;

string getNewConstr()
{
    return "C" + to_string(constrCnt++);
}

Solution gurobiSolver(vector<Edge>& avEdges, vb& fixedEdge)
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
        print(e);
        printf("Fixed = %d\n", fixedEdge[cnt/2] ? 1 : 0);
        mEdges[cnt] = {e.u, e.v, e.len, cnt};
        mEdges[cnt+1] = {e.v, e.u, e.len, cnt+1};
        Nmin[e.v].push_back(cnt);
        Nmin[e.u].push_back(cnt+1);
        d += e.len;
        cnt += 2;
    }
    try 
    {
         // Create an environment
        GRBEnv env = GRBEnv(true);
        //env.set("LogFile", "mip.log");
        env.set("OutputFlag", "0");
        env.start();
        

        // Create an empty model
        GRBModel model = GRBModel(env);

        if(computeValues)
        {
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
            computeValues = false;
        }

        // Create binary variables
        GRBVar* x = model.addVars(m, GRB_BINARY);
        GRBVar* y = model.addVars(n*n, GRB_BINARY);
        GRBVar* z = model.addVars(n*n*n, GRB_BINARY);

        // Create integer variables
        GRBVar* delta = model.addVars(n, GRB_INTEGER);
        GRBVar* eta = model.addVars(n, GRB_INTEGER);
        GRBVar* rho = model.addVars(n*n*n, GRB_INTEGER);

        

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
        cnt = 0;
        for(Edge& e : mEdges)
        {           
            model.addConstr(delta[e.v] >= delta[e.u] + x[e.id]*e.len - (1 - x[e.id])*d, getNewConstr());            // (06)
            model.addConstr(delta[e.v] <= delta[e.u] + x[e.id]*e.len + (1 - x[e.id])*d, getNewConstr());            // (07)
            model.addConstr(eta[e.v] >= eta[e.u] + x[e.id] - (1 - x[e.id])*n, getNewConstr());                      // (10)
            model.addConstr(eta[e.v] <= eta[e.u] + x[e.id] + (1 - x[e.id])*n, getNewConstr());                      // (11)
            model.addConstr(y[getIdx[0][e.u][e.v]] >= x[e.id], getNewConstr());                                     // (14)
            // Constraint to ensure the fixed edges are in solution
            if(fixedEdge[cnt])
            {
                model.addConstr(x[2*cnt]+x[2*cnt+1] == 1, getNewConstr());
            }
            cnt++;
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
        return sol;
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
        model.reset();    

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

bool inline leq(double a, double b)
{
    return a < b || abs(a-b) < EPS;
}

// evolutionary/genetic algorithm
struct Evolutionary
{
    vector<Solution> solutions;
    vector<int> parents;
    vector<double> fitness;
    int popSize, numPar;
    Evolutionary(int popSize, int numPar)
    {
        solutions.resize(popSize);
        parents.resize(numPar);
        fitness.resize(popSize);
        this->popSize = popSize;
        this->numPar = numPar;
    }

    Solution run()
    {
        genInitialPop(popSize);
        int gen = 1;
        double maxObj, minObj;
        Solution best;
        best.objective = DBL_MAX;
        double fitSum;
        double rngDbl;
        double accVal;
        int rngInt;

        //Mersenne Twister: Good quality random number generator
        std::mt19937 rng; 
        //Initialize with non-deterministic seeds
        rng.seed(seed); 

        while(gen <= 5)
        {
            minObj = DBL_MAX;
            maxObj = 0;
            // Evaluate solutions
            for(Solution& sol : solutions)
            {
                sol.computeObjectiveFun();
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
                    fitness[i] = 1.0 - (solutions[i].objective - minObj)/(maxObj - minObj);
                fitSum += fitness[i];
            }
            // selecting numPar parents
            // Never select worst solution found? (fitness = 0)
            std::uniform_real_distribution<double> distrib(0.0, fitSum);
            for(int i = 0; i < numPar; ++i)
            {
                rngDbl = distrib(rng);
                accVal = 0.0;
                for(int j = 0; j < popSize; ++j)
                {
                    if(leq(rngDbl, accVal + fitness[j]))    // solution chosen
                    {
                        parents[i] = j;
                        break;
                    }
                    accVal += fitness[j];
                }
            }
            // Crossover between parents
            vector<Solution> offspring;
            for(int i = 0; i < numPar; ++i)
                for(int j = i+1; j < numPar; ++j)
                    offspring.push_back(crossover(solutions[parents[i]], solutions[parents[j]]));
                    //offspring.push_back(solutions[parents[i]]);
            return best;
            for(Solution& sol : offspring)
            {
                rngInt = rand()%2;
                if(rngInt)
                {
                    sol.mutateInserting();
                }
                else
                {
                    sol.mutateRemoving();
                }
                if(sol.objective < best.objective)      // update if solution is better
                {
                    best = sol;
                }
            }

            gen++;
        }
        return best;
    }

    Solution crossover(Solution& s1, Solution& s2)
    {
        printf("Crossovering solutions...\n");
        print(s1);
        print(s2);
        bool equal = true;
        vector<Edge> avEdges;
        vb fixedEdge;
        int nFixedEdges = 0;
        for(int i = 0; i < m; ++i)
        {
            if((!s1.usedEdge[i]) && (!s2.usedEdge[i]))
                continue;
            // Edge used in at least one tree
            avEdges.push_back(edges[i]);
            if(s1.usedEdge[i] && s2.usedEdge[i])
            {
                fixedEdge.push_back(true);
                nFixedEdges++;
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
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            sol = gurobiSolver(avEdges, fixedEdge);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            //std::cout << "Solver time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms with " << (int) avEdges.size() << " edges which " << nFixedEdges << " are already set\n";
            sol.computeObjectiveFun();
        }
        return sol;
    }

    /* Generate popSize initial solutions (trees) by shuffling the edges
	   and inserting the edges like Kruskal Algorithm */
    void genInitialPop(int popSize)
    {
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
            solutions[i] = sol;
        }
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
    req.resize(n, vector<double>(n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            cin >> req[i][j];
            req[j][i] = req[i][j];
        }
    }
    /*vb fixedEdge(m, false);
    fixedEdge[0] = fixedEdge[3] = true;
    Solution s = gurobiSolver(edges, fixedEdge);
    s.computeObjectiveFun();
    print(s);
    return 0;*/
    Evolutionary ev(5, 3);
    Solution best = ev.run();
    print(best);
    return 0;
}