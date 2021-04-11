#include <bits/stdc++.h>

using namespace std;

#define mp(a, b) make_pair(a, b)
#define vb vector<bool>
#define vi vector<int>
#define ii pair<int, int>
#define EPS 1e-3
#define TIMEOUT 1200

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

chrono::steady_clock::time_point a, b, c;

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

vector<int> penalty;
const double alpha = 0.3;
double ft1, lambda;

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

    // Mutate when considering to remove a random edge - O(m*n^2)
    void mutateRemoving()
    {
        this->fillDist();
        // selecting edge to remove
        vi possibleEdges(n-1);
        int idx = 0;
        for(auto it = usedEdges.begin(); it != usedEdges.end(); ++it)
        {
            possibleEdges[idx++] = *it;
        }
        assert(idx == n-1);
        while(true)
        {
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
            // Find the best edge to add when removing the chosen edge
            double minObj, curObj;
            int addEdge;
            minObj = curObj = this->objective;
            addEdge = rdInt;
            vector<int> possibleEdges;
            for(int i = 0; i < m; ++i)
            {
                // XOR is used to ensure the edge is connecting the two disconnected sets A and B
                if(!hasEdge(i) && (inA[edges[i].u]^inA[edges[i].v]))
                {
                    possibleEdges.push_back(i);
                }            
            }
            shuffle(possibleEdges.begin(), possibleEdges.end(), default_random_engine(seed));
            int cnt = 0;
            for(int& i : possibleEdges)
            {
                cnt++;
                curNode = edges[i].u;
                list<int> nodesToUpdate;
                for(int j = 0; j < n; ++j)
                {
                    if(inA[curNode]^inA[j]) // need to be updated
                    {
                        curObj -= dist[curNode][j]*req[curNode][j];
                        dist[curNode][j] = dist[j][curNode] = edges[i].len + dist[edges[i].v][j];
                        curObj += dist[curNode][j]*req[curNode][j];
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
                            curObj -= dist[ainfo.v][j]*req[ainfo.v][j];
                            dist[ainfo.v][j] = dist[j][ainfo.v] = ainfo.len + dist[curNode][j];
                            curObj += dist[ainfo.v][j]*req[ainfo.v][j];
                        }
                    }
                }
                if(curObj < minObj)
                {
                    minObj = curObj;
                    addEdge = i;
                }
                if(cnt >= 100)
                    break;
            }
            // Insert the best edge in solution
            Edge bestEdge = edges[addEdge];
            //this->usedEdge[addEdge] = true;
            usedEdges.insert(addEdge);
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
                    dist[curNode][i] = dist[i][curNode] = bestEdge.len + dist[neighbor][i];
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
                            dist[ainfo.v][i] = dist[i][ainfo.v] = ainfo.len + dist[curNode][i];
                        }
                    }
                }
            }
            // assertion to check if distances are correctly updated
            double tmp = 0;
            for(int i = 0; i < n; ++i)
                for(int j = i+1; j < n; ++j)
                    tmp += dist[i][j]*req[i][j];
            assert(eq(tmp, this->objective));
            break;
        }
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

    // Neighbors of current solution
    void mutateRemoving(Solution& best, int& notImproving)
    {
        this->fillDist();
        // selecting edge to remove
        vi possibleEdges(n-1);
        int idx = 0;
        for(auto it = usedEdges.begin(); it != usedEdges.end(); ++it)
        {
            possibleEdges[idx++] = *it;
        }
        assert(idx == n-1);
        while(true)
        {
            int rdInt = possibleEdges[rand()%(int) possibleEdges.size()];
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
            // Find the best edge to add when removing the chosen edge
            double minObj, curObj;
            int addEdge;
            minObj = curObj = this->objective;
            addEdge = rdInt;
            vector<int> possibleEdges;
            for(int i = 0; i < m; ++i)
            {
                // XOR is used to ensure the edge is connecting the two disconnected sets A and B
                if(!hasEdge(i) && (inA[edges[i].u]^inA[edges[i].v]))
                {
                    possibleEdges.push_back(i);
                }            
            }
            shuffle(possibleEdges.begin(), possibleEdges.end(), default_random_engine(seed));
            int cnt = 0;
            for(int& i : possibleEdges)
            {
                cnt++;
                Solution tmp = *this;
                tmp.usedEdges.insert(i);
                tmp.adj[edges[i].u].push_back(AdjInfo(edges[i].v, edges[i].len, edges[i].id));
                tmp.adj[edges[i].v].push_back(AdjInfo(edges[i].u, edges[i].len, edges[i].id));
                tmp.computeObjectiveFun();
                if(lt(tmp.objective, best.objective))
                {
                    best = tmp;
                    printf("Objective = %.10f\n", best.objective);
                    notImproving = 0;
                }
                if(tmp.objective < minObj)
                {
                    minObj = tmp.objective;
                    addEdge = i;
                }
                if(cnt >= 100)
                    break;
            }
            // Insert the best edge in solution
            Edge bestEdge = edges[addEdge];
            //this->usedEdge[addEdge] = true;
            usedEdges.insert(addEdge);
            this->adj[bestEdge.u].push_back(AdjInfo(bestEdge.v, bestEdge.len, bestEdge.id));
            this->adj[bestEdge.v].push_back(AdjInfo(bestEdge.u, bestEdge.len, bestEdge.id));
            this->computeObjectiveFun();
            break;
        }
    }

    void perturbate()
    {
        vector<int> uEdges(n-1);
        int cnt = 0;
        for(auto it = usedEdges.begin(); it != usedEdges.end(); ++it)
        {
            uEdges[cnt++] = *it;
        }
        shuffle(uEdges.begin(), uEdges.end(), default_random_engine(seed));
        Solution tmp;
        Edge* e;
        UnionFind uf(n);
        for(int i = 0; i < (n-1)/2; ++i)
        {
            e = &edges[uEdges[i]];
            tmp.usedEdges.insert(e->id);
            tmp.adj[e->u].push_back(AdjInfo(e->v, e->len, e->id));
            tmp.adj[e->v].push_back(AdjInfo(e->u, e->len, e->id));
            uf.unionSet(e->u, e->v);
        }
        vector<pair<double, int>> KruskalRST;
        for(int i = 0; i < m; ++i)
        {
            if(tmp.usedEdges.find(i) == tmp.usedEdges.end())
            {
                KruskalRST.push_back(mp(edges[i].len, edges[i].id));
            }
        }
        sort(KruskalRST.begin(), KruskalRST.end());
        for(pair<double, int> pp : KruskalRST)
        {
            e = &edges[pp.second];
            if(!uf.isSameSet(e->u, e->v))
            {
                uf.unionSet(e->u, e->v);
                tmp.usedEdges.insert(e->id);
                tmp.adj[e->u].push_back(AdjInfo(e->v, e->len, e->id));
                tmp.adj[e->v].push_back(AdjInfo(e->u, e->len, e->id));
                if((int) tmp.usedEdges.size() == n-1)
                    break;
            }
        }
        tmp.computeObjectiveFun();
        *this = tmp;
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


void localSearch(Solution& best)
{
    Solution bestIt, imp, tmp;
    int ellapsed;
    printf("Objective = %.10f\n", best.objective);
    bestIt = best;
    int notImproving;
    notImproving = 0;
    do
    {
        notImproving++;
        imp = bestIt;
        c = chrono::steady_clock::now();
        ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
        if(ellapsed >= TIMEOUT)
        {
            break;
        }
        tmp = bestIt;
        tmp.mutateRemoving(best, notImproving);
        if(lt(tmp.objective, imp.objective))
        {
            imp = tmp;
            notImproving = 0;
        }
        if(lt(tmp.objective, best.objective))
        {
            best = tmp;
        }
        bestIt = imp;
    } while(notImproving < 25);
}


Solution GLS(Solution& start)
{
    fill(penalty.begin(), penalty.end(), 0);
    ft1 = 0.0;
    lambda = 0.0;
    Solution best, bestIt, imp, tmp, newSol;
    best = start;
    int ellapsed;
    printf("Objective = %.10f\n", best.objective);
    bestIt = best;
    int notImproving;
    notImproving = 0;
    do
    {
        notImproving++;
        imp = bestIt;
        c = chrono::steady_clock::now();
        ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
        if(ellapsed >= TIMEOUT)
        {
            return best;
        }
        tmp = bestIt;
        tmp.mutateRemoving(best, notImproving);
        if(lt(tmp.objective, imp.objective))
        {
            notImproving = 0;
            imp = tmp;
        }
        if(lt(tmp.objective, best.objective))
        {
            best = tmp;
        }
        bestIt = imp;
    } while(notImproving < 25);
    ft1 = bestIt.objective;
    lambda = alpha*ft1/(n-1);
    std::uniform_real_distribution<double> distrib(0.0, 1.0);
    double rngVal;
    int bestImprovement = 0;
    do
    {
        bestImprovement++;
        double mn, mx, util;
        mn = DBL_MAX;
        mx = 0;
        for(auto it = bestIt.usedEdges.begin(); it != bestIt.usedEdges.end(); ++it)
        {
            mn = min(mn, edges[*it].len);
            mx = max(mx, edges[*it].len);
        }
        for(auto it = bestIt.usedEdges.begin(); it != bestIt.usedEdges.end(); ++it)
        {
            if(mn == mx) 
                util = 0.1;
            else
                util = (edges[*it].len-mn)/(mx-mn);
            util = min(util, 0.8);
            rngVal = distrib(rng);
            if(leq(rngVal, util))
            {
                penalty[*it]++;
            }
        }
        notImproving = 0;
        double bestObj = best.objective;
        do
        {
            notImproving++;
            imp = bestIt;
            c = chrono::steady_clock::now();
            ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
            if(ellapsed >= TIMEOUT)
            {
                return best;
            }
            tmp = bestIt;
            tmp.mutateRemoving(best, notImproving);
            if(lt(tmp.objective, imp.objective))
            {  
                notImproving = 0;
                imp = tmp;
            }
            if(lt(tmp.objective, best.objective))
            {
                best = tmp;
            }
            bestIt = imp;
        } while (notImproving < 25);
        
        if(lt(best.objective, bestObj))
        {
            bestObj = best.objective;
            bestImprovement = 0;
        }

    } while (bestImprovement >= 25);
    
    return best;
}

// evolutionary/genetic algorithm
struct Evolutionary
{
    vector<Solution> solutions;
    vector<Solution> offspring;
    int popSize;                // initial population size
    int numGen;                 // number of generations
    int numCrossover;           // number of crossovers
    int offspringSize;
    Evolutionary(int popSize, int numGen, int numCrossover)
    {
        this->popSize = popSize;
        this->numGen = numGen;
        this->numCrossover = numCrossover;
        this->offspringSize = numCrossover+popSize;
        solutions.resize(popSize);
        offspring.resize(offspringSize);
    }

    Solution run()
    {
        return ILS();
    }

    void SimulatedAnnealing(Solution& best)
    {
        int numSol = 100;
        vector<Solution> iniSol(numSol);
        genRandom(iniSol, numSol);
        double sum, avg, sd;
        sum = 0.0;
        for(int i = 0; i < numSol; ++i)
        {
            sum += iniSol[i].objective;
        }
        avg = (double) sum/numSol;
        sd = 0;
        for(int i = 0; i < numSol; ++i)
        {
            sd += (iniSol[i].objective-avg)*(iniSol[i].objective-avg);
        }
        sd = sqrt((double) sd/numSol);
        double T = 2*sd;
        int iter, lastImprove;
        iter = lastImprove = 0;
        Solution tmp, curIt;
        curIt = best;
        std::uniform_real_distribution<double> distrib(0.0, 1.0);
        double rngVal;
        int ellapsed;
        while(iter < 200000 && lastImprove < 10000)
        {
            lastImprove++;
            iter++;
            c = chrono::steady_clock::now();
            ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
            if(ellapsed >= TIMEOUT)
            {
                break;
            }
            tmp = curIt;
            tmp.mutateRemoving();
            if(lt(tmp.objective, best.objective))
            {
                lastImprove = 0;
                best = tmp;
                printf("Obj = %.10f\n", best.objective);
            }
            if(lt(tmp.objective, curIt.objective))
            {
                curIt = tmp;
            }
            else
            {
                double prob = exp((curIt.objective-tmp.objective)/(T+EPS));
                prob = min(prob, 1.0);
                prob = max(prob, 0.0);
                assert(leq(0.0, prob) && leq(prob, 1.0));
                rngVal = distrib(rng);
                if(lt(rngVal, prob))
                {
                    curIt = tmp;
                }
            }
            T *= 0.99;
        }
    }

    Solution ILS()
    {
        Solution sol;
        genMinPath(sol);
        Solution best = sol;
        int ellapsed;
        int notImproving = 0;
        do
        {
            notImproving++;
            c = chrono::steady_clock::now();
            ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
            if(ellapsed >= TIMEOUT)
            {
                break;
            }
            localSearch(sol);
            if(lt(sol.objective, best.objective))
            {
                best = sol;
                notImproving = 0;
            }
            sol.perturbate();
        } while(notImproving < 10);
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
            // Calling a greedy shortest path tree from a random node
            buildMinPathSolution(avEdges, sol);
            sol.computeObjectiveFun();
        }
        return sol;
    }

    void tournamentSelection(vector<Solution>& offspring)
    {
        int a, b;
        vector<int> nTimes(offspringSize, 0);
        set<int> s;
        for(int i = 0; i < offspringSize; ++i)
            s.insert(i);
        int idx = 0;
        while(s.size())
        {
            auto it1 = s.begin();
            advance(it1, rand()%(s.size()));
            nTimes[*it1]++;
            a = *it1;
            s.erase(it1);
            auto it2 = s.begin();
            advance(it2, rand()%(s.size()));
            nTimes[*it2]++;
            b = *it2;
            s.erase(it2);
            if(lt(offspring[b].objective, offspring[a].objective))
                a = b;
            solutions[idx++] = offspring[a];
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

    void genMinPath(Solution& sol)
    {
        printf("MinPathPop\n");
        // generate adjacency list to perform Dijkstra
        vector<AdjInfo> adj[n];
        for(Edge& e : edges)
        {
            adj[e.u].push_back(AdjInfo(e.v, e.len, e.id));
            adj[e.v].push_back(AdjInfo(e.u, e.len, e.id));
        }
        int cur;
        vector<double> dist(n, DBL_MAX);
        vector<int> uEdge(n, -1);
        cur = rand()%n;
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
    }

    void genRandom(vector<Solution>& iniSol, int numSol)
    {
        vector<Edge> cpy = edges;
        int numForests;
        for(int i = 0; i < numSol; ++i)
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
            iniSol[i] = sol;
        }
    }

};

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        printf("usage: ./simpleEvo popSize numGen numCrossovers < inputFile\n");
        return -1;
    }
    cin >> n >> m;
    edges.resize(m);
    for(int i = 0; i < m; ++i)
    {
        cin >> edges[i].u >> edges[i].v >> edges[i].len;
        edges[i].id = i;
    }
    seen.resize(n);
    penalty.resize(m);
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
    for(int seedid = 0; seedid < 10; ++seedid)
    {
        seed = seedVector[seedid];
        printf("seed = %u\n", seed);
        //Initialize seeds
        srand(seed);
        rng.seed(seed);
        Evolutionary ev(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
        a = chrono::steady_clock::now();
        Solution best = ev.run();
        printf("Best Value Found = %.10f\n", best.objective);
        b = chrono::steady_clock::now();
        cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(b - a).count() << endl;
        log << best.objective << "," <<  std::chrono::duration_cast<std::chrono::seconds>(b - a).count() << endl;
        double tmp = best.objective;
        best.computeObjectiveFun();
        assert(eq(best.objective, tmp));
    }
    log.close();
    return 0;
}
