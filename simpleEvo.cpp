#include <bits/stdc++.h>

using namespace std;

#define mp(a, b) make_pair(a, b)
#define vb vector<bool>
#define vi vector<int>
#define ii pair<int, int>
#define EPS 1e-3
#define TIMEOUT 3600

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
struct Solution
{
    vector<vector<AdjInfo>> adj;
    vb usedEdge;
    double objective;
    
    Solution()
    {
        adj.resize(n, vector<AdjInfo>());
        usedEdge.resize(m, false);
        objective = 0;
    }
    
    void clear()
    {
        for(int i = 0; i < n; ++i)
        {
            adj[i].clear();
        }
        fill(usedEdge.begin(), usedEdge.end(), false);
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
        return usedEdge[idx];
    }

    // Add Edge and Fix Randomly
    void addEdgeAndFix(int id)
    {
        int rdInt = id;
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
        usedEdge[rdInt] = true;
        this->adj[edge.u].push_back(AdjInfo(edge.v, edge.len, edge.id));
        this->adj[edge.v].push_back(AdjInfo(edge.u, edge.len, edge.id));
        cur = edge.v;
        int e2, rmEdge;
        e2 = rmEdge = edge.id;         // e2 represents the removed edge
        // traverse all edges in cycle
        vector<int> possibleEdges;
        while(cur != edge.u)
        {
            e2 = eIdx[cur];
            possibleEdges.push_back(e2);
            cur = dad[cur];
        }
        rmEdge = possibleEdges[rand()%((int) possibleEdges.size())];
        // remove the edge s.t. objective function is minimum
        edge = edges[rmEdge];
        this->removeEdge(edge);
        this->computeObjectiveFun();
    }

    // Remove Edge and Fix Randomly
    void removeEdgeAndFix(int id)
    {
        int rdInt = id;
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
        int addEdge;
        vector<int> possibleEdges;
        for(int i = 0; i < m; ++i)
        {
            // XOR is used to ensure the edge is connecting the two disconnected sets A and B
            if(!hasEdge(i) && (inA[edges[i].u]^inA[edges[i].v]))
            {
                possibleEdges.push_back(i);
            }            
        }
        addEdge = possibleEdges[rand()%(int) possibleEdges.size()];
        // Insert the best edge in solution
        Edge bestEdge = edges[addEdge];
        this->usedEdge[addEdge] = true;
        this->adj[bestEdge.u].push_back(AdjInfo(bestEdge.v, bestEdge.len, bestEdge.id));
        this->adj[bestEdge.v].push_back(AdjInfo(bestEdge.u, bestEdge.len, bestEdge.id));
        this->computeObjectiveFun();
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
            print(edges[i]);
    }
    putchar('\n');
    printf("Objective value = %.2f\n", s.objective);
    putchar('\n');
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
        genRandomPop();
        for(Solution& sol : solutions)
        {
            assert(lt(0, sol.objective) && lt(sol.objective, DBL_MAX));
        }
        int gen = 1;
        double maxObj, minObj;
        Solution best;
        best.objective = DBL_MAX;
        int ellapsed;
        double curBestVal = DBL_MAX;
        Solution* tmpBest;
        while(gen <= numGen)
        {
            c = chrono::steady_clock::now();
            ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
            if(ellapsed >= TIMEOUT/2.0)
            {
                break;
            }
            int nEdgesUsed = 0;
            vb edgesUsed(m, false);
            for(Solution& sol : solutions)
            {
                for(int j = 0; j < m; ++j)
                {
                    if(sol.usedEdge[j])
                    {
                        if(!edgesUsed[j])
                        {
                            edgesUsed[j] = true;
                            nEdgesUsed++;
                        }
                    }
                }
            }
            if(nEdgesUsed == n-1)
                break;
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
            // Crossover between parents
            int id1, id2;
            for(int i = 0; i < numCrossover; ++i)
            {
                id1 = rand()%(popSize);
                id2 = rand()%(popSize);
                offspring[i] = crossover(solutions[id1], solutions[id2]);
            }
            int idx = numCrossover;
            for(int i = 0; i < popSize; ++i)
            {
                offspring[idx++] = solutions[i];
            }
            tournamentSelection(offspring);
            if(lt(best.objective, curBestVal))
            {
                curBestVal = best.objective;
            }
            gen++;
            printf("Best so far = %.10f\n", best.objective);
        }
        SimulatedAnnealing(best);
        return best;
    }

    void SimulatedAnnealing(Solution& best)
    {
        int numSol = 1000;
        vector<Solution> iniSol(numSol);
        genRandom(iniSol, numSol);
        c = chrono::steady_clock::now();
        int ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
        if(ellapsed >= TIMEOUT)
        {
            return ;
        }
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
        int rdInt;
        while(iter < 20000 && lastImprove < 2000)
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
            rdInt = rand()%m;
            if(tmp.usedEdge[rdInt])
            {
                tmp.removeEdgeAndFix(rdInt);
            }
            else
            {
                tmp.addEdgeAndFix(rdInt);
            }   
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

    Solution crossover(Solution& s1, Solution& s2)
    {
        vector<Edge> avEdges;
        vb fixedEdge;
        vector<bool> setEdges(m, false);
        Solution sol;
        vector<pair<double, int>> KruskalRST;
        UnionFind uf(n);
        for(int i = 0; i < m; ++i)
        {
            if(s1.usedEdge[i] && s2.usedEdge[i])
            {
                sol.usedEdge[i] = true;
                sol.adj[edges[i].u].push_back(AdjInfo(edges[i].v, edges[i].len, edges[i].id));
                sol.adj[edges[i].v].push_back(AdjInfo(edges[i].u, edges[i].len, edges[i].id));
                uf.unionSet(edges[i].u, edges[i].v);
            }
            else if(s1.usedEdge[i] || s2.usedEdge[i])
            {
                KruskalRST.push_back(mp(edges[i].len, edges[i].id));
            }
        }
        sort(KruskalRST.begin(), KruskalRST.end());
        Edge* e;
        for(pair<double, int> pp : KruskalRST)
        {
            e = &edges[pp.second];
            if(!uf.isSameSet(e->u, e->v))
            {
                uf.unionSet(e->u, e->v);
                sol.usedEdge[e->id] = true;
                sol.adj[e->u].push_back(AdjInfo(e->v, e->len, e->id));
                sol.adj[e->v].push_back(AdjInfo(e->u, e->len, e->id));
            }
        }
        sol.computeObjectiveFun();
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
            sol.computeObjectiveFun();
            solutions[i] = sol;
        }
    }

    void genRandom(vector<Solution>& iniSol, int numSol)
    {
        vector<Edge> cpy = edges;
        int numForests;
        int ellapsed;
        for(int i = 0; i < numSol; ++i)
        {
            c = chrono::steady_clock::now();
            ellapsed = std::chrono::duration_cast<std::chrono::seconds>(c - a).count();
            if(ellapsed >= TIMEOUT)
            {
                break;
            }
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
    for(int seedid = 0; seedid < 5; ++seedid)
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
