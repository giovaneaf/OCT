#include <bits/stdc++.h>

using namespace std;

#define mp(a, b) make_pair(a, b)
#define vb vector<bool>
#define vi vector<int>
#define ii pair<int, int>
#define EPS 1e-3
#define TIMEOUT 300

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
vector<vector<AdjInfo>> adjList;

// seed used to generate random numbers
unsigned seed;
// seeds used for testing
unsigned seedVector[] = {280192806, 871237442, 2540188929, 107472404, 3957311442, 316851227, 619606212, 1078082709, 916212990, 698598169};

//Mersenne Twister: Good quality random number generator
std::mt19937 rng;

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
    double h;
    
    Solution()
    {
        adj.resize(n, vector<AdjInfo>());
        objective = 0;
        h = 0;
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
        h = objective;
        double acc = 0.0;
        for(auto it = usedEdges.begin(); it != usedEdges.end(); ++it)
        {
            acc += penalty[*it]*edges[*it].len;
        }
        h += lambda*acc;
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
            minObj = curObj = this->h;
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
                    notImproving = 0;
                    printf("Objective = %.10f\n", best.objective);
                }
                if(tmp.h < minObj)
                {
                    minObj = tmp.h;
                    addEdge = i;
                    break;
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

chrono::steady_clock::time_point a, b, c;

void genRandomSol(Solution& sol)
{
    vector<Edge> cpy = edges;
    int numForests;
    shuffle(begin(cpy), end(cpy), default_random_engine(seed));
    UnionFind uf(n);
    numForests = n;
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
}

Solution GLS()
{
    fill(penalty.begin(), penalty.end(), 0);
    ft1 = 0.0;
    lambda = 0.0;
    Solution best, bestIt, imp, tmp, newSol;
    genRandomSol(best);
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
        if(lt(tmp.h, imp.h))
        {
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
            if(lt(tmp.h, imp.h))
            {
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


int main()
{
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
        a = chrono::steady_clock::now();
        Solution best = GLS();
        printf("Best Value Found = %.10f\n", best.objective);
        b = chrono::steady_clock::now();
        cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(b - a).count() << endl;
        log << best.objective << "," <<  std::chrono::duration_cast<std::chrono::seconds>(b - a).count() << endl;
    }
    log.close();
    return 0;
}