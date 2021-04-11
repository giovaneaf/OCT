#include <bits/stdc++.h>
#include <gurobi_c++.h>

using namespace std;

struct Edge 
{
    int u, v;
    double len;
    int id;
    Edge(int u = 0, int v = 0, double len = 0.0, int id = 0) : u(u), v(v), len(len), id(id) {}
};

vector<vector<double>> req;

// Vectors used for gurobi
vector<vector<double>> reqMin2;
vector<vector<int>> Nmin;
vector<vector<vector<int>>> getIdx;
vector<Edge> mEdges;

string getNewConstr()
{
    static int constrCnt = 0;
    return "C" + to_string(constrCnt++);
}
inline void print(Edge& e)
{
    printf("(%d, %d, %.2f, %d)\n", e.u, e.v, e.len, e.id);
}


int n, m;

void newFormulation()
{
    double d = 0;
    cin >> n >> m;
    m *= 2;
    mEdges.resize(m);
    Nmin.resize(n, vector<int>());
    for(int i = 0; i < m; i += 2)
    {
        cin >> mEdges[i].u >> mEdges[i].v >> mEdges[i].len;
        mEdges[i].id = i;
        mEdges[i+1] = {mEdges[i].v, mEdges[i].u, mEdges[i].len, i+1};
        Nmin[mEdges[i].v].push_back(i);
        Nmin[mEdges[i].u].push_back(i+1);
        d += mEdges[i].len;
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
    GRBEnv env = GRBEnv(true);
    try 
    {
        env.set("OutputFlag", "1");
        env.set("TimeLimit", "1800");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);
        
        static bool computeValues = true;

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
            int cnt = 0;
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
            model.addConstr(delta[e.v] >= delta[e.u] + x[e.id]*e.len - (1 - x[e.id])*d, getNewConstr());    // (06)
            model.addConstr(delta[e.v] <= delta[e.u] + x[e.id]*e.len + (1 - x[e.id])*d, getNewConstr());    // (07)
            model.addConstr(eta[e.v] >= eta[e.u] + x[e.id] - (1 - x[e.id])*n, getNewConstr());              // (10)
            model.addConstr(eta[e.v] <= eta[e.u] + x[e.id] + (1 - x[e.id])*n, getNewConstr());              // (11)
            model.addConstr(y[getIdx[0][e.u][e.v]] >= x[e.id], getNewConstr());                                     // (14)
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

        /*for(int i = 0; i < m; ++i)
        {
            cout << x[i].get(GRB_StringAttr_VarName) << " "
                << x[i].get(GRB_DoubleAttr_X) << '\n';
        }*/

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << '\n';


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
}

void flowFormulation()
{
    cin >> n >> m;
    vector<Edge> avEdges(m);
    for(int i = 0; i < m; i++)
    {
        cin >> avEdges[i].u >> avEdges[i].v >> avEdges[i].len;
        avEdges[i].id = i;
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
    // Variables used for solver
    GRBEnv env = GRBEnv(true);
    vector<vector<int>> getIdxFlow;
    vector<double> O;
    vector<vector<int>> Nmin, Nplus;
    vector<Edge> mEdges;
    double d = 0;
    getIdxFlow.resize(n, vector<int>(2*m));
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
        env.set("OutputFlag", "1");
        env.set("TimeLimit", "1800");
        env.start();
        O.resize(n, 0.0);
        for(int u = 0; u < n; ++u)
        {
            for(int v = u+1; v < n; ++v)
            {
                O[u] += req[u][v];
            }
        }
        int cnt = 0;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < m; ++v)
            {
                getIdxFlow[u][v] = cnt++;
            }
        }
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
        // Optimize model
        model.optimize();
        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << '\n';
        for(int i = 0; i < m/2; ++i)
        {
            if(x[i].get(GRB_DoubleAttr_X) > 0.99)
            {
                printf("(%d, %d, %f)\n", avEdges[i].u, avEdges[i].v, avEdges[i].len);
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
}

void rootedFormulation()
{
    cin >> n >> m;
    vector<Edge> avEdges(m);
    double bigM = 0.0;
    vector<vector<double>> cost(n, vector<double>(n, 1000000000));
    for(int i = 0; i < m; i++)
    {
        cin >> avEdges[i].u >> avEdges[i].v >> avEdges[i].len;
        avEdges[i].id = i;
        cost[avEdges[i].u][avEdges[i].v] = avEdges[i].len; 
        cost[avEdges[i].v][avEdges[i].u] = avEdges[i].len;
        bigM = max(bigM, avEdges[i].len);
    }
    bigM *= (n-1);
    req.resize(n, vector<double>(n));
    for(int i = 0; i < n; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            cin >> req[i][j];
            req[j][i] = req[i][j];
        }
    }
    // Variables used for solver
    GRBEnv env = GRBEnv(true);
    vector<vector<int>> getIdxFlow;
    vector<double> O;
    vector<vector<int>> Nmin, Nplus;
    vector<Edge> mEdges;
    getIdxFlow.resize(n, vector<int>(2*m));
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
        cnt += 2;
    }
    try 
    {
        env.set("OutputFlag", "1");
        env.set("TimeLimit", "1800");
        env.start();
        int cnt = 0;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < n; ++v)
            {
                if(u == v) continue;
                getIdxFlow[u][v] = cnt++;
            }
        }
        // Create an empty model
        GRBModel model = GRBModel(env);
        // Create binary variables
        GRBVar* x = model.addVars(n*(n-1), GRB_BINARY);
        GRBVar* p = model.addVars(n*(n-1), GRB_BINARY);
        // Create continuous variables
        GRBVar* d = model.addVars(n*n, GRB_CONTINUOUS);
        GRBLinExpr obj, linexpr, linexpr2;;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < n; ++v)
            {
                if(u == v) continue;
                obj.addTerms(&req[u][v], &d[getIdxFlow[u][v]], 1);
            }
        }
        model.setObjective(obj, GRB_MINIMIZE);              // (a)
        const double one = 1.0;
        for(int j = 0; j < n; ++j)
        {
            linexpr.clear();
            for(int i = 0; i < n; ++i)
            {
                if(i == j) continue;
                linexpr.addTerms(&one, &x[getIdxFlow[i][j]], 1);
            }
            model.addConstr(linexpr == 1 - (j == 0 ? 1 : 0));                   // (b)
        }
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; ++j)
            {
                if(i == j) continue;
                model.addConstr(x[getIdxFlow[i][j]] <= p[getIdxFlow[i][j]]);        // (c)
            }
        }

        for(int i = 0; i < n; i++)
        {
            for(int j = i+1; j < n; ++j)
            {
                model.addConstr(p[getIdxFlow[i][j]] + p[getIdxFlow[j][i]] <= 1);        // (d)
            }
        }

        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                if(i == j) continue;
                for(int k = 0; k < n; ++k)
                {
                    if(k == i || k == j) continue;
                    model.addConstr(p[getIdxFlow[i][j]] + x[getIdxFlow[j][k]] <= 1 + p[getIdxFlow[i][k]]); // (e)
                    model.addConstr(p[getIdxFlow[i][k]] + x[getIdxFlow[j][k]] <= 1 + p[getIdxFlow[i][j]]); // (f)
                    model.addConstr(d[getIdxFlow[i][j]] >= d[getIdxFlow[i][k]] + cost[k][j] - bigM*(2-x[getIdxFlow[k][j]] - p[getIdxFlow[i][k]])); // (g)
                    model.addConstr(d[getIdxFlow[i][j]] >= d[getIdxFlow[i][k]] + cost[k][j] - bigM*(1-x[getIdxFlow[k][j]] + p[getIdxFlow[i][j]] + p[getIdxFlow[j][i]])); // (h)
                }
            }
        }

        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                //printf("%d %d %f\n", i, j, cost[i][j]);
                if(i == j) continue;
                model.addConstr(d[getIdxFlow[i][j]] >= cost[i][j]*x[getIdxFlow[i][j]]);                    // (k)
            }
        }
        // Optimize model
        model.optimize();
        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << '\n';
        /*for(int i = 0; i < n; ++i)
        {
            for(int j = i+1; j < n; ++j)
            {
                cout << req[i][j] << ' ';
                cout << d[getIdxFlow[i][j]].get(GRB_StringAttr_VarName) << " "
                    << d[getIdxFlow[i][j]].get(GRB_DoubleAttr_X) << '\n';
            }
        }*/
        /*for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                if(i == j) continue;
                if(x[getIdxFlow[i][j]].get(GRB_DoubleAttr_X) > 0.99)
                {
                    printf("(%d, %d, %f)\n", i, j, cost[i][j]);
                }
            }
        }
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                if(i == j) continue;
                printf("%d %d %f\n", i, j, d[getIdxFlow[i][j]].get(GRB_DoubleAttr_X));
            }
        }
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                if(i == j) continue;
                printf("%d %d %f\n", i, j, p[getIdxFlow[i][j]].get(GRB_DoubleAttr_X));
            }
        }*/
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
}

int main()
{
    cout << fixed << setprecision(10);
    //newFormulation();
    //flowFormulation();
    rootedFormulation();
    return 0;
}