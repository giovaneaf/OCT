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

int main()
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
    try 
    {
         // Create an environment
        GRBEnv env = GRBEnv(true);
        //env.set("LogFile", "mip.log");
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
    return 0;
}

// Knapsack example (Gurobi)
/*
int main()
{
    try 
    {
         // Create an environment
        GRBEnv env = GRBEnv(true);
        //env.set("LogFile", "mip.log");
        env.start();

         // Create an empty model
        GRBModel model = GRBModel(env);

        double W[] = {1, 2, 3};
        double V[] = {2, 3, 4};
        double C = 4;
        int sz = sizeof(W)/sizeof(double);

        // Create variables
        GRBVar* e = model.addVars(sz, GRB_BINARY);

        GRBLinExpr obj;

        obj.addTerms(V, e, sz);

        // Set objective: maximize x + y + 2 z
        model.setObjective(obj, GRB_MAXIMIZE);

        GRBLinExpr cons;

        cons.addTerms(W, e, sz);

        // Add constraint: x + 2 y + 3 z <= 4
        //model.addConstr(x + 2 * y + 3 * z <= C, "c0");
        model.addConstr(cons, GRB_LESS_EQUAL, C, "c0");

        // Optimize model
        model.optimize();

        for(int i = 0; i < sz; ++i)
        {
            cout << e[i].get(GRB_StringAttr_VarName) << " "
                << e[i].get(GRB_DoubleAttr_X) << endl;
        }

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
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
    return 0;
}
*/
