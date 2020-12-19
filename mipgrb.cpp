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

vector<Edge> mEdges;
vector<vector<double>> req;
vector<vector<int>> Nmin;

// Vectors used for gurobi
vector<vector<double>> reqMin2;

string getNewConstr()
{
    static int cnt = 0;
    return "C" + to_string(cnt++);
}

int n, m;

inline int get(int u, int v)
{
    return u*n+v;
}
inline int get(int u, int v, int w)
{
    return u*n*n+v*n+w;
}


int main()
{
    int d = 0;
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
            reqMin2.resize(n, vector<double>(n));
            // Computed needed values for formulation
            for(int u = 0; u < n; ++u)
            {
                for(int v = u+1; v < n; ++v)
                {
                    reqMin2[u][v] = -2*req[u][v];
                }
            }
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
                    obj.addTerms(&reqMin2[u][v], &rho[get(w,u,v)], 1);
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
            model.addConstr(y[get(e.u, e.v)] >= x[e.id], getNewConstr());                                     // (14)
        }

        for(int v = 0; v < n; ++v)
        {
            model.addConstr(y[get(v, v)] == 1, getNewConstr());   //(13)
            for(int u = 0; u < n; ++u)
            {
                linexpr.addTerms(&one, &y[get(u, v)], 1);
            }
            model.addConstr(linexpr == eta[v]+1);               //(12)
            linexpr.clear();
        }

        int uv, uw, vw, uvw, wuv, wu, wv;
        for(int u = 0; u < n; ++u)
        {
            for(int v = 0; v < n; ++v)
            {
                if(u == v) continue;
                uv = get(u, v);
                for(int w = 0; w < n; ++w)
                {
                    uw = get(u, w);
                    vw = get(v, w);
                    uvw = get(u, v, w);
                    wuv = get(w, u, v);
                    wu = get(w, u);
                    wv = get(w, v);
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
        }

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << '\n';*/


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
