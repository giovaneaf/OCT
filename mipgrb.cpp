#include <bits/stdc++.h>
#include <gurobi_c++.h>

using namespace std;

struct Edge 
{
    int u, v, len, id;
    Edge(int u = 0, int v = 0, int len = 0, int id = 0) : u(u), v(v), len(len), id(id) {}
};

vector<Edge> mEdges;
vector<vector<double>> req;
vector<vector<int>> Nmin;

string getNewConstr()
{
    static int cnt = 0;
    return "C" + to_string(cnt++);
}

int main()
{
    int d;
    int n, m;
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

        // Create binary variables
        GRBVar* x = model.addVars(m, GRB_BINARY);
        GRBVar* y = model.addVars(n*n, GRB_BINARY);
        GRBVar* z = model.addVars(n*n*n, GRB_BINARY);

        // Create integer variables
        GRBVar* delta = model.addVars(n, GRB_INTEGER);
        GRBVar* eta = model.addVars(n, GRB_INTEGER);
        GRBVar* rho = model.addVars(n*n*n, GRB_INTEGER);

        GRBLinExpr obj;
        int idx = 0;
        for(int u = 0; u < n; ++u)
        {
            for(int v = u+1; v < n; ++v)
            {
                idx = u*n*n+v*n;
                obj.addTerms(&req[u][v], &delta[u], 1);
                obj.addTerms(&req[u][v], &delta[v], 1);
                const double tmp = -2*req[u][v];
                for(int w = 0; w < n; ++w)
                {
                    obj.addTerms(&tmp, &rho[idx+w], 1);
                }
            }
        }

        model.setObjective(obj, GRB_MINIMIZE);

        GRBLinExpr linexpr;

        int root = 0;
        const double one = 1.0;
        for(int& inEdge : Nmin[root])
        {
            linexpr.addTerms(&one, &x[inEdge], 1);
        }

        model.addConstr(linexpr == 0, getNewConstr());
        model.addConstr(delta[root] == 0, getNewConstr());
        model.addConstr(eta[root] == 0, getNewConstr());

        for(int u = 1; u < n; ++u)
        {
            linexpr.clear();
            for(int& inEdge : Nmin[u])
            {
                linexpr.addTerms(&one, &x[inEdge], 1);
            }
            model.addConstr(linexpr == 1, getNewConstr());
            model.addConstr(delta[u] <= 0, getNewConstr());
            model.addConstr(eta[u] <= 0, getNewConstr());
        }

        //model.addConstr(delta[0])




        /*cons.addTerms(W, e, sz);

        // Add constraint: x + 2 y + 3 z <= 4
        //model.addConstr(x + 2 * y + 3 * z <= C, "c0");
        model.addConstr(cons, GRB_LESS_EQUAL, C, "c0");

        // Optimize model
        model.optimize();*/

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
