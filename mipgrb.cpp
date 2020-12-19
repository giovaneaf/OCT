#include <bits/stdc++.h>
#include <gurobi_c++.h>

using namespace std;

int main()
{
    
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
