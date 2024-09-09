## "python 'C:\Program Files\IBM\ILOG\CPLEX_Studio2211\python\setup.py' install"

import inst
import sys
print(sys.version)
import cplex

read = inst.Read('ot_comb/instancias/1987/4/scp49.txt')

nrows,ncolumns, cj, E = read.read_inst()

# Step 3: Create a CPLEX model
model = cplex.Cplex()

# Step 4: Define the decision variables
x = [f'x{i}' for i in range(1,ncolumns+1)]
lb = [0 for i in range(ncolumns)]
ub = [1 for i in range(ncolumns)]
xtypes = ['C' for _ in range(ncolumns)]

A = [[f'a{i}{j}' for j in range(1,ncolumns+1)] for i in range(nrows)]
alb = [[0 for _ in range(ncolumns)] for _ in range(nrows)]
aub = [[1 if j in E[i] else 0 for j in range(ncolumns) ] for i in range(1,nrows+1)]
Atypes = [['B' for _ in range(ncolumns)] for _ in range(nrows)]

model.variables.add(names=x, lb=lb, ub=ub, types=xtypes)

for i in range(nrows):
    model.variables.add(names=A[i], lb=alb[i], ub=aub[i], types=Atypes[i])

# Step 5: Set the objective function
model.objective.set_linear(list(zip(x, cj)))
model.objective.set_sense(model.objective.sense.minimize)

# for i in range(nrows):
for i in range(1,nrows+1):
     model.linear_constraints.add(
                                    lin_expr=[cplex.SparsePair(ind=[x[j-1] for j in E[i]], val=[1]*len(E[i]))],
                                    senses=['G'],
                                    rhs=[1],
                                    names=[f'c{i}']
                                )

constraint_names = model.linear_constraints.get_names()

# Step 7: Solve the model
model.solve()

# Step 8: Retrieve and display the results
solution = model.solution
print("Objective value:", solution.get_objective_value())
print("Values of decision variables:")
for var_name in x:
    if solution.get_values(var_name) > 0:
        print(f"{var_name} = {solution.get_values(var_name)}")

