import inst
import cplex
import numpy as np
from cplex import infinity

def ub_linear(c,S):

    model = cplex.Cplex()

    ncolumns = len(S[0])
    nrows = len(S)
    
    # Step 4: Define the decision variables
    x = [f'x{i}' for i in range(0,ncolumns)]
    xlb = [0 for i in range(ncolumns)]
    xlb[-1] =1
    xlb[2]=0.3
    xub = [1 for i in range(ncolumns)]
    xtypes = ['C' for _ in range(ncolumns)]

    model.variables.add(names=x, lb=xlb, ub=xub, types=xtypes)
    model.objective.set_linear(list(zip(x, c)))
    model.objective.set_sense(model.objective.sense.minimize)
   
    for i in range(len(S)):

        model.linear_constraints.add(
                                        lin_expr=[cplex.SparsePair(ind=[x[j] for j in range(len(S[i])) if S[i][j] != 0], val=[1]*sum(S[i]))],
                                        senses=['G'],
                                        rhs=[.8],
                                        names=[f'c{i}']
                                    )
    model.solve()

    solution = model.solution
    xk = np.zeros(len(S[0]))
    for var_name in x:
        i = int(var_name.replace('x',''))
        xk[i] = solution.get_values(var_name)
        # if solution.get_values(var_name) > 0:
        print(f"{var_name} = {solution.get_values(var_name)}")
    return solution.get_objective_value(), xk

def subgradient(c,E, nrows, ncolumns, max_iter):

    I = {}
    for i in E:
        for j in E[i]:
            if j not in I:I[j] = []
            I[j].append(i)

    uk = np.zeros(nrows)
    subgradk = np.zeros(nrows)
   
    eps = 0.02
    pi = 2
    k=0

    # Step 3: Create a CPLEX model
    model = cplex.Cplex()
    
    # Step 4: Define the decision variables
    x = [f'x{i}' for i in range(1,ncolumns+1)]
    xtypes = ['C' for _ in range(ncolumns)]
    ub = sum(c)

    model.objective.set_sense(model.objective.sense.minimize)
    while k < max_iter:

        cl = [c[j] - sum([uk[i-1] for i in I[j+1]]) for j in range(ncolumns)]
        xlb = [0 if cl[j] > 0 else 1 for j in range(ncolumns)]
        xub = [1 for _ in range(ncolumns)]
        
        model.variables.add(names=x, lb=xlb, ub=xub, types=xtypes)
        model.objective.set_linear(list(zip(x, cl)))
        model.objective.set_offset(-sum(uk))
        
        # Step 7: Solve the model
        model.solve()

        solution = model.solution
        zk = solution.get_objective_value()
        print("Values of decision variables:")
        xk = np.zeros(ncolumns)
        for var_name in x:
            i = int(var_name.replace('x',''))
            xk[i-1] = solution.get_values(var_name)
            # if solution.get_values(var_name) > 0:
            # print(f"{var_name} = {solution.get_values(var_name)}")
        print(f'objective value = {zk}')

        for i in range(nrows):
            subgradk[i] = 1-sum([xk[j-1] for j in E[i+1]])
        
        T = pi * (ub - zk) / sum(subgradk[i]**2 for i in range(nrows))

        for i in range(nrows):
            uk[i] = max(0,uk[i] + (1+eps)*T*subgradk[i])
        k += 1

def main():
    read = inst.Read('ot_comb/instancias/1987/4/scp49.txt')
    nrows,ncolumns, cj, E = read.read_inst()
    subgradient(cj,E, nrows, ncolumns, 10)



if __name__ == '__main__':
    main()