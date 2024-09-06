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
    print(I)
    uk = np.zeros(nrows)
    subgradk = np.zeros(nrows)
   
    eps = 0.02
    pi = 2
    k=0

    # Step 3: Create a CPLEX model
    model = cplex.Cplex()
    
    # Step 4: Define the decision variables
    x = [f'x{i}' for i in range(0,ncolumns)]
    xtypes = ['C' for _ in range(ncolumns)]


    model.objective.set_sense(model.objective.sense.minimize)
    while k < max_iter:

        cl = [c[j] - sum([uk[i] for i in I[j]]) for j in range(ncolumns)]
        xlb = [0 if cl[j] > 0 else 1 for j in range(ncolumns)]
        xub = [1 for j in range(ncolumns)]
    
        model.variables.add(names=x, lb=xlb, ub=xub, types=xtypes)
        model.objective.set_linear(list(zip(x, cl)))
        model.objective.set_offset(-sum(uk))
        constraint_names = model.linear_constraints.get_names()

        # Step 7: Solve the model
        model.solve()

        solution = model.solution
        zk = solution.get_objective_value()
        print("Values of decision variables:")
        xk=[]
        for var_name in x:
            i = int(var_name.replace('x',''))
            xk[i] = solution.get_values(var_name)
            # if solution.get_values(var_name) > 0:
            print(f"{var_name} = {solution.get_values(var_name)}")
        print(f'objective value = {zk}')

        for i in range(nrows):
            subgradk[i] = 1-sum([xk[j] for j in E[i]])
            print(subgradk[i])
        
        T = pi * (ub - zk) / sum(subgradk[i]**2 for i in range(nrows))

        for i in range(nrows):
            uk[i] = max(0,uk[i] + (1+eps)*T*subgradk[i])
        k += 1
        print(uk)

def subgradient2(c,S,max_iter):

    ncolumns = len(S[0])
    nrows = len(S)

    Ji = []
    for i in range(len(S)):
        J = []
        for j in range(len(S[i])):
            if S[i][j] > 0: J.append(j) 
        Ji.append(J)


    Ij = []
    for j in range(ncolumns):
        I = []
        for i in range(len(S)):
            if S[i][j] > 0: I.append(i) 
        Ij.append(I)

    # lagrangean multiplier 
    uk = np.zeros(nrows)
    subgradk = np.zeros(nrows)
   
    eps = 0.02
    pi = 2
    k=0
    x = np.zeros(ncolumns)
    ub, xk = ub_linear(c,S)
    
    for i in range(nrows):
        subgradk[i] = 1-sum([xk[j] for j in Ji[i]])
        
    T = pi * (sum(c) - ub) / sum(subgradk[i]**2 for i in range(nrows))

    for i in range(nrows):
        uk[i] = max(0,uk[i] + (1+eps)*T*subgradk[i])

    # Step 3: Create a CPLEX model
    model = cplex.Cplex()

    ncolumns = len(S[0])
    
    # Step 4: Define the decision variables
    x = [f'x{i}' for i in range(0,ncolumns)]
    xtypes = ['C' for _ in range(ncolumns)]
    model.objective.set_sense(model.objective.sense.minimize)
   
    while k < max_iter:

        cl = [c[j] - sum([uk[i] for i in Ij[j]]) for j in range(ncolumns)]
        xlb = [0 if cl[j] > 0 else 1 for j in range(ncolumns)]
        xub = [1 for j in range(ncolumns)]
 
        model.variables.add(names=x, lb=xlb, ub=xub, types=xtypes)
        model.objective.set_linear(list(zip(x, cl)))
        model.objective.set_offset(-sum(uk))

        constraint_names = model.linear_constraints.get_names()

        # Step 7: Solve the model
        model.solve()

        solution = model.solution
        zk = solution.get_objective_value()
        print("Values of decision variables:")
        for var_name in x:
            i = int(var_name.replace('x',''))
            xk[i] = solution.get_values(var_name)
            # if solution.get_values(var_name) > 0:
            print(f"{var_name} = {solution.get_values(var_name)}")
        print(f'objective value = {zk}')

        for i in range(nrows):
            subgradk[i] = 1-sum([xk[j] for j in Ji[i]])
            print(subgradk[i])
        
        T = pi * (ub - zk) / sum(subgradk[i]**2 for i in range(nrows))

        for i in range(nrows):
            uk[i] = max(0,uk[i] + (1+eps)*T*subgradk[i])
        k += 1
        print(uk)
    # print(list(zip(x, cj)))
    # # Step 5: Set the objective function
    
    # 
    # # for i in range(nrows):
    # for i in range(len(S)):

    #     model.linear_constraints.add(
    #                                     lin_expr=[cplex.SparsePair(ind=[x[j] for j in range(len(S[i])) if S[i][j] != 0], val=[1]*sum(S[i]))],
    #                                     senses=['G'],
    #                                     rhs=[1],
    #                                     names=[f'c{i}']
    #                                 )

    # constraint_names = model.linear_constraints.get_names()

    # # Step 7: Solve the model
    # model.solve()

    # # Step 8: Retrieve and display the results
    # solution = model.solution
    # print("Objective value:", solution.get_objective_value())
    # print("Values of decision variables:")
    # for var_name in x:
    #     if solution.get_values(var_name) > 0:
    #         print(f"{var_name} = {solution.get_values(var_name)}")


def main():
    read = inst.Read('instancias/1987/4/scp49.txt')
    nrows,ncolumns, cj, E = read.read_inst()
    subgradient(cj,E, nrows, ncolumns, 2)



if __name__ == '__main__':
    main()