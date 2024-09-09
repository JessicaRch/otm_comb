from guloso_pcc import Greedy
import inst
import cplex
import numpy as np

def lagrangean(c,Ji,nrows,ncolumns, max_iter=20):

    # initial params 
    pi = 2
    eps = 0.02
    k = 0
    zk = 0
    half = np.floor(max_iter/20)
    counter = 0

    # initialization 
    greedy = Greedy(c,Ji, nrows, ncolumns)
    Ij = greedy.Si
    Jstar = greedy.pcc()

    ub = 0 
    for i in Jstar: ub+= c[i]

    u = {i:0 for i in Ji}
    new_c = {j:c[j] - sum([u[i] for i in Ij[j] ]) for j in c}

    x = [f'x{i}' for i in c]
    xtypes = ['C' for _ in c]

    # cplex params
    model = cplex.Cplex()
    log_file = open('cplex_log.txt', 'w') # to save log in a separate file

    model.objective.set_sense(model.objective.sense.minimize)
    xlb = [0 for _ in range(ncolumns)]
    xub = [1 for _ in range(ncolumns)]
    model.variables.add(names=x, lb=xlb, ub=xub, types=xtypes)

    steps = []
    Zmax = 0
    while k < max_iter and (ub-Zmax) > 1e-3 and (pi >  5e-4):
        uk = 0
        for i in Ji: uk+=u[i]

        model.objective.set_linear([(x[i-1], new_c[i]) for i in new_c])
        model.objective.set_offset(uk)

        # Step 7: Solve the model
        model.solve()
        
        model.set_results_stream(log_file)

        solution = model.solution
        zk = solution.get_objective_value()

        xj = {j:1 if new_c[j] < 0 else 0 for j in new_c}
        subgrad = {i:1-sum([xj[j] for j in Ji[i]]) for i in Ji}

        # finding new upper bound 

        Jk = []
        for i in xj: 
            if xj[i]:  Jk.append(i)  
        Jstar = greedy.pcc(Jk)

        ub = 0 
        for i in Jstar: ub+= c[i]
        
        subT = 0
        for i in subgrad: subT += subgrad[i]**2
        if subT > 0:
            T = pi* (ub-zk) / subT

        for i in u: u[i] = max(0,u[i]+(1+eps)*T*subgrad[i])
        for j in new_c: new_c[j] = c[j] - sum([u[i] for i in Ij[j]])

        # fixando variáveis
        for i in xj: 
            if (not xj[i]) and (new_c[i] > 0) and ((zk+new_c[i]) > ub):
                xub[i-1] = 0
            if (xj[i]) and ((zk-new_c[i]) > ub):
                xlb[i-1]=1
        steps.append(T)
        k+=1
        counter+=1
        
        if Zmax < zk:
            counter = 0
            Zmax = zk
            xk = np.zeros(ncolumns)
            for var_name in x:
                i = int(var_name.replace('x',''))
                xk[i-1] = solution.get_values(var_name)
                # if solution.get_values(var_name):
                # print(f"{var_name} = {solution.get_values(var_name)}")

        if counter  > half:
            pi = pi/2
            counter = 0
    z = 0 
    for j in xj: 
        if xj[j]: z+=c[j]   

    print(f'steps={steps}')
    print(f'pi={pi},Zmax={Zmax}, z={zk}, ub={ub}, xlb={xlb} \n')
    log_file.close()

def main():
    read = inst.Read('instancias/1987/4/scp41.txt')
    nrows,ncolumns, c, E = read.read_inst()
    # c = {1:2,2:3,3:4,4:5}
    # E = {1:[1,3],2:[1,4],3:[2,3,4]}
    # ncolumns = 4
    # nrows = 3
    lagrangean(c,E,nrows, ncolumns,500)


if __name__ == '__main__':
    main()