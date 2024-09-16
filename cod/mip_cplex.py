## "python 'C:\Program Files\IBM\ILOG\CPLEX_Studio2211\python\setup.py' install"
import cplex
import sys

package_directory = '/home/jessicarichards/Documents/mat_dout/ot_comb/codigos/trab_final/cod'
sys.path.append(package_directory)

from pckgs import inst


def mip_cplex(nrows, ncolumns, cj, E, vartype='B',show=False):
    # Step 3: Create a CPLEX model
    model = cplex.Cplex()

    # Step 4: Define the decision variables
    x = [f'x{i}' for i in range(1,ncolumns+1)]
    lb = [0 for i in range(ncolumns)]
    ub = [1 for i in range(ncolumns)]
    xtypes = [vartype for _ in range(ncolumns)]

    model.variables.add(names=x, lb=lb, ub=ub, types=xtypes)


    # Step 5: Set the objective function
    model.objective.set_linear([(x[i-1], cj[i]) for i in cj])
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
    if show:
        print("Objective value:", solution.get_objective_value())
        print("Values of decision variables:")
        for var_name in x:
            if solution.get_values(var_name) > 0:
                print(f"{var_name} = {solution.get_values(var_name)}")

    return solution.get_objective_value()

def main():

    read = inst.Read('trab_final/instancias/1987/A/scpa3.txt')
    nrows,ncolumns, cj, E = read.read_inst()
    zcont = mip_cplex(nrows, ncolumns, cj, E,'C')
    zbin = mip_cplex(nrows, ncolumns, cj, E,'B')
    print(f'binary:{zbin}')
    print(f'continous:{zcont}')

if __name__ =='__main__':
    main()