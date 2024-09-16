import sys 
import cplex
import numpy as np
import time

import os
import pandas as pd

package_directory = '/home/jessicarichards/Documents/mat_dout/ot_comb/codigos/trab_final/cod'
sys.path.append(package_directory)
from pckgs import inst
from pckgs.guloso_pcc_2 import Greedy
from mip_cplex import mip_cplex

class Otm:

    def __init__(self) -> None:

        pass



    def lagrangean(self,c,Ji,nrows,ncolumns, max_iter=20, fix=False, ub_otm=0, ub_relx=0, info=False):

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
        Jmin = Jstar

        ub = 0 
        for i in Jstar: ub+= c[i]

        u = {i:0 for i in Ji}
        new_c = {j:c[j] - sum([u[i] for i in Ij[j] ]) for j in c}

        x = [f'x_{i}' for i in c]
        xtypes = ['C' for _ in c]

        # cplex params   
        log_file = open('cplex_log.txt', 'w') # to save log in a separate file

        xlb = [0 for _ in range(ncolumns)]
        xub = [1 for _ in range(ncolumns)]

        steps = []
        Zmax = 0

        tmp_heur = 0
        tmp_lag = 0
        
        allZ = []
        allUb = []
        fixed=True
        while k < max_iter and (ub-Zmax) > 1 and (pi >  5e-4) and abs(ub-ub_otm)>1e-5 :
            
            uk = 0
            for i in Ji: uk+=u[i]

            s = time.time()
            # Set the display level

            model = cplex.Cplex()   

            model.parameters.mip.display.set(0)
            model.variables.add(names=x, lb=xlb, ub=xub, types=xtypes)

            model.objective.set_sense(model.objective.sense.minimize)

            model.objective.set_linear([(x[i-1], new_c[i]) for i in new_c])
            model.objective.set_offset(uk)

            # Step 7: Solve the model
            model.solve()
        
            model.set_results_stream(log_file)
            model.set_results_stream(log_file)

            solution = model.solution
            zk = solution.get_objective_value()
            # if k > 2:
            #     print(zk)
            #     break
            # else:print(f'k={k},zk={zk}, uk={uk}, {new_c}')
            f = time.time()
            tmp_lag += (f-s)

            # heuristica
            s = time.time()
            xj = {j:1 if new_c[j] < 0 else 0 for j in new_c}
            subgrad = {i:1-sum([xj[j] for j in Ji[i]]) for i in Ji}

            if Zmax < zk:
                
                counter = 0
                Zmax = zk
                # finding new upper bound 

                Jk = []
                for i in xj: 
                    if xj[i]:  Jk.append(i)  

                Jstar = greedy.pcc(Jk)
                ubk = 0 
                for i in Jstar: ubk+= c[i]

                if ubk < ub:
                    ub = ubk
                    Jmin = Jstar

            subT = 0
            for i in subgrad: subT += subgrad[i]**2
            if subT > 0:
                T = pi* (ub-zk) / subT
            
            for i in u: u[i] = max(0,u[i]+(1+eps)*T*subgrad[i])
            for j in new_c: new_c[j] = c[j] - sum([u[i] for i in Ij[j]])

            # fixando variáveis
            fixed=False
            if fix:
                for i in xj: 
                    if (xub[i-1] > 0) and (not xj[i]) and (new_c[i] > 0) and ((zk+new_c[i]) > ub):
                        xub[i-1] = 0
                        greedy.sets.remove(i)
                        fixed = True
                    if (xj[i]) and ((zk-new_c[i]) > ub):
                        xlb[i-1] = 1
                        fixed = True
            f = time.time()
            tmp_heur += (f-s)
            # for var_name in x:
            #     if solution.get_values(var_name):
            #         print(f"{var_name} = {solution.get_values(var_name)}")

            steps.append(T)
            k+=1
            counter+=1
            allZ.append(zk)
            allUb.append(ub)

            if counter  > half:
                pi = pi/2
                counter = 0
        # print(f'steps={steps}')
        log_file.close()
        if not info:
            return Zmax, ub, k, tmp_heur, tmp_lag, pi
    
        return Zmax, ub, k, tmp_heur, tmp_lag, pi, allZ, allUb
# def remove_variable(model, data):
def main2():
    read = inst.Read('trab_final/instancias/1987/4/scp41.txt')
    nrows,ncolumns, c, E = read.read_inst()
    otm = Otm()
    s = time.time()
    # c = {1:2,2:3,3:4,4:5}
    # E = {1:[1,3],2:[1,4],3:[2,3,4]}
    # ncolumns = 4
    # nrows = 3
    zcont = mip_cplex(nrows, ncolumns, c, E,'C')
    zbin = mip_cplex(nrows, ncolumns, c, E,'B')
    Zmax, ub, k, tmp_heur, tmp_lag, pi = otm.lagrangean(c,E,nrows, ncolumns,1e4, True)   
    f = time.time()
    print(f'tmp={f-s}')



def main():

    results = 'results.csv'
    if not os.path.exists(results):
        columns = ['fname','Zmax', 'ub','niter','tmp_heur','tmp_lag','fix_var','pi', 'best_relax','best_mip']
        df = pd.DataFrame(columns=columns)
        df.to_csv(results)
    else:
        df = pd.read_csv(results)
        columns = df.columns
        
    otm = Otm()

    main1987 = 'trab_final/instancias/1987'
    dir1987 = os.listdir(main1987)

    for fix_var in [True, False]:
        for sdir in dir1987[::-1]:
            dirname = f'{main1987}/{sdir}'
            files = os.listdir(dirname)
            if sdir in ['4','5','6','A']: continue
           
            for f in files:
                # if f !='scp41.txt': continue
                fpath = f'{main1987}/{sdir}/{f}'
                read = inst.Read(fpath)
                nrows,ncolumns, c, E = read.read_inst()
                zcont = mip_cplex(nrows, ncolumns, c, E,'C')
                zbin = mip_cplex(nrows, ncolumns, c, E,'B')
                Zmax, ub, k, tmp_heur, tmp_lag, pi = otm.lagrangean(c,E,nrows, ncolumns,1e4, fix_var,zbin)

                data =  {
                            'fname':f,
                            'Zmax':Zmax, 
                            'ub':ub,
                            'niter':k,
                            'tmp_heur':tmp_heur,
                            'tmp_lag':tmp_lag,
                            'fix_var':fix_var,
                            'pi':pi,
                            'best_relax':zcont,
                            'best_mip':zbin
                        }
                # Make data frame of above data
                df = pd.DataFrame([data])
                # append data frame to CSV file
                df.to_csv(results, mode='a', index=False, header=False)
        
if __name__ == '__main__':
    main()