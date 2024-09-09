import inst

class Greedy:

    def __init__(self, c:list, E:dict, nrows:int, ncolumns:int):
        self.c = c
        self.E = E
        self.nrows = nrows
        self.ncolumns = ncolumns
        self.Si = self.create_S(E)
        
    def create_S(self,E:dict):
        ''' for a given dict of si and it's feasible sets return feasible sets of S'''
        Si = {}
        for i in E:
            for j in E[i]:
                if j not in Si:
                    Si[j] = [i]
                else:
                    Si[j].append(i)
        return Si
    
    def greedy_pcc(self, available_sets, chosen_s=[]):
        '''
            greedy heuristic for chossing the initial set of the set covering problem
        '''
        nsi = 0
        cp_max = max(self.c)
        Jstar = []
        found=True
        while nsi < self.nrows and found:
            
            cp_min = cp_max
            found = False
            min_set = 0
            for j in available_sets:
                Sp = self.Si[j]
                nSp = 0
                for s in Sp:
                    if s not in chosen_s:
                        nSp +=1

                if nSp > 0:
                    found = True
                    if self.c[j]/nSp < cp_min:
                        cp_min = self.c[j]/nSp
                        min_set = j

            if found:
                Jstar.append(min_set)
                available_sets.remove(min_set)

                for si in self.Si[min_set]:
                    if si not in chosen_s:
                        chosen_s.append(si)
                        nsi+=1

        return Jstar, nsi, chosen_s
    
    def pcc(self, Jk=[]):

        if len(Jk) > 0:
            sets = [i for i in Jk]
            remaining_sets = [i  for i in range(1, self.ncolumns+1) if i not in Jk]
        else:
            sets = [i for i in range(1, self.ncolumns+1)]
        
        Jstar, nsi, chosen_s = self.greedy_pcc(sets,[])
        if nsi == self.nrows:
            return Jstar
        
        nJ, nsremain, cremain = self.greedy_pcc(remaining_sets, chosen_s)
        return Jstar+nJ
        
def main():
    read = inst.Read('instancias/1987/4/scp49.txt')
    nrows,ncolumns, c, E = read.read_inst()
    c = {1:2,2:3,3:4,4:5}
    E = {1:[1,3],2:[1,4],3:[2,3,4]}
    ncolumns = 4
    nrows = 3
    greedy = Greedy(c,E, nrows, ncolumns)
    Jstar = greedy.pcc()


if __name__ == '__main__':
    main()