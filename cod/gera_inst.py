# This file will create instances for queen graph coloring problem

def create_graph(n:int, filename:str=''):
    '''
        return (V,E) queen graph with n² vertices

        n: size of the board 
    '''

    V = set(range(0,n*n))
    E = {}
    # row edges
    # for k in range(n):
    #     for i in range(k*n, (k+1)*n): 
    #         E[i] = []
    #         for j in range(i+1,n*(k+1)):
    #             e = (i,j)
    #             E[i].append(e)

    # column e
    for k in range(n):
        for i in range(n*k, (k+1)*n):
            E[i] = []    
            for j in range(k*n+1,(k+1)*n):
                E[i].append((i,n*j+i)) 
    return V,E


if __name__ == '__main__':
    print(create_graph(4,'o'))
