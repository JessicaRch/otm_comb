import os
import numpy as np
from math import inf

class Read:

    def __init__(self, filename):
        self.filename = filename
    
    def read_inst(self):
        '''
            return empty if file does not exist and instaces if it does
        '''

        try:
            with open(self.filename,'r') as f:
                read = f.read()
                try:
                    return self.str2list(read)
                except:
                    return self.realworld(read)
        except:
            print(f'File {self.filename} does not exists')
            return []
        
    def realworld(self,txt):
        '''
            return number of rows, number of columns, costs, edges
            read instances of real world application 
        '''
        lines = txt.split('\n')
        nrows,ncolumns = self.strlist2numlist(lines[0])
        cj = np.zeros(ncolumns)
        j,i=0,0
        E = {}
        while j < ncolumns: # read costs
            j+=1
            info = self.strlist2numlist(lines[j])
            cj[j-1] = info[0]
            E[j] = np.zeros(info[1], dtype=int)
            for k in range(info[1]):
                E[j][k] = info[k+2]
           

        # read Edges
        return nrows, ncolumns, cj, E

    
    def str2list(self, txt):
        '''
            retun instances
            given a text return a list of instances were the first
        '''
        lines = txt.split('\n')
        nrows,ncolumns = self.strlist2numlist(lines[0])
        i = 0
        j=0
        cj = {}
        while i < ncolumns:
            j+=1
            for num in self.strlist2numlist(lines[j]):
                cj[i+1] = num
                i+=1

        E = {}
        for i in range(1,nrows+1):
            j += 1
            l = self.strlist2numlist(lines[j])
            E[i] = np.zeros(l[0],dtype=int)     
            k = 0
            while k < l[0]:
                j += 1
                for num in self.strlist2numlist(lines[j]):
                    E[i][k] = num
                    k+=1 
        return nrows, ncolumns, cj, E

    def strlist2numlist(self, array):
        ''' given a array of strings eliminate empty values and transform to num'''
        return [self.str2num(a)  for a in array.split(' ') if a != '']

    def str2num(self, num):
        '''given a str return its float or int values'''
        try:
            return int(num)
        except:
            try:
                return float(num)
            except:
                print(f'The given string {num} is not a number')
                return -inf
            

if __name__ == '__main__':
    filename = 'instancias/real_world/rail507.txt'
    read = Read(filename)
    print(read.read_inst())