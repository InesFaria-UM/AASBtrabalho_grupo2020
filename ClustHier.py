    
from BinaryTree import BinaryTree
from MatrixNum import MatrixNum

class ClustHier:

    def __init__(self, matdists):
        self.matdists = matdists
        
    def distance (self, tree1, tree2):
        c1=tree1.getCluster()
        c2=tree2.getCluster()
        sd=0.0
        for i in range(len(c1)):
            for j in range(len(c2)):
                sd+= self.matdists.getValue(c1[i],c2[j])
        return sd/len(c1)*len(c2)
    
    def executeClustering(self):
        trees = []
        tableDist = self.matdists.copy()
        for i in range(self.matdists.numRows()):
            t = BinaryTree(i)
            trees.append(t) # adds leaves to the list
        for k in range(self.matdists.numRows(), 1, -1):
            # choose minimum value from matrix not in diagonal
            mins = tableDist.minDistIndexes()
            i = mins[0]
            j = mins[1]
            n = BinaryTree(-1, tableDist.getValue(i, j)/2.0, trees[i], trees[j])
            if k>2:
                trees.pop(i)
                trees.pop(j)
                tableDist.removeRow(i)
                tableDist.removeRow(j)
                tableDist.removeCol(i)
                tableDist.removeCol(j)
                dists = []
                for x in range(len(trees)):
                    dists.append(self.distance(n, trees[x]))
                tableDist.addRow(dists)
                cdists = []
                for y in range(len(dists)):
                    cdists.append(dists[y])
                cdists.append(0.0)
                tableDist.addCol(cdists)
                trees.append(n)
            else: return n
            tableDist.printmat()

def test():
    m = MatrixNum(5,5)
    m.setValue(1, 0, 2)
    m.setValue(0, 1, 2)
    m.setValue(2, 0, 5)
    m.setValue(0, 2, 5)
    m.setValue(3, 0, 7)
    m.setValue(0, 3, 7)
    m.setValue(4, 0, 9)
    m.setValue(0, 4, 9)
    m.setValue(2, 1, 4)
    m.setValue(1, 2, 4)
    m.setValue(3, 1, 6)
    m.setValue(1, 3, 6)
    m.setValue(4, 1, 7)
    m.setValue(1, 4, 7)
    m.setValue(3, 2, 4)
    m.setValue(2, 3, 4)
    m.setValue(4, 2, 6)
    m.setValue(2, 4, 6)
    m.setValue(4, 3, 3)
    m.setValue(3, 4, 3)
    hc = ClustHier(m)
    arv = hc.executeClustering()
    arv.printtree()
    
if __name__ == '__main__': 
    test()
    
# Result for test
# Root:[1, 0, 2, 4, 3] Dist.: 3.25
#     Left:[1, 0, 2] Dist.: 2.25
#         Left:[1, 0] Dist.: 1.0
#             Left:[1] Dist.: 0
#             Right:[0] Dist.: 0
#         Right:[2] Dist.: 0
#     Right:[4, 3] Dist.: 1.5
#         Left:[4] Dist.: 0
#         Right:[3] Dist.: 0
