
import sys

class BinaryTree:

    def __init__(self, val, dist=0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right
        
    def getCluster(self):
        res = []
        if self.value >= 0:
            res.append(self.value)
        else:
            if self.left != None:
                res.extend(self.left.getCluster())
            if self.right!= None:
                res.extend(self.right.getCluster())
        return res
    
    def printtree(self):
        self.printtreerec(0, "Root")
    
    def printtreerec (self, level, side):
        for i in range(level): sys.stdout.write("\t")
        al = self.getCluster();
        sys.stdout.write(side + ":" + str(al)+ " Dist.: " + str(self.distance) + "\n")
        if self.value < 0:
            if (self.left != None): 
                self.left.printtreerec(level+1, "Left")
            else: 
                sys.stdout.write("Null")
            if (self.right != None): 
                self.right.printtreerec(level+1, "Right")
            else: 
                sys.stdout.write("Null\n")


def test():              
    a=BinaryTree(2)
    b=BinaryTree(3)
    c=BinaryTree(4)
    d=BinaryTree(1)
   
    e=BinaryTree(-1,2.0,a,b)
    f=BinaryTree(-1,1.5,c,d)
    g=BinaryTree(-1,4.5,e,f)
    g.printtree()
    
if __name__ == '__main__':
    test()
    
# Result for test
# Root:[1, 2, 3, 4] Dist.: 4.0
#     Left:[1, 2] Dist.: 2.0
#         Left:[1] Dist.: 0
#         Right:[2] Dist.: 0
#     Right:[3, 4] Dist.: 3.0
#         Left:[3] Dist.: 0
#         Right:[4] Dist.: 0

