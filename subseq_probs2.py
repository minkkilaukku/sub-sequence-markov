class Node:
    def __init__(self):
        self.children = {}
        self.word = ""
        self.index = None
        self.isWord = False

    def addChild(self, symb):
        toAdd = Node()
        toAdd.word = self.word + symb
        self.children[symb] = toAdd
        return toAdd
        
    def giveIndices(self, currInd=0):
        self.index = currInd
        i = currInd+1
        for c in self.children:
            i = self.children[c].giveIndices(i)
        return i

    def findNodeWith(self, s):
        if len(s)==0: return self
        if s[0] in self.children:
            return self.children[s[0]].findNodeWith(s[1:])
        return None

    def toString(self, pad=''):
        ret = ""
        for c in self.children:
            ender = '_' if self.children[c].isWord else ''
            ret += pad+'|'+c+ender+"\n"+self.children[c].toString(pad+' ')
        return ret


def getAllNodes(tree):
    yield tree
    for c in tree.children:
        for cNode in getAllNodes(tree.children[c]):
            yield cNode
    

def makeTreeOfRuns(runs):
    root = Node()
    for r in runs:
        for i, c in enumerate(r):
            p = root.findNodeWith(r[:i+1])
            endChild = p
            if p == None:
                p = root.findNodeWith(r[:i])
                endChild = p.addChild(c)
            if i==len(r)-1:
                endChild.isWord = True
    return root


def getMat(symbs, seqs):
    """ make the transition matrix for the subseqs """
    S = len(symbs)
    t = makeTreeOfRuns(seqs)
    n = t.giveIndices()
    a = [[0]*n for _ in range(n)]
    enderInds = []
    for node in getAllNodes(t):
        i1 = node.index
        if node.isWord:
            a[i1][i1] = 1
            enderInds.append(i1)
            continue
        for symb in symbs:
            s = node.word + symb
            node2 = None
            while not node2: #at least the root will be found with s==''
                node2 = t.findNodeWith(s)
                if not node2:
                    s = s[1:]
            a[i1][node2.index] += 1/S
    A = Matrix(QQ, a)
    
    #remove the final states
    nonFinInds = [x.index for x in getAllNodes(t) if not x.isWord]
    Q = A[nonFinInds, nonFinInds]
    dimQ = Q.dimensions()[0]
    #and add a single one final state to the end
    b = [list(Q[i])+[1-sum(Q[i])] for i in range(dimQ)] + [[0]*dimQ+[1]]
    
    return Matrix(CDF, b)

def getProb(symbols, seqs, n):
    Q = getMat(symbols, seqs)
    return (Q**n)[0][-1]



def getApproxFormula(symbs, seqs):
    Q = getMat(symbs, seqs)
    dimQ = Q.dimensions()[0]
    D, P = Q.eigenmatrix_right()
    
    #eigen vectors are the columns of P i.e rows of P.T
    eigs = sorted([(D[k][k], P.T[k]) for k in range(D.dimensions()[0])],
             key = lambda x: abs(x[0]), reverse=True)

    #reordering w.r.t magnitude of eigvals
    D = matrix([[0 if j!=k else eigs[k][0] for j in range(dimQ)] for k in range(dimQ)])
    P = matrix([eigs[k][1] for k in range(dimQ)])
    
    lastBaseV = vector([0 if k<dimQ-1 else 1 for k in range(dimQ)])
    invLastCol = P.T.solve_right(lastBaseV)
    pFirstRow = P.T[0]
    aCoeffs = [pFirstRow[k]*invLastCol[k] for k in range(dimQ)]
    eigVals = [eigs[k][0] for k in range(dimQ)]
    
    return (aCoeffs, eigVals)


runs = ["123456", "112233", "445566", "111222", "333444", "555666",
        "111111", "222222", "333333", "444444", "555555", "666666"]

symbols = set(c for c in "".join(runs))

aCoeffs, eigVals = getApproxFormula(symbols, runs)
tarkkuus = 4
approxFormula = lambda n: sum(aCoeffs[k]*eigVals[k]**n for k in range(tarkkuus)).real_part()

print ("Approximate formula")
formulaTerms = []
for k in range(tarkkuus):
    formulaTerms.append("(%s)*(%s)^n" %(str(aCoeffs[k]), eigVals[k]))
print ("\n+ ".join(formulaTerms))

print ("\nValues approx formula and true value")
for j in range(30):
    print ("%d: \t %s \t %f" %(j, approxFormula(j), getProb(symbols, runs, j)))


#Q = getMat(symbols, runs)
#print ((Q**20)[0][-1])
