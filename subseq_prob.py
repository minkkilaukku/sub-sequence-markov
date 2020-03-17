def getProb(symbols, seqs, n):
    Q = getMat(symbols, seqs)
    M = Q.dimensions()[0]
    return (Q**n)[0][M-1]



def getMat(symbols, seqs):

    def cleanSeqs(arr):
        ret = list(arr)
        changed = True
        while changed:
            changed = False
            for x in ret:
                hadADelete = False
                for y in ret:
                    if x!=y and y in x:
                        ret.remove(x)
                        hadADelete = True
                        changed = True
                        break
                if hadADelete:
                    break
        return ret
                    
    
    seqs = cleanSeqs(seqs)

    symbolsLen = len(symbols)
    prefixes = set([""]) #first state is '' empty start state

    for seq in seqs:
        for i in range(1, len(seq)):
            prefixes.add(seq[0:i])

    prefixes = sorted(list(prefixes))
    M = len(prefixes)+1 #last state (index M-1) is the absorbing one

    def findFallsBackTo(appendedRun):
        for i in range(0,len(appendedRun)):
            if appendedRun[i:] in prefixes:
                return appendedRun[i:]
        return ""

    def incrState(stateDict, state):
        if state in stateDict:
            stateDict[state] += 1
        else:
            stateDict[state] = 1

    def someEndInSeqs(s):
        for i in range(len(s)):
            if s[i:] in seqs: return True
        return False

    def getToStateVals(stateStr):
        ret = {}
        for i in symbols:
            stPlusI = stateStr+str(i)
            if someEndInSeqs(stPlusI):
                incrState(ret, M-1)
            elif stPlusI in prefixes:
                incrState( ret, prefixes.index(stPlusI) )
            else:
                fallsToPrefix = findFallsBackTo(stPlusI)
                incrState(ret, prefixes.index(fallsToPrefix))
        return ret    


    a = []
    aInts = []
    for i in range(M):
        a += [[]]
        aInts += [[]]
        for j in range(M):
            a[i] += [0]
            aInts[i] += [0]

    for i in range(M-1):
        ps = getToStateVals(prefixes[i])
        for toI in ps:
            a[toI][i] = ps[toI]/symbolsLen
            aInts[toI][i] = ps[toI]

    a[M-1][M-1] = 1
    aInts[M-1][M-1] = symbolsLen
    return Matrix(CDF, a).T


symbs = list(str(c) for c in range(1,7))
finals = ['123456',
'112233',
'445566',
'111222',
'333444',
'555666',
'111111',
'222222',
'333333',
'444444',
'555555',
'666666']


Q = getMat(symbs, finals)
dimQ = Q.dimensions()[0]
D, P = Q.eigenmatrix_right()


#check how the eigenvectors go
#for k in range(dimQ):
#    print (Q*P.T[k]-D[k][k]*P.T[k])
    
#print ("\n")

#put columns (i.e. rows of P.T)
eigs = sorted([(D[k][k], P.T[k]) for k in range(D.dimensions()[0])],
             key = lambda x: abs(x[0]), reverse=True)

#reordering w.r.t order of eigvals
D = matrix([[0 if j!=k else eigs[k][0] for j in range(dimQ)] for k in range(dimQ)])
P = matrix([eigs[k][1] for k in range(dimQ)])


lastBaseV = vector([0 if k<dimQ-1 else 1 for k in range(dimQ)])
invLastCol = P.T.solve_right(lastBaseV)
pFirstRow = P.T[0]
aCoeffs = [pFirstRow[k]*invLastCol[k] for k in range(dimQ)]

tarkkuus = 3
approxFormula = lambda n: sum(aCoeffs[k]*D[k][k]**n for k in range(tarkkuus)).real_part()
print ("Approximate formula")
formulaTerms = []
for k in range(tarkkuus):
    formulaTerms.append("(%s)*(%s)^n" %(str(aCoeffs[k]), D[k][k]))
print ("\n+ ".join(formulaTerms))

print ("\nValues approx formula and true value")
for j in range(30):
    print ("%d: \t %s \t %f" %(j, approxFormula(j), getProb(symbs, finals, j)))



#print ("last column of inverse of P")
#print (invLastRow)
#print ("eigs:")
#print ("\n".join(str(x) for x in eigs[:3]))

#check the eigen vectors went right, remember now the rows are the eig.vecs
#for k in range(dimQ):
#    print (Q*P[k]-D[k][k]*P[k])

    
#print (", ".join(str(D[k][k]) for k in range(D.dimensions()[0])))
#show (Q)
#show ([P.T, D, P.T.inverse()])
#print (P.T*D - Q*P.T)


solFStr = """
1 -1.0011546714408708*0.9997646566235854^n
- (0.012616475011663993 + 0.020083245284924466*I)*(0.07445021698430404 - 0.18248047052695532*I)^n
- (0.012616475011663509 - 0.020083245284924618*I)*(0.07445021698430393 + 0.18248047052695507*I)^n
"""
