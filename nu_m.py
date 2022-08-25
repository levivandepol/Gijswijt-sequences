from mpmath import *
#warning: 0.1 and mpf(1/10) are incorrect; use mpf(1)/mpf(10) instead
mp.dps=30


def krul(L):
    l=len(L)
    candidate=1
    period=0
    for x in range(1,int((l+2)/2)):
        M=L[l-x:l]
        N=L[l-2*x:l-x]
        i=1
        while(M==N and i*x<l):
            i=i+1
            N=L[l-(i+1)*x:l-i*x]
        if i>candidate:
            candidate=i
    return candidate

def V_and_nu(maximum):   
  V_strings=[]
  nu_values=[]
  for m in range(1,maximum+1):
  
    A=[] #A^{(m+2)}
    number=0 #number of Ttjes encountered
    set=[] #V_{m+1}
    for i in range(500):
        k=krul(A)
        A.append(max(m+2,k))
        if k<m+1:
            number+=1
            set.append(number-1)
        if k==m+1:
            number+=1
    value=mpf(0)
    for s in set:
        value=value+mpf(1)/mpf((m+1)**s)
    value=value*(mpf(m)/mpf(m+1))
    nu_values.append(value)
    V_strings.append(set[0:20])  
  for m in range(2,maximum+1):
        print('nu_',m,' : ',value,sep='')
  #it turns out that the algorithm that calculates the V strings
  #has output R when m=1 is inserted.
  #however, it does not output nu1
  print("R : ",V_strings[0])
  for m in range(2,maximum+1):
    print("V_",m," : ",V_strings[m-1],sep='')

def test(one,two):#whenever we perform this test, `one' is a prefix of `two'.
    l=len(one)
    m=len(two)
    return(one[2*l-m:l]==two[l:m])

def QR(x):
  x=79       # we will calculate iota1^{-1}(2^x)  
  A=[]       # A^{(3)}
  number=0   # number of T strings encountered
  Q_x=[]     # V_{2}= set of y\leq x-1 such that krul(T_{y+1}^{(2)})=1
  V=[]       # set of T_y^{(2)} for y=1,2,...,x
  while(number <= x):
    k=krul(A)
    A.append(max(3,k))
    if k<2:
        number+=1
        Q_x.append(number-1)
        del A[-1]
        V.append(A.copy())
        A.append(3)
    if k==2:
        number+=1
        del A[-1]
        V.append(A.copy())
        A.append(3)
  R_x=[] 
  for i in Q_x:
      if i==x:
          break
      j=i+1
      while(j<x):
          if not test(V[i],V[j]):
              R_x.append([i,j])
          j=j+1
  result=1+len(Q_x)
  
# uncomment the following line to print Q_79
#  print("Q_",x,": ",Q_x,sep="")
#  uncomment the following line to print R_79
#  print('R_',x,': ',R_x,sep="")
  nu1=mpf(0)
  for [a,b] in R_x:
      result=result+2**(x-b-1)
  print('nu1 lower bound:',mpf(result-81)/mpf(2**79))
  print('nu1 upper bound:',mpf(result-1)/mpf(2**79))
  print('value of iota1^-1(2^',x,'): ', result,sep='')


V_and_nu(10)
QR(79)

