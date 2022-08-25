import math
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
 
def epsilon(m):
  A=[m+1] #A^{(m+1)}
  number=1 #number of S strings encountered
  position=0 #position of end of last S
  value=mpf(1) #approximation for epsilonm
#the value 2000/m is chosen to get enough accurracy for at least 20 decimals
  for i in range(1,int(2000/m)):
      k=krul(A)
      A.append(max(m+1,k))
      if k<m+1:
 #         print(number,i-position)
          value=value+mpf(i-position)/mpf((m+1)**number)
          position=mpf(i)
          number+=1
  print('epsilon_',m,': ',value,sep='')


def Glue_strings(number):
  print('First',number,'glue strings:')
  A2=[2]
  S=[2]
  i=0
  while(True):
      c=krul(A2)
      if c==1:
        print(S)
        i=i+1
        if i==number:
            break
        S=[]
      A2.append(max(c,2))
      S.append(max(c,2))

def term(m,i):
    return (m+2)*((2*i-1)**(math.log(m+2)/math.log(m+1)))/((m+1)**i)

def Sum(m,t):
    result=0
    for i in range((m+1)*(t-1),(m+1)*(t-1)+100):
        result += term(m,i)
    return result

def upper_bound(m,t):
    return 3.5*((m+1)**(t-1)-1)*Sum(m,t)

for m in range(1,11):
    epsilon(m)

# Uncomment the following line to print the first 100 first order glue strings
#Glue_strings(100)

# Uncomment the following lines to verify the estimates made in Proposition 7.9
#print(upper_bound(1,80))
#print(upper_bound(2,80))
#print(upper_bound(3,5))