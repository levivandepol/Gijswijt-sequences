import math

# Replace the value of Limit by 513 to obtain the value of beta_{356}^{(1)}.
Limit = 80 # number of entries that we will print

def ord(a,m): # the number of factors m in a
    e = 0
    while (a % (m**(e+1)) == 0):
        e += 1
    return e

def numberToBase(n, b): # https://stackoverflow.com/questions/2267362/how-to-convert-an-integer-to-a-string-in-any-base
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits

def first_nonzero_bits(expansion): # returns the location of the first two nonzero bits of expansion
    length = len(expansion)        # if they exist
    t = 0
    while (expansion[t]==0):
        t += 1
    u = t+1
    while (u<length and expansion[u] == 0):
        u += 1
    if u < length:
        return (t+1,u+1)
    else:
        return [t+1]
    
def check_S(m,a,b): # check whether (a,b) is in S_m
    expansion_a = numberToBase(a,m)
    expansion_b = numberToBase(b,m)
    B = len(expansion_b)
    A = len(expansion_a)
    for i in range(B-A):
        expansion_a.append(0)
    i = 0
    while(i<B and expansion_a[i] == expansion_b[i]):
        i += 1
    if i == B:
        return True
    if 2 * expansion_a[i] < expansion_b[i]:
        return False
    i += 1
    while(i<B and expansion_a[i] == expansion_b[i]):
        i += 1
    if i<B:
        return False
    return True
    
def check_image_iota(m,a):  # checks whether a is in the image of iota_m
     global stop            # here m>1
     l=m+1                
     a_l=a                  
     while a_l>0:
         expansion = numberToBase(a_l,l)
         if max(expansion)>m:
             return False
         a_l = iota_lists[l][ord(a_l,l)+1]
         l += 1         
     if l >= 6:             # this will be the last iteration for this m.
        stop=True
     return True

def check_image_iota_1(a):  # checks whether a is in the image of iota_1
     global stop           
     l=2
     a_l=a                  
     while a_l>0:
         expansion = numberToBase(a_l,l)
         if max(expansion)>1:
             return False
         c=first_nonzero_bits(expansion)
         if len(c)>1:
              (a,b) = c
              (a,b) = (iota_lists[l][a],iota_lists[l][b])
              if check_S(l+1,a,b):
                  return False         
         a_l = iota_lists[l][ord(a_l,l)+1]
         l += 1         
     if l >= 6:             # this will be the last iteration for this m.
        stop=True
     return True

def construct_image_iota(m): #constructs the beginning of the function iota_m
    global stop
    global Limit
    stop = False
    a = 0
    while(True):
         if a == Limit+1 and m == 1: # removing this line would cause the algorithm to run for decades
              break
         if m == 1:
             result = check_image_iota_1(a)
         else:
             result = check_image_iota(m,a)
         if result == True:
             iota_lists[m].append(a)
         if stop:
             return
         a += 1
       
def tau(m,a):
    return len_pi(m+1,chi(m+1,iota_lists[m][a]))

def chi(m,b):
    expansion = numberToBase(b,m)
    expansion = ''.join(str(v) for v in expansion)
    return int(expansion[::-1],m+1)

def len_pi(m,c):
    if c==0:
        return 0
    result=0
    for p in range(int(math.log(c,m+1))+1):
#        print(m,c,p)
        result += (int(c/(m+1)**p) - int(c/(m+1)**(p+1))) * (1+tau_lists[m][p+1])
    return result
    
def construct_tau(m): #constructs the beginning of the function tau^{(m)}
    global stop
    global Limit
    stop=False
    L=len(iota_lists[m])
    for a in range(1,L):
       tau_lists[m].append(tau(m,a))
  


# we will define iota_lists such that iota_m[p] equals iota_m(p) from the article
iota_lists=[-1,[-1],[-1],[-1],[-1],[-1],[-1]] # since iota_m(0) is undefined, we will never use the first values of each sublist.
# we also define tau_lists such that tau_lists[m][p] equals tau^{(m)}(p) from the article
tau_lists=[-1,[-1],[-1],[-1],[-1],[-1],[-1]]

for m in range(5,0,-1):
    construct_image_iota(m)
    print('iota_',m,': ',iota_lists[m][1:],sep='')
    
for m in range(5,0,-1):
    construct_tau(m)
    if m == 1: # to prevent printing too large lists, we take the first 80 elements 
        tau_lists[1] = tau_lists[1][0:Limit+1]
    print('tau_',m,': ',tau_lists[m][1:],sep='')

sigma_1=[]
beta_1=[1]
tau_1=tau_lists[1]
L=len(tau_1)
for t in range(1,L-1):
    S=tau_1[t+1]-tau_1[t]
    sigma_1.append(S)
    beta_1.append(2*beta_1[-1]+S)    

print('sigma_1:',sigma_1)
print('beta_1:',beta_1)

print('tau_1 is Entry A092432 in the OEIS')
print('sigma_1 is Entry A091579 in the OEIS.')
print('beta_1 is Entry A091411 in the OEIS.')
print('We print the values of iota_1 which are at most ',Limit,'.',sep='')
print('This code does not calculate curling numbers of strings.')
print('Instead, it uses the iota functions to construct the sequences.')