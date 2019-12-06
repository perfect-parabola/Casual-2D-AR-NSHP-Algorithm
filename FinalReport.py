import math, numpy as np

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

class H:
    def __init__(self, theta):
        self.theta = theta

    def GetRho(self, j, k):
        if (j+k)%2 ==1 and j >=0 and k >=0:
            return 0

        elif (j + k) % 2 == 0 and j >= 0 and k >= j:
            return ((1-math.sqrt(1-4*(self.theta**2)))/(2*self.theta))**k

        elif (j + k) % 2 == 0 and k >= 0 and j > k:
            response =  ((1-math.sqrt(1-4*(self.theta**2)))/(2*self.theta))**k
            max_a = int((j-k)/2-1)
            for a in range(0, max_a+1):
                response += -math.sqrt(1-4*(self.theta**2))*nCr(2*a+k, a+k)*self.theta**(2*a+k)
            return response
        elif j<0:
            return self.GetRho(-j, k)
        elif k <0:
            return self.GetRho(j, -k)
        else:
            print("GetRho Error")

    def GetRhos(self, j, k, n):
        if n == 0:
            return np.array([self.GetRho(j, k)])
        elif n >0:
            list1 = [self.GetRho(j-i, k+n) for i in range(1, n+1)]
            list1 = np.array(list1)
            list2 = [self.GetRho(j - n, k + n -i) for i in range(1, n+1)]
            list2 = np.array(list2)
            list3 = [self.GetRho(j - n, k - i) for i in range(1, n+1)]
            list3 = np.array(list3)
            list4 = [self.GetRho(j - n + i, k - n) for i in range(1, n+1)]
            list4 = np.array(list4)
            return np.concatenate([list1, list2, list3, list4])
        else:
            print("GetRhos Error")

    def GetR(self, a, b):
        if a == 0 and b == 0:
            return np.array([1])
        elif a == 0:
            return self.GetRhos(0,0,b)
        else:
            R1 = np.stack([self.GetRhos(i, -a, b) for i in range(1, a+1)], axis=0)
            R2 = np.stack([self.GetRhos(a, i-a, b) for i in range(1, a+1)], axis=0)
            R3 = np.stack([self.GetRhos(a, i, b) for i in range(1, a+1)], axis=0)
            R4 = np.stack([self.GetRhos(a-i, a, b) for i in range(1, a+1)], axis=0)
            return np.vstack((R1, R2, R3, R4))

    def GetPhi(self, a, b, c):
        if a == 2 and b == 1 and c == 0:
            return np.matmul(np.linalg.inv(self.GetR(1, 1)), self.GetR(1, 2))
        elif a == 2 and b == 1:
            return np.matmul(np.linalg.inv(self.GetR(1, 1)), self.GetR(1, c+2))
        elif a - b  == 1:
            x = self.GetR(b, b+c+1)
            for i in range(1, b):
                x += -np.matmul(self.GetR(b, i), self.GetPhi(b, i, c+1))
            return np.matmul(np.linalg.inv(self.GetTHETA(b)), x)
        else:
            return self.GetPhi(a-1, b, c+1) - np.matmul(self.GetPhi(a-1, b, 0), self.GetPhi(a, a-1, c))

    def GetTHETA(self, m):
        response = self.GetR(m, m)
        for a in range(1, m):
            response += -np.matmul(self.GetR(m, a), self.GetPhi(m, a, 0))
        return response

    def Geth(self, m):
        response = self.GetR(m, 0)
        for a in range(1, m):
            response += - np.matmul(self.GetR(m, a), self.GetBeta(m-1, a))
        return response

    def GetBeta(self, m, a):
        if m == 1 and a == 1:
            return np.matmul(np.linalg.inv(self.GetR(1, 1)), self.GetR(1, 0))
        elif m == a:
            return np.matmul(np.linalg.inv(self.GetTHETA(m)), self.Geth(m))
        else:
            return self.GetBeta(m-1, a) - np.matmul(self.GetPhi(m, a, 0), self.GetBeta(m, m))

    def GetLambda(self, m):
        if m == 2:
            return 1 - np.matmul(self.GetR(0, 1), self.GetBeta(1, 1))
        return self.GetLambda(m-1) - np.matmul(np.transpose(self.Geth(m)), self.GetBeta(m, m))

    def GetSigmaSquare(self, m):
        return self.GetLambda(m+1)


