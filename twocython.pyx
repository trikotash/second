import cython
cimport cython
import numpy as np
cimport numpy as np

cdef class diff():
    
    cpdef explicit_method(self,init,n):
        cdef:
            np.ndarray U
            double K
            int i

        U = np.copy(init)
        K = 1.0 
        for i in range(1,n):
            U[1:n-1,i] = U[1:n-1,i-1] + K*(U[2:n,i-1] + U[0:n-2,i-1] - 2*U[1:n-1,i-1])
        
        return U


    cpdef implicit_method(self,init,n):
        cdef:
            double [:,:] m_tech
            double [:] A,B,C,D,alpha,beta,y
            double K
            int i
        
        m_tech = memoryview(init)
        K = 1.0
        y = np.zeros(n)
        alpha = np.zeros(n)
        beta = np.zeros(n)
        A = np.zeros(n) - 1
        B = np.zeros(n) + 2 + 2/K
        C = np.zeros(n) - 1
        D = np.zeros(n) 
        for j in range(1,n):
            D[1] = m_tech[1,j-1]*(2/K - 2) + m_tech[2,j-1] + m_tech[0,j-1] + m_tech[0,j]
            y[1] = B[1]
            alpha[1] = -C[1]/y[1]
            beta[1] = D[1]/y[1]
            for i in range(2,n-2):
                D[i] = m_tech[i,j-1]*(2/K - 2) + m_tech[i+1,j-1] + m_tech[i-1,j-1]
                y[i] = B[i] + A[i]*alpha[i-1]
                alpha[i] = -C[i]/y[i]
                beta[i] = (D[i] - A[i]*beta[i-1])/y[i]
            D[n-2] = m_tech[n-2,j-1]*(2/K - 2) + m_tech[n-1,j-1] + m_tech[n-3,j-1] + m_tech[n-1,j]
            y[n-2] = B[n-2] + A[n-2]*alpha[n-3]
            beta[n-2] = (D[n-2] - A[n-2]*beta[n-3])/y[n-2]

            m_tech[n-2,j] = beta[n-2]
            for i in range(n-3,0,-1):
                m_tech[i,j] = alpha[i]*m_tech[i+1,j] + beta[i]
        
        return init
                
    
    
    cpdef implicit_newton_cooling(self,init,n):
        cdef:
            double [:,:] m_tech
            double [:] A,B,C,D,alpha,beta,y
            double K,Ue,h
            int i
        
        K = 1.0

        Ue = 60
        h = 0.5
    
        m_tech = memoryview(init)
        y = np.zeros(n)
        alpha = np.zeros(n)
        beta = np.zeros(n)
        A = np.zeros(n) - 1
        B = np.zeros(n) + 2 + 2/K
        C = np.zeros(n) - 1
        D = np.zeros(n) 
        for l in range(1,n):
            D[1] = m_tech[1,l-1]*(2/K - 2) + m_tech[2,l-1] + m_tech[0,l-1] + m_tech[0,l] - h*(m_tech[1,l-1] - Ue)*2/K
            y[1] = B[1]
            alpha[1] = -C[1]/y[1]
            beta[1] = D[1]/y[1]
            for i in range(2,n-2):
                D[i] = m_tech[i,l-1]*(2/K - 2) + m_tech[i+1,l-1] + m_tech[i-1,l-1] - h*(m_tech[i,l-1] - Ue)*2/K
                y[i] = B[i] + A[i]*alpha[i-1]
                alpha[i] = -C[i]/y[i]
                beta[i] = (D[i] - A[i]*beta[i-1])/y[i]
            D[n-2] = m_tech[n-2,l-1]*(2/K - 2) + m_tech[n-1,l-1] + m_tech[n-3,l-1] + m_tech[n-1,l] - h*(m_tech[n-2,l-1] - Ue)*2/K
            y[n-2] = B[n-2] + A[n-2]*alpha[n-3]
            beta[n-2] = (D[n-2] - A[n-2]*beta[n-3])/y[n-2]
            m_tech[n-2,l] = beta[n-2]
            for j in range(n-3,0,-1):
                m_tech[j,l] = alpha[j]*m_tech[j+1,l] + beta[j]
        
        return init
    
    


    
    
    