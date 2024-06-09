""" - Get time response of vibrational systems

    - Get natural frequencies of systems
    
    - Work with 1 DOF systems and MDOF systems
    
    - Systems with free and force vibration
    
    Classes
    -------
    
     - DOF1
     
       - For systems of 1 DOF
       
     - MDOF
     
       - For MDOF systems
       
     Functions
     ---------
     
      - DOF1
      
        - Creates a 1 DOF system
        
      - MDOF
      
        - Creates a MDOF system
        
      - A
      
        - Return the space-state matrix A
        
      - B
      
        - Return the space-state matrix B
        
      - plot_responses
      
        - Plot the displacement response of the system
        
     Attributes
     ----------
     
      - wn
      
        - Return the natural frequencies of the system
        
      - wd
      
        - Return the damped frequencies of the system
        
      - Damping_Ratio
      
        - Return the damping ratios of the system
        
        """


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.linalg import eigvals


class DOF1:
    
    def __init__(self, M, K, C):
        
        self.M = M
        self.K = K
        self.C = C
        
        self.wn = None
        self.wd = None
        self.cc = None
        self.z = None
        self.r = None
        self.y = None
    
    def DOF1(self, CI, t, F0=None, w=None):

        """ Response of a 1 DOF Free Vibration System
        
            Parameters
            ------------------------------------------
            
            M - Mass (kg)
            
            K - Stiffiness (N*m)
            
            C - Damping (N*s/m)
            
            CI - Initial Conditions
            
            t - Time Span (s)
            
            For Harmonic excited systems, enter the following values:
                
            F - Force Magnitude (N)
            
            w - Excitation Frequency (rad/s)
            
            """



        self.wn = np.sqrt(self.K/self.M)
        self.cc = 2*np.sqrt(self.K*self.M)
        self.z = self.C/self.cc

        if self.z > 1:
            self.wd = self.wn*(np.sqrt((1-(self.z**2))+0j))
        else:
            self.wd = self.wn*(np.sqrt(1-(self.z**2)))
        
        if w==None:
            self.r = None
        
        else:
            self.r = w/self.wn
            
        if F0==None:

            if self.z<1:

                A = np.sqrt(((CI[0])**2) + ((((CI[1]) + (self.z*self.wn*CI[0]))/self.wd)**2))
                psi = np.arctan((self.wd*CI[0])/((CI[1]) + (self.z*self.wn*CI[0])))

                x = A*np.exp(-self.z*self.wn*t)*np.sin((self.wd*t)+psi)

            elif self.z == 1:

                x = (CI[0]*np.exp(-self.wn*t)) + ((CI[1] + (self.wn*CI[0]))*t*np.exp(-self.wn*t))

            else:

                S1 = (-self.z*self.wn) + (self.wn*np.sqrt((self.z**2)-1))
                S2 = (-self.z*self.wn) - (self.wn*np.sqrt((self.z**2)-1))

                A1 = ((S2*CI[0]) - CI[1])/(S2-S1)
                B1 = ((-S1*CI[0]) + CI[1])/(S2-S1)

                x = (A1*np.exp(S1*t)) + (B1*np.exp(S2*t))

            self.y = [x, t]

        else:
            
            m = self.M
            k = self.K
            c = self.C

            def dof1(S, t):
                x, v = S
                return [v, 
                        ((F0*np.sin(w*t))/m) - ((k/m)*x) -((c/m)*v)]

            y = odeint(dof1, CI, t)
            self.y = y.T
        
    def plot_responses(self, t):
        plt.title('Displacement Response')
        plt.plot(t, self.y[0])
        plt.xlabel('Time (s)')
        plt.ylabel('Displacement (m)')
        
        return plt.show()
            
                
class MDOF:
           
    def __init__(self, M, K, C):
        
        self.M = M
        self.K = K
        self.C = C
        self.n = len(M)
        self.resp = None
        self.Damping_Ratio = None
        self.wn = None
        self.wd = None
        
    def A(self):
        I = np.eye(self.n)
        Ze = np.zeros_like(I)
        A = np.vstack((np.hstack((Ze, I)), np.hstack((-(np.matmul(np.linalg.inv(self.M),self.K)), -(np.matmul(np.linalg.inv(self.M),self.C))))))
        return A

    
    def B(self):
        Ze = np.zeros_like(self.M)
        B =  np.vstack((Ze, np.linalg.inv(self.M)))
        return B
       
    def MDOF(self, CI, t, F=None, w = None, name=''):

        """ Response of MDOF Vibration Systems
        
            Parameters
            ----------
            
            M - Mass (kg) [Array]
            
            K - Stiffiness (N*m) [Array]
            
            C - Damping (N*s/m) [Array]
            
            CI - Initial Conditions [Array]
            
            t - Time Span (s)
            
            F - Sinusoidal Harmonic Forces (N) [Array]
            
              If excitation force is sinusoidal, type 'sin' in input
              
              If excitation force is a cosine wave, type 'cos' in input
              
            w - Excitation Frequencies (rad/s) [Array]
            
            """

        def sys(x, t):
            if name == 'cos':
                u = F*np.cos(w*t)
            else:
                u = F*np.sin(w*t)
            ddx = np.dot(self.A(), x) + np.dot(self.B(), u)
            
            return ddx
        
        sysn = odeint(sys, CI, t)
        self.resp = sysn.T
        
        self.wn = np.unique(np.absolute(eigvals(self.A())))[:self.n]
        
        self.wd = np.unique(np.absolute(np.imag(eigvals(self.A()))))[:self.n]
        
        self.Damping_Ratio = np.sqrt(-(((self.wd**2)/((self.wn**2)))) + np.repeat(1, len(self.M)))
        
        return sysn

    def plot_responses(self, t):
        
        """Plot Responses of MDOF Vibration Systems"""
        
        ax =[]
        
        for i in range(self.n):
            ax.append(str(i+1))
        
        f, ax = plt.subplots(self.n, 1, sharex=True)
        f.figure.set_figwidth(10)
        f.set_figheight(8)
        
        f.suptitle('Displacement Response')
        ax[self.n-1].set_xlabel('Time (s)')
        
        for i in range(self.n):
            
            desloc = np.array(self.resp[i])
            text = f'Displacement {i+1} (m)'
            
            ax[i].set_ylabel(text)
            ax[i].plot(t, self.resp[i])
            ax[i].set_xlim(min(t), max(t))
            ax[i].set_ylim(-abs(desloc.min()) - (abs(desloc.min())/10), abs(desloc.max())+(abs(desloc.max())/10))
            plt.tight_layout()
            
        return plt.show






