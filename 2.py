import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from twocython import diff


class plot():
    def __init__(self,U,x,t):
        self.U = U
        self.x = x
        self.t = t
        
    def plot_x(self):
        plt.plot(self.x,self.U[:,0*self.N_t//8], color='b', label='fot t=0.0')
        plt.plot(self.x,self.U[:,1*self.N_t//8], label='fot t=1.25')
        plt.plot(self.x,self.U[:,2*self.N_t//8], label='fot t=2.5')
        plt.plot(self.x,self.U[:,3*self.N_t//8], label='fot t=3.75')
        plt.plot(self.x,self.U[:,4*self.N_t//8], label='fot t=5.0')
        plt.plot(self.x,self.U[:,5*self.N_t//8], label='fot t=6.25')
        plt.plot(self.x,self.U[:,6*self.N_t//8], label='fot t=7.5')
        plt.plot(self.x,self.U[:,7*self.N_t//8], label='fot t=8.75')
        plt.plot(self.x,self.U[:,8*self.N_t//8 -1], label='fot t=10.0')
        plt.xlabel('x')
        plt.ylabel('температура')
        plt.title('Решение для разных времен')
        plt.legend(loc='right')
        plt.show()


       
    def plot_t(self):
        plt.plot(self.t, self.U[0,:], color='b', label='x = 0.0')
        plt.plot(self.t, self.U[1*self.N_x//10,:],  label='x = 1')
        plt.plot(self.t, self.U[2*self.N_x//10,:],  label='x = 2')
        plt.plot(self.t, self.U[4*self.N_x//10,:],  label='x = 4')
        plt.plot(self.t, self.U[5*self.N_x//10,:],  label='x = 5')
        plt.xlabel('время')
        plt.ylabel('температора')
        plt.title('Решения для разных координат')
        plt.legend(loc='upper right')
        plt.show()
   
    def plot_3dim(self):
        s = go.Surface(z=self.U, x=self.t, y=self.x, colorscale='inferno', opacity=0.9, hidesurface=False)
        s.contours.x = dict(size=0.5, show=True, start=0, end=self.T)
        s.contours.y = dict(size=0.5, show=True, start=0, end=self.L)
        layout = go.Layout(title = '1D heat equation', xaxis = dict(title = 'x coordinat'), yaxis = dict(title = 'time'))
        fig = go.Figure(data=[s], layout=layout)
        fig.update_layout(title='1D heat equation', xaxis = dict(title = 'x coordinat'), yaxis = dict(title = 'time'), autosize=False, width=1500, height=1000, scene = { "xaxis": {"nticks": 20}, "zaxis": {"nticks": 4}, 'camera_eye': {"x": 1, "y": 1, "z": 0.5}, "aspectratio": {"x": 1, "y": 1, "z": 0.5} })
        fig.update_traces(contours_z=dict(show=True, usecolormap=False, project_z=True))      
        fig.write_html('second.html')
        fig.show()        
        

n = 1000
init = np.zeros((n,n))
x = np.linspace(0,10,n)
t = np.linspace(0,10,n)
#начальные условия:
init[0,:] = 0
init[n-1,:] = 0
init[:,0] = 100
sol = diff().implicit_method(init,n)
graph = plot(sol,x,t)
graph.plot_x()
graph.plot_t()
graph.plot_3dim()