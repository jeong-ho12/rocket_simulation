from scipy.integrate import odeint
import numpy as np
from transform import Transformer
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

class Rocket:
    def __init__(
        self,
        mass_struct,
        mass_pro,
        t_b,
        position,
        velocity
    ):
        self.mass_pro = mass_pro                    ## propellant mass
        self.mass_struct = mass_struct              ## structure mass
        self.t_b = t_b                              ## burnout time
        self.position = position                    ## position
        self.velocity = velocity                    ## velocity
        self.mass = mass_pro + mass_struct          ## mass
        self.zeroparam = np.hstack((self.velocity, self.position, self.mass)) ## zeroparam
        self.force_effect = []

    ## differential equation: 
    ## input: velocity, position, mass, drag coefficient, TVC angle, heading angle
    def force_burn(self, zeroparam, t, Cd, psi, theta, phi, alpha, beta):
        vx, vy, vz, px, py, pz, m = zeroparam

        T = 400/self.t_b*np.sin(np.pi/self.t_b*t)
        M = Transformer().body_to_earth(np.array([psi, theta, phi]))
        con_vec = Transformer().TVC_convential_vector(alpha, beta)

        v = np.array([vx,  vy,  vz])
        p = np.array([px,  py,  pz])
        g = np.array([ 0,   0, 9.8])

        rho = 1.225*(1-2.256e-5*p[2])**5.256
        d = 0.16
        S = np.pi*d*d/4
        k = 0.5*rho*S*Cd*np.linalg.norm(v)

        m_dot = -self.mass_pro/self.t_b
        vx_dot, vy_dot, vz_dot = (M@con_vec*T -k*v)/m -g
        px_dot, py_dot, pz_dot = v
        
        return np.array([vx_dot, vy_dot, vz_dot, px_dot, py_dot, pz_dot, m_dot])
        
    ## differential equation: free fall
    def force_free(self, zeroparam, t, Cd, psi, theta, phi, alpha, beta):
        vx, vy, vz, px, py, pz, m = zeroparam

        M = Transformer().body_to_earth(np.array([psi, theta, phi]))

        v = np.array([vx,  vy,  vz])
        p = np.array([px,  py,  pz])
        g = np.array([ 0,   0, 9.8])

        rho = 1.225*(1-2.256e-5*p[2])**5.256
        d = 0.16
        S = np.pi*d*d/4
        K = 0.5*rho*S*Cd*np.linalg.norm(v)**2

        if v[2] >= 0:
            m_dot = 0
            vx_dot, vy_dot, vz_dot = -K*v/np.linalg.norm(v)/m -g
            px_dot, py_dot, pz_dot = v
        else:
            m_dot = 0
            vx_dot, vy_dot, vz_dot = -K*v/np.linalg.norm(v)/m -g
            px_dot, py_dot, pz_dot = v

        return np.array([vx_dot, vy_dot, vz_dot, px_dot, py_dot, pz_dot, m_dot])

    def calcul_force_effect(self, external_element):
        ## 이륙 ~ 연소 종료까지
        if realTime <= self.t_b:
            t1 = np.linspace(0,0.1,11)
            self.force_effect = odeint(self.force_burn, self.zeroparam, t1, args=tuple(external_element))

        ## 연소종료 ~ 자유낙하
        else :
            t2 = np.linspace(0,0.1,11)
            self.force_effect = odeint(self.force_free, self.zeroparam, t2, args=tuple(external_element))


if __name__ == '__main__':
    t = np.linspace(0,0.1,11)   # calcul 0.1s
    drag_coeff = 0.2            # 항력 계수
    mass_struct = 2            # 구조체 질량 2kg
    mass_pro = 0.15             # 추진제 질량 0.15kg
    burnTime = 2                # 연소 시간 2s # 연소 시간을 바꿔도 총 충격량은 200ns로 고정. 충격량 변화를 위해선, 미분방정식 속 추력 함수 참고.
    realTime = 0
    v = np.array([0, 0, 100])
    p = np.array([0, 0, 0])
    TVC_angle = np.array([0, 0])
    heading_angle = np.array([0, 0, 1])
    external_element = np.hstack((drag_coeff, TVC_angle, heading_angle))

    rocket = Rocket(mass_struct, mass_pro, burnTime, p, v)

    ax = plt.axes(projection="3d")

    plt.gca().set_xlim(-15,15)
    plt.gca().set_ylim(-15,15)
    plt.gca().set_zlim(-10,10)  
    plt.quiver(0,0,0,1,0,0,length=5,color='red')
    plt.quiver(0,0,0,0,1,0,length=5,color='green')
    plt.quiver(0,0,0,0,0,1,length=5,color='blue')

    while realTime <= 20:
        rocket.calcul_force_effect(external_element)
        rocket.zeroparam = rocket.force_effect[-1,:]
        realTime += 0.1

        for i in range(len(rocket.force_effect[:,0])):
            ax.cla()
            ax.set_xlim(0, 500)
            ax.set_ylim(0, 500)
            ax.set_zlim(0, 500)
          
            ax.quiver(rocket.force_effect[i,3], rocket.force_effect[i,4], rocket.force_effect[i,5], \
                rocket.force_effect[i,0], rocket.force_effect[i,1], rocket.force_effect[i,2], \
                linewidth=2.0, arrow_length_ratio=0.0, color='black')
            
            
        plt.pause(0.1)
    
    plt.show()
