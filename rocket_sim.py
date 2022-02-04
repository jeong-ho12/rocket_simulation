import math

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d

class Rocket:
    def __init__(
        self,
        mass_struct,
        mass_pro,
        t_b,
        drag_coeff,
        position = np.array([0,0,0]),
        velocity = np.array([0,0,0]),
    ):
        self.mass_pro = mass_pro                    ## 연료 질량
        self.mass_struct = mass_struct              ## 구조체 질량
        self.t_b = t_b                              ## 연소 시간
        self.drag = drag_coeff                      ## 항력 계수
        self.position = position                    ## 로켓 위치
        self.velocity = velocity                    ## 로켓 속도
        self.mass = mass_pro+mass_struct            ## 로켓 총 무게
        self.zeroparam = [0,self.mass,0]
        self.force_effect = []

    ## 연소 중 미분방정식
    def force_burn(self,zeroparam,t,b):
        T = 400/self.t_b*math.sin(math.pi/self.t_b*realTime)       ## 임펄스 200인 sin 함수 모양의 추력
        v ,m, h = zeroparam
        g = 9.8
        b = b
        # rho = 1.225*(1-2.256e-5*h)**5.256
        d = 0.16
        # S = np.pi*d*d/4

        
        dm_dt = -self.mass_pro/self.t_b
        dv_dt = T/m -g - b*v/m
        dh_dt = v
        return [dv_dt,dm_dt, dh_dt]
        
    ## 자유낙하 미분방정식
    def force_free(self,zeroparam,t,b):
        v, m, h= zeroparam
        g = 9.8
        b = b
        d = 0.16


        dm_dt = 0
        dv_dt = -g - b*v/m
        dh_dt = v
        return [dv_dt,dm_dt,dh_dt]

    ## 외력을 계산하기 위한 함수
    def calcul_force_effect(self):
        ## 이륙 ~ 연소 종료까지
        if realTime <= self.t_b:
            t1 = np.linspace(0,0.1,11)
            self.force_effect = odeint(self.force_burn,self.zeroparam,t1,args=(self.drag,))

        ## 연소종료 ~ 자유낙하
        else :
            t2 = np.linspace(0,0.1,11)
            self.force_effect = odeint(self.force_free,self.zeroparam,t2,args=(self.drag,))
            print(self.force_effect)




if __name__ == '__main__':
    t = np.linspace(0,0.1,11)   # calcul 0.1s
    drag_coeff = 0.2            # 항력 계수
    mass_struct = 2            # 구조체 질량 2kg
    mass_pro = 0.15             # 추진제 질량 0.15kg
    burnTime = 2                # 연소 시간 2s # 연소 시간을 바꿔도 총 충격량은 200ns로 고정. 충격량 변화를 위해선, 미분방정식 속 추력 함수 참고.
    realTime = 0
    rocket = Rocket(mass_struct,mass_pro,burnTime,drag_coeff,velocity=[0,0,0])   # 초기 속도 0,0,0
    while realTime<=20:
        rocket.calcul_force_effect()
        rocket.zeroparam = rocket.force_effect[-1]
        realTime += 0.1
