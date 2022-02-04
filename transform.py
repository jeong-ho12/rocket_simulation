import numpy as np
import matplotlib.pyplot as plt

class Transformer():
    def __init___(self):
        pass
    
    def body_to_earth(self, rocAtt):
        
        psi = rocAtt[0] * np.pi/180
        theta = rocAtt[1] * np.pi/180
        phi = rocAtt[2] * np.pi/180
        
        y = np.array([[np.cos(phi), -np.sin(phi), 0],
                      [np.sin(phi),  np.cos(phi), 0],
                      [          0,            0, 1]])

        p = np.array([[ np.cos(theta), 0, np.sin(theta)],
                      [             0, 1,             0],
                      [-np.sin(theta), 0, np.cos(theta)]])

        r = np.array([[1,           0,            0],
                      [0, np.cos(psi), -np.sin(psi)],
                      [0, np.sin(psi),  np.cos(psi)]])
                     
        PP = np.array([[0, 0, -1],
                       [0, 1,  0],
                       [1, 0,  0]])
             
        transformation_mtx = y@p@r
                    
        return transformation_mtx
    
    def TVC_convential_vector(self, theta1, theta2):
        convential_vector = np.array([np.cos(theta1)*np.cos(theta2),
                                      np.sin(theta2),
                                      -np.cos(theta2)*np.sin(theta1)])
        return convential_vector
