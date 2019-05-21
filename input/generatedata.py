import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt
import json


class ThetaLogistic:

  def __init__(self, b0=0.15, b2=-0.125, b3=0.1, Q=0.47**2, R=0.39**2, L0=5):
    self.b0 = b0
    self.b2 = b2
    self.b3 = b3
    self.Q = Q
    self.R = R
    self.L0 = L0
    print(L0)

  def simulate(self, T):
    Y = np.zeros(T)
    L = np.zeros(T)
    L[0] = self.L0
    Y[0] = L[0] + np.sqrt(self.R)*randn()
    for t in range(1,T):
      L[t] = L[t-1] + self.b0 + self.b2*np.exp(self.b3*L[t-1]) + np.sqrt(self.Q)*randn()
      Y[t] = L[t] + np.sqrt(self.R)*randn()

    return (L,Y)







def main():
  mdl = ThetaLogistic()

  T = 150
  L,Y = mdl.simulate(T)

  with open("simulation_data.json", "w") as f:
    json.dump({"L": L.tolist(), "y":Y.tolist()}, f)


  fig, (ax1, ax2) = plt.subplots(2,1)
  ax1.plot(range(T), L)
  ax1.plot(range(T), Y)

  plt.show()



if __name__ == "__main__":
  main()

