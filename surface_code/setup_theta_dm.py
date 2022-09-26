import numpy as np

def main():
    # Want theta from 1e-7 to pi/4, with 20 points per decade
    # pi/4 ~0.785xxx, so stop there
    # Want deltaM from 1e-9 to 1e4, with 20 points per decade
    
    m_exp  = np.arange(-9, 4, 1)
    m_mant = np.arange(0.5, 10.5, 0.5)
    t_exp  = np.arange(-8, 0, 1)
    t_mant = np.copy(m_mant)

    with open("deltaMs", "w") as m:
        for i in m_exp:
            for j in m_mant:
                m.write(f"{j}E{i}\n")
        m.write("1.0E4\n")
    
    with open("thetas", "w") as t:
        for i in t_exp:
            for j in t_mant:
                if i != -1:
                    t.write(f"{j}E{i}\n")
                elif (i == -1 and j <= 7.5):
                    t.write(f"{j}E{i}\n")
                elif (i == -1 and j > 7.5):
                    t.write("7.85E-1\n")
                    break
                    
#    with open("thetas", "w") as t:
#        for i in t_exp:
#            for j in t_mant:
#                if i != -1:
#                    t.write(f"{j}E{i}\n")
#                elif i == -1 and j <= 7.5:
#                    t.write(f"{j}E{i}\n")
#                elif i == -1 and j > 7.5:
#                    t.write("7.85E-1\n")
#                    break

if __name__ == "__main__":
    main()
