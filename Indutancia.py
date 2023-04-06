import numpy as np
import matplotlib.pyplot as plt
def redKron5(a):
    #Entrada uma matrix 4x4
    #Será feito a redução de kron
    f = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            f[i, j] = a[i, j] - a[i, 4]*a[4, j]/a[4, 4]
        
    return f

def redKron4(a):
    #Entrada uma matrix 4x4
    #Será feito a redução de kron
    f = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            f[i, j] = a[i, j] - a[i, 3]*a[3, j]/a[3, 3]
        
    return f

class Wire:
    def __init__(self, h, r, f, flecha=0, Ds=None):
        self.h = h
        self.r = r
        self.rfict=0.778*r
        self.f = f
        self.Ds = Ds
        self.flecha=flecha
        self.hreal=self.h-0.7*self.flecha

class WireMatrix(Wire):
    def __init__(self, wires, dist_matrix):
        self.wires = wires
        self.dist_matrix = dist_matrix

    def GroundDistance(self):
        n_wires = len(self.wires)
        D = np.zeros((n_wires, n_wires))
        for i in range(n_wires):
            for j in range(i+1, n_wires):
                D[i,j] = np.sqrt(self.dist_matrix[i,j]**2 + 4*self.wires[i].hreal*self.wires[j].hreal)
                D[j,i] = D[i,j]

        return D
    
    def MonophasicInductance(self, i=0, j=1):
        d12 = self.dist_matrix[i,j]
        D12 = self.GroundDistance()[i,j]
        h1 = self.wires[i].h
        rfict_1 = self.wires[i].rfict
        L = 2*10**(-4) * (np.log(d12/rfict_1) + np.log(2*h1/D12))
        freq=self.wires[0].f
        X = 2*np.pi*freq*L
        return X
    
    def ThreeConductorsMonophaseInductance(self):
        #Deve-se inserir na ordem a0, a'1, b2, b'3,c4 
        #d1 simboliza aa' e d2 simboliza bb'
        d0=self.dist_matrix[0,1]
        d03=self.dist_matrix[0,3]
        d21=self.dist_matrix[2,1]
        d2=self.dist_matrix[2,3]
        d41=self.dist_matrix[4,1]
        d43=self.dist_matrix[4,3]
        freq=self.wires[0].f
        dm=(d0*d03*d21*d2*d41*d43)**(1/6)

        Ds_1=self.wires[0].Ds
        Ds_2=self.wires[1].Ds

        d02=self.dist_matrix[0,2]
        d04=self.dist_matrix[0,4]
        d13=self.dist_matrix[1,3]
        DSI=(Ds_1**3 * d02**4 *d04**2)**(1/9)
        DSR=(Ds_2**2 * d13**2)**(1/4)
        print(DSR)
        Li=2*10**(-4)*np.log(dm/DSI)
        Lr=2*10**(-4)*np.log(dm/DSR)

        [Xi, Xr]=[2*np.pi*freq*Li, 2*np.pi*freq*Lr]
        return Xi, Xr
    def ThreePhaseInductanceWithoutGround(self):
        n_wires = len(self.wires)
        f = np.zeros((n_wires, n_wires))
        for i in range(n_wires):
            for j in range(n_wires):
                if i == j:
                    f[i,j] = 2*10**(-4)*np.log(1/self.wires[i].Ds)
                else:
                    f[i,j] = 2*10**(-4)*np.log(1/self.dist_matrix[i,j])
        
        L0 = f[0,0] - 1/2*(f[0,1] + f[0,2])
        L1 = f[1,1] - 1/2*(f[1,2] + f[1,0])
        L2 = f[2,2] - 1/2*(f[2,0] + f[2,1])
        freq=self.wires[0].f

        Ls=(L0+L1+L2)/3
        
        [X0, X1, X2, Xs]=[2*np.pi*freq*L0, 2*np.pi*freq*L1, 2*np.pi*freq*L2, 2*np.pi*freq*Ls]

        return X0, X1, X2, Xs
    
    def ThreePhaseInductanceWithoutGroundTrancado(self):
        Dm=1
        n_wires = len(self.wires)
        for i in range(n_wires):
            for j in range(n_wires):
                if i!=j and i>j:
                    Dm*=self.dist_matrix[i,j]

        Dm=Dm**(1/n_wires)
        freq=self.wires[0].f
        Xs=2*np.pi*freq*2*10**(-4)*np.log(Dm/self.wires[0].Ds)
      
        return Xs
    
    def ThreePhaseInductandeGround(self):
        n_wires = len(self.wires)
        f = np.zeros((n_wires, n_wires))
        D=self.GroundDistance()
        for i in range(n_wires):
            for j in range(n_wires):
                if i == j:
                    f[i,j] = 2*10**(-4)*np.log(2*self.wires[i].h/self.wires[i].Ds)
                else:
                    f[i,j] = 2*10**(-4)*np.log(D[i,j]/self.dist_matrix[i,j])
        
        L0 = f[0,0] - 1/2*(f[0,1] + f[0,2])
        L1 = f[1,1] - 1/2*(f[1,2] + f[1,0])
        L2 = f[2,2] - 1/2*(f[2,0] + f[2,1])
        print(L0)
        print(L1)
        print(L2)
        Ls=(L0+L1+L2)/3
        [x0,x1,x2,xs]=[2*np.pi*60*L0, 2*np.pi*60*L1, 2*np.pi*60*L2, 2*np.pi*60*Ls]

        return x0, x1, x2, xs
    
    def CalcAuxVariables(self):
        Dmi=1
        Dm=1
        hm=1
        D=self.GroundDistance()
        n_wires = len(self.wires)
        for i in range(n_wires):
            for j in range(n_wires):
                if i!=j and i>j:
                    Dm*=self.dist_matrix[i,j]
                    Dmi*=D[i,j]
                if i==j:
                    hm*=self.wires[i].hreal
        
        Dm=Dm**(1/n_wires)
        Dmi=Dmi**(1/n_wires)
        hm=hm**(1/n_wires)
        
        return Dm, Dmi, hm
    def ThreePhaseInductanceGroundTrancado(self):
        Dm, Dmi, hm = self.CalcAuxVariables()
        Xs=2*np.pi*self.f*2*10**(-4)*(np.log(2*hm/self.wires[0].Ds)-np.log(Dmi/Dm))
        return Xs
    
    def ThreePhaseInductanceforLightingTransposto(self):
        n_wires = len(self.wires)
        hm=1
        Dm=1
        Dmi=1
        d=self.dist_matrix
        D=self.GroundDistance()
        for i in range(n_wires):
            for j in range(n_wires):
                if i!=j and i>j:
                    Dm*=self.dist_matrix[i,j]
                    Dmi*=D[i,j]
                if i==j:
                    hm*=self.wires[i].hreal
        n_wires = len(self.wires)
        hm=hm**(1/n_wires)
        Dm=Dm**(1/n_wires)
        Dmi=Dmi**(1/n_wires)

        freq=self.wires[0].f
        f_ = 2e-4*np.log(2*hm/self.wires[0].Ds)
        f_r = 2e-4*np.log(2*self.wires[3].h/self.wires[1].Ds)
        Dmi = (D[0,1]*D[1,2]*D[2,0])**(1/3)
        Dm = (d[0,1]*d[1,2]*d[2,0])**(1/3)
        Dmir = (D[0,3]*D[1,3]*D[2,3])**(1/3)
        Dmr = (d[0,3]*d[1,3]*d[2,3])**(1/3)
        f = 2*10**(-4)*np.log(Dmi/Dm)
        fr = 2*10**(-4)*np.log(Dmir/Dmr)
        F = np.array([[f_, f, f, fr],
                      [f, f_, f, fr],
                      [f, f, f_, fr],
                      [fr, fr, fr, f_r]])
        x = F*2*np.pi*freq
        y = x[0,0] - x[0,1]
        
        return y

    def ThreePhaseInductanceforLighting(self):
        f=np.zeros((4,4))
        D=self.GroundDistance()
        print(D)
        d=self.dist_matrix
        freq=self.wires[0].f

        for i in range(4):
            for j in range(4):
                if i==j:
                    f[i,i]=2*10**(-4) *np.log(2*self.wires[i].h/self.wires[i].Ds)
                else:
                    f[i, j] = 2*10**(-4)*np.log(D[i,j]/d[i,j])
        x = f*2*np.pi*freq
        print("Matriz f \n", f)
        
        r = redKron4(f)
        print("Matriz reduzida \n", r)
        y = (r[0,0]+r[1,1]+r[2,2])/3 - (r[0,1]+r[0,2]+r[1,2])/3
        return y
    
    
    def ThreePhaseInductanceforLightingWithoutGround(self):
        f=np.zeros((4,4))
        D=self.GroundDistance()
        print(D)
        d=self.dist_matrix
        freq=self.wires[0].f

        for i in range(4):
            for j in range(4):
                if i==j:
                    f[i,i]=2*10**(-4) *np.log(1/self.wires[i].Ds)
                else:
                    f[i, j] = 2*10**(-4)*np.log(1/d[i,j])
        x = f*2*np.pi*freq
        print("Matriz f \n", f)
        
        r = redKron4(f)
        print("Matriz reduzida \n", r)
        y = (r[0,0]+r[1,1]+r[2,2])/3 - (r[0,1]+r[0,2]+r[1,2])/3
        return y
    
    
    def PlotBizurado(self, delta):
        n_wires = len(self.wires)
        f = np.zeros((n_wires, n_wires))
        VecDs=np.linspace(self.wires[0].Ds*(1-delta), self.wires[0].Ds*(1+delta), 500)
        X=np.zeros(500)
        P=np.zeros(500)
        i=0
        print(VecDs)
        for Ds in VecDs:
            
            for i in range(n_wires):
                for j in range(n_wires):
                    if i == j:
                        f[i,j] = 2*10**(-4)*np.log(1/Ds)
                    else:
                        f[i,j] = 2*10**(-4)*np.log(1/self.dist_matrix[i,j])
             
            L0 = f[0,0] - 1/2*(f[0,1] + f[0,2])
            L1 = f[1,1] - 1/2*(f[1,2] + f[1,0])
            L2 = f[2,2] - 1/2*(f[2,0] + f[2,1])
            freq=self.wires[0].f

            Ls=(L0+L1+L2)/3

            [X0, X1, X2, Xs]=[2*np.pi*freq*L0, 2*np.pi*freq*L1, 2*np.pi*freq*L2, 2*np.pi*freq*Ls]
            P[i]=1/Xs
            i+=1
            
                
        plt.plot(VecDs, P)
        plt.show()
        
    def PlotBizuradoGuarni(self, delta):
        Ds_1 = self.wires[0].Ds
        n_wires = len(self.wires)
        VecDs = np.linspace(Ds_1*(1-delta), Ds_1*(1+delta), 500)
        Xs = np.zeros(500)
        P = np.zeros(500)
        for k in range(len(VecDs)):
            f = np.zeros((n_wires, n_wires))
            for i in range(n_wires):
                for j in range(n_wires):
                    if i == j:
                        f[i,j] = 2*10**(-4)*np.log(1/VecDs[k])
                    else:
                        f[i,j] = 2*10**(-4)*np.log(1/self.dist_matrix[i,j])
            
            L0 = f[0,0] - 1/2*(f[0,1] + f[0,2])
            L1 = f[1,1] - 1/2*(f[1,2] + f[1,0])
            L2 = f[2,2] - 1/2*(f[2,0] + f[2,1])
            freq=self.wires[0].f

            Ls=(L0+L1+L2)/3
            
            Xs[k] = 2*np.pi*freq*Ls
            P[k] = 1/Xs[k]
        plt.plot(VecDs, P)
        plt.show()
        
    def ThreePhaseInductance2forLighting(self):
        f=np.zeros((5,5))
        D=self.GroundDistance()
        print("Matrid distancia: \n", D)
        d=self.dist_matrix
        freq=self.wires[0].f

        for i in range(5):
            for j in range(5):
                if i==j:
                    f[i,i]=18 * 1e6 *np.log(2*self.wires[i].hreal/self.wires[i].r)
                else:
                    f[i, j] = 18 * 1e6 *np.log(D[i,j]/d[i,j])
        print("Matriz de campo elétrico \n", f)
        #x = f*2*np.pi*freq
        raux = redKron5(f)
        print("Reduzida: \n", raux)
        r=redKron4(raux)
        print("Reduzida final: \n", r)
        y = (r[0,0]+r[1,1]+r[2,2])/3 - (r[0,1]+r[0,2]+r[1,2])/3
        return y
        
