



import numpy as np

'''シュミレーション設定'''

grid_X=0
grid_Y=0
grid_Z=0


'''分割数(x*y*x)'''
SimuTime=0
'''シュミレーションする現実時間[s]'''
mleng_X=0
mleng_Y=0
mleng_Z=0
'''各軸方向のメッシュ間の距離[cm]'''
DT=0
'''微小時間[s]'''




'''各種物性値'''

E_com_In=0  
''' Inの結合エネルギー [eV]'''
SpeHeat_InN=0
SpeHeat_GaN=0
''' InN,GaNの比熱 [J/kg*K)'''
ThermalCon_GaN=0
''' InN,GaNの熱伝導率[w/mK]'''
Density_InN=0
Density_GaN=0
''' InN,GaNの密度[g/cm3]'''




'''計算用変数'''



log=np.ones([grid_X,grid_Y,grid_Z,int(SimuTime/DT)+1,1,1,1,1])
'''記録用
    座標(x,y,z):
    ステップ(時間の指標):
    温度[K]:
    熱量輸送速度[1/J*s]:
    Inの濃度[N/cm3]:
    In輸送速度[1/N*s]'''




'''一時保存変数
       代入による計算ロスを無くすため、後で消す'''

DInM_XY=np.ones([grid_X-2,grid_Y-2,grid_Z-2,1])
'''In拡散係数(xy面上)'''
DInM_Z=np.ones([grid_X-2,grid_Y-2,grid_Z-2,1])
'''In拡散係数(z方向)'''
Dq=np.ones([grid_X-2,grid_Y-2,grid_Z-2,1])
'''熱拡散係数(xyz方向)'''
A_ND_X=DT/(mleng_X*mleng_X)
'''Δt/Δx^2'''
A_ND_Y=DT/(mleng_Y*mleng_Y)
'''Δt/Δx^2'''
A_ND_Z=DT/(mleng_Z*mleng_Z)
'''Δt/Δx^2'''



'''計算過程'''


for Step in range(int(SimuTime/DT)):
    

    
    '''拡散係数の更新'''


    Dq=Dq=np.ones([grid_X-2,grid_Y-2,grid_Z-2,1])
    '''熱拡散係数(xyz方向)'''
    
    DInM_XY=np.ones([grid_X-2,grid_Y-2,grid_Z-2,1])
    '''In拡散係数(xy面上)'''
    DInM_Z=np.ones([grid_X-2,grid_Y-2,grid_Z-2,1])
    '''In拡散係数(z方向)'''
    





    '''熱移流拡散について'''
   

    log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step+1,1,0,0,0] = (2-2*A_ND_X*Dq)*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,1,0,0,0]
    + A_ND_X*Dq*( log[1:grid_X-2,2:grid_Y-1,2:grid_Z-1,Step,1,0,0,0] + log[3:grid_X,2:grid_Y-1,2:grid_Z-1,Step,1,0,0,0])
    '''熱拡散(X)'''
    -(DT/(2*mleng_X))*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,1,0,0]* (log[1:grid_X-2,2:grid_Y-1,2:grid_Z-1,Step,1,0,0,0]
       - log[3:grid_X,2:grid_Y-1,2:grid_Z-1,Step+1,1,0,0,0])
    '''熱移流(X)'''

    log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step+1,1,0,0,0] = (2-2 *A_ND_Y *Dq) *log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,1,0,0,0]
    + A_ND_Y*Dq*( log[2:grid_X-1,1:grid_Y-2,2:grid_Z-1,Step,1,0,0,0] + log[2:grid_X-1,3:grid_Y,2:grid_Z-1,Step,1,0,0,0])
    '''熱拡散(Y)'''
    -(DT/(2*mleng_X))*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,1,0,0]* (log[2:grid_X-1,1:grid_Y-2,2:grid_Z-1,Step,1,0,0,0]
       - log[2:grid_X-1,3:grid_Y,2:grid_Z-1,Step,1,0,0,0])
    '''熱移流(Y)'''

    log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step+1,1,0,0,0] = (2-2 *A_ND_Z *Dq) *log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,1,0,0,0]
    + A_ND_Z*Dq*( log[2:grid_X-1,2:grid_Y-1,1:grid_Z-2,Step+1,1,0,0,0] + log[2:grid_X-1,2:grid_Y-1,3:grid_Z,Step,1,0,0,0])
    '''熱拡散(Z)'''
    -(DT/(2*mleng_X))*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,1,0,0]* (log[2:grid_X-1,1:grid_Y-2,2:grid_Z-1,Step,1,0,0,0]
       - log[2:grid_X-1,2:grid_Y-1,3:grid_Z,Step,1,0,0,0])
    '''熱移流(Z)'''


    '''物質移流拡散について'''

    log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step+1,0,0,1,0] = (2-2 *A_ND_X *DInM_XY) *log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,0,1,0]
    + A_ND_X*DInM_XY*( log[1:grid_X-2,2:grid_Y-1,2:grid_Z-1,Step,0,0,1,0] + log[3:grid_X,2:grid_Y-1,2:grid_Z-1,Step,0,0,1,0])
    '''物質拡散(X)'''
    -(DT/(2*mleng_X))*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,0,0,1]* (log[1:grid_X-2,2:grid_Y-1,2:grid_Z-1,Step,0,0,1,0]
       - log[3:grid_X,2:grid_Y-1,2:grid_Z-1,Step+1,0,0,1,0])
    '''物質移流(X)'''

    log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step+1,0,0,1,0] = (2-2 *A_ND_X *DInM_XY) *log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,0,1,0]
    + A_ND_X*DInM_XY*( log[2:grid_X-1,1:grid_Y-2,2:grid_Z-1,Step,0,0,1,0] + log[2:grid_X-1,3:grid_Y,2:grid_Z-1,Step,0,0,1,0])
    '''物質拡散(Y)'''
    -(DT/(2*mleng_X))*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,0,0,1]* (log[2:grid_X-1,1:grid_Y-2,2:grid_Z-1,Step,0,0,1,0]
       - log[2:grid_X-1,3:grid_Y,2:grid_Z-1,Step+1,0,0,1,0])
    '''物質移流(Y)'''

    log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step+1,0,0,1,0] = (2-2 *A_ND_X *DInM_XY) *log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,0,1,0]
    + A_ND_X*DInM_XY*( log[2:grid_X-1,2:grid_Y-1,1:grid_Z-2,Step,0,0,1,0] + log[2:grid_X-1,2:grid_Y-1,3:grid_Z,Step,0,0,1,0])
    '''物質拡散(Z)'''
    -(DT/(2*mleng_X))*log[2:grid_X-1,2:grid_Y-1,2:grid_Z-1,Step,0,0,0,1]* (log[2:grid_X-1,2:grid_Y-1,1:grid_Z-2,Step,0,0,1,0]
       - log[2:grid_X-1,2:grid_Y-1,3:grid_Z,Step+1,0,0,1,0])
    '''物質移流(Z)'''


   