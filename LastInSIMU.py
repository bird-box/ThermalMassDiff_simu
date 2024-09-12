import numpty as np


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



'''シュミレーション設定'''

grid_n=0
'''分割数(n*n*n)'''
SimuTime=0
'''シュミレーションする現実時間[s]'''
Dleng_X=0
Dleng_Y=0
Dleng_Z=0
'''各軸方向のメッシュ間の距離[cm]'''
DT=0
'''微小時間[s]'''



'''計算用変数'''

Tem_B=np.ones([grid_n,3])
InCon_B=np.ones([grid_n,3])
'''記録用：手前時刻における三次元の温度[K]/Inの濃度[N/cm3]'''

Tem_A=np.ones([grid_n,3])
InCon_A=np.ones([grid_n,3])
'''記録用：未来時刻における三次元の温度[K]/Inの濃度[N/cm3]'''



sTem_B_XB=np.ones([grid_n-2,3])
sTem_B_XA=np.ones([grid_n-2,3])

sInCon_B_XB=np.ones([grid_n-2,3])
sInCon_B_XA=np.ones([grid_n-2,3])
'''計算用：手前時間における X軸方向 1戻り/進み 温度[K] /Inの濃度[N/cm3]'''

sTem_B_YB=np.ones([grid_n-2,3])
sTem_B_YA=np.ones([grid_n-2,3])

sInCon_B_YB=np.ones([grid_n-2,3])
sInCon_B_YA=np.ones([grid_n-2,3])
'''計算用：手前時間における Y軸方向 1戻り/進み 温度[K] /Inの濃度[N/cm3]'''

sTem_B_ZB=np.ones([grid_n-2,3])
sTem_B_ZA=np.ones([grid_n-2,3])

sInCon_B_ZB=np.ones([grid_n-2,3])
sInCon_B_ZA=np.ones([grid_n-2,3])
'''計算用：手前時間における Z軸方向 1戻り/進み 温度[K] /Inの濃度'''



SaveTem=([grid_n,3,SimuTime])
SaveInCon=([grid_n,3,SimuTime])
'''保存用 温度/Inの濃度'''