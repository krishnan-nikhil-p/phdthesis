import numpy as np
import glob
import os 
#import cv2

def full_out_traj(file,tTot,atoms):
    i=0
    ln =0
    cell_str = []
    cell_arr = np.zeros((tTot,atoms,4))
    with open (file) as myfile:

        for line in myfile:

            if "Atoms. Timestep:" in line:
                #print(ln)
                i+=1
                cell_str = []
                ln=0
            if "Atoms. Timestep:" not in line and str(atoms) not in line:
                cell_str.append(line)
                ln+=1


            if ln==atoms:
                cell_arr[i-1] = np.array([np.array(cell_str[i].split( )[:]).astype(float) for i in range(len(cell_str))])

    return cell_arr



def find_tTot_traj(file):
    i=0
    with open (file) as myfile:

        for line in myfile:
            if "Atoms. Timestep:" in line:
                i+=1

            
    return i 


def find_tot_atoms_traj(file):
    with open(file) as f:
        first_line = f.readline()


    return int(first_line)
    
    
    
totT =find_tTot_traj('traj.xyz')
totA =  find_tot_atoms_traj('traj.xyz')
arrT = full_out_traj('traj.xyz', totT, totA)


scales = np.array([10,25,50,100,150,200,300,500])

Hf_10 = np.zeros((len(arrT[::5]),80,10))
Hf_diff_10 = np.zeros((len(arrT[::5]),80,10))

Hf_25 = np.zeros((len(arrT[::5]),200,25))
Hf_diff_25 = np.zeros((len(arrT[::5]),200,25))

Hf_50 = np.zeros((len(arrT[::5]),400,50))
Hf_diff_50 = np.zeros((len(arrT[::5]),400,50))


Hf_100 = np.zeros((len(arrT[::5]),800,100))
Hf_diff_100 = np.zeros((len(arrT[::5]),800,100))

Hf_150 = np.zeros((len(arrT[::5]),1200,150))
Hf_diff_150 = np.zeros((len(arrT[::5]),1200,150))

Hf_200 = np.zeros((len(arrT[::5]),1600,200))
Hf_diff_200 = np.zeros((len(arrT[::5]),1600,200))

Hf_300 = np.zeros((len(arrT[::5]),2400,300))
Hf_diff_300 = np.zeros((len(arrT[::5]),2400,300))

Hf_500 = np.zeros((len(arrT[::5]),4000,500))
Hf_diff_500 = np.zeros((len(arrT[::5]),4000,500))




for a,arr in enumerate(arrT[::5]):
    xa = arr[arr[:,0]==1,1]
    ya = arr[arr[:,0]==1,2]

    xb = arr[arr[:,0]==2,1]
    yb = arr[arr[:,0]==2,2]

    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [80,10])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [80,10])

    Hf_10[a] =Ha/(Ha+Hb)
    Hf_diff_10[a] = np.abs(Ha-Hb) / (Ha+Hb)


    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [200,25])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [200,25])
    #Hf = np.zeros(Ha.shape)
    Hf_25[a] =Ha/(Ha+Hb)
    Hf_diff_25[a] = np.abs(Ha-Hb) / (Ha+Hb)


    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [400,50])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [400,50])
    Hf_50[a] =Ha/(Ha+Hb)
    Hf_diff_50[a] = np.abs(Ha-Hb) / (Ha+Hb)


    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [800,100])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [800,100])
    Hf_100[a] =Ha/(Ha+Hb)
    Hf_diff_100[a] = np.abs(Ha-Hb) / (Ha+Hb)


    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [1200,150])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [1200,150])
    Hf_150[a] =Ha/(Ha+Hb)
    Hf_diff_150[a] = np.abs(Ha-Hb) / (Ha+Hb)



    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [1600,200])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [1600,200])
    Hf_200[a] =Ha/(Ha+Hb)
    Hf_diff_200[a] = np.abs(Ha-Hb) / (Ha+Hb)

  


    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [2400,300])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [2400,300])
    Hf_300[a] =Ha/(Ha+Hb)
    Hf_diff_300[a] = np.abs(Ha-Hb) / (Ha+Hb)


    Ha, xedgesA, yedgesA = np.histogram2d(xa,ya,bins= [4000,500])
    Hb, xedgesB, yedgesB = np.histogram2d(xb,yb,bins= [4000,500])
    Hf_500[a] =Ha/(Ha+Hb)
    Hf_diff_500[a] = np.abs(Ha-Hb) / (Ha+Hb)


np.save('binArray10.npy',Hf_10)
np.save('binArrayDiff10.npy',Hf_diff_10)

np.save('binArray25.npy',Hf_25)
np.save('binArrayDiff25.npy',Hf_diff_25)

np.save('binArray50.npy',Hf_50)
np.save('binArrayDiff50.npy',Hf_diff_50)

np.save('binArray100.npy',Hf_100)
np.save('binArrayDiff100.npy',Hf_diff_100)

np.save('binArray150.npy',Hf_150)
np.save('binArrayDiff150.npy',Hf_diff_150)

np.save('binArray200.npy',Hf_200)
np.save('binArrayDiff200.npy',Hf_diff_200)

np.save('binArray300.npy',Hf_300)
np.save('binArrayDiff300.npy',Hf_diff_300)

np.save('binArray500.npy',Hf_500)
np.save('binArrayDiff500.npy',Hf_diff_500)




