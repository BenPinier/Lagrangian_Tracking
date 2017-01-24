# -*-coding:Latin-1 -* 
#usr/bash/python
import numpy as np


from decimal import Decimal
import math as lm
import LOADFLOW


import time
 

from scipy.interpolate import RegularGridInterpolator
import os
import shutil
import sys
import random 
from scipy.interpolate import interp1d

cimport numpy as np
cimport cython

@cython.boundscheck(False)


cpdef Correlation_Vel_lag(long Nb_particles,double tf,int y0_init,int N_t,int N_y, str DIR_Vel_out, str DIR_Particles_out,str DIR_RES):
	cdef int i,j,a,t,count,cu_index_y
	cdef str line, ftype,line_p,line_v
	cdef double t0 = 0
	cdef double dt = 0.025
	cdef double dt_data = 0.05
	cdef double pos,v0,vt,cu_mean

	cdef double c_u = 0
	cdef double tmpu0,tmpv0,tmpw0
	cdef double tmpu1,tmpv1,tmpw1

	cdef double kappa = 0.4


	#T_0 = [0,400,800,1200]
	#cdef int size_T0 = 4

	T_0 = [12500]
	cdef int size_T0 = 1

#	T_0 = [0,400,800,1200]#
#	cdef int size_T0 = 4
	cdef int Nb_part_tot = size_T0*Nb_particles

	cdef list_line
#
#	cdef int N_t = 400
#	cdef int N_y
	cdef int N_yint 

#	N_y = 192/8+2

	cdef np.ndarray[double, ndim =1] y_values = np.zeros(N_y)
	cdef np.ndarray[double, ndim =1] nb_part_y= np.zeros(N_y)



	N_yint = 360


	cdef np.ndarray[double, ndim =1] y_interp = np.zeros(N_yint)


#	if N_y == 129:#
#		for i in range(N_yint):
#			y_interp[i] = (i*1.)/180.
			

#	if N_y == 192:
#		for i in range(N_yint):
#			y_interp[i] = (i*1.)/3600



	for i in range(N_yint):
		y_interp[i] = (i*1.)/N_yint
			

	cdef np.ndarray[double, ndim =2] vu = np.zeros((Nb_part_tot,N_t))
	cdef np.ndarray[double, ndim =2] vv = np.zeros((Nb_part_tot,N_t))
	cdef np.ndarray[double, ndim =2] vw = np.zeros((Nb_part_tot,N_t))

						
	cdef np.ndarray[double, ndim =1] MoyenneU = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] MoyenneV = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] MoyenneW = np.zeros(N_t)


	cdef np.ndarray[double, ndim =1] Ruu_num = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] Rvv_num = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] Rww_num = np.zeros(N_t) 

	cdef np.ndarray[double, ndim =1] Ruu_denum = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] Rvv_denum = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] Rww_denum = np.zeros(N_t)


	cdef np.ndarray[double, ndim =1] R_uu= np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] R_vv = np.zeros(N_t)
	cdef np.ndarray[double, ndim =1] R_ww = np.zeros(N_t)
	
	cdef object FILE_V,FILE_P,file_y 
#	cdef str DIR_Vel_out = 'Results/Velocities/'
#	cdef str DIR_Particles_out = 'Results/Tracking/'

	cdef double epsilon = 1e-15



	cdef np.ndarray[double, ndim = 1] c_line = np.zeros(3)
	cdef np.ndarray[double, ndim = 1] c_line2 = np.zeros(3)
	cdef int index_y=y0_init

	cdef int find,k

	cdef double ustar  = 0.035

	#file_y = open("yp.dat",'r')
	file_y = open("yp4.dat",'r')
#	file_y = open("yp1000.dat",'r')
	for i in range(N_y):
		y_values[i] = float(file_y.readline())
		print i,y_values[i]


# --------  Evaluation de la RMS eulerienne ------------

	cdef np.ndarray[double, ndim=1] sigma_u1 = np.zeros(N_y)	
	cdef np.ndarray[double, ndim=1] sigma_v1 = np.zeros(N_y)	
	cdef np.ndarray[double, ndim=1] sigma_w1 = np.zeros(N_y)	
	cdef np.ndarray[double, ndim=1] U_mean1 = np.zeros(N_y)	

	cdef np.ndarray[double, ndim=1] sigma_u = np.zeros(N_yint)	
	cdef np.ndarray[double, ndim=1] sigma_v = np.zeros(N_yint)	
	cdef np.ndarray[double, ndim=1] sigma_w = np.zeros(N_yint)	
	cdef np.ndarray[double, ndim=1] U_mean = np.zeros(N_yint)	

	if 0:
		FILE_SIG = open('upup.dat')

		for i in range(65):
			line = FILE_SIG.readline()
			sigma_u1[i] = np.sqrt(Decimal(line))

		for i in range(65):
			sigma_u1[65+i] = sigma_u1[63-i]
		
		FILE_SIG.close()

		FILE_SIG = open('vpvp.dat')

		for i in range(65):
			line = FILE_SIG.readline()
			sigma_v1[i] = np.sqrt(Decimal(line))

		for i in range(65):
			sigma_v1[65+i] = sigma_v1[63-i]

		FILE_SIG.close()

		FILE_SIG = open('wpwp.dat')

		for i in range(65):
			line = FILE_SIG.readline()
			sigma_w1[i] = np.sqrt(Decimal(line))

		for i in range(65):
			sigma_w1[65+i] = sigma_w1[63-i]

		FILE_SIG.close()


		# Evaluation de la vitesse moyenne

		FILE_SIG = open('uplus.dat')

		for i in range(65):
			line = FILE_SIG.readline()
			U_mean1[i] = Decimal(line.split()[2])

		for i in range(65):
			U_mean1[65+i] = U_mean1[63-i]

		FILE_SIG.close()

	cdef int Re = 1000
	cdef double yplus
	for i in range(N_y):
		yplus = y_values[i]*Re

		yplus = Re*(min(y_values[i],2.-y_values[i]))

		if yplus < 30:
			U_mean1[i] = ustar*(2.5*lm.log(1.0+0.4*yplus)+7.8*(1.0-lm.exp(-yplus/11.0) - yplus/11.0*lm.exp(-0.33*yplus)));
		else:
			U_mean1[i] = ustar*(1.0/kappa*lm.log(yplus)+5.5);


# Interpolation des moyennes

#	f = interp1d(y_values,sigma_u1,kind='linear')
#	sigma_u = f(y_interp)
	
#	f = interp1d(y_values,sigma_v1,kind='linear')
#	sigma_v = f(y_interp)

#	f = interp1d(y_values,sigma_w1,kind='linear')
#	sigma_w = f(y_interp)

	f = interp1d(y_values,U_mean1,kind='linear')
	U_mean = f(y_interp)



	if 1:
		for i in range(N_y):
			print i,y_values[i],U_mean1[i]
		for i in range(N_yint):
			print i,y_interp[i],U_mean[i]



# --------------------  Evaluation des moyennes ------------------------

#	sys.exit()
	count = 0


	for i in range(Nb_part_tot-1):

		if np.mod(i,100) ==0:
			print i	
	#		print vu[i-1,0] 
		FILE_V = open(DIR_Vel_out+"y_"+y0_init.__str__()+"_P"+(np.mod(i,Nb_particles)+1).__str__().zfill(4)+"_t0_"+T_0[(int)(i/Nb_particles)].__str__()+".dat")
		FILE_P = open(DIR_Particles_out+"y_"+y0_init.__str__()+"_P"+(np.mod(i,Nb_particles)+1).__str__().zfill(4)+"_t0_"+T_0[(int)(i/Nb_particles)].__str__()+".dat")

		t=0
		for t in range(N_t):
			line_p = (FILE_P.readline())
			line_v  = (FILE_V.readline())
#			print line_p,t

			pos = Decimal(line_p.split()[1])

			cu_index_y = int(np.ceil(pos*1000.))
			cu_mean = U_mean[cu_index_y]

			vu[i,t] = Decimal(line_v.split()[0])
#			vu[i,t] -= cu_mean
#			vu[i,t] /= sigma_u[cu_index_y]
			vv[i,t] = Decimal(line_v.split()[1])
#			vv[i,t] /= sigma_v[cu_index_y]
			vw[i,t] = Decimal(line_v.split()[2])
#			vw[i,t] /= sigma_w[cu_index_y]
			
#			t = t + 1

		FILE_V.close()
		FILE_P.close()

	print "Fin des lectures"

	for t in range(N_t):
#		print vu[0,0];
		MoyenneU[t] = np.average(vu[0:1024,t])
		MoyenneV[t] = np.average(vv[0:1024,t])
		MoyenneW[t] = np.average(vw[0:1024,t])

#		print MoyenneU[t],t

#	for t in range(N_t,0,-1):
#		MoyenneU[t-1] = np.average(MoyenneU[0:t])
#		MoyenneV[t-1] = np.average(MoyenneV[0:t])
#		MoyenneW[t-1] = np.average(MoyenneW[0:t])
#		print MoyenneU[t-1],t-1
		


# --------------------  Evaluation de la Correlation ------------------------

#	sys.exit()
	count = 0
	t=0

	for i in range(0,Nb_part_tot-1):

		if np.mod(i,100) ==0:
			print i	

		tmpu0 = vu[i,0]-MoyenneU[0]
		tmpv0 = vv[i,0]-MoyenneV[0]
		tmpw0 = vw[i,0]-MoyenneW[0]

		

		Ruu_denum[0] +=  tmpu0*tmpu0 
		Rvv_denum[0] +=  tmpv0*tmpv0
		Rww_denum[0] +=  tmpw0*tmpw0


		Ruu_num[0] +=  tmpu0*tmpu0 
		Rvv_num[0] +=  tmpv0*tmpv0
		Rww_num[0] +=  tmpw0*tmpw0

	#	print tmpv0

		t = 1
		j=1
		for t in range(1,N_t):


			tmpu1 = vu[i,t]-MoyenneU[t]
			tmpv1 = vv[i,t]-MoyenneV[t]
			tmpw1 = vw[i,t]-MoyenneW[t]

			Ruu_denum[t] +=  tmpu1*tmpu1
			Rvv_denum[t] +=  tmpv1*tmpv1
			Rww_denum[t] +=  tmpw1*tmpw1

			Ruu_num[t] +=  tmpu0*tmpu1
			Rvv_num[t] +=  tmpv0*tmpv1
			Rww_num[t] +=  tmpw0*tmpw1	

			#	print 'ok'
			j+=1
				
			t = t + 1
			count += 1

		#print 'fin'

		print i
	#	FILE_P.close()

	
	Ruu_denum[:] = np.sqrt(Ruu_denum[:]/count)
	Ruu_num[:] = (Ruu_num[:])/count

	Rvv_denum[:] = np.sqrt(Rvv_denum[:]/count)
	Rvv_num[:] = (Rvv_num[:])/count

	Rww_denum[:] = np.sqrt(Rww_denum[:]/count)
	Rww_num[:] = Rww_num[:]/count



		
	for t in range(N_t):

		R_uu[t] =  Ruu_num[t]/(Ruu_denum[0]*Ruu_denum[t]+epsilon)
		R_vv[t] =  Rvv_num[t]/(Rvv_denum[0]*Rvv_denum[t]+epsilon)
		R_ww[t] =  Rww_num[t]/(Rww_denum[0]*Rww_denum[t]+epsilon)

	#	print R_vv[index_y,0] 


	#for index_y in range(N_y):

	FILE_RES = open(DIR_RES+"Correlation_Lagrangienne_y_"+y0_init.__str__()+".dat","w")
	for t in range(N_t):
		print t
		if np.mod(t,1) == 0:
			FILE_RES.write(R_uu[t].__str__()+" "+R_vv[t].__str__()+" "+R_ww[t].__str__()+"\n")
	FILE_RES.close()

		
