
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

cimport numpy as np
cimport cython


@cython.cdivision(True) 

@cython.boundscheck(False)



cdef class Var_Lag:
	cdef np.double_t[6] Lx
	cdef np.double_t[6] Ly
	cdef np.double_t[6] Lz

	

cdef void  Velocity_to_grid(np.ndarray Grid,np.ndarray U,np.ndarray V,np.ndarray W, int n_x,int n_y,int n_z):
	cdef unsigned int i,j,k
	for i in range(n_x):
		for j in range(n_y):
			for k in range(n_z):
				Grid[i,j,k,0] = U[i,j,k]
				Grid[i,j,k,1] = V[i,j,k]
				Grid[i,j,k,2] = W[i,j,k]

	for j in range(n_y):
		for k in range(n_z):
			Grid[n_x,j,k,0] = U[0,j,k]
			Grid[n_x,j,k,1] = V[0,j,k]
			Grid[n_x,j,k,2] = W[0,j,k]


	for j in range(n_y):
		for k in range(n_x):
			Grid[k,j,n_z,0] = U[k,j,0]
			Grid[k,j,n_z,1] = V[k,j,0]
			Grid[k,j,n_z,2] = W[k,j,0]

	for j in range(n_y):
			Grid[n_x,j,n_z,0] = U[0,j,0]
			Grid[n_x,j,n_z,1] = V[0,j,0]
			Grid[n_x,j,n_z,2] = W[0,j,0]




cdef void Update_Vel(np.ndarray[double, ndim=3]  U,np.ndarray[double, ndim=3]  V,np.ndarray[double, ndim=3]  W,int case,int n_x,int n_y,int n_z,double dx,double dz,np.ndarray[double,ndim=1]	y_v,double t):			

	if case == 1:
		for i in range(n_x):
			for j in range(n_y):
				for k in range(n_z):	
					U[i,j,k] = 1. - (y_v[j]-1.)*	(y_v[j]-1.)
					V[i,j,k] = 0.
					W[i,j,k] = 0.5*lm.cos(0.25*i*dx)






				
cdef void init_Particles(int n_p,np.ndarray  X,np.ndarray  y_v,int n_x,int n_y,int n_z,int y0_init,double lx,double lz,double part_per_plane):
	cdef int c = 0
	cdef int i,j,k
	#while c < n_p:
	i=0
	j=0
	k=0
	c=0


	k = y0_init
	while  c < n_p:
		i=0
		while i <  part_per_plane and  c < n_p:
			j=0
			while j <  part_per_plane and  c < n_p:
				X[c,1,0] = y_v[k]
				X[c,0,0] = lx/( part_per_plane *1.1)*(i+1)
				X[c,2,0] = lz/( part_per_plane *1.1)*(j+1)
				c=c+1
				j=j+1
			i=i+1
		k=k+1


cdef double PolyLagrange(int i,double xhi) nogil:
	cdef double res = 0.



	if i== 1:
		res += 6.*xhi
		xhi *= xhi
		res += -5.0*xhi
		xhi *= xhi
		res += -5.0*xhi
		xhi *= xhi
		res += 5.0*xhi
		xhi *= xhi
		res += -xhi

		res = res / 120.


	if i == 2:
		res += -6.*xhi
		xhi *= xhi
		res += 16.*xhi
		xhi *= xhi
		res += -xhi
		xhi *= xhi
		res += -4.0*xhi
		xhi *= xhi
		res += xhi

		res = res / 24.

	if i == 3:
		res = 12.
		res += -4.*xhi
		xhi *= xhi
		res += -15.*xhi
		xhi *= xhi
		res += 5.*xhi
		xhi *= xhi
		res += 3.*xhi
		xhi *= xhi
		res += -xhi

		res = res / 12.

	if i == 4:
		res += 12.*xhi
		xhi *= xhi
		res += 8.*xhi
		xhi *= xhi
		res += -7.0*xhi
		xhi *= xhi
		res += -2.0*xhi
		xhi *= xhi
		res += 2.0*xhi

		res = res / 12.

	if i == 5:
		res += -6.*xhi
		xhi *= xhi
		res += -1.*xhi
		xhi *= xhi
		res += 7.*xhi
		xhi *= xhi
		res += xhi
		xhi *= xhi
		res += -1.0*xhi

		res = res / 24.

	if i == 6:
		res += 4.*xhi
		xhi *= xhi*xhi
		res += -5.*xhi
		xhi *= xhi*xhi
		res += xhi


		res = res / 120.

	else:
		res =0.

	return res
	

cdef double LagrangianInterp(double X0,double X1,double X2,np.ndarray[double,ndim=4]  Grid,double dx,np.ndarray[double,ndim=1]   y_v,double dz,int nx, int ny,int nz,int index):

	cdef double res = 0.

	cdef np.ndarray[double,ndim=1] Lx =  np.zeros(6)
	cdef np.ndarray[double,ndim=1] Ly =  np.zeros(6)
	cdef np.ndarray[double,ndim=1] Lz =  np.zeros(6)

	cdef double y_1,y_2,y_3,y_4,y_5,y_6

	cdef double xhi_x = (X0-np.floor(X0/dx))/dx
	cdef double xhi_z = (X2-np.floor(X2/dz))/dz

	cdef int x_index,z_index

	
	

	cdef int find = 0
	cdef int k = 0
	if X1 >= y_v[ny-1]:
		k = ny
	else:		
		while find == 0:
			if X1 >= y_v[k] and  X1 <= y_v[k+1]:
				find += 1
			k += 1

	# Si il y a assez de point pour gerer le mur
	if k >= 2:

		y_1 = y_v[k-2]
		y_2 = y_v[k-1]
		y_3 = y_v[k]
		y_4 = y_v[k+1]
		y_5 = y_v[k+2]
		y_6 = y_v[k+3]


		Lx[0] = PolyLagrange(0,xhi_x)
		Lx[1] = PolyLagrange(1,xhi_x)
		Lx[2] = PolyLagrange(2,xhi_x)
		Lx[3] = PolyLagrange(3,xhi_x)
		Lx[4] = PolyLagrange(4,xhi_x)
		Lx[5] = PolyLagrange(5,xhi_x)


		Lz[0] = PolyLagrange(0,xhi_z)
		Lz[1] = PolyLagrange(1,xhi_z)
		Lz[2] = PolyLagrange(2,xhi_z)
		Lz[3] = PolyLagrange(3,xhi_z)
		Lz[4] = PolyLagrange(4,xhi_z)
		Lz[5] = PolyLagrange(5,xhi_z)

		Ly[0] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_1-y_2)*(y_1-y_3)*(y_1-y_4)*(y_1-y_5)*(y_1-y_6))
		Ly[1] = (X1-y_1)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_2-y_1)*(y_2-y_3)*(y_2-y_4)*(y_2-y_5)*(y_2-y_6))
		Ly[2] = (X1-y_2)*(X1-y_1)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_3-y_2)*(y_3-y_1)*(y_3-y_4)*(y_3-y_5)*(y_3-y_6))
		Ly[3] = (X1-y_2)*(X1-y_3)*(X1-y_1)*(X1-y_5)*(X1-y_6)/((y_4-y_2)*(y_4-y_3)*(y_4-y_1)*(y_4-y_5)*(y_4-y_6))
		Ly[4] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_1)*(X1-y_6)/((y_5-y_2)*(y_5-y_3)*(y_5-y_4)*(y_5-y_1)*(y_5-y_6))
		Ly[5] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_1)/((y_6-y_2)*(y_6-y_3)*(y_6-y_4)*(y_6-y_5)*(y_6-y_1))

		for i in xrange(6):
			x_index = np.mod(np.floor(X0/dx),nx+1)
			for j in xrange(6):
				z_index = np.mod(np.floor(X2/dz),nz+1)
				for l in xrange(6):
					res +=  Grid[x_index,k+l-2,z_index,index]*Lx[i]*Ly[l]*Lz[j]

	else:
		res = TrilinearInterp(X0,X1,X2,Grid,dx,y_v,dz,nx,ny,nx,index)


	return res


cdef double TrilinearInterp(double X0,double X1,double X2,np.ndarray[double,ndim=4]  Grid,double dx,np.ndarray[double,ndim=1]   y_v,double dz,int nx, int ny,int nz,int index):
	cdef double  res =0

	# Reperer la particule

	cdef int x_0 = np.floor(X0/dx)
	cdef int x_1 = (np.floor(X0/dx))+1


	cdef int z_0 = np.floor(X2/dz)
	cdef int z_1 = (np.floor(X2/dz))+1

	cdef int find = 0
	cdef int k = 0
	if X1 >= y_v[y_v.size-1]:
		k = y_v.size
	else:		
		while find == 0:
			if X1 >= y_v[k] and  X1 <= y_v[k+1]:
				find += 1
			k += 1

	cdef int y_1 =  k 
	cdef int y_0 =  k-1

	#print x_0,x_1,z_0,z_1,y_0,y_1,X1

	cdef double xd = (X0-x_0*dx)/(x_1*dx-x_0*dx)
	cdef double yd = (X1-y_v[y_0])/(y_v[y_1]-y_v[y_0])
	cdef double zd = (X2-z_0*dz)/(z_1*dz-z_0*dz)

	cdef double c00 = Grid[x_0,y_0,z_0,index]*(1.-xd) + Grid[x_1,y_0,z_0, index]*(xd) 
	cdef double  c01 = Grid[x_0,y_0,z_1,index]*(1.-xd) + Grid[x_1,y_0,z_1, index]*(xd) 
	cdef double c10 = Grid[x_0,y_1,z_0,index]*(1.-xd) + Grid[x_1,y_1,z_0, index]*(xd) 
	cdef double  c11 = Grid[x_0,y_1,z_1,index]*(1.-xd) + Grid[x_1,y_1,z_1, index]*(xd) 

	cdef double c0 = c00*(1.-yd) + c10*yd
	cdef double c1 = c01*(1.-yd) + c11*yd

	res = c0*(1.-zd) + c1*zd

	

	#print Grid[x_0,y_0,z_0,index],Grid[x_1,y_0,z_0,index],Grid[x_0,y_1,z_0,index],Grid[x_1,y_0,z_0,index],res
	return res

cpdef void Lagrangian(FILE_P):

	file_param = open(FILE_P,'r')

	line = file_param.readline()
	cdef long Nb_particles = int(line)
	print Nb_particles
	line = file_param.readline()
	cdef long Nb_part_per_plane = int(line)
	print Nb_part_per_plane
	line = file_param.readline()
	cdef double tf = Decimal(line)
	print tf
	line = file_param.readline()
	cdef int y0_init = int(line)
	print y0_init
	line = file_param.readline()
	cdef double dt = Decimal(line)
	print dt
	line = file_param.readline()
	cdef double dt_data = Decimal(line)
	print dt_data
	line = file_param.readline()
	cdef double nb_iter_win = Decimal(line)
	print nb_iter_win
	line = file_param.readline()
	cdef double lx = Decimal(line)
	print lx
	line = file_param.readline()
	cdef double ly = Decimal(line)
	print ly
	line = file_param.readline()
	cdef double lz = Decimal(line)
	print lz

	line = file_param.readline()
	cdef int n_x = int(line)
	print n_x
	line = file_param.readline()
	cdef int n_y = int(line)
	print n_y
	line = file_param.readline()
	cdef int n_z = Decimal(line)
	print n_z
	line = file_param.readline()
	cdef FILE_Vitesses = line[:-1]

	print FILE_Vitesses

	line = file_param.readline()
	cdef DIR_Particle_out = line[:-1]

	line = file_param.readline()
	cdef DIR_Vel_out = line[:-1]
	
	cdef double t0 = 0



	test = 1

	cdef long decalage = 300
#	#cdef double tf = 1
#	cdef double dt = 0.05
#	cdef double dt_data = 0.05

#	cdef double nb_iter_win = 1

#	cdef double lx = 4.*lm.pi
#	cdef double ly = 2.
#	cdef double lz = 4.*lm.pi/3.

#	cdef int n_x = 128
#	cdef int n_y = 129
#	cdef int n_z = 128

	cdef double epsilon = 0.024

	cdef int i = 0


	cdef np.ndarray[double, ndim=4]  Grid_t =  np.zeros((n_x+1,n_y,n_z+1,3))
	cdef np.ndarray[double, ndim=4] Grid_t1 =  np.zeros((n_x+1,n_y,n_z+1,3))
	cdef np.ndarray[double, ndim=4] Grid_t2 =  np.zeros((n_x+1,n_y,n_z+1,3))


	cdef np.ndarray[double, ndim=4] Grid_current_t = np.zeros((n_x+1,n_y,n_z+1,3))
	cdef np.ndarray[double, ndim=4] Grid_next_t = np.zeros((n_x+1,n_y,n_z+1,3))


	FILE_Uinst  = FILE_Vitesses+'/ux'
	FILE_Vinst  = FILE_Vitesses+'/uy'
	FILE_Winst  = FILE_Vitesses+'/uz'

	cdef np.double_t[6][6] XXX

	print FILE_Uinst

	FILE_Umean  = '/Volumes/MyBook/Channel_flow/channel/umean.dat'
	FILE_Vmean  = '/Volumes/MyBook/Channel_flow/channel/vmean.dat'
	FILE_Wmean  = '/Volumes/MyBook/Channel_flow/channel/wmean.dat'


	cdef np.ndarray[double, ndim=3] History_Particles = np.zeros((Nb_particles,3,(tf-t0)/dt+2))
	cdef np.ndarray[double, ndim=3] History_Vel_Particles = np.zeros((Nb_particles,3,(tf-t0)/dt+2))

	cdef double dx = lx/(n_x)
	cdef double dz = lz/(n_z)

	cdef np.ndarray[double, ndim=1] x_values = np.linspace(0,lx,n_x+1)
	cdef np.ndarray[double, ndim=1]  z_values = np.linspace(0,lz,n_z+1)
	cdef np.ndarray[double, ndim=1]  y_values= np.zeros(n_y)


	

	# Lecture des points


	file_y = open(FILE_Vitesses+"/yp.dat",'r')
	for i in range(n_y):
		y_values[i] = float(file_y.readline())

	# Grille de calcul

	Grid = np.meshgrid(x_values, y_values, z_values, indexing='ij')
	# Initialisation des premiers champs de vitesse
	cdef np.ndarray[double, ndim=3]  U = np.zeros((n_x,n_y,n_z))
	cdef np.ndarray[double, ndim=3]  V= np.zeros((n_x,n_y,n_z))
	cdef np.ndarray[double, ndim=3]  W= np.zeros((n_x,n_y,n_z))
	cdef int current_index_t = 0

	print int(lm.sqrt(Nb_part_per_plane))

	init_Particles(Nb_particles,History_Particles,y_values,n_x,n_y,n_z,y0_init,lx,lz,int(lm.sqrt(Nb_part_per_plane)))


	cdef double current_time = t0
	cdef double ratio_time = 0
	cdef np.ndarray[double, ndim=1]  k1 = np.zeros(3)
	cdef np.ndarray[double, ndim=1]  k2 = np.zeros(3)
	cdef np.ndarray[double, ndim=1]  k3 = np.zeros(3)
	cdef np.ndarray[double, ndim=1]  k4 = np.zeros(3)
	cdef np.ndarray[double, ndim=1]  yn = np.zeros(3)
	cdef double iter_t = 0

	cdef double dt2 = 0.5*dt
	cdef double dt6 = dt/6.0

	cdef double time_clock = time.time()
	cdef double t2

	cdef double y_00,y_01,y_02


	current_index_t += 1

	if test == 0:

		U = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Uinst+str(current_index_t+decalage).zfill(5))
		V = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Vinst+str(current_index_t+decalage).zfill(5))
		W = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Winst+str(current_index_t+decalage).zfill(5))
	else:
		Update_Vel( U,V,W,1,n_x,n_y,n_z,dx,dz,y_values,0)		

			
	Velocity_to_grid(Grid_current_t,U,V,W,n_x,n_y,n_z)

	if test == 0:
		U = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Uinst+str(current_index_t+decalage+1).zfill(5))
		V = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Vinst+str(current_index_t+decalage+1).zfill(5))
		W = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Winst+str(current_index_t+decalage+1).zfill(5))
	else:
		Update_Vel( U,V,W,1,n_x,n_y,n_z,dx,dz,y_values,dt)

	Velocity_to_grid(Grid_next_t,U,V,W,n_x,n_y,n_z)

	print "Nb_p: "+Nb_particles.__str__()
	while current_time < tf-dt:
		t2 = time.time() - time_clock 
		print "Iteration ",iter_t.__str__()," t= ",current_time.__str__()," Temps de calcul: ",t2.__str__(),"s",np.mod(current_time+epsilon,dt_data) , np.mod(current_time,dt_data) < dt


		
		ratio_time_0 = np.mod(iter_t,nb_iter_win)/nb_iter_win


		print ratio_time_0
		ratio_time_1 = np.mod(iter_t+0.5,nb_iter_win)/nb_iter_win


		Grid_t = (1.-ratio_time_0)*Grid_current_t + ratio_time_0*Grid_next_t
		Grid_t1 = (1.-ratio_time_1)*Grid_current_t + ratio_time_1*Grid_next_t

		
		# Next temporal window
		if np.mod(iter_t+1,nb_iter_win) == 0:
		
			current_index_t=  current_index_t + 1
			print "Nouvelle fenetre temporelle ",current_index_t.__str__()
			Grid_current_t[:,:,:,:] = Grid_next_t[:,:,:,:]

			if test == 0:

				U = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Uinst+str(current_index_t+decalage+1).zfill(5))
				V = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Vinst+str(current_index_t+decalage+1).zfill(5))
				W = LOADFLOW.loadflow(n_x,n_y,n_z,FILE_Winst+str(current_index_t+decalage+1).zfill(5))
			else:
				Update_Vel( U,V,W,1,n_x, n_y,n_z,dx,dz,y_values,current_time)

			Velocity_to_grid(Grid_next_t,U,V,W,n_x,n_y,n_z)


		ratio_time_2 = np.mod(iter_t+1,nb_iter_win)/nb_iter_win

		print ratio_time_0,ratio_time_1,ratio_time_2


		

		Grid_t2 = (1.-ratio_time_2)*Grid_current_t + ratio_time_2*Grid_next_t

		# Iteration RK-4 pour chaue points

		for i in range(Nb_particles):


			y_00  = History_Particles[i,0,iter_t]
			y_01 = History_Particles[i,1,iter_t]
			y_02 = History_Particles[i,2,iter_t]


			k1[0] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t,dx,y_values,dz,n_x,n_y,n_z,0)
			k1[1] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t,dx,y_values,dz,n_x,n_y,n_z,1)
			k1[2] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t,dx,y_values,dz,n_x,n_y,n_z,2)



			k2[0] =  LagrangianInterp((y_00+dt2*k1[0])%lx,max(0,y_01+dt2*k1[1]),(y_02+dt2*k1[2])%lz,Grid_t1,dx,y_values,dz,n_x,n_y,n_z,0)
			k2[1] =  LagrangianInterp((y_00+dt2*k1[0])%lx,max(0,y_01+dt2*k1[1]),(y_02+dt2*k1[2])%lz,Grid_t1,dx,y_values,dz,n_x,n_y,n_z,1)
			k2[2] =  LagrangianInterp((y_00+dt2*k1[0])%lx,max(0,y_01+dt2*k1[1]),(y_02+dt2*k1[2])%lz,Grid_t1,dx,y_values,dz,n_x,n_y,n_z,2)

			k3[0] =  LagrangianInterp((y_00+dt2*k2[0])%lx,max(0,y_01+dt2*k2[1]),(y_02+dt2*k2[2])%lz,Grid_t1,dx,y_values,dz,n_x,n_y,n_z,0)
			k3[1] =  LagrangianInterp((y_00+dt2*k2[0])%lx,max(0,y_01+dt2*k2[1]),(y_02+dt2*k2[2])%lz,Grid_t1,dx,y_values,dz,n_x,n_y,n_z,1)
			k3[2] =  LagrangianInterp((y_00+dt2*k2[0])%lx,max(0,y_01+dt2*k2[1]),(y_02+dt2*k2[2])%lz,Grid_t1,dx,y_values,dz,n_x,n_y,n_z,2)



			k4[0] =  LagrangianInterp((y_00+dt*k3[0])%lx,max(0,y_01+dt*k3[1]),(y_02+dt*k3[2])%lz,Grid_t2,dx,y_values,dz,n_x,n_y,n_z,0)
			k4[1] =  LagrangianInterp((y_00+dt*k3[0])%lx,max(0,y_01+dt*k3[1]),(y_02+dt*k3[2])%lz,Grid_t2,dx,y_values,dz,n_x,n_y,n_z,1)
			k4[2] =  LagrangianInterp((y_00+dt*k3[0])%lx,max(0,y_01+dt*k3[1]),(y_02+dt*k3[2])%lz,Grid_t2,dx,y_values,dz,n_x,n_y,n_z,2)

		



			History_Vel_Particles[i,0,iter_t] =(k1[0])
			History_Vel_Particles[i,1,iter_t] =(k1[1])
			History_Vel_Particles[i,2,iter_t] =(k1[2])
	
			History_Particles[i,0,iter_t+1] = y_00 + dt6*(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0])
			History_Particles[i,1,iter_t+1] = y_01 + dt6*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])
			History_Particles[i,2,iter_t+1] = y_02 + dt6*(k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2])

			if History_Particles[i,1,iter_t+1] < 0:
				History_Particles[i,1,iter_t+1] = 0



			if History_Particles[i,0,iter_t+1] >= lx or History_Particles[i,0,iter_t+1] < 0: 
				History_Particles[i,0,iter_t+1] = np.mod(History_Particles[i,0,iter_t+1],lx)

			if History_Particles[i,2,iter_t+1] > lz or History_Particles[i,2,iter_t+1] < 0: 
				History_Particles[i,2,iter_t+1] = np.mod(History_Particles[i,2,iter_t+1],lz)


		# Fin de la modification spatiale

		current_time += dt
		iter_t += 1

		
	# Fin de la dynamique

		y_00  = History_Particles[i,0,iter_t]
		y_01 = History_Particles[i,1,iter_t]
		y_02 = History_Particles[i,2,iter_t]

			

		k1[0] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t,dx,y_values,dz,n_x,n_y,n_z,0)
		k1[1] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t,dx,y_values,dz,n_x,n_y,n_z,1)
		k1[2] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t,dx,y_values,dz,n_x,n_y,n_z,2)




		History_Vel_Particles[i,0,iter_t] = k1[0]
		History_Vel_Particles[i,1,iter_t] = k1[1]
		History_Vel_Particles[i,2,iter_t] = k1[2]



	# Fichiers de sorties



	k = y0_init
	FILE_RES = open(DIR_Particle_out+"/y_"+k.__str__()+"_P"+".dat","w")
	for i in range(Nb_particles):
		if np.mod(i,Nb_part_per_plane) == 0 and i > 0:
			FILE_RES.close()
			k = k+1
			print np.mod(i,Nb_part_per_plane)

		FILE_RES = open(DIR_Particle_out+"/y_"+(k).__str__()+"_P"+(np.mod(i,Nb_part_per_plane)+1).__str__().zfill(4)+".dat","w")
		for t in range(int(iter_t)):	
			FILE_RES.write(History_Particles[i,0,t].__str__()+" ")
			FILE_RES.write(History_Particles[i,1,t].__str__()+" ")
			FILE_RES.write(History_Particles[i,2,t].__str__()+" ")
			FILE_RES.write("\n")
		
	k = y0_init
	FILE_RES = open(DIR_Vel_out+"/y_"+k.__str__()+"_P"+".dat","w")
	for i in range(Nb_particles):
		if np.mod(i,Nb_part_per_plane) == 0 and i > 0:
			FILE_RES.close()
			k = k+1
			print np.mod(i,Nb_part_per_plane)

		FILE_RES = open(DIR_Vel_out+"/y_"+(k).__str__()+"_P"+(np.mod(i,Nb_part_per_plane)+1).__str__().zfill(4)+".dat","w")
		for t in range(int(iter_t)):	
			FILE_RES.write(History_Vel_Particles[i,0,t].__str__()+" ")
			FILE_RES.write(History_Vel_Particles[i,1,t].__str__()+" ")
			FILE_RES.write(History_Vel_Particles[i,2,t].__str__()+" ")
			FILE_RES.write("\n")

	

