#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <ctime>
#include <string> 
// VARIABLES GLOBALES


#define ReTau 390
#define SCHM 5

double Lx[6];
double Ly[6];
double Lz[6];


double y_1,y_2,y_3,y_4,y_5,y_6;
double x_0,x_1,x_2,x_3,x_4,x_5,z_0,z_1,z_2,z_3,z_4,z_5;
int x0,x1,z0,z1,y_a,y_b,find;
int ix_1,iz_1;

double xhi_x;
double xhi_z;

double xd,yd,zd,c00,c01,c10,c11,c0,c1;

unsigned int i_index,j_index,ij_index;
unsigned int i_index2,j_index2,ij_index2;
int nynz,it ;



#include <sstream>

 template < typename T, typename S > inline T positive_modulo(T i, S n) {
    return fmod(fmod(i,n)+n,n);
}

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
	std::ostringstream stm ;
	stm << n ;
	return stm.str() ;
    }
}

// Lecture du fichier de donnees

void Chargement_parametres(std::string& FILE_arg,int* Nb_particules,int* Nb_particules_par_ligne,double* tf,int* y0_init,double* dt,double* dt_data,int* nb_iter_win,double* lx,double* ly,double* lz,int* nx,int* ny,int* nz,std::string& FILE_VITESSE,std::string& DIR_Particules_output,std::string& DIR_Vitesses_output){

	std::ifstream file(FILE_arg.c_str());
	
	std::string str;

	std::getline(file, str);
	*Nb_particules = atoi(str.c_str());

	std::getline(file, str);
	*Nb_particules_par_ligne = atoi(str.c_str());

	std::getline(file, str);
	*tf = atof(str.c_str());

	std::getline(file, str);
	*y0_init = atoi(str.c_str());

	std::getline(file, str);
	*dt = atof(str.c_str());

	std::getline(file, str);
	*dt_data = atof(str.c_str());

	std::getline(file, str);
	*nb_iter_win = atoi(str.c_str());

	std::getline(file, str);
	*lx = atof(str.c_str());

	std::getline(file, str);
	*ly = atof(str.c_str());

	std::getline(file, str);
	*lz = atof(str.c_str());

	std::getline(file, str);
	*nx = atoi(str.c_str());

	std::getline(file, str);
	*ny = atoi(str.c_str());

	std::getline(file, str);
	*nz = atoi(str.c_str());

	std::getline(file, FILE_VITESSE);

	std::getline(file, DIR_Particules_output);

	std::getline(file, DIR_Vitesses_output);

}



// Lecture de y

void Lecture_Y(double*Y,int ny,std::string& file_str){
	std::ifstream file(file_str.c_str());
	std::string str;
	for(int i=0;i<=ny;i++){
		std::getline(file, str);
		Y[i]=atof(str.c_str());
		//std::cout << str << " " << Y[i] <<std::endl;
	}

}


void Affiche_Param(int* Nb_particules,int* Nb_particules_par_ligne,double* tf,int* y0_init,double* dt,double* dt_data,int* nb_iter_win,double* lx,double* ly,double* lz,int* nx,int* ny,int* nz,std::string& FILE_VITESSE,std::string& DIR_Particules_output,std::string& DIR_Vitesses_output){

	std::cout << "Affichage des parametres" << std::endl;
	
	std::cout << *Nb_particules << std::endl;
	std::cout << *Nb_particules_par_ligne  << std::endl;
	std::cout << *tf << std::endl;
	std::cout << *y0_init  << std::endl;
	std::cout << *dt << std::endl;
	std::cout << *dt_data << std::endl;
	std::cout << *nb_iter_win << std::endl;
	std::cout << *lx << std::endl;
	std::cout << *ly << std::endl;
	std::cout << *lz << std::endl;
	std::cout << *nx << std::endl;
	std::cout << *ny << std::endl;
	std::cout << *nz << std::endl;

	std::cout << FILE_VITESSE << std::endl;
	std::cout << DIR_Particules_output << std::endl;
	std::cout << DIR_Vitesses_output << std::endl;

}

// Lecture des champs de vitesses depuis f90

void loadflow_1000(double* U,double* V,double* W,std::string file_str,int nx,int ny,int nz){

	std::ifstream ifs;

	ifs.open(file_str.c_str(),std::ios::in);

	
	for(int k=0;k<nx;k++){
			for(int j=0;j<ny;j++){
				for(int i=0;i<nz;i++){
					ifs >> U[k*ny*nz+j*nz+i]; 
					ifs >> V[k*ny*nz+j*nz+i]; 
					ifs >> W[k*ny*nz+j*nz+i]; 

				}
			}
	}


}

void loadflow(double* U,std::string file_str,int nx,int ny,int nz){

	std::ifstream ifs;
	ifs.open(file_str.c_str(),std::ios::in|std::ios::binary);


		
	std::streampos size = ifs.tellg();

	char* buff;
	buff = new char[size];
	

		int c=0;
		for(int k=0;k<nz;k++){
			for(int j=0;j<ny;j++){
				for(int i=0;i<nx;i++){
					ifs.read((char*) &U[i*ny*nz+j*nz+k],8);
			//		std::cout << U[i*ny*nz+j*nz+k] << " " << std::endl;
					c++;
				}
			}
		 }

	
	free(buff);


}


void Lire_Flow(double* U,double* V,double* W,int nx,int ny,int nz,std::string FILE_U,std::string FILE_V, std::string FILE_W,std::string FILE_VITESSE,int current_index_t){

	if(ReTau < 1000){
		if(current_index_t < 10){
		    loadflow(U,FILE_U+"00000"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(V,FILE_V+"00000"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(W,FILE_W+"00000"+patch::to_string(current_index_t),nx,ny,nz);
		}
		else if(current_index_t < 100){
		    loadflow(U,FILE_U+"0000"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(V,FILE_V+"0000"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(W,FILE_W+"0000"+patch::to_string(current_index_t),nx,ny,nz);
		}
		else if(current_index_t < 1000){
		    loadflow(U,FILE_U+"000"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(V,FILE_V+"000"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(W,FILE_W+"000"+patch::to_string(current_index_t),nx,ny,nz);
		}
		else if(current_index_t< 10000){
		    loadflow(U,FILE_U+"00"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(V,FILE_V+"00"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(W,FILE_W+"00"+patch::to_string(current_index_t),nx,ny,nz);
		}
		else{
		    loadflow(U,FILE_U+"0"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(V,FILE_V+"0"+patch::to_string(current_index_t),nx,ny,nz);
		    loadflow(W,FILE_W+"0"+patch::to_string(current_index_t),nx,ny,nz);
		}	
	}
	else{

		std::string file_str = FILE_VITESSE + "/Data"+patch::to_string(current_index_t)+".dat";
		loadflow_1000(U,V,W,file_str,nx,ny, nz);


	}




}




// Gestion des fichiers de sorties


void Enregistrement_champs(std::string& FILE_RES_str,double* H_U,double* H_V,double* H_W,int y0_init,int Nb_particules,int Nb_part_per_plane, int iter_t,int nt,int t0){

	std::ofstream ofs;
	std::string FILE_RES;
	int k = y0_init;



	for(int i=0;i<Nb_particules;i++){
		if(fmod(i,Nb_part_per_plane) == 0 && i > 0){
			ofs.close();
			k = k+1;
			//print np.mod(i,Nb_part_per_plane);
		}


		if(fmod(i,Nb_part_per_plane)+1 < 10)	
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P000"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+"_t0_"+patch::to_string(t0)+".dat";
		else if(fmod(i,Nb_part_per_plane)+1 < 100)	
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P00"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+"_t0_"+patch::to_string(t0)+".dat";
		else if(fmod(i,Nb_part_per_plane)+1 < 1000)	
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P0"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+"_t0_"+patch::to_string(t0)+".dat";
		else
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+"_t0_"+patch::to_string(t0)+".dat";

			
		ofs.open(FILE_RES.c_str());

		for(int j=0;j<iter_t;j++){
			ofs << patch::to_string(H_U[i*nt+j]);
			ofs << " ";
			ofs << patch::to_string(H_V[i*nt+j]);
			ofs << " ";
			ofs << patch::to_string(H_W[i*nt+j]);
			ofs << std::endl;
		}

		ofs.close();
	}

}


// Initialisation des vitesses


void init_Particles(int n_p,double* X_1,double* X_2,double*X_3,double*  y_v,int nx,int ny,int nz,int nt,int y0_init,double lx,double lz,double part_per_plane){
	int c,i,j,k;
	i=0;
	j=0;
	k=0;
	c=0;
	for(int i=0;i<ny;i++)
		std::cout << y_v[i] << " " << i << std::endl;


	k = y0_init;
	while(c < n_p){
		i=0;
		while (i <  (int)sqrt(part_per_plane) &&  c < n_p){
			j=0;
			while (j <  (int)sqrt(part_per_plane) &&  c < n_p){
			//	std::cout << c << " " << k << " " << y_v[k] <<std::endl;
				X_2[c*nt] = y_v[k];
				//X_1[c*nt] = 10*lx/nx+0.01;
				X_1[c*nt] = lx/( sqrt(part_per_plane) *1.2)*(i+1);
                X_3[c*nt]= lz/(sqrt(part_per_plane)*1.1)*(j+1);
				c=c+1;
				j=j+1;
				
			}
			i=i+1;
		}
		k=k+1;
	}
}

// Mise a jour des vitesses
void Update_Vel(double*  U,double*  V,double*  W,int cas,int nx,int ny,int nz,double dx,double dz,double* y_v,double t){		


	int c = 0;

	if(cas == 1){
		for(int i=0;i<nx;i++){
			i_index = i*ny*(nz);
			for(int j=0;j<ny;j++){
				ij_index = i_index+j*(nz);
				for(int k=0;k<nz;k++){
				        U[c] = y_v[j]*(2-y_v[j]);
					V[c] = 0.;
					W[c] = 0.8*sin((double)i*dx*2);	
					c++;
					}
			}
		}
	}
    
    if(cas == 2){
        for(int i=0;i<nx;i++){
            i_index = i*ny*(nz);
            for(int j=0;j<ny;j++){
                ij_index = i_index+j*(nz);
                for(int k=0;k<nz;k++){
                    U[c] = 1.0+cos(t);
                    V[c] = 0.1*sin(t);
                    W[c] = t*cos(t)/25;
                    c++;
                }
            }
        }
    }
}



// Vitesse vers grille

void  Velocity_to_grid(double* Grid,double* U, int nx,int ny,int nz){
	unsigned int i,j,k;

	nx=nx-1;
	nz=nz-1;


	int cU = 0;
	int cG = 0;
	
	for(i=0;i<nx;i++){
	 	i_index = i*ny*(nz);
	 	i_index2= i*ny*(nz-1);
		for(j=0;j<ny;j++){
			ij_index = i_index+j*(nz);
			ij_index2 = i_index2+j*(nz-1);
			for(k=0;k<nz;k++){
				Grid[ij_index+k] = U[ij_index2+k];
			}
		}
	}

	for(j=0;j<ny;j++){
		j_index = j*(nz);
		j_index2 = j*(nz-1);
		for( k=0;k<nz;k++){
			Grid[(nx)*ny*(nz)+j_index+k] = U[j_index2+k];
		}
	}


	for(i=0;i<nx-1;i++){
		i_index = i*ny*(nz);
		i_index2= i*ny*(nz-1);
		for(j=0;j<ny;j++){
			ij_index = i_index+j*(nz);
			ij_index2 = i_index2+j*(nz-1);
			Grid[ij_index+nz-1] = U[ij_index2];
		}
	}

	i_index = (nx)*ny*(nz);
	i_index2= (nx-1)*ny*(nz-1);
	for(j=0;j<ny;j++){
		ij_index = i_index+j*(nz);
		ij_index2 = i_index2+j*(nz-1);
		Grid[ij_index+nz] = U[ij_index2];
	}
}


// ------ Methodes numeriques ------

// Lineaire


double TrilinearInterp(double X0,double X1,double X2,double* Grid,double dx,double*  y_v,double dz,int nx, int ny,int nz,int index){
	double  res =0;

//	std::cout << "Lineaire "   << X0 << " "<< X1 << " " << X2 << std::endl;

	// Reperer la particule

	x0 = floor(X0/dx);
	x1 = floor(X0/dx)+1;
	ix_1 = fmod(x_1,nx);

	/*if(X0>12)
		std::00 << x_1 <<" "<<x_0 << " "<< X0<< " "<< std::endl;
*/

	z0 = floor(X2/dz);
	z1 = floor(X2/dz)+1;
	iz_1 = fmod(z_1,nz);

	if(X2< 0)
		std::cout << x1 <<" "<<x0 << " "<< X0<< " " << iz_1 <<" "<<z0 << " "<< X2<< fmod(X2,nz) << std::endl;

	find = 0;
	int k = 0;
	
	if (X1 >= y_v[ny-1])
		k = ny;
	else{
		while (find == 0){
			if( X1 >= y_v[k] &&  X1 <= y_v[k+1])
				find += 1;
			k += 1;
		}
	}
	
	y_b =  k ;
	y_a =  k-1;

//	std::cout << x_0 << " " << z_0 << " " << y_a << std::endl;

	xd = (X0-x0*dx)/(x1*dx-x0*dx);
	yd = (X1-y_v[y_a])/(y_v[y_b]-y_v[y_a]);
	zd = (X2-z0*dz)/(z1*dz-z0*dz);


//	std::cout << xd << " " << yd << " " <<zd << std::endl;

	c00 = Grid[x0*ny*nz+y_a*nz+z0]*(1.-xd) + Grid[ix_1*ny*nz+y_a*nz+z0]*(xd) ;
	c01 = Grid[x0*ny*nz+y_a*nz+iz_1]*(1.-xd) + Grid[ix_1*ny*nz+y_a*nz+iz_1]*(xd) ;
	c10 = Grid[x0*ny*nz+y_b*nz+z0]*(1.-xd) + Grid[ix_1*ny*nz+y_b*nz+z0]*(xd) ;
	c11 = Grid[x0*ny*nz+y_b*nz+iz_1]*(1.-xd) + Grid[ix_1*ny*nz+y_b*nz+iz_1]*(xd) ;

	c0 = c00*(1.-yd) + c10*yd;
	c1 = c01*(1.-yd) + c11*yd;

	res = c0*(1.-zd) + c1*zd;

	

//	std::cout << res << std::endl;
	return res;

}

// Interpolation Lagrangienne

double PolyLagrange(int i,double xhi){
	double res = 0.;
	i=i+1;

	if (i== 1){
		res = 6.*xhi;
		xhi *= xhi;
		res += -5.0*xhi;
		xhi *= xhi;
		res += -5.0*xhi;
		xhi *= xhi;
		res += 5.0*xhi;
		xhi *= xhi;
		res += -xhi;

		res = res / 120.;
	}

	else if (i == 2){
		res = 12.*xhi;
		xhi *= xhi;
		res += 16.*xhi;
		xhi *= xhi;
		res += xhi;
		xhi *= xhi;
		res += -4.0*xhi;
		xhi *= xhi;
		res += -xhi;

		res = res / 24.;
	}
	else if (i == 3){
		res = 12.;
		res += -4.*xhi;
		xhi *= xhi;
		res += -15.*xhi;
		xhi *= xhi;
		res += 5.*xhi;
		xhi *= xhi;
		res += 3.*xhi;
		xhi *= xhi;
		res += -xhi;

		res = res / 12.;
	}
	else if(i == 4){
		res = 12.*xhi;
		xhi *= xhi;
		res += 8.*xhi;
		xhi *= xhi;
		res += -7.0*xhi;
		xhi *= xhi;
		res += -2.0*xhi;
		xhi *= xhi;
		res += 2.0*xhi;

		res = res / 12.;
	}
	else if(i == 5){
		res = -6.*xhi;
		xhi *= xhi;
		res += -1.*xhi;
		xhi *= xhi;
		res += 7.*xhi;
		xhi *= xhi;
		res += xhi;
		xhi *= xhi;
		res += -1.0*xhi;

		res = res / 24.;
	}
	else if (i == 6){
		res = 4.*xhi;
		xhi *= xhi;
		xhi *= xhi;
		res += -5.*xhi;
		xhi *= xhi;
		xhi *= xhi;
		res += xhi;

		res = res / 120.;
	}
	else
		res =0.;

	return res;
}	




double LagrangianInterp(double X0,double X1,double X2,double*  Grid,double dx,double*  y_v,double dz,int nx, int ny,int nz,int index){

	double res = 0.;


	xhi_x = (X0 - floor(X0/dx)*dx)/dx;
	xhi_z = (X2 - floor(X2/dz)*dz)/dz;

	//	std::cout << X0 << " " <<  xhi_x << " "<<X2 << " " << floor(X2/dz) <<" " <<  xhi_z << std::endl;
	//std::cout << floor(X0/dx)*dx << " " <<dx*(floor(X0/dx)+1) << std::endl;

	
	int x_index,z_index;
	double coeff = 0;
	int i,j,l;


	x_2 = floor(X0/dx)*dx;
	x_0 = x_2 - 2*dx;
	x_1 = x_2 -dx;
	x_3 = x_2 +dx;
	x_4 = x_2 +2.*dx;
	x_5 = x_2+3.*dx;

	z_2 = floor(X2/dz)*dz;
	z_0 = z_2 -2.*dz;
	z_1 = z_2 - dz;
	z_3 = z_2 + dz;
	z_4 = z_2 + 2*dz;
	z_5 = z_2 + 3*dz;

	
	//std::cout << z_0 << " " << z_1 << " " << z_2 << " " << z_3 << " " << z_4 << " " << z_5 << std::endl;

//	std::cout  << "Lagrange " << X0 << " "<< X1 << " " << X2 << std::endl;

	find = 0;
	int k = 0;
	if (X1 >= y_v[ny-1])
		k = ny;
	else{
		while (find == 0){
			if( X1 >= y_v[k] &&  X1 <= y_v[k+1])
				find += 1;
			k += 1;
		}
	}

	k--;

	// Si il y a assez de point pour gerer le mur
	if(k >= 2){

		y_1 = y_v[k-2];
		y_2 = y_v[k-1];
		y_3 = y_v[k];
		y_4 = y_v[k+1];
		y_5 = y_v[k+2];
		y_6 = y_v[k+3];


		Lx[0] = -(xhi_x+1.)*xhi_x*(xhi_x-1.)*(xhi_x-2.)*(xhi_x-3.)/120.;
		Lx[1] = (xhi_x+2.)*xhi_x*(xhi_x-1.)*(xhi_x-2.)*(xhi_x-3.)/24.;
		Lx[2] = -(xhi_x+2.)*(xhi_x+1.)*(xhi_x-1.)*(xhi_x-2.)*(xhi_x-3.)/12.;
		Lx[3] = (xhi_x+2.)*(xhi_x+1.)*xhi_x*(xhi_x-2.)*(xhi_x-3.)/12.;
		Lx[4] = -(xhi_x+2.)*(xhi_x+1.)*xhi_x*(xhi_x-1.)*(xhi_x-3.)/24.;
		Lx[5] = (xhi_x+2.)*(xhi_x+1.)*xhi_x*(xhi_x-1.)*(xhi_x-2.)/120.;


		//Lz[0] = (X2-z_1)*(X2-z_2)*(X3-z_3)*
		Lz[0] = -(xhi_z+1)*xhi_z*(xhi_z-1.)*(xhi_z-2.)*(xhi_z-3.)/120.;
		Lz[1] = (xhi_z+2.)*xhi_z*(xhi_z-1.)*(xhi_z-2.)*(xhi_z-3.)/24.;
		Lz[2] = -(xhi_z+2.)*(xhi_z+1.)*(xhi_z-1.)*(xhi_z-2.)*(xhi_z-3.)/12.;
		Lz[3] = (xhi_z+2.)*(xhi_z+1.)*xhi_z*(xhi_z-2.)*(xhi_z-3.)/12.;
		Lz[4] = -(xhi_z+2.)*(xhi_z+1.)*xhi_z*(xhi_z-1.)*(xhi_z-3.)/24.;
		Lz[5] = (xhi_z+2.)*(xhi_z+1.)*xhi_z*(xhi_z-1.)*(xhi_z-2.)/120.;

		Ly[0] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_1-y_2)*(y_1-y_3)*(y_1-y_4)*(y_1-y_5)*(y_1-y_6));
		Ly[1] = (X1-y_1)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_2-y_1)*(y_2-y_3)*(y_2-y_4)*(y_2-y_5)*(y_2-y_6));
		Ly[2] = (X1-y_2)*(X1-y_1)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_3-y_2)*(y_3-y_1)*(y_3-y_4)*(y_3-y_5)*(y_3-y_6));
		Ly[3] = (X1-y_2)*(X1-y_3)*(X1-y_1)*(X1-y_5)*(X1-y_6)/((y_4-y_2)*(y_4-y_3)*(y_4-y_1)*(y_4-y_5)*(y_4-y_6));
		Ly[4] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_1)*(X1-y_6)/((y_5-y_2)*(y_5-y_3)*(y_5-y_4)*(y_5-y_1)*(y_5-y_6));
		Ly[5] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_1)/((y_6-y_2)*(y_6-y_3)*(y_6-y_4)*(y_6-y_5)*(y_6-y_1));




		

		/*	std::cout << Ly[0] << " "<< Ly[1] << " "<< Ly[2] << " "<< Ly[3] << " "<< Ly[4]<< " "<<Ly[5]<<" " <<Ly[0]+Ly[1]+Ly[2]+Ly[3]+Ly[4]+Ly[5] << " " <<std::endl;
		std::cout << Lx[0] << " "<< Lx[1] << " "<< Lx[2] << " "<< Lx[3] << " "<< Lx[4]<< " "<<Lx[5] << " " <<xhi_x << " "<<  Lx[0]+Lx[1]+Lx[2]+Lx[3]+Lx[4]+Lx[5]<< std::endl;
		std::cout << Lz[0] << " "<< Lz[1] << " "<< Lz[2] << " "<< Lz[3] << " "<< Lz[4]<< " "<<Lz[5] <<" " <<xhi_z << " "<<  Lz[0]+Lz[1]+Lz[2]+Lz[3]+Lz[4]+Lz[5]<<  std::endl;
		*/

		x_index = (int)positive_modulo(floor(X0/dx)-2,nx);
		z_index = (int)positive_modulo(floor(X2/dz)-2,(nz));

		for(i=0;i<6;i++){
		  z_index = (int)positive_modulo(floor(X2/dz)-2,(nz));
			for(j=0;j<6;j++){
				for(l=0;l<6;l++){
				  
					res +=  Grid[x_index*ny*nz+(k+l-2)*nz+z_index]*Lx[i]*Ly[l]*Lz[j];	
					coeff += Lx[i]*Ly[l]*Lz[j];
					
					if(l==100)
					  std::cout << res << " " << Grid[x_index*ny*nz+(k+l-2)*nz+z_index] << " " << z_index<<std::endl;
				}
				z_index=(z_index+1)%nz;
			}
			x_index=(x_index+1)%nx;
		}
	}
	else
		res = TrilinearInterp(X0,X1,X2,Grid,dx,y_v,dz,nx,ny,nz,index);

	//	std::cout << res << " " << sin(X0) << " " << floor(X0/dx) << " " << floor(X2/dz) << " " << y_v[k]  << " "<<0.1*sin(X2)<< " "<<X1<<" " << coeff <<std::endl<<std::endl;



	return res;

}



// Allocation de la memoire
void Memoire_Grille_3D(double*** Grid,int nx,int ny,int nz){


	Grid = (double***)malloc(nx*sizeof(**Grid));

	std::cout << nx << " " << ny<< " " << nz << std::endl;

	for(int i=0;i<nx;i++){
		Grid[i] = (double**)malloc(ny*sizeof(*Grid));
		for(int j=0;j<ny;j++){
			Grid[i][j] = (double*)malloc(nz*sizeof(Grid));
		}
	}	

}


void free_memoire_3D(double*** Grid,int nx,int ny,int nz){
	for(int i=nx;i<-1;i--){
		for(int j=ny;j<-1;j--){
			free(Grid[i][j]);
		}
		free(Grid[i]);
	}	

}

void Memoire_Grille_4D(double**** Grid,int nx,int ny,int nz,int dim){

std::cout << nx << " " << ny<< " " << nz <<" " << dim<< std::endl;
	Grid = (double****)malloc(nx*sizeof(***Grid));

	for(int i=0;i<nx;i++){
		Grid[i] = (double***)malloc(ny*sizeof(**Grid));
		for(int j=0;j<ny;j++){
			Grid[i][j] = (double**)malloc(nz*sizeof(*Grid));
			for(int k=0;k<nz;k++){
				Grid[i][j][k] = (double*)malloc(dim*sizeof(Grid));
			}
		}
	}	

}


void free_memoire_4D(double**** Grid,int nx,int ny,int nz,int dim){
	for(int i=nx;i<-1;i--){
		for(int j=ny;j<-1;j--){
			for(int k=nz;k<-1;k--){
				free(Grid[i][j][k]);
			}
			free(Grid[i][j]);
		}
		free(Grid[i]);
	}	

}


using namespace::std;

int main(int argc, char** argv){
	int i,j,k,l;
	// Variables de la simulation Lagrangienne

	std::string FILE_arg = argv[1];
    
    int t0_win = atoi(argv[2]);

	std::cout << FILE_arg << std::endl;

	 int Nb_particules;
	 int Nb_particules_par_ligne;
	 double tf; // temps final
	 int y0_init;
	 double dt;
	 double dt_data;
	 int nb_iter_win;

	 double lx,ly,lz;
	 int nx,ny,nz;

	 double dx,dz;

	int test = 1;


	 std::string FILE_VITESSE, DIR_Particules_output, DIR_Vitesses_output;


	 std::cout << "Chargement des parametrs .... ";
	Chargement_parametres( FILE_arg,&Nb_particules,&	Nb_particules_par_ligne,&tf,&y0_init,&dt,&dt_data,&nb_iter_win,&lx,&ly,&lz,&nx,&ny,&nz,FILE_VITESSE,DIR_Particules_output,DIR_Vitesses_output);
	std::cout << "....OK " << std::endl;
	Affiche_Param(&Nb_particules,&Nb_particules_par_ligne,&tf,&y0_init,&dt,&dt_data,&nb_iter_win,&lx,&ly,&lz,&nx,&ny,&nz,FILE_VITESSE,DIR_Particules_output,DIR_Vitesses_output);


	std::string Y_file = FILE_VITESSE+"/yp4.dat";

	if(ReTau < 1000){
		 dx = lx/nx;
		 dz = lz/nz;
	 }
	 else{
		 dx = lx/(nx-1);
		 dz = lz/(nz-1);
	 }

	 double nt = floor(tf/dt+1);


	double* History_Particles_U;
	double* History_Particles_V;
	double* History_Particles_W;
	
	double* History_Vel_Particles_U;
	double* History_Vel_Particles_V;
	double* History_Vel_Particles_W;

	double* Grid_current_t_U;
	double* Grid_current_t_V;
	double* Grid_current_t_W;

	double* Grid_next_t_U;
	double* Grid_next_t_V;
	double* Grid_next_t_W;

	double*  Grid_t_U;
	double*  Grid_t_V;
	double*  Grid_t_W;

	double*  Grid_t1_U;
	double*  Grid_t1_V;
	double*  Grid_t1_W;

	double*  Grid_t2_U;
	double*  Grid_t2_V;
	double*  Grid_t2_W;

	double*  Grid_t3_U;
	double*  Grid_t3_V;
	double*  Grid_t3_W;


	double*  Grid_t4_U;
	double*  Grid_t4_V;
	double*  Grid_t4_W;

	double*  Grid_t5_U;
	double*  Grid_t5_V;
	double*  Grid_t5_W;




	std::cout <<  Nb_particules << std::endl;


	std::cout <<floor(tf/dt+2) << std::endl;

	History_Particles_U = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Particles_V = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Particles_W = (double*)malloc(Nb_particules*nt*sizeof(double));

	History_Vel_Particles_U = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Vel_Particles_V = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Vel_Particles_W = (double*)malloc(Nb_particules*nt*sizeof(double));

	if(dt_data > 1){
	
		Grid_current_t_U= (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_current_t_V= (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_current_t_W= (double*)malloc((nx)*ny*(nz)*sizeof(double));

		Grid_next_t_U= (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_next_t_V= (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_next_t_W= (double*)malloc((nx)*ny*(nz)*sizeof(double));

	}
	else{

		Grid_t_U = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t_V = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t_W = (double*)malloc((nx)*ny*(nz)*sizeof(double));

		Grid_t1_U = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t1_V = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t1_W = (double*)malloc((nx)*ny*(nz)*sizeof(double));

		Grid_t2_U = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t2_V = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t2_W = (double*)malloc((nx)*ny*(nz)*sizeof(double));

		Grid_t3_U = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t3_V = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t3_W = (double*)malloc((nx)*ny*(nz)*sizeof(double));

		Grid_t4_U = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t4_V = (double*)malloc((nx)*ny*(nz)*sizeof(double));
		Grid_t4_W = (double*)malloc((nx)*ny*(nz)*sizeof(double));
	}




	double* x_values;
	double* z_values;
	double* y_values;

	x_values = (double*)malloc((nx+1)*sizeof(double));
	z_values = (double*)malloc((nz+1)*sizeof(double));
	y_values = (double*)malloc(ny*sizeof(double));

	
	for(int i=0;i<=nx;i++){
	  x_values[i] = i*dx;
	  std::cout << x_values[i]<<std::endl;}

	for(int i=0;i<=nz;i++){
		z_values[i] = i*dz;
		std::cout << z_values[i]<<std::endl;
}
	Lecture_Y(y_values,ny,Y_file);


	std::string FILE_U = FILE_VITESSE+"/ux";
	std::string FILE_V = FILE_VITESSE+"/uy";
	std::string FILE_W = FILE_VITESSE+"/uz";

	double* U;
	double* V;
	double* W;

	U = (double*)malloc((nx)*ny*(nz)*sizeof(double));
	V = (double*)malloc((nx)*ny*(nz)*sizeof(double));
	W = (double*)malloc((nx)*ny*(nz)*sizeof(double));
	
	int current_index_t = 0;
	double current_time = 0;
	double ratio_time = 0;


	double  k1[3];
	double  k2[3];
	double  k3[3];
	double  k4[3];
	double 	k5[3];
	double  yn[3];




	double* k1_AMC_U;
	double* k1_AMC_V;
	double* k1_AMC_W;

	double* k2_AMC_U;
	double* k2_AMC_V;
	double* k2_AMC_W;

	double* k3_AMC_U;
	double* k3_AMC_V;
	double* k3_AMC_W;

	double* k4_AMC_U;
	double* k4_AMC_V;
	double* k4_AMC_W;

	double* k5_AMC_U;
	double* k5_AMC_V;
	double* k5_AMC_W;

	k1_AMC_U = (double*)malloc(Nb_particules*sizeof(double));
	k1_AMC_V = (double*)malloc(Nb_particules*sizeof(double));
	k1_AMC_W = (double*)malloc(Nb_particules*sizeof(double));

	k2_AMC_U = (double*)malloc(Nb_particules*sizeof(double));
	k2_AMC_V = (double*)malloc(Nb_particules*sizeof(double));
	k2_AMC_W = (double*)malloc(Nb_particules*sizeof(double));

	k3_AMC_U = (double*)malloc(Nb_particules*sizeof(double));
	k3_AMC_V = (double*)malloc(Nb_particules*sizeof(double));
	k3_AMC_W = (double*)malloc(Nb_particules*sizeof(double));

	k4_AMC_U = (double*)malloc(Nb_particules*sizeof(double));
	k4_AMC_V = (double*)malloc(Nb_particules*sizeof(double));
	k4_AMC_W = (double*)malloc(Nb_particules*sizeof(double));

	k5_AMC_U = (double*)malloc(Nb_particules*sizeof(double));
	k5_AMC_V = (double*)malloc(Nb_particules*sizeof(double));
	k5_AMC_W = (double*)malloc(Nb_particules*sizeof(double));
	
	int iter_t = 0;

	double dt2 = 0.5*dt;
	double dt6 = dt/6.0;

	time_t time_clock;
	time_t cu_time;
	time(&time_clock);
	double t2;

	double y_00,y_01,y_02;
	std:: cout << "Itilialisation des particules ... " << std::endl;
	init_Particles(Nb_particules,History_Particles_U,History_Particles_V,History_Particles_W,y_values,nx,ny,nz,nt,y0_init,lx,lz,Nb_particules_par_ligne);
	std:: cout << "Itilialisation des particules ok " << std::endl;

	current_index_t = t0_win+1;


	std::cout << current_index_t << endl; 

	if(dt_data > 1){

		Lire_Flow(Grid_current_t_U,Grid_current_t_V,Grid_current_t_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t);	
		Lire_Flow(Grid_next_t_U,Grid_next_t_V,Grid_next_t_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+1);
	}
	else{
		Lire_Flow(Grid_t_U,Grid_t_V,Grid_t_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t);	
		Lire_Flow(Grid_t1_U,Grid_t1_V,Grid_t1_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+1);	
		Lire_Flow(Grid_t2_U,Grid_t2_V,Grid_t2_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+2);	
		Lire_Flow(Grid_t3_U,Grid_t3_V,Grid_t3_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+3);	
		Lire_Flow(Grid_t4_U,Grid_t4_V,Grid_t4_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+4);	
	}





	// Initialisation pour la methode d'Adams Moulton

	if(SCHM == 5){

		cout << "Initialisation pour les shemas AMC " << endl;

		// Methode d'Euler pour y1
		cout << "Ordre 1  " << endl;

		for(i=0;i<Nb_particules;i++){


			it = i*nt;

			k1_AMC_U[i] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_U,dx,y_values,dz,nx,ny,nz,0);
			k1_AMC_V[i] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_V,dx,y_values,dz,nx,ny,nz,1);
			k1_AMC_W[i] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_W,dx,y_values,dz,nx,ny,nz,2);

			History_Particles_U[it+1] = y_00 + dt*k1_AMC_U[i];
			History_Particles_V[it+1] = y_01 + dt*k1_AMC_V[i];
			History_Particles_W[it+1] = y_02 + dt*k1_AMC_W[i];

			if(History_Particles_V[it+1] < 0)
				History_Particles_V[it+1] = 0;

			if(History_Particles_U[it+1] >= lx || History_Particles_U[it+1] < 0)
				History_Particles_U[it+1] = positive_modulo(History_Particles_U[it+1],lx);

			if (History_Particles_W[it+1] > lz || History_Particles_W[it+1] < 0) 
				History_Particles_W[it+1] = positive_modulo(History_Particles_W[it+1],lz);
		}

		// Methode d'ordre 2 pour y2
		cout << "Ordre 2  " << endl;
		
		for(i=0;i<Nb_particules;i++){

			it = i*nt;

			k2_AMC_U[i] =  LagrangianInterp(History_Particles_U[it+1],History_Particles_V[it+1],History_Particles_W[it+1],Grid_t1_U,dx,y_values,dz,nx,ny,nz,0);
			k2_AMC_V[i] =  LagrangianInterp(History_Particles_V[it+1],History_Particles_V[it+1],History_Particles_W[it+1],Grid_t1_V,dx,y_values,dz,nx,ny,nz,1);
			k2_AMC_W[i] =  LagrangianInterp(History_Particles_W[it+1],History_Particles_V[it+1],History_Particles_W[it+1],Grid_t1_W,dx,y_values,dz,nx,ny,nz,2);

			History_Particles_U[it+2] = History_Particles_U[it+1] + dt*(1.50*k2_AMC_U[i]-0.5*k1_AMC_U[i]);
			History_Particles_V[it+2] = History_Particles_V[it+1] + dt*(1.50*k2_AMC_V[i]-0.5*k1_AMC_V[i]);
			History_Particles_W[it+2] = History_Particles_W[it+1] + dt*(1.50*k2_AMC_W[i]-0.5*k1_AMC_W[i]);

			if(History_Particles_V[it+2] < 0)
				History_Particles_V[it+2] = 0;

			if(History_Particles_U[it+2] >= lx || History_Particles_U[it+2] < 0)
				History_Particles_U[it+2] = positive_modulo(History_Particles_U[it+2],lx);

			if (History_Particles_W[it+2] > lz || History_Particles_W[it+2] < 0) 
				History_Particles_W[it+2] = positive_modulo(History_Particles_W[it+2],lz);
		}

				// Methode d'ordre 3 pour y3
		cout << "Ordre 3  " << endl;
		
		for(i=0;i<Nb_particules;i++){

			it = i*nt;

			k3_AMC_U[i] =  LagrangianInterp(History_Particles_U[it+2],History_Particles_V[it+2],History_Particles_W[it+2],Grid_t2_U,dx,y_values,dz,nx,ny,nz,0);
			k3_AMC_V[i] =  LagrangianInterp(History_Particles_U[it+2],History_Particles_V[it+2],History_Particles_W[it+2],Grid_t2_V,dx,y_values,dz,nx,ny,nz,1);
			k3_AMC_W[i] =  LagrangianInterp(History_Particles_U[it+2],History_Particles_V[it+2],History_Particles_W[it+2],Grid_t2_W,dx,y_values,dz,nx,ny,nz,2);


			History_Particles_U[it+3] = History_Particles_U[it+2] + dt*(23.0/11.0*k3_AMC_U[i]-4.0/3.0*k2_AMC_U[i]+5.0/12.0*k1_AMC_U[i]);
			History_Particles_V[it+3] = History_Particles_V[it+2] + dt*(23.0/11.0*k3_AMC_V[i]-4.0/3.0*k2_AMC_V[i]+5.0/12.0*k1_AMC_V[i]);
			History_Particles_W[it+3] = History_Particles_W[it+2] + dt*(23.0/11.0*k3_AMC_W[i]-4.0/3.0*k2_AMC_W[i]+5.0/12.0*k1_AMC_W[i]);


			if(History_Particles_V[it+3] < 0)
				History_Particles_V[it+3] = 0;

			if(History_Particles_U[it+3] >= lx || History_Particles_U[it+3] < 0)
				History_Particles_U[it+3] = positive_modulo(History_Particles_U[it+3],lx);

			if (History_Particles_W[it+3] > lz || History_Particles_W[it+3] < 0) 
				History_Particles_W[it+3] = positive_modulo(History_Particles_W[it+3],lz);

		}
		
	
	// Methode d'ordre 4 pour y4
		cout << "Ordre 4  " << endl;
		
		for(i=0;i<Nb_particules;i++){
			it = i*nt;

			k4_AMC_U[i] =  LagrangianInterp(History_Particles_U[it+3],History_Particles_V[it+3],History_Particles_W[it+3],Grid_t3_U,dx,y_values,dz,nx,ny,nz,0);
			k4_AMC_V[i] =  LagrangianInterp(History_Particles_U[it+3],History_Particles_V[it+3],History_Particles_W[it+3],Grid_t3_V,dx,y_values,dz,nx,ny,nz,1);
			k4_AMC_W[i] =  LagrangianInterp(History_Particles_U[it+3],History_Particles_V[it+3],History_Particles_W[it+3],Grid_t3_W,dx,y_values,dz,nx,ny,nz,2);

			History_Particles_U[it+4] = History_Particles_U[it+3] + dt*(55.0/24.0*k4_AMC_U[i]-59.0/24.0*k3_AMC_U[i]+37.0/24.0*k2_AMC_U[i]-3./8.*k1_AMC_U[i]);
			History_Particles_V[it+4] = History_Particles_V[it+3] + dt*(55.0/24.0*k4_AMC_V[i]-59.0/24.0*k3_AMC_V[i]+37.0/24.0*k2_AMC_V[i]-3./8.*k1_AMC_V[i]);
			History_Particles_W[it+4] = History_Particles_W[it+3] + dt*(55.0/24.0*k4_AMC_W[i]-59.0/24.0*k3_AMC_W[i]+37.0/24.0*k2_AMC_W[i]-3./8.*k1_AMC_W[i]);

			if(History_Particles_V[it+4] < 0)
				History_Particles_V[it+4] = 0;

			if(History_Particles_U[it+4] >= lx || History_Particles_U[it+4] < 0)
				History_Particles_U[it+4] = positive_modulo(History_Particles_U[it+4],lx);

			if (History_Particles_W[it+4] > lz || History_Particles_W[it+4] < 0) 
				History_Particles_W[it+4] = positive_modulo(History_Particles_W[it+4],lz);
		}

		// Methode d'ordre 5 pour y5
		cout << "Ordre 5  " << endl;
		
		for(i=0;i<Nb_particules;i++){
			it = i*nt;

			k5_AMC_U[i] =  LagrangianInterp(History_Particles_U[it+4],History_Particles_V[it+4],History_Particles_W[it+4],Grid_t4_U,dx,y_values,dz,nx,ny,nz,0);
			k5_AMC_V[i] =  LagrangianInterp(History_Particles_U[it+4],History_Particles_V[it+4],History_Particles_W[it+4],Grid_t4_V,dx,y_values,dz,nx,ny,nz,1);
			k5_AMC_W[i] =  LagrangianInterp(History_Particles_U[it+4],History_Particles_V[it+4],History_Particles_W[it+4],Grid_t4_W,dx,y_values,dz,nx,ny,nz,2);

			History_Particles_U[it+5] = History_Particles_U[it+4] + dt*(1901.0/720.0*k5_AMC_U[i]-1387.0/360.0*k4_AMC_U[i]+109.0/30.0*k3_AMC_U[i]-637./360.*k2_AMC_U[i]+251./720.*k1_AMC_U[i]);
			History_Particles_V[it+5] = History_Particles_V[it+4] + dt*(1901.0/720.0*k5_AMC_V[i]-1387.0/360.0*k4_AMC_V[i]+109.0/30.0*k3_AMC_V[i]-637./360.*k2_AMC_V[i]+251./720.*k1_AMC_V[i]);
			History_Particles_W[it+5] = History_Particles_W[it+4] + dt*(1901.0/720.0*k5_AMC_W[i]-1387.0/360.0*k4_AMC_W[i]+109.0/30.0*k3_AMC_W[i]-637./360.*k2_AMC_W[i]+251./720.*k1_AMC_W[i]);

			if(History_Particles_V[it+5] < 0)
				History_Particles_V[it+5] = 0;

			if(History_Particles_U[it+5] >= lx || History_Particles_U[it+5] < 0)
				History_Particles_U[it+5] = positive_modulo(History_Particles_U[it+4],lx);

			if (History_Particles_W[it+5] > lz || History_Particles_W[it+5] < 0) 
				History_Particles_W[it+5] = positive_modulo(History_Particles_W[it+5],lz);
		}

		current_index_t  += 4;

		cout << current_index_t << endl;
		current_time = 4*dt;
		iter_t = 4;
	}


	double ratio_time_0,ratio_time_1,ratio_time_2;

	int nynz = ny*(nz);
	int dnt = 3*nt;

	int i_index,j_index,ij_index;


	if(SCHM == 4){

		while(current_time < tf-dt){

			time(&cu_time);
			std::cout << "Iteration: " << iter_t <<  " file "<<  current_index_t <<  " Temps: " << current_time <<" Temps de calcul " << cu_time - time_clock<<  std::endl;



			if(dt_data > 1){

				ratio_time_0 = fmod(iter_t,nb_iter_win)/nb_iter_win;
				ratio_time_1 = fmod(iter_t+0.5,nb_iter_win)/nb_iter_win;

				i_index = 0;
				for(i=0;i<nx;i++){
					for(j=0;j<ny;j++){
						for(k=0;k<nz;k++){
								Grid_t_U[i_index] =  (1.-ratio_time_0)*Grid_current_t_U[i_index] + ratio_time_0*Grid_next_t_U[i_index];
								Grid_t_V[i_index] =  (1.-ratio_time_0)*Grid_current_t_V[i_index] + ratio_time_0*Grid_next_t_V[i_index];
								Grid_t_W[i_index] =  (1.-ratio_time_0)*Grid_current_t_W[i_index] + ratio_time_0*Grid_next_t_W[i_index];


																		
								Grid_t1_U[i_index]  =  (1.-ratio_time_1)*Grid_current_t_U[i_index]  + ratio_time_1*Grid_next_t_U[i_index];
								Grid_t1_V[i_index]  =  (1.-ratio_time_1)*Grid_current_t_V[i_index]  + ratio_time_1*Grid_next_t_V[i_index];
								Grid_t1_W[i_index]  =  (1.-ratio_time_1)*Grid_current_t_W[i_index]  + ratio_time_1*Grid_next_t_W[i_index];

								i_index++;
					
						}
					}
				}


		
				//Next temporal window
		
				if(fmod(iter_t+1,nb_iter_win) == 0 and iter_t != 0){
					std::cout << "Nouvelle fenetre: " <<current_index_t+1<< std::endl;
					current_index_t=  current_index_t + 2;


					i_index = 0;
					for(i=0;i<nx;i++){
						for(j=0;j<ny;j++){
							for(k=0;k<nz;k++){
								Grid_current_t_U[i_index] = Grid_next_t_U[i_index] ;
								Grid_current_t_V[i_index] = Grid_next_t_V[i_index] ;
								Grid_current_t_W[i_index] = Grid_next_t_W[i_index];
								i_index++;
							}
						}
					}


					Lire_Flow(Grid_next_t_U,Grid_next_t_V,Grid_next_t_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+1);

				}


				ratio_time_2 = fmod(iter_t+1,nb_iter_win)/nb_iter_win;

				i_index = 0;
				for(i=0;i<nx;i++){
					for(j=0;j<ny;j++){
						for(k=0;k<nz;k++){
							Grid_t2_U[i_index] =  (1.-ratio_time_2)*Grid_current_t_U[i_index] + ratio_time_2*Grid_next_t_U[i_index];
							Grid_t2_V[i_index] =  (1.-ratio_time_2)*Grid_current_t_V[i_index] + ratio_time_2*Grid_next_t_V[i_index];
							Grid_t2_W[i_index] =  (1.-ratio_time_2)*Grid_current_t_W[i_index] + ratio_time_2*Grid_next_t_W[i_index];
							i_index++;
						}
					}
				}

			}
			else{
				i_index = 0;
				if(iter_t != 0){
					current_index_t=  current_index_t + 2;
					for(i=0;i<nx;i++){
						for(j=0;j<ny;j++){
							for(k=0;k<nz;k++){
								Grid_t_U[i_index] = Grid_t2_U[i_index] ;
								Grid_t_V[i_index] = Grid_t2_V[i_index] ;
								Grid_t_W[i_index] = Grid_t2_W[i_index];
								i_index++;
							}
						}
					}

					Lire_Flow(Grid_t1_U,Grid_t1_V,Grid_t1_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+1);	
					Lire_Flow(Grid_t2_U,Grid_t2_V,Grid_t2_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t+2);	
				}
			}
	




			//#pragma omp parallel for private (i,k1,k2,k3,k4)
			for(i=0;i<Nb_particules;i++){

				it = i*nt;

				y_00 = History_Particles_U[it+iter_t];
				y_01 = History_Particles_V[it+iter_t];
				y_02 = History_Particles_W[it+iter_t];

			//	std::cout << y_00 << " "  << y_01<< " " << y_02<< std::endl;


				k1[0] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_U,dx,y_values,dz,nx,ny,nz,0);
				k1[1] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_V,dx,y_values,dz,nx,ny,nz,1);
				k1[2] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_W,dx,y_values,dz,nx,ny,nz,2);

	//			std::cout << "Fin des interpolations 1" << std::endl;

	//			std::cout << fmod(y_00+dt2*k1[0]+lx,lx) << " " << lx << std::endl;

				k2[0] =  LagrangianInterp(positive_modulo(y_00+dt2*k1[0],lx),fmax(0,y_01+dt2*k1[1]),positive_modulo(y_02+dt2*k1[2],lz),Grid_t1_U,dx,y_values,dz,nx,ny,nz,0);
				k2[1] =  LagrangianInterp(positive_modulo(y_00+dt2*k1[0],lx),fmax(0,y_01+dt2*k1[1]),positive_modulo(y_02+dt2*k1[2],lz),Grid_t1_V,dx,y_values,dz,nx,ny,nz,1);
				k2[2] =  LagrangianInterp(positive_modulo(y_00+dt2*k1[0],lx),fmax(0,y_01+dt2*k1[1]),positive_modulo(y_02+dt2*k1[2],lz),Grid_t1_W,dx,y_values,dz,nx,ny,nz,2);

	//			std::cout << "Fin des interpolations 2" << std::endl;

				k3[0] =  LagrangianInterp(positive_modulo(y_00+dt2*k2[0],lx),fmax(0,y_01+dt2*k2[1]),positive_modulo(y_02+dt2*k2[2],lz),Grid_t1_U,dx,y_values,dz,nx,ny,nz,0);
				k3[1] =  LagrangianInterp(positive_modulo(y_00+dt2*k2[0],lx),fmax(0,y_01+dt2*k2[1]),positive_modulo(y_02+dt2*k2[2],lz),Grid_t1_V,dx,y_values,dz,nx,ny,nz,1);
				k3[2] =  LagrangianInterp(positive_modulo(y_00+dt2*k2[0],lx),fmax(0,y_01+dt2*k2[1]),positive_modulo(y_02+dt2*k2[2],lz),Grid_t1_W,dx,y_values,dz,nx,ny,nz,2);

	//			std::cout << "Fin des interpolations 3" << std::endl;

				k4[0] = LagrangianInterp(positive_modulo(y_00+dt*k3[0],lx),fmax(0,y_01+dt*k3[1]),positive_modulo(y_02+dt*k3[2],lz),Grid_t2_U,dx,y_values,dz,nx,ny,nz,0);
				k4[1] = LagrangianInterp(positive_modulo(y_00+dt*k3[0],lx),fmax(0,y_01+dt*k3[1]),positive_modulo(y_02+dt*k3[2],lz),Grid_t2_V,dx,y_values,dz,nx,ny,nz,1);
				k4[2] = LagrangianInterp(positive_modulo(y_00+dt*k3[0],lx),fmax(0,y_01+dt*k3[1]),positive_modulo(y_02+dt*k3[2],lz),Grid_t2_W,dx,y_values,dz,nx,ny,nz,2);


	//			std::cout << "Fin des interpolations 4" << std::endl;

				History_Vel_Particles_U[it+iter_t] =k1[0];
				History_Vel_Particles_V[it+iter_t] =k1[1];
				History_Vel_Particles_W[it+iter_t] =k1[2];

			/*	if(iter_t < 10){
					std::cout << k1[0] << " "<<  y_00 << " " << y_01 << " " << y_02 << std::endl;
				}
	*/

				
				History_Particles_U[it+iter_t+1] = y_00 + dt6*(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0]);
				History_Particles_V[it+iter_t+1] = y_01 + dt6*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1]);
				History_Particles_W[it+iter_t+1] = y_02 + dt6*(k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2]);

				if(History_Particles_V[it+iter_t+1] < 0)
					History_Particles_V[it+iter_t+1] = 0;



				if(History_Particles_U[it+iter_t+1] >= lx || History_Particles_U[it+iter_t+1] < 0)
					History_Particles_U[it+iter_t+1] = positive_modulo(History_Particles_U[it+iter_t+1],lx);

				if (History_Particles_W[it+iter_t+1] > lz || History_Particles_W[it+iter_t+1] < 0) 
					History_Particles_W[it+iter_t+1] = positive_modulo(History_Particles_W[it+iter_t+1],lz);

				//	std::cout << History_Particles_U[it+iter_t+1]  << " "  << History_Particles_V[it+iter_t+1]<< " " <<History_Particles_W[it+iter_t+1] << std::endl;
			}

			current_time += dt;
			iter_t += 1;
		}
	}
	else if(SCHM == 5){

		while(current_time < tf-dt){

			time(&cu_time);
			std::cout << "Iteration: " << iter_t <<  " file "<<  current_index_t <<  " Temps: " << current_time <<" Temps de calcul " << cu_time - time_clock<<  std::endl;

			current_index_t =  current_index_t + 1;
			iter_t++;
			
			for(i=0;i<Nb_particules;i++){

				k1_AMC_U[i] =  k2_AMC_U[i];
				k1_AMC_V[i] =  k2_AMC_U[i];
				k1_AMC_W[i] =  k2_AMC_U[i];

				k2_AMC_U[i] =  k3_AMC_U[i];
				k2_AMC_V[i] =  k3_AMC_U[i];
				k2_AMC_W[i] =  k3_AMC_U[i];

				k3_AMC_U[i] =  k4_AMC_U[i];
				k3_AMC_V[i] =  k4_AMC_U[i];
				k3_AMC_W[i] =  k4_AMC_U[i];

				k4_AMC_U[i] =  k5_AMC_U[i];
				k4_AMC_V[i] =  k5_AMC_U[i];
				k4_AMC_W[i] =  k5_AMC_U[i];

			}

			std::cout << "Lecture du champs: " << current_index_t << endl;

			Lire_Flow(Grid_t_U,Grid_t_V,Grid_t_W,nx,ny,nz,FILE_U,FILE_V,FILE_W,FILE_VITESSE,current_index_t);	

			for(i=0;i<Nb_particules;i++){

				it = i*nt;
				k5_AMC_U[i] =  LagrangianInterp(History_Particles_U[it+iter_t],History_Particles_V[it+iter_t] ,History_Particles_W[it+iter_t],Grid_t_U,dx,y_values,dz,nx,ny,nz,0);
				k5_AMC_V[i] =  LagrangianInterp(History_Particles_U[it+iter_t],History_Particles_V[it+iter_t] ,History_Particles_W[it+iter_t],Grid_t_V,dx,y_values,dz,nx,ny,nz,1);
				k5_AMC_W[i] =  LagrangianInterp(History_Particles_U[it+iter_t],History_Particles_V[it+iter_t] ,History_Particles_W[it+iter_t],Grid_t_W,dx,y_values,dz,nx,ny,nz,2);

				History_Particles_U[it+iter_t+1] = History_Particles_U[it+iter_t] + dt*(1901.0/720.0*k5_AMC_U[i]-1387.0/360.0*k4_AMC_U[i]+109.0/30.0*k3_AMC_U[i]-637./360.*k2_AMC_U[i]+251./720.*k1_AMC_U[i]);
				History_Particles_V[it+iter_t+1] = History_Particles_V[it+iter_t] + dt*(1901.0/720.0*k5_AMC_V[i]-1387.0/360.0*k4_AMC_V[i]+109.0/30.0*k3_AMC_V[i]-637./360.*k2_AMC_V[i]+251./720.*k1_AMC_V[i]);
				History_Particles_W[it+iter_t+1] = History_Particles_W[it+iter_t] + dt*(1901.0/720.0*k5_AMC_W[i]-1387.0/360.0*k4_AMC_W[i]+109.0/30.0*k3_AMC_W[i]-637./360.*k2_AMC_W[i]+251./720.*k1_AMC_W[i]);
			

				if(History_Particles_V[it+iter_t+1] < 0)
					History_Particles_V[it+iter_t+1] = 0;

				if(History_Particles_U[it+iter_t+1] >= lx || History_Particles_U[it+iter_t+1] < 0)
					History_Particles_U[it+iter_t+1] = positive_modulo(History_Particles_U[it+iter_t+1],lx);

				if (History_Particles_W[it+iter_t+1] > lz || History_Particles_W[it+iter_t+1] < 0) 
					History_Particles_W[it+iter_t+1] = positive_modulo(History_Particles_W[it+iter_t+1],lz);
			}

				//	std::cout << History_Particles_U[it+iter_t+1]  << " "  << History_Particles_V[it+iter_t+1]<< " " <<History_Particles_W[it+iter_t+1] << std::endl;

		current_time += dt;
		}
	}


	Enregistrement_champs(DIR_Particules_output,History_Particles_U,History_Particles_V,History_Particles_W,y0_init,Nb_particules,Nb_particules_par_ligne,iter_t,nt,t0_win);
	Enregistrement_champs(DIR_Vitesses_output,History_Vel_Particles_U,History_Vel_Particles_V,History_Vel_Particles_W,y0_init,Nb_particules, Nb_particules_par_ligne, iter_t,nt,t0_win);

	// Liberation de la memoire



	free(U);
	free(V);
	free(W);
	
	free(History_Particles_U);
	free(History_Particles_V);
	free(History_Particles_W);
	
	free(History_Vel_Particles_U);
	free(History_Vel_Particles_V);
	free(History_Vel_Particles_W);

	free(Grid_t_U);
	free(Grid_t_V);
	free(Grid_t_W);

	free(Grid_t1_U);
	free(Grid_t1_V);
	free(Grid_t1_W);


	free(Grid_t2_U);
	free(Grid_t2_V);
	free(Grid_t2_W);

	free(Grid_t3_U);
	free(Grid_t3_V);
	free(Grid_t3_W);
	
	free(Grid_t4_U);
	free(Grid_t4_V);
	free(Grid_t4_W);
	
	free(Grid_t5_U);
	free(Grid_t5_V);
	free(Grid_t5_W);
	
	free(x_values);
	free(z_values);
	free(y_values);

	free(Grid_current_t_U);
	free(Grid_current_t_V);
	free(Grid_current_t_W);

	free(Grid_next_t_U);
	free(Grid_next_t_V);
	free(Grid_next_t_W);

	free(k1_AMC_U);
	free(k1_AMC_V);
	free(k1_AMC_W);
	
	free(k2_AMC_U);
	free(k2_AMC_V);
	free(k2_AMC_W);
	
	free(k3_AMC_U);
	free(k3_AMC_V);
	free(k3_AMC_W);
	
	free(k4_AMC_U);
	free(k5_AMC_V);
	free(k4_AMC_W);
	
	free(k5_AMC_U);
	free(k5_AMC_V);
	free(k5_AMC_W);

return 0;
}
