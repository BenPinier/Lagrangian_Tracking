#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string> 
// VARIABLES GLOBALES

double Lx[6];
double Ly[6];
double Lz[6];


double y_1,y_2,y_3,y_4,y_5,y_6;
int x_0,x_1,z_0,z_1,y_a,y_b,find;

double xhi_x;
double xhi_z;

double xd,yd,zd,c00,c01,c10,c11,c0,c1;

unsigned int i_index,j_index,ij_index;
unsigned int i_index2,j_index2,ij_index2;
int nynz,it ;
#include <sstream>

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

// Gestion des fichiers de sorties


void Enregistrement_champs(std::string& FILE_RES_str,double* H_U,double* H_V,double* H_W,int y0_init,int Nb_particules,int Nb_part_per_plane, int iter_t,int nt){

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
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P000"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+".dat";
		else if(fmod(i,Nb_part_per_plane)+1 < 100)	
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P00"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+".dat";
		else if(fmod(i,Nb_part_per_plane)+1 < 1000)	
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P0"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+".dat";
		else
			FILE_RES = FILE_RES_str+"/y_"+patch::to_string(k)+"_P"+patch::to_string(fmod(i,Nb_part_per_plane)+1)+".dat";

			
		ofs.open(FILE_RES.c_str());

		for(int j=0;j<iter_t;j++){
			ofs << patch::to_string(H_U[i*nt+j]);
			ofs << " ";
			ofs << patch::to_string(H_V[i*nt+j]);
			ofs << " ";
			ofs << patch::to_string(H_W[i*nt+j]);
			ofs << std::endl;
		}
	}

}


// Initialisation des vitesses


void init_Particles(int n_p,double* X_1,double* X_2,double*X_3,double*  y_v,int nx,int ny,int nz,int nt,int y0_init,double lx,double lz,double part_per_plane){
	int c,i,j,k;
	i=0;
	j=0;
	k=0;
	c=0;
	for(int i=0;i<65;i++)
		std::cout << y_v[i] << std::endl;


	k = y0_init;
	while(c < n_p){
		i=0;
		while (i <  part_per_plane &&  c < n_p){
			j=0;
			while (j <  part_per_plane &&  c < n_p){
			//	std::cout << c << " " << k << " " << y_v[k] <<std::endl;
				X_2[c*nt] = y_v[k];
				X_1[c*nt] = lx/( part_per_plane *1.1)*(i+1);
				X_3[c*nt]= lz/( part_per_plane *1.1)*(j+1);
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
	if(cas == 1){
		for(int i=0;i<nx;i++){
			i_index = i*ny*(nz);
			for(int j=0;j<ny;j++){
				ij_index = i_index+j*(nz);
				for(int k=0;k<nz;k++){
					U[ij_index+k] = 1. - (y_v[j]-1.)*	(y_v[j]-1.);
					V[ij_index+k] = 0.;
					W[ij_index+k] = 0.5*cos(0.25*i*dx);
					//W[ij_index+k] = k*dz;
				}
			}
		}
	}
}



// Vitesse vers grille

void  Velocity_to_grid(double* Grid,double* U, int nx,int ny,int nz){
	unsigned int i,j,k;

	
	for(i=0;i<nx-1;i++){
	 	i_index = i*ny*(nz);
	 	i_index2= i*ny*(nz-1);
		for(j=0;j<ny;j++){
			ij_index = i_index+j*(nz);
			ij_index2 = i_index2+j*(nz-1);
			for(k=0;k<nz-1;k++){
				Grid[ij_index+k] = U[ij_index2+k];
			}
		}
	}

	for(j=0;j<ny;j++){
		j_index = j*(nz);
		j_index2 = j*(nz-1);
		for( k=0;k<nz-1;k++){
			Grid[nx*ny*(nz)+j_index+k] = U[j_index2+k];
		}
	}


	for(i=0;i<nx-1;i++){
		i_index = i*ny*(nz);
		i_index2= i*ny*(nz-1);
		for(j=0;j<ny;j++){
			ij_index = i_index+j*(nz);
			ij_index2 = i_index2+j*(nz-1);
			Grid[ij_index+nz] = U[ij_index2];
		}
	}

	i_index = (nx-1)*ny*(nz);
	i_index2= (nx-2)*ny*(nz-1);
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

	x_0 = floor(X0/dx);
	x_1 = floor(X0/dx)+1;


	z_0 = floor(X2/dz);
	z_1 = floor(X2/dz)+1;



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

	xd = (X0-x_0*dx)/(x_1*dx-x_0*dx);
	yd = (X1-y_v[y_a])/(y_v[y_b]-y_v[y_a]);
	zd = (X2-z_0*dz)/(z_1*dz-z_0*dz);


//	std::cout << xd << " " << yd << " " <<zd << std::endl;

	c00 = Grid[x_0*ny*nz+y_a*nz+z_0]*(1.-xd) + Grid[x_1*ny*nz+y_a*nz+z_0]*(xd) ;
	c01 = Grid[x_0*ny*nz+y_a*nz+z_1]*(1.-xd) + Grid[x_1*ny*nz+y_a*nz+z_1]*(xd) ;
	c10 = Grid[x_0*ny*nz+y_b*nz+z_0]*(1.-xd) + Grid[x_1*ny*nz+y_b*nz+z_0]*(xd) ;
	c11 = Grid[x_0*ny*nz+y_b*nz+z_1]*(1.-xd) + Grid[x_1*ny*nz+y_b*nz+z_1]*(xd) ;

	c0 = c00*(1.-yd) + c10*yd;
	c1 = c01*(1.-yd) + c11*yd;

	res = c0*(1.-zd) + c1*zd;

	

//	std::cout << res << std::endl;
	return res;

}

// Interpolation Lagrangienne

double PolyLagrange(int i,double xhi){
	double res = 0.;

	if (i== 1){
		res += 6.*xhi;
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
		res += -6.*xhi;
		xhi *= xhi;
		res += 16.*xhi;
		xhi *= xhi;
		res += -xhi;
		xhi *= xhi;
		res += -4.0*xhi;
		xhi *= xhi;
		res += xhi;

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
		res += 12.*xhi;
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
		res += -6.*xhi;
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
		res += 4.*xhi;
		xhi *= xhi*xhi;
		res += -5.*xhi;
		xhi *= xhi*xhi;
		res += xhi;

		res = res / 120.;
	}
	else
		res =0.;

	return res;
}	




double LagrangianInterp(double X0,double X1,double X2,double*  Grid,double dx,double*  y_v,double dz,int nx, int ny,int nz,int index){

	double res = 0.;


	xhi_x = (X0-floor(X0/dx))/dx;
	xhi_z = (X2-floor(X2/dz))/dz;

	int x_index,z_index;

	int i,j,l;

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

	// Si il y a assez de point pour gerer le mur
	if(k >= 2){

		y_1 = y_v[k-2];
		y_2 = y_v[k-1];
		y_3 = y_v[k];
		y_4 = y_v[k+1];
		y_5 = y_v[k+2];
		y_6 = y_v[k+3];


		Lx[0] = PolyLagrange(0,xhi_x);
		Lx[1] = PolyLagrange(1,xhi_x);
		Lx[2] = PolyLagrange(2,xhi_x);
		Lx[3] = PolyLagrange(3,xhi_x);
		Lx[4] = PolyLagrange(4,xhi_x);
		Lx[5] = PolyLagrange(5,xhi_x);


		Lz[0] = PolyLagrange(0,xhi_z);
		Lz[1] = PolyLagrange(1,xhi_z);
		Lz[2] = PolyLagrange(2,xhi_z);
		Lz[3] = PolyLagrange(3,xhi_z);
		Lz[4] = PolyLagrange(4,xhi_z);
		Lz[5] = PolyLagrange(5,xhi_z);

		Ly[0] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_1-y_2)*(y_1-y_3)*(y_1-y_4)*(y_1-y_5)*(y_1-y_6));
		Ly[1] = (X1-y_1)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_2-y_1)*(y_2-y_3)*(y_2-y_4)*(y_2-y_5)*(y_2-y_6));
		Ly[2] = (X1-y_2)*(X1-y_1)*(X1-y_4)*(X1-y_5)*(X1-y_6)/((y_3-y_2)*(y_3-y_1)*(y_3-y_4)*(y_3-y_5)*(y_3-y_6));
		Ly[3] = (X1-y_2)*(X1-y_3)*(X1-y_1)*(X1-y_5)*(X1-y_6)/((y_4-y_2)*(y_4-y_3)*(y_4-y_1)*(y_4-y_5)*(y_4-y_6));
		Ly[4] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_1)*(X1-y_6)/((y_5-y_2)*(y_5-y_3)*(y_5-y_4)*(y_5-y_1)*(y_5-y_6));
		Ly[5] = (X1-y_2)*(X1-y_3)*(X1-y_4)*(X1-y_5)*(X1-y_1)/((y_6-y_2)*(y_6-y_3)*(y_6-y_4)*(y_6-y_5)*(y_6-y_1));

		for(i=0;i<6;i++){
			x_index = fmod(floor(X0/dx),(nx));
			for(j=0;j<6;j++){
				z_index = fmod(floor(X2/dz),(nz));
				for(l=0;l<6;l++){
					res +=  Grid[x_index*nynz+(k+l-2)*nz+z_index]*Lx[i]*Ly[l]*Lz[j];
				}
			}
		}
	}
	else
		res = TrilinearInterp(X0,X1,X2,Grid,dx,y_v,dz,nx,ny,nz,index);

//	std::cout << res << " " << 1.-(y_v[1]-1.)*(y_v[1]-1.) << " " << floor(X0/dx) << " " << floor(X2/dz) << " " << y_v[k]  << " "<< Grid[int(floor(X0/dx))*ny*(nz)+1*(nz)+int(floor(X2/dz))] <<std::endl<<std::endl;
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
	 
	Chargement_parametres( FILE_arg,&Nb_particules,&	Nb_particules_par_ligne,&tf,&y0_init,&dt,&dt_data,&nb_iter_win,&lx,&ly,&lz,&nx,&ny,&nz,FILE_VITESSE,DIR_Particules_output,DIR_Vitesses_output);
	Affiche_Param(&Nb_particules,&Nb_particules_par_ligne,&tf,&y0_init,&dt,&dt_data,&nb_iter_win,&lx,&ly,&lz,&nx,&ny,&nz,FILE_VITESSE,DIR_Particules_output,DIR_Vitesses_output);


	std::string Y_file = FILE_VITESSE+"/yp.dat";

	 dx = lx/nx;
	 dz = lz/nz;

	 double nt = floor(tf/dt);


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



	std::cout <<  Nb_particules << std::endl;


	std::cout <<floor(tf/dt+2) << std::endl;

	History_Particles_U = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Particles_V = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Particles_W = (double*)malloc(Nb_particules*nt*sizeof(double));

	History_Vel_Particles_U = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Vel_Particles_V = (double*)malloc(Nb_particules*nt*sizeof(double));
	History_Vel_Particles_W = (double*)malloc(Nb_particules*nt*sizeof(double));


	
	Grid_current_t_U= (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_current_t_V= (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_current_t_W= (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));

	Grid_next_t_U= (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_next_t_V= (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_next_t_W= (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	
	/*Memoire_Grille_3D(History_Particles, Nb_particules,3,floor(tf/dt+2));
	Memoire_Grille_3D(History_Vel_Particles, Nb_particules,3,floor(tf/dt+2));*/

	Grid_t_U = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_t_V = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_t_W = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));

	Grid_t1_U = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_t1_V = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_t1_W = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));

	Grid_t2_U = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_t2_V = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));
	Grid_t2_W = (double*)malloc((nx+1)*ny*(nz+1)*sizeof(double));

/*	Memoire_Grille_4D(Grid_t,nx+1,ny,nz+1,3);
	Memoire_Grille_4D(Grid_t1,nx+1,ny,nz+1,3);
	Memoire_Grille_4D(Grid_t2,nx+1,ny,nz+1,3);*/


		
/*	Memoire_Grille_4D(Grid_current_t,nx+1,ny,nz+1,3);
	Memoire_Grille_4D(Grid_next_t,nx+1,ny,nz+1,3);
*/	

	double* x_values;
	double* z_values;
	double* y_values;

	x_values = (double*)malloc((nx+1)*sizeof(double));
	z_values = (double*)malloc((nz+1)*sizeof(double));
	y_values = (double*)malloc(ny*sizeof(double));

	
	for(int i=0;i<=nx;i++)
		x_values[i] = i*dx;


	for(int i=0;i<=nz;i++)
		z_values[i] = i*dz;

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
	
/*	Memoire_Grille_3D(U,nx+1,ny,nz+1);
	Memoire_Grille_3D(V,nx+1,ny,nz+1);
	Memoire_Grille_3D(W,nx+1,ny,nz+1);
*/
	int current_index_t = 0;
	double current_time = 0;
	double ratio_time = 0;
	double  k1[3];
	double  k2[3];
	double  k3[3];
	double  k4[3];
	double  yn[3];
	int iter_t = 0;

	double dt2 = 0.5*dt;
	double dt6 = dt/6.0;

	time_t time_clock;
	time(&time_clock);
	double t2;

	double y_00,y_01,y_02;
	std:: cout << "Itilialisation des particules ... " << std::endl;
	init_Particles(Nb_particules,History_Particles_U,History_Particles_V,History_Particles_W,y_values,nx,ny,nz,nt,y0_init,lx,lz,Nb_particules_par_ligne);
	std:: cout << "Itilialisation des particules ok " << std::endl;

	current_index_t += 1;


	Update_Vel( U,V,W,1,nx,ny,nz,dx,dz,y_values,0);		
	Velocity_to_grid(Grid_current_t_U,U,nx+1,ny,nz+1);
	Velocity_to_grid(Grid_current_t_V,V,nx+1,ny,nz+1);
	Velocity_to_grid(Grid_current_t_W,W,nx+1,ny,nz+1);

	Update_Vel( U,V,W,1,nx,ny,nz,dx,dz,y_values,dt);
	Velocity_to_grid(Grid_next_t_U,U,nx+1,ny,nz+1);
	Velocity_to_grid(Grid_next_t_V,V,nx+1,ny,nz+1);
	Velocity_to_grid(Grid_next_t_W,W,nx+1,ny,nz+1);


	double ratio_time_0,ratio_time_1,ratio_time_2;

	int nynz = ny*(nz+1);
	int dnt = 3*nt;

	int i_index,j_index,ij_index;

	while(current_time < tf-dt){


		std::cout << "Iteration: " << iter_t << " Temps: " << current_time << std::endl;

		ratio_time_0 = fmod(iter_t,nb_iter_win)/nb_iter_win;
		ratio_time_1 = fmod(iter_t+0.5,nb_iter_win)/nb_iter_win;


		for(i=0;i<=nx;i++){
			i_index = i*nynz;
			for(j=0;j<ny;j++){
				ij_index = i*nynz+j*(nz+1);
				for(k=0;k<=nz;k++){
						Grid_t_U[ij_index+k] =  (1.-ratio_time_0)*Grid_current_t_U[ij_index+k] + ratio_time_0*Grid_next_t_U[ij_index+k];
						Grid_t_V[ij_index+k] =  (1.-ratio_time_0)*Grid_current_t_V[ij_index+k] + ratio_time_0*Grid_next_t_V[ij_index+k];
						Grid_t_W[ij_index+k] =  (1.-ratio_time_0)*Grid_current_t_W[ij_index+k] + ratio_time_0*Grid_next_t_W[ij_index+k];

/*if(k==0)
							std::cout << Grid_t_W[ij_index+k]  << std::endl;*/
																		
						Grid_t1_U[ij_index+k]  =  (1.-ratio_time_1)*Grid_current_t_U[ij_index+k]  + ratio_time_1*Grid_next_t_U[ij_index+k];
						Grid_t1_V[ij_index+k]  =  (1.-ratio_time_1)*Grid_current_t_V[ij_index+k]  + ratio_time_1*Grid_next_t_V[ij_index+k];
						Grid_t1_W[ij_index+k]  =  (1.-ratio_time_1)*Grid_current_t_W[ij_index+k]  + ratio_time_1*Grid_next_t_W[ij_index+k];
					
				}
			}
		}


		//Next temporal window
		if(fmod(iter_t+1,nb_iter_win) == 0){
	//		std::cout << "Nouvelle fenetre: " <<current_index_t+1<< std::endl;
			current_index_t=  current_index_t + 1;

			for(i=0;i<=nx;i++){
				i_index = i*nynz;
				for(j=0;j<ny;j++){
					ij_index = i*nynz+j*(nz+1);
					for(k=0;k<=nz;k++){
						Grid_current_t_U[ij_index+k] = Grid_next_t_U[ij_index+k] ;
						Grid_current_t_V[ij_index+k] = Grid_next_t_V[ij_index+k] ;
						Grid_current_t_W[ij_index+k] = Grid_next_t_W[ij_index+k];
					}
				}
			}

			Update_Vel( U,V,W,1,nx,ny,nz,dx,dz,y_values,current_time);
			Velocity_to_grid(Grid_next_t_U,U,nx+1,ny,nz+1);
			Velocity_to_grid(Grid_next_t_V,V,nx+1,ny,nz+1);
			Velocity_to_grid(Grid_next_t_W,W,nx+1,ny,nz+1);
		}


		ratio_time_2 = fmod(iter_t+1,nb_iter_win)/nb_iter_win;


		for(i=0;i<=nx;i++){
			i_index = i*nynz;
			for(j=0;j<ny;j++){
				ij_index = i*nynz+j*(nz+1);
				for(k=0;k<=nz;k++){
					Grid_t2_U[ij_index+k] =  (1.-ratio_time_2)*Grid_current_t_U[ij_index+k] + ratio_time_2*Grid_next_t_U[ij_index+k];
					Grid_t2_V[ij_index+k] =  (1.-ratio_time_2)*Grid_current_t_V[ij_index+k] + ratio_time_2*Grid_next_t_V[ij_index+k];
					Grid_t2_W[ij_index+k] =  (1.-ratio_time_2)*Grid_current_t_W[ij_index+k] + ratio_time_2*Grid_next_t_W[ij_index+k];
				}
			}
		}


		for(i=0;i<Nb_particules;i++){

			it = i*nt;

			y_00 = History_Particles_U[it+iter_t];
			y_01 = History_Particles_V[it+iter_t];
			y_02 = History_Particles_W[it+iter_t];

		//	std::cout << y_00 << " "  << y_01<< " " << y_02<< std::endl;


			k1[0] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_U,dx,y_values,dz,nx+1,ny,nz+1,0);
			k1[1] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_V,dx,y_values,dz,nx+1,ny,nz+1,1);
			k1[2] =  LagrangianInterp(y_00,y_01 ,y_02 ,Grid_t_W,dx,y_values,dz,nx+1,ny,nz+1,2);

//			std::cout << "Fin des interpolations 1" << std::endl;

			k2[0] =  LagrangianInterp(fmod(y_00+dt2*k1[0],lx),fmax(0,y_01+dt2*k1[1]),fmod(y_02+dt2*k1[2],lz),Grid_t1_U,dx,y_values,dz,nx+1,ny,nz+1,0);
			k2[1] =  LagrangianInterp(fmod(y_00+dt2*k1[0],lx),fmax(0,y_01+dt2*k1[1]),fmod(y_02+dt2*k1[2],lz),Grid_t1_V,dx,y_values,dz,nx+1,ny,nz+1,1);
			k2[2] =  LagrangianInterp(fmod(y_00+dt2*k1[0],lx),fmax(0,y_01+dt2*k1[1]),fmod(y_02+dt2*k1[2],lz),Grid_t1_W,dx,y_values,dz,nx+1,ny,nz+1,2);

//			std::cout << "Fin des interpolations 2" << std::endl;

			k3[0] =  LagrangianInterp(fmod(y_00+dt2*k2[0],lx),fmax(0,y_01+dt2*k2[1]),fmod(y_02+dt2*k2[2],lz),Grid_t1_U,dx,y_values,dz,nx+1,ny,nz+1,0);
			k3[1] =  LagrangianInterp(fmod(y_00+dt2*k2[0],lx),fmax(0,y_01+dt2*k2[1]),fmod(y_02+dt2*k2[2],lz),Grid_t1_V,dx,y_values,dz,nx+1,ny,nz+1,1);
			k3[2] =  LagrangianInterp(fmod(y_00+dt2*k2[0],lx),fmax(0,y_01+dt2*k2[1]),fmod(y_02+dt2*k2[2],lz),Grid_t1_W,dx,y_values,dz,nx+1,ny,nz+1,2);

//			std::cout << "Fin des interpolations 3" << std::endl;

			k4[0] =  LagrangianInterp(fmod(y_00+dt*k3[0],lx),fmax(0,y_01+dt*k3[1]),fmod(y_02+dt*k3[2],lz),Grid_t2_U,dx,y_values,dz,nx+1,ny,nz+1,0);
			k4[1] =  LagrangianInterp(fmod(y_00+dt*k3[0],lx),fmax(0,y_01+dt*k3[1]),fmod(y_02+dt*k3[2],lz),Grid_t2_V,dx,y_values,dz,nx+1,ny,nz+1,1);
			k4[2] =  LagrangianInterp(fmod(y_00+dt*k3[0],lx),fmax(0,y_01+dt*k3[1]),fmod(y_02+dt*k3[2],lz),Grid_t2_W,dx,y_values,dz,nx+1,ny,nz+1,2);


//			std::cout << "Fin des interpolations 4" << std::endl;

			History_Vel_Particles_U[it+iter_t] =(k1[0]);
			History_Vel_Particles_V[it+iter_t] =(k1[1]);
			History_Vel_Particles_W[it+iter_t] =(k1[2]);



				
			History_Particles_U[it+iter_t+1] = y_00 + dt6*(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0]);
			History_Particles_V[it+iter_t+1] = y_01 + dt6*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1]);
			History_Particles_W[it+iter_t+1] = y_02 + dt6*(k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2]);

			if(History_Particles_V[it+iter_t+1] < 0)
				History_Particles_V[it+iter_t+1] = 0;



			if(History_Particles_U[it+iter_t+1] >= lx || History_Particles_U[it+iter_t+1] < 0)
				History_Particles_U[it+iter_t+1] = fmod(History_Particles_U[it+iter_t+1],lx);

			if (History_Particles_W[it+iter_t+1] > lz || History_Particles_W[it+iter_t+1] < 0) 
				History_Particles_W[it+iter_t+1] = fmod(History_Particles_W[it+iter_t+1],lz);

			//	std::cout << History_Particles_U[it+iter_t+1]  << " "  << History_Particles_V[it+iter_t+1]<< " " <<History_Particles_W[it+iter_t+1] << std::endl;
		}

		current_time += dt;
		iter_t += 1;
	}


	Enregistrement_champs(DIR_Particules_output,History_Particles_U,History_Particles_V,History_Particles_W,y0_init,Nb_particules,Nb_particules_par_ligne,iter_t,nt);
	Enregistrement_champs(DIR_Vitesses_output,History_Vel_Particles_U,History_Vel_Particles_V,History_Vel_Particles_W,y0_init,Nb_particules, Nb_particules_par_ligne, iter_t,nt);

	// Liberation de la memoire



	
/*	free_memoire_3D(U,nx,ny,nz);
	free_memoire_3D(V,nx,ny,nz);
	free_memoire_3D(W,nx,ny,nz);
	free_memoire_3D(History_Particles, Nb_particules,3,floor((tf)/dt+2));
	free_memoire_3D(History_Vel_Particles, Nb_particules,3,floor((tf)/dt+2));
	free_memoire_4D(Grid_current_t,nx,ny,nz,3);
	free_memoire_4D(Grid_next_t,nx,ny,nz,3);*/



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
	
	free(x_values);
	free(z_values);
	free(y_values);
/*
	free(Grid_current_t_U);
	free(Grid_current_t_V);
	free(Grid_current_t_W);

	free(Grid_next_t_U);
	free(Grid_next_t_V);
	free(Grid_next_t_W);
*/
return 0;
}
