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
#include <vector>

#include <sstream>

template < typename T> 
T average(T A,int c){

	T res = 0;
	for(int i=0;i<c;i++)
		res += A[i]/c;

return res;

} 

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

std::vector<double> split_to_double(const std::string& s, char seperator)
{
   std::vector<double> output;

    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );

        output.push_back(atof(substring.c_str()));

        prev_pos = ++pos;
    }

    output.push_back(atof(s.substr(prev_pos, pos-prev_pos).c_str())); // Last word

    return output;
}

void Chargement_parametres(std::string& FILE_arg,int* Nb_particules,int* Nb_particules_par_temps,double* tf,int* y0_init,int* ny,double* dt_data,std::string& FILE_POS,std::string& DIR_Particules_output){

	std::ifstream file(FILE_arg.c_str());
	
	std::string str;

	std::getline(file, str);
	*Nb_particules = atoi(str.c_str());

	std::getline(file, str);
	*Nb_particules_par_temps = atoi(str.c_str());

	std::getline(file, str);
	*tf = atof(str.c_str());

	std::getline(file, str);
	*y0_init = atoi(str.c_str());

	std::getline(file, str);
	*ny = atoi(str.c_str());

	std::getline(file, str);
	*dt_data = atof(str.c_str());



	
	std::getline(file, FILE_POS);

//	std::getline(file, DIR_Particules_output);

}





using namespace::std;


int main(int argc,char** argv){



	int i,j,k,l;

	std::string FILE_arg = argv[1];
	std::cout << FILE_arg << std::endl;

	 int Nb_particules;
	 int Nb_particules_par_ligne;
	 double tf; // temps final
	 int y0_init;
	 double dt;
	 double dt_data;
	 int nb_iter_win;
	 int ny;

	std::string FILE_VITESSE,DIR_Vitesses_output;

	std::ifstream cu_file;

	std::ofstream cu_file_out; 

	std::string cu_file_str;


	Chargement_parametres( FILE_arg,&Nb_particules,&	Nb_particules_par_ligne,&tf,&y0_init,&ny,&dt_data,FILE_VITESSE,DIR_Vitesses_output);


	std::cout << Nb_particules << " " << Nb_particules_par_ligne << " " << tf << " "<< y0_init << " " << ny << " " << dt_data << " " << FILE_VITESSE << std::endl;

	double* TensorAxx;
	double* TensorAyy;
	double* TensorAzz;
	double* TensorAxy;
	double* TensorAxz;
	double* TensorAyz;

	std::string cu_FILE;
					

	TensorAxx = (double*)malloc(ny*sizeof(double));
	TensorAyy = (double*)malloc(ny*sizeof(double));
	TensorAzz = (double*)malloc(ny*sizeof(double));
	TensorAxy = (double*)malloc(ny*sizeof(double));
	TensorAxz = (double*)malloc(ny*sizeof(double));
	TensorAyz = (double*)malloc(ny*sizeof(double));

	std::vector<double> cu_Pos,Pos0;
	double PrevPos[3];
	std::string last_line,line;

	int Array_T[4];

	int count;

	Array_T[0] = 0;
	Array_T[1] = 400;
	Array_T[2] = 800;
	Array_T[3] = 1200;


	for(i=1;i<ny;i++){
		TensorAxx[i] = 0.;
		TensorAyy[i] = 0.;
		TensorAzz[i] = 0.;
		TensorAxy[i] = 0.;
		TensorAxz[i] = 0.;
		TensorAyz[i] = 0.;
	}

	std::cout << ny << " " << std::endl;

	double U_avg;
	double* cU;
	double* cU0;
	double* cW;
	double* cW0;
	double* cV;
	double* cV0;
	double driftU,driftV,driftW;

	double lx = 12.566370614359172;
	double lz = 4.1887902047863905;

	double decalage_x,decalage_z;
	
	cU = (double*)malloc(Nb_particules_par_ligne*4*sizeof(double));
	cU0 = (double*)malloc(Nb_particules_par_ligne*4*sizeof(double));

	cW = (double*)malloc(Nb_particules_par_ligne*4*sizeof(double));
	cW0 = (double*)malloc(Nb_particules_par_ligne*4*sizeof(double));

	cV = (double*)malloc(Nb_particules_par_ligne*4*sizeof(double));
	cV0 = (double*)malloc(Nb_particules_par_ligne*4*sizeof(double));


	int lmax = 1500;

	for(i=1;i<ny;i++){
		std::cout<< i << std::endl;
		count = 0;
		driftU = 0;
		driftV = 0;
		driftW = 0;
		for(j=0;j<4;j++){

		
			for(k=1;k<=Nb_particules_par_ligne;k+=10){
				if(k<10)
					cu_file_str = FILE_VITESSE+"/y_"+patch::to_string(i)+"_P000"+patch::to_string(k)+"_t0_"+patch::to_string(Array_T[j])+".dat";
				else if(k<100)
					cu_file_str = FILE_VITESSE+"/y_"+patch::to_string(i)+"_P00"+patch::to_string(k)+"_t0_"+patch::to_string(Array_T[j])+".dat";
				else if(k<1000)
					cu_file_str = FILE_VITESSE+"/y_"+patch::to_string(i)+"_P0"+patch::to_string(k)+"_t0_"+patch::to_string(Array_T[j])+".dat";
				else if(k<10000)
					cu_file_str = FILE_VITESSE+"/y_"+patch::to_string(i)+"_P"+patch::to_string(k)+"_t0_"+patch::to_string(Array_T[j])+".dat";

				cu_file.open(cu_file_str.c_str());

				l=0;
				getline (cu_file,line);
				Pos0 = split_to_double(line,' ');
				cU0[count] = Pos0[0];
				cW0[count] = Pos0[1];
				cV0[count] = Pos0[2];
				decalage_x = 0.;
				decalage_z = 0.;

				while (getline (cu_file,line) && l<=lmax)
				{
						cu_Pos = split_to_double(line,' ');


										
					 	TensorAxx[i] += (-Pos0[0]+cu_Pos[0])*(-Pos0[0]+cu_Pos[0]);
						TensorAzz[i] += (-Pos0[1]+cu_Pos[1])*(-Pos0[1]+cu_Pos[1]);
						TensorAyy[i] += (-Pos0[2]+cu_Pos[2]+decalage_z)*(-Pos0[2]+cu_Pos[2]+decalage_z);
						TensorAxz[i] += (-Pos0[0]+cu_Pos[0])*(-Pos0[1]+cu_Pos[1]);
						TensorAxy[i] += (-Pos0[0]+cu_Pos[0])*(-Pos0[2]+cu_Pos[2]);
						TensorAyz[i] += (-Pos0[1]+cu_Pos[1])*(-Pos0[2]+cu_Pos[2]);


						Pos0[0] = cu_Pos[0];
						Pos0[1] = cu_Pos[1];
						Pos0[2] = cu_Pos[2];
					
					l++;
				}
				cu_file.close();
				count++;
				Pos0.clear();
				cu_Pos.clear();

				
			}
		}
     	std::cout << driftU/count << " " <<driftV/count << " "<<driftW /count << endl;
		driftU /= count;
		driftW /=count;
		driftV /=count;
		driftW;

		
/*
		for(k=0;k<count;k++){
			TensorAxx[i] += (-cU0[k]+cU[k]-driftU)*(-cU0[k]+cU[k]-driftU)/count;
			TensorAzz[i] += (-cW0[k]+cW[k]-driftW)*(-cW0[k]+cW[k]-driftW)/count;
			TensorAyy[i] += (-cV0[k]+cV[k]-driftV)*(-cV0[k]+cV[k]-driftV)/count;

			TensorAxz[i] += (-cU0[k]+cU[k]-driftU)*(-cW0[k]+cW[k]-driftW)/count;
			TensorAxy[i] += (-cU0[k]+cU[k]-driftU)*(-cV0[k]+cV[k]-driftV)/count;
			TensorAyz[i] += (-cW0[k]+cW[k]-driftW)*(-cV0[k]+cV[k]-driftV)/count;
		}*/
	}


	for(i=1;i<ny;i++){
		TensorAxx[i] /= (double)count;
		TensorAyy[i] /= (double)count;
		TensorAzz[i] /= (double)count;
		TensorAxy[i] /= (double)count;
		TensorAxz[i] /= (double)count;
		TensorAyz[i] /= (double)count;
	}



	cu_file_out.open("Results/TensorA_Vel2_Re180.dat", std::ofstream::out );

	 for(i=0;i<ny;i++){
		cu_file_out << TensorAxx[i]<< " ";
		cu_file_out << TensorAyy[i]<< " ";
		cu_file_out << TensorAzz[i]<< " ";
		cu_file_out << TensorAxy[i]<< " ";
		cu_file_out << TensorAxz[i]<< " ";
		cu_file_out << TensorAyz[i]<< " ";
		cu_file_out << std::endl;
	}

	cu_file_out.close();


	free(TensorAxx);
	free(TensorAyy);
	free(TensorAzz);
	free(TensorAxz);
	free(TensorAxy);
	free(TensorAyz);
						

return 0;

}
